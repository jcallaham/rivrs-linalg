# Data Model: Workspace Reuse & Per-Supernode Allocation Optimization

**Feature**: 021-workspace-reuse
**Date**: 2026-02-22

## New Entities

### FactorizationWorkspace

Pre-allocated reusable buffers for the per-supernode factorization loop. Sized to `max_front_size` from symbolic analysis. One instance per sequential factorization, or one per rayon worker thread for parallel factorization.

**Fields:**
- `frontal_data: Mat<f64>` — Dense buffer, max_front x max_front. Reused as each supernode's frontal matrix.
- `frontal_row_indices: Vec<usize>` — Buffer for frontal matrix row indices, capacity max_front.
- `delayed_cols_buf: Vec<usize>` — Buffer for collecting delayed columns from children, capacity bounded by max fully-summed columns.
- `global_to_local: Vec<usize>` — Global-to-local row index mapping, length n (matrix dimension). Folded in from the separate `shared_g2l` / `G2L_BUF` variables per plan decision D4.

**Lifecycle:**
1. Allocated once before the factorization loop (or lazily on first use for thread-local)
2. Before each supernode: zero the m x m lower triangle subregion of `frontal_data`, clear `frontal_row_indices` and `delayed_cols_buf`
3. During supernode processing: used as workspace, contents are temporary
4. After supernode: outputs (FrontFactors, ContributionBlock) are extracted into owned allocations; workspace buffers are NOT freed
5. After factorization loop completes: workspace is dropped (sequential) or returned to thread-local storage (parallel)

**Invariant:** The workspace is never shared between concurrent supernodes. In sequential mode, it is passed as `&mut`. In parallel mode, each rayon worker has its own instance via `Cell`-based move semantics.

### AptpKernelWorkspace

Pre-allocated reusable buffers for the BLAS-3 inner loop inside `factor_inner` and `two_level_factor`. Sized to `inner_block_size` (default 32) and `max_front_size`.

**Fields:**
- `backup_data: Mat<f64>` — Block backup for restore-on-failure, max_front x inner_block_size
- `l11_temp: Mat<f64>` — Copy of L11 block for TRSM aliasing avoidance, inner_block_size x inner_block_size
- `ld_workspace: Mat<f64>` — L*D product workspace for update_trailing/cross_terms, max_front x inner_block_size
- `copy_workspace: Mat<f64>` — Copy buffer for L21/L_panel aliasing avoidance, max_front x inner_block_size

**Lifecycle:**
1. Allocated once at the start of `aptp_factor_in_place` (or passed in from FactorizationWorkspace)
2. Reused across all block iterations within a single kernel call
3. Contents overwritten each block iteration — no zeroing needed between blocks
4. Freed (or returned to workspace) when kernel call completes

**Relationship to FactorizationWorkspace:** AptpKernelWorkspace can be a sub-struct of FactorizationWorkspace or allocated independently within `aptp_factor_in_place`. The key constraint is that dimensions depend on `max_front_size` (known from symbolic) and `inner_block_size` (a configuration constant).

## Modified Entities

### FrontalMatrix (existing)

**Current:** Owns `data: Mat<f64>` and `row_indices: Vec<usize>`, allocated fresh per supernode.

**After change:** Fields become references into (or are populated from) FactorizationWorkspace buffers. The FrontalMatrix struct may become a lightweight view/handle rather than an owning container.

**Key constraint:** `scatter_original_entries_multi`, `extend_add`, and `aptp_factor_in_place` all take `&mut FrontalMatrix` or `MatMut` — they operate on the workspace buffer in place.

### ContributionBlock (existing)

**Current:** Owns `data: Mat<f64>` (copied from frontal trailing submatrix) and `row_indices: Vec<usize>`.

**After change (US2):** The copy in `extract_contribution` may be eliminated or optimized. Since the frontal workspace is reused, the contribution block must still own its data — but the extraction can potentially be done via ownership transfer of a subregion rather than element-by-element copy.

### factor_single_supernode (existing function)

**Current signature (conceptual):**
```
fn factor_single_supernode(
    sn_idx, supernode, child_contributions,
    matrix, perm_fwd, perm_inv, options, scaling,
    global_to_local: &mut Vec<usize>,
) -> Result<SupernodeResult>
```

**After change:** Gains a `workspace: &mut FactorizationWorkspace` parameter.

### aptp_factor_in_place (existing function)

**Current:** Allocates col_order, d, pivot_log, panel_perm_buf, row_perm_buf, col_order_buf, and per-block backup/copies internally.

**After change:** May accept an `&mut AptpKernelWorkspace` parameter for the per-block temporary buffers (backup, l11_copy, ld workspace). The output allocations (col_order, d, pivot_log) remain per-call since they are returned.

## Unchanged Entities

- **FrontFactors**: No changes. Per-supernode output stored in AptpNumeric.
- **AptpSymbolic**: No changes. Provides max_front_size for workspace sizing.
- **AptpNumeric**: No changes to stored data. May change how factorization loop is orchestrated.
- **SparseLDLT**: No changes to public API. Workspace is internal.
- **MixedDiagonal**: No changes. Per-supernode output (d11 in FrontFactors).

## Entity Relationships

```
AptpSymbolic ---provides max_front_size---> FactorizationWorkspace (sizing)
                                               |
                                               | contains
                                               v
                                          AptpKernelWorkspace
                                               |
factor_single_supernode ---borrows &mut---> FactorizationWorkspace
    |                                          |
    | calls                                    | provides buffers to
    v                                          v
aptp_factor_in_place ---borrows &mut---> AptpKernelWorkspace (or internal)
    |
    | returns owned
    v
AptpFactorResult (col_order, d, pivot_log) ---> FrontFactors (per-supernode output)
```

## Memory Budget

For a matrix with `max_front_size = M` and `inner_block_size = IB` (default 32):

| Buffer | Size (bytes) | Typical (M=1000) |
|--------|-------------|-------------------|
| frontal_data | 8 * M^2 | 8 MB |
| frontal_row_indices | 8 * M | 8 KB |
| delayed_cols_buf | 8 * M | 8 KB |
| backup_data | 8 * M * IB | 256 KB |
| l11_temp | 8 * IB^2 | 8 KB |
| ld_workspace | 8 * M * IB | 256 KB |
| copy_workspace | 8 * M * IB | 256 KB |
| **Total** | **~8 * M^2 + overhead** | **~8.8 MB** |

For parallel factorization with T threads: total workspace is T * (8 * M^2 + overhead). With M=1000 and T=8, this is ~70 MB — acceptable given that the factor storage itself is typically much larger.
