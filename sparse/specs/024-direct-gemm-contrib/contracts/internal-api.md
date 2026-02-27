# Internal API Contracts: Direct GEMM into Contribution Buffer

**Feature**: 024-direct-gemm-contrib
**Date**: 2026-02-24 (revised)

This feature modifies internal (`pub(crate)`) APIs only. No public API changes.

## Modified Contracts

### C1: FactorizationWorkspace::new

**Current**: `fn new(max_front: usize, n: usize) -> Self`
**Modified**: Same signature. Now also allocates `contrib_buffer: Mat<f64>` sized to `max_front × max_front`.

**Postcondition added**: `self.contrib_buffer.nrows() >= max_front && self.contrib_buffer.ncols() >= max_front`

### C2: FactorizationWorkspace::ensure_capacity

**Current**: `fn ensure_capacity(&mut self, max_front: usize, n: usize)`
**Modified**: Also grows `contrib_buffer` if needed.

**Postcondition added**: `self.contrib_buffer.nrows() >= max_front`

### C3: update_trailing (restricted region)

**Current**: Computes symmetric rank-k update on full trailing submatrix `A[ts..m, ts..m]` using lower-triangular GEMM.

**Modified**: Restricted to FS region + cross-terms only:
- Lower-triangular GEMM on `A[ts..p, ts..p]` (FS×FS, region 1)
- Rectangular GEMM on `A[p..m, ts..p]` (NFS×FS cross-term, region 2)
- **Skips** `A[p..m, p..m]` (NFS×NFS, region 3)

**New parameter**: `num_fully_summed` (the boundary `p`) to define the FS/NFS split.

**Precondition**: `ts <= p <= m`

**Postcondition**: `A[p..m, p..m]` is unchanged (retains original assembled values, no Schur complement updates applied).

### C4: Deferred contribution GEMM (new function)

**New**: Computes the NFS×NFS Schur complement in a single GEMM after the blocking loop.

```
fn compute_contribution_gemm(
    frontal_data: &Mat<f64>,        // workspace with L21, assembled NFS×NFS
    num_fully_summed: usize,        // p: boundary between FS and NFS
    num_eliminated: usize,          // ne: columns successfully eliminated
    m: usize,                       // frontal matrix dimension
    d: &MixedDiagonal,             // diagonal from blocking loop
    contrib_buffer: MatMut<f64>,    // output: pre-allocated contribution buffer
    par: Par,                       // parallelism control
)
```

**Precondition**:
- `frontal_data[p..m, 0..ne]` contains L21 for NFS rows (from per-block TRSM)
- `frontal_data[p..m, p..m]` contains original assembled values (untouched by blocking loop)
- `contrib_buffer` is sized >= `(m - p) × (m - p)`
- `contrib_buffer` contents on entry are irrelevant (will be overwritten)

**Behavior**:
1. Copies assembled NFS×NFS values from `frontal_data[p..m, p..m]` into `contrib_buffer[0..nfs, 0..nfs]` (lower triangle)
2. Applies rank-`ne` symmetric update in-place: `contrib_buffer -= L21_NFS * D * L21_NFS^T` (using `alpha=-1, beta=1`)

**Postcondition**:
- `contrib_buffer[0..nfs, 0..nfs]` contains `assembled_NFS_NFS - L21_NFS * D * L21_NFS^T` (lower triangle)
- where `nfs = m - p`

### C5: extract_contribution (index-only + delayed copy)

**Current**: Allocates `Mat::zeros(size, size)`, copies entire trailing submatrix from workspace, builds row indices.

**Modified**: The NFS×NFS portion is already in `contrib_buffer` (from C4). This function now:
1. Copies the small delayed portion (rows `ne..p`) and cross-terms from the workspace into the buffer (only when `ne < p`, i.e., delayed columns exist)
2. Builds `row_indices` and `num_delayed`
3. Moves `contrib_buffer` into the returned `ContributionBlock`

**New parameter**: `contrib_buffer: Mat<f64>` (moved in, already containing NFS×NFS data)

**Postcondition**: No per-supernode allocation. Returned `ContributionBlock.data` is the moved-in buffer.

### C6: extend_add (sequential path — new return value)

**Current**: `fn extend_add(parent: &mut FrontalMatrix<'_>, child: &ContributionBlock, global_to_local: &[usize])`

**Modified**: Takes ownership of child contribution and returns the consumed buffer for recycling:
```
fn extend_add(
    parent: &mut FrontalMatrix<'_>,
    child: ContributionBlock,          // Takes ownership (was &)
    global_to_local: &[usize],
) -> Mat<f64>                          // Returns consumed data buffer for recycling
```

**Postcondition**: Returned `Mat<f64>` can be reused as `workspace.contrib_buffer`.

### C7: extend_add_mapped (sequential path — new return value)

**Current**: `fn extend_add_mapped(parent: &mut FrontalMatrix<'_>, child: &ContributionBlock, ea_row_map: &[u32])`

**Modified**: Same ownership transfer as C6.

### C8: factor_single_supernode (deferred GEMM integration)

**Current**: Calls `aptp_factor_in_place` then `extract_contribution`.

**Modified**: Between factorization and extraction, calls the deferred contribution GEMM (C4) to write NFS×NFS into `workspace.contrib_buffer`. Then calls the simplified `extract_contribution` (C5) which adds delayed-column data and moves the buffer into `ContributionBlock`.

## Unchanged Contracts

- `aptp_factor_in_place` — dispatches to `factor_inner` / `two_level_factor` / `tpp_factor_as_primary`. Returns `AptpFactorResult` with `num_eliminated`. The deferred GEMM is called after this returns, not inside it.
- `update_cross_terms` — handles failed-column updates within the FS region. Unaffected by NFS deferral.
- `scatter_original_entries_multi` — scatters original matrix entries into frontal workspace. Unaffected.
- `AptpNumeric::factor` — orchestration changes (swap logic), no signature changes
- `SparseLDLT` public API — no changes
- `aptp_solve` — no changes
