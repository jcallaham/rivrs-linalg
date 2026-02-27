# Data Model: Contribution Workspace Reuse

**Feature**: 023-contrib-workspace-reuse
**Date**: 2026-02-23

## Modified Entities

### FactorizationWorkspace (extended)

Existing pre-allocated reusable buffer for the per-supernode factorization loop.

**Existing fields** (unchanged):
- `frontal_data: Mat<f64>` — dense frontal matrix buffer (max_front × max_front)
- `frontal_row_indices: Vec<usize>` — row indices for current front
- `delayed_cols_buf: Vec<usize>` — delayed columns from children
- `global_to_local: Vec<usize>` — global→local index mapping (length n)

**New fields**:
- `frontal_data_alt: Mat<f64>` — second frontal buffer for dual-buffer ping-pong (max_front × max_front). Used by parent while child occupies `frontal_data`.
- `frontal_row_indices_alt: Vec<usize>` — row indices for the alternate buffer.
- `contribution_pool: Vec<Mat<f64>>` — free-list of reusable dense matrices for contribution block extraction. Buffers are taken before extraction and returned after extend-add consumption. Sized lazily (first use allocates, subsequent uses recycle).

**Invariants**:
- At most one of `frontal_data` / `frontal_data_alt` is actively being factored at any time
- `contribution_pool` buffers may be any size ≥ the contribution they're used for
- In parallel mode, each thread owns its own workspace instance (no sharing)

**Lifecycle**:
- Created once at factorization start (`FactorizationWorkspace::new`)
- Reused across all supernodes in the factorization
- Frontal buffers swap roles at each parent/child transition
- Pool buffers cycle: taken → written → stored in ContributionBlock → consumed → returned

### ContributionBlock (modified behavior, same structure)

Schur complement and delayed columns from a factored supernode.

**Fields** (unchanged):
- `data: Mat<f64>` — dense trailing submatrix
- `row_indices: Vec<usize>` — global permuted indices
- `num_delayed: usize` — delayed column count

**Changed behavior**:
- `data` may now originate from the contribution pool (a recycled `Mat<f64>`) rather than a fresh `Mat::zeros`
- When consumed by extend-add, the `Mat<f64>` should be returned to the pool rather than dropped
- In the direct extend-add case (last child, sequential path), no `ContributionBlock` is created at all — the contribution is read directly from the frontal workspace

### FrontalMatrix (unchanged)

Borrowed view into the frontal workspace. No structural changes.

**New usage**: In the direct extend-add path, the trailing submatrix of a factored `FrontalMatrix` is read directly for extend-add, using the contribution's row indices computed from the APTP result permutation and the frontal_row_indices.

### AssemblyMaps (potentially extended)

Precomputed scatter index structure.

**Potential new field**:
- Direct extend-add row mappings: for the zero-delay case, mapping from frontal workspace positions `[ne..m]` to parent frontal positions, allowing direct extend-add without `global_to_local` lookups.

**Decision**: Whether to precompute these mappings or compute them on-the-fly depends on profiling. The existing `ea_map` already provides child contribution → parent row mappings; direct extend-add may be able to reuse these with an offset adjustment for `ne`.

## New Internal Concepts

### Buffer Role Assignment

In the dual-buffer sequential path, each frontal buffer has a role:
- **Active**: Currently being factored (holds this supernode's frontal matrix)
- **Accumulating**: Receiving extend-add contributions for a parent supernode

Roles swap when transitioning from child processing to parent processing.

### Contribution Pool Protocol

The contribution pool follows a simple take/return discipline:
1. **Take**: Before `extract_contribution`, take a buffer from the pool. If pool is empty, allocate a new `Mat<f64>`.
2. **Size check**: If the taken buffer is too small for the contribution, replace it with a new allocation (return the small one to the pool or discard).
3. **Zero**: Zero the triangle of the taken buffer that will be written.
4. **Write**: `extract_contribution` writes into the taken buffer.
5. **Store**: The buffer becomes `ContributionBlock.data`.
6. **Consume**: After extend-add, the `Mat<f64>` is detached from the consumed `ContributionBlock`.
7. **Return**: The detached `Mat<f64>` is returned to the pool for reuse.

## State Transitions

### Supernode Processing (Sequential DFS Path)

```
IDLE → CHILD_FACTORING → CHILD_EXTEND_ADD → [repeat for more children] → PARENT_SCATTER → PARENT_FACTORING → PARENT_EXTRACT → IDLE
```

- `CHILD_FACTORING`: Child's frontal matrix assembled and factored in active buffer
- `CHILD_EXTEND_ADD`: Child's contribution extended-add'd into accumulating buffer (either direct from active buffer for last child, or from pooled ContributionBlock for earlier children)
- `PARENT_SCATTER`: Original CSC entries scattered into accumulating buffer
- `PARENT_FACTORING`: Accumulating buffer becomes active buffer, parent factored
- `PARENT_EXTRACT`: Parent's contribution extracted (pool or direct for parent's parent)
