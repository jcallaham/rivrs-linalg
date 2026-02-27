# Research: Triangular Solve & Solver API

**Feature**: 016-triangular-solve-api
**Date**: 2026-02-16

## R1: faer Dense Triangular Solve API

**Decision**: Use `solve_unit_lower_triangular_in_place_with_conj` for forward L11 solve and `solve_unit_upper_triangular_in_place_with_conj` with `l11.transpose()` for backward L11^T solve.

**Rationale**: faer does not have a direct "solve lower triangular transpose" function. The standard pattern (used in faer's own supernodal LBLT solver) is to transpose L11 into an upper triangular matrix and call the upper-triangular solver. For real f64, conjugation is a no-op so `Conj::No` is used throughout.

**Alternatives considered**:
- Manual triangular solve loop: rejected — faer's implementation is BLAS-optimized and handles blocking/SIMD.
- `solve_lower_triangular_in_place` (non-unit): rejected — our L11 has unit diagonal by construction.

**Key signatures**:
```rust
use faer::linalg::triangular_solve::{
    solve_unit_lower_triangular_in_place_with_conj,
    solve_unit_upper_triangular_in_place_with_conj,
};
// Forward: solve_unit_lower_triangular_in_place_with_conj(l11.as_ref(), Conj::No, local.as_mat_mut(), Par::Seq)
// Backward: solve_unit_upper_triangular_in_place_with_conj(l11.as_ref().transpose(), Conj::No, local.as_mat_mut(), Par::Seq)
```

**Gotcha**: These functions take `MatRef`/`MatMut` (2D), not column vectors. A `&mut [f64]` local buffer must be wrapped as a column matrix view via `faer::col::from_slice_mut` or `Mat` construction.

## R2: faer MemStack / StackReq Pattern

**Decision**: Use `StackReq::new::<f64>(max_front_size)` for workspace sizing. Allocate `MemBuffer::new(req)` once, create `MemStack::new(&mut buf)`, and use `temp_mat_zeroed` to allocate per-supernode temporary buffers from the stack.

**Rationale**: This follows faer's convention exactly. The workspace is reused across supernodes (each supernode's allocation is released when moving to the next), so the total allocation is just the maximum per-supernode need.

**Key pattern**:
```rust
use faer::dyn_stack::{MemStack, StackReq, MemBuffer};
use faer::linalg::{temp_mat_zeroed, temp_mat_scratch};

// Compute requirement (from AptpNumeric stats)
let max_size = numeric.stats().max_front_size;
let req = temp_mat_scratch::<f64>(max_size, 1);

// User allocates once:
let mut mem = MemBuffer::new(req);
let stack = MemStack::new(&mut mem);

// Inside solve, per-supernode:
let (mut tmp, _remaining) = temp_mat_zeroed::<f64, _, _>(r, 1, stack);
```

**Alternatives considered**:
- Heap-allocate Vec per supernode: rejected — violates "no heap in solve hot path" requirement.
- Pre-allocate max-sized Vec outside stack: simpler but doesn't compose with faer's stack convention.

**Gotcha**: `temp_mat_uninit` is unsafe (uninitialized memory). Use `temp_mat_zeroed` for safety unless profiling shows it's a bottleneck. faer's internal solvers use `temp_mat_uninit` for performance.

**Gotcha**: `MemBuffer` may need direct import from `dyn_stack` crate if not re-exported by faer. Verify at compile time.

## R3: Dense Matrix-Vector Multiply for Scatter/Gather

**Decision**: Use `matmul_with_conj` from `faer::linalg::matmul` for both the forward scatter (`tmp = L21 * local`) and backward gather-update (`local -= L21^T * tmp`).

**Rationale**: faer's matmul dispatches to optimized BLAS-2 kernels when one dimension is 1 (GEMV case). Using the same API for both directions keeps the code consistent with faer's supernodal solve implementation.

**Key patterns**:
```rust
use faer::linalg::matmul::matmul_with_conj;
use faer::{Accum, Conj, Par};

// Forward scatter: tmp = L21 * local
matmul_with_conj(tmp.rb_mut(), Accum::Replace, l21.as_ref(), Conj::No, local.as_ref(), Conj::No, 1.0, Par::Seq);

// Backward update: local -= L21^T * tmp
matmul_with_conj(local.rb_mut(), Accum::Add, l21.as_ref().transpose(), Conj::No, tmp.as_ref(), Conj::No, -1.0, Par::Seq);
```

**Alternatives considered**:
- Manual dot-product loops: rejected — less readable and doesn't benefit from faer's SIMD.
- `Col`-based operations: rejected — matmul works on `MatRef`/`MatMut`; column views convert via `.as_mat()`.

## R4: Sparse Backward Error Computation

**Decision**: Use `sparse_dense_matmul` from `faer::sparse::linalg::matmul` to compute `A * x` without dense conversion, then compute `||A*x - b|| / (||A||_F * ||x|| + ||b||)`.

**Rationale**: The existing `validate::backward_error()` calls `a.to_dense()` which is O(n^2) memory — impractical for SuiteSparse matrices with n > 5000. `sparse_dense_matmul` operates directly on `SparseColMatRef` with O(nnz) work.

**Key pattern**:
```rust
use faer::sparse::linalg::matmul::sparse_dense_matmul;

// For symmetric A stored as lower triangle CSC, sparse_dense_matmul
// computes only A_lower * x, NOT the full symmetric product. To get
// the correct Ax for backward error, use two passes:
//   ax  = A_lower * x           (standard sparse matmul)
//   ax += A_lower^T * x         (transpose matmul)
//   ax -= diag(A) * x           (subtract diagonal counted twice)
// This is O(nnz) work and O(n) extra memory (for diagonal extraction).
let mut ax = Mat::<f64>::zeros(n, 1);
sparse_dense_matmul(ax.as_mut(), Accum::Replace, a.as_ref(), x.as_mat(), 1.0, Par::Seq);
sparse_dense_matmul(ax.as_mut(), Accum::Add, a.as_ref().transpose(), x.as_mat(), 1.0, Par::Seq);
// Subtract diagonal contribution (counted twice above):
for j in 0..n {
    let col_start = a.col_ptr()[j];
    let col_end = a.col_ptr()[j + 1];
    for idx in col_start..col_end {
        if a.row_idx()[idx] == j {
            ax[(j, 0)] -= a.val()[idx] * x[j];
        }
    }
}
// residual = ax - b, then compute norms
```

**Note**: For `||A||_F` computation on sparse matrices, we can iterate over `a.val()` directly: `sqrt(sum of val[i]^2 * multiplier)` where the multiplier is 2 for off-diagonal entries and 1 for diagonal (since we store lower triangle of symmetric A).

**Gotcha**: `sparse_dense_matmul` currently ignores the `par` parameter for CSC matrices (sequential internally). Not a concern for Phase 7.

**Gotcha**: `sparse_dense_matmul` on lower-triangle CSC computes `A_lower * x`, NOT the full symmetric `A * x`. For symmetric backward error, the two-pass approach above (A_lower + A_lower^T - diag) is required. This is a common pitfall when testing with matrices stored as lower triangle only.

## R5: Scaling Coordinate System

**Decision**: Store MC64 scaling factors in elimination order on `SparseLDLT`. Permute from original order during `analyze()`.

**Rationale**: SPRAL stores `fkeep%scaling` in elimination order: `fkeep%scaling(i) = scaling(akeep%invp(i))`. Scaling is applied after permutation in both factor (scatter step) and solve (pre/post scaling). Storing in elimination order means scaling application is a simple element-wise multiply without index remapping at solve time.

**Transform**: Given `orig_scaling` from `match_order_metis()` and `perm_fwd` from symbolic analysis:
```
elim_scaling[i] = orig_scaling[perm_fwd[i]]  for i = 0..n
```

**Application in factor** (scatter step): `scaled_val = scaling[perm_inv[orig_row]] * val * scaling[perm_inv[orig_col]]`. Since `perm_inv[orig_row]` maps original row → elimination index, and scaling is in elimination order, this is just `scaling[elim_row] * val * scaling[elim_col]`.

**Application in solve**: Before forward solve: `rhs_perm[i] *= scaling[i]`. After backward solve: `rhs_perm[i] *= scaling[i]`. Symmetric scaling (same S on both sides).

## R6: col_indices Already Encodes local_perm

**Decision**: No runtime `local_perm` application needed during solve. The gather/scatter uses `col_indices` directly.

**Rationale**: `extract_front_factors` in `numeric.rs` constructs `col_indices` as:
```rust
let col_indices: Vec<usize> = local_perm[..ne]
    .iter()
    .map(|&lp| frontal.row_indices[lp])
    .collect();
```
This means `col_indices[i]` = global permuted index of the column at **factored position** `i`. Since L11 and D11 are stored in factored order, gathering `rhs[col_indices[i]]` produces values already aligned with L11's column ordering.

The plan document (written before implementation) describes "apply local_perm to reorder" as a separate step, but this was absorbed into `col_indices` construction during Phase 6 implementation. No plan update needed — the spec is satisfied without explicit runtime permutation.

**Verification**: Unit tests will confirm gather produces correct values by comparing against analytically computed results on hand-constructed matrices.

## R7: Zero Pivot Handling in MixedDiagonal

**Decision**: Modify `MixedDiagonal::solve_in_place()` to handle zero pivots gracefully instead of panicking.

**Rationale**: SPRAL's `action=true` (default) behavior: zero pivots are stored as `d[2*i] = 0.0` and the D-solve naturally computes `x[i] *= 0.0`, zeroing the corresponding solution component. Our `MixedDiagonal` currently has `debug_assert!(d != 0.0)` which panics in debug builds.

**Changes**:
- For 1x1 pivot with `d == 0.0`: set `x[col] = 0.0` (instead of dividing by zero)
- For 2x2 block with `det == 0.0`: set `x[col] = 0.0` and `x[partner] = 0.0`
- Remove `debug_assert!(d != 0.0)` and `debug_assert!(det != 0.0)`
- The `debug_assert!(self.num_delayed() == 0)` remains — delayed columns in D is a programming error

**Note**: Phase 6 already records `zero_pivots` in `FactorizationStats`. The solve just needs to handle them without panicking.

## R8: faer's LBLT Solve as Reference Implementation

**Decision**: Model our solve implementation on faer's `SupernodalIntranodeLbltRef::solve_in_place_no_numeric_permute_with_conj` (cholesky.rs lines 2261-2373), adapted for our `FrontFactors` structure.

**Rationale**: faer's LBLT solve handles 1x1+2x2 pivots with per-supernode gather/scatter — the exact same pattern we need. Key differences from our implementation:
- faer's supernodes are contiguous (can slice global vector directly); ours may not be (must gather via `col_indices`)
- faer uses a `subdiag` array to detect 2x2 blocks; we use `MixedDiagonal` with `PivotType` enum
- faer's pivot permutation is per-global-index; ours is absorbed into `col_indices`

The overall structure (forward/D/backward with gather/matmul/scatter) is identical.
