# Research: Multifrontal Numeric Factorization

**Date**: 2026-02-15
**Feature**: 015-multifrontal-factorization

## Decision 1: Frontal Matrix Factorization Strategy

**Decision**: Pass the entire frontal matrix (F11 + F21 + F22) to `aptp_factor_in_place()` as a single dense matrix, rather than factoring F11 separately and then computing L21/Schur complement via BLAS-3 (TRSM + GEMM).

**Rationale**: The Phase 5 kernel's `update_schur_1x1` and `update_schur_2x2` functions use `m = a.nrows()` (the full matrix dimension), so Schur complement updates from eliminated columns propagate to ALL trailing rows, including those beyond `num_fully_summed`. After `aptp_factor_in_place()` returns:

- **L11** is in the lower triangle of columns `0..num_eliminated`
- **L21** is already computed at rows `[k..m, 0..num_eliminated]` (via column-by-column rank-1/rank-2 updates)
- **F22** (contribution block) is already updated at `[k..m, k..m]`

This eliminates the need for separate TRSM (L21 computation) and GEMM (Schur complement update) steps in Phase 6.

**Alternatives considered**:
- **Separate BLAS-3 (TRSM + GEMM)**: Factor only F11, then compute `L21 = F21 * L11^{-T} * D11^{-1}` via `solve_lower_triangular_in_place` and `F22 -= L21 * D11 * L21^T` via `matmul`. This gives better cache utilization for large fronts (BLAS-3 vs BLAS-2) but adds implementation complexity. Deferred to Phase 8.1 as a performance optimization.
- **Hybrid**: Use BLAS-2 (pass entire front) for small fronts, BLAS-3 for large fronts. Adds complexity without benefit until parallelism is added.

**Impact on spec**: FR-005 and FR-006 (L21 computation and Schur complement) are satisfied implicitly by the Phase 5 kernel when the entire frontal matrix is passed. They describe _what_ must be true of the output, not _how_ to compute it.

**Caveat**: The `swap_symmetric` function within APTP permutes F21 columns when swapping columns within the fully-summed region. The local permutation from `AptpFactorResult.perm` must be used to map L21 columns back to their original global indices.

## Decision 2: Unified Supernode Abstraction

**Decision**: Create a `SupernodeInfo` struct that provides a uniform interface over both supernodal and simplicial `AptpSymbolic` decompositions. For simplicial, each column is treated as a trivial 1-column supernode.

**Rationale**: AptpSymbolic's supernodal accessors (`supernode_begin()`, `supernode_pattern()`, etc.) return `Option<_>`, returning `None` for simplicial decompositions. Rather than branching on `is_supernodal()` throughout the factorization code, a precomputed `Vec<SupernodeInfo>` provides:
- Column range `[begin..end)` per supernode
- Off-diagonal row pattern
- Parent supernode index (or None for root)
- Children list

For simplicial: each column j produces a SupernodeInfo with `begin=j, end=j+1`, pattern derived from the simplicial CSC structure (`col_ptr`/`row_idx`), and parent derived from `etree()`.

**Alternatives considered**:
- **Direct dispatch on `is_supernodal()`**: Branches at every access point. Error-prone and verbose.
- **Always force supernodal in faer**: Not supported by faer's API — the choice is internal.

## Decision 3: Contribution Block Representation

**Decision**: The contribution block from a child supernode is the trailing `(m - ne) x (m - ne)` submatrix of the factored frontal matrix, where `ne = num_eliminated`. It includes both delayed columns (positions `ne..k`, which become additional fully-summed columns at the parent) and non-fully-summed rows/columns (positions `k..m`).

**Rationale**: After `aptp_factor_in_place()`:
- Positions `0..ne`: factored (stored as FrontFactors)
- Positions `ne..k`: delayed columns with partially-updated values
- Positions `k..m`: non-fully-summed with Schur complement updates applied

The entire `[ne..m, ne..m]` submatrix constitutes what gets extend-added into the parent. The first `(k - ne)` rows/columns are the delayed columns that become fully-summed at the parent; the remaining `(m - k)` are non-fully-summed (partially-summed at parent).

Row/column indices for the contribution block:
- Delayed positions: global indices from `AptpFactorResult.delayed_cols`
- Non-fully-summed positions: original `frontal_row_indices[k..m]` (unchanged by swaps)

**Alternatives considered**:
- **Copy contribution block into separate allocation**: Extra allocation and copy. Instead, we extract the contribution submatrix from the frontal matrix before deallocating it.
- **Store contribution as sparse**: Contribution blocks are dense (standard in multifrontal methods).

## Decision 4: Assembly Tree Traversal

**Decision**: Iterate supernodes in index order `0..n_supernodes` for bottom-up (postorder) traversal. faer guarantees that parent indices are always greater than child indices.

**Rationale**: From the analysis of faer's `SymbolicSupernodalCholesky` and `AptpSymbolic::supernode_parent()`: the parent supernode index is always strictly greater than the child. This means iterating `s = 0, 1, ..., n_supernodes - 1` naturally processes all children before their parents. No explicit topological sort is needed.

For simplicial supernodes derived from `etree()`, the same property holds (elimination tree parent indices are always greater than child indices by construction).

**Precomputation**: Build a `children: Vec<Vec<usize>>` map at the start of factorization by iterating all supernodes and recording each supernode's children. This enables efficient iteration over children during extend-add.

## Decision 5: Dense BLAS Operations via faer

**Decision**: Use faer's existing dense APIs. No new external dependencies needed.

**Key APIs** (verified in faer 0.22 source):
- **Triangular solve**: `faer::linalg::triangular_solve::solve_lower_triangular_in_place(L, rhs, Par::Seq)` — used if we ever need explicit TRSM (Phase 9 BLAS-3 optimization)
- **Matrix multiply**: `faer::linalg::matmul::matmul(dst, Accum::Add, lhs, rhs, alpha, Par::Seq)` — used for reconstruction validation
- **Submatrix views**: `MatMut::submatrix_mut(row_start, col_start, nrows, ncols)` — zero-cost views into frontal matrix partitions

**Already used in codebase**: `matmul` is used in `src/validate.rs` for reconstruction tests. Same patterns apply to Phase 6.

## Decision 6: Reconstruction Test for Distributed Factors

**Decision**: Implement a test utility that reassembles the global L and D matrices from per-supernode `FrontFactors`, then uses the existing `reconstruction_error()` function for validation.

**Rationale**: The existing `reconstruction_error()` in `src/validate.rs` operates on dense L and D matrices. Phase 6 produces per-supernode factors distributed across `Vec<FrontFactors>`. A utility function `reassemble_global_factors(numeric: &AptpNumeric, symbolic: &AptpSymbolic) -> (Mat<f64>, MixedDiagonal)` will:

1. Allocate dense L (n x n) initialized to identity
2. Allocate global MixedDiagonal (n entries)
3. For each supernode s with FrontFactors:
   - Place L11 entries at global positions using `col_indices`
   - Place L21 entries at global positions using `row_indices` and `col_indices`
   - Copy D11 entries into global MixedDiagonal at corresponding positions
4. Apply permutations (fill-reducing + local APTP permutations)

This utility is for testing only (O(n^2) memory), not for the solve path. Phase 7 uses the distributed factors directly.

## Decision 7: Handling Delayed Columns at Parent

**Decision**: When a child supernode has delayed columns, they are embedded in the contribution block and become additional fully-summed columns at the parent. The parent's frontal matrix is dynamically sized to accommodate the symbolic structure plus any delayed columns from all children.

**Rationale**: The contribution block naturally includes delayed columns (at positions `ne..k` in the factored frontal matrix). During extend-add at the parent:
1. Map contribution block's row indices to parent's local indices
2. Add contribution block entries into the parent's frontal matrix
3. The delayed columns' entries end up in the parent's fully-summed region (since they map to columns that belong to this or an ancestor supernode)

The parent's frontal matrix size = (k_parent + d_total + r_parent) where:
- k_parent = symbolic supernode size
- d_total = total delayed columns from all children
- r_parent = off-diagonal pattern size

This requires dynamic sizing at each supernode, but the symbolic analysis provides the base size and the `pivot_buffer_estimates()` provide a hint for expected delays.
