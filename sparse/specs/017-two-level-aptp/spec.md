# Feature Specification: Two-Level APTP Factorization

**Feature Branch**: `017-two-level-aptp`
**Created**: 2026-02-16
**Status**: Draft
**Input**: User description: "Implement Phase 8.1 — two-level APTP factorization with BLAS-3 blocking for large frontal matrices"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Large Indefinite System Factorization (Priority: P1)

A solver user provides a large sparse symmetric indefinite matrix (e.g., from a finite-element optimization problem) where some frontal matrices exceed 256 rows during multifrontal factorization. The solver automatically applies two-level blocked APTP within those large fronts, using BLAS-3 operations (TRSM/GEMM) to factor the frontal matrix efficiently while maintaining the same numerical stability guarantees as the existing single-level kernel.

**Why this priority**: This is the core deliverable of Phase 8.1 — without two-level blocking, large fronts use O(n^2) rank-1/rank-2 updates per column instead of O(1) BLAS-3 calls per block, dominating total factorization time (e.g., ncvxqp3 with max_front 2447 takes 53s single-level). 91% of supernodal SuiteSparse matrices have max_front > 256.

**Independent Test**: Can be fully tested by factoring a dense matrix of size 512+ using the two-level kernel in isolation and verifying reconstruction error < 1e-12, then verifying backward error < 5e-11 on SuiteSparse CI matrices through the full SparseLDLT pipeline.

**Acceptance Scenarios**:

1. **Given** a frontal matrix with dimension > 256, **When** the APTP kernel is invoked, **Then** the two-level algorithm is used (outer blocks of nb columns, inner sub-blocks of ib columns) and the resulting factorization has reconstruction error < 1e-12.
2. **Given** the 9 CI SuiteSparse matrices with MatchOrderMetis ordering, **When** the full SparseLDLT pipeline is run with two-level APTP, **Then** backward error is no worse than the Phase 7 single-level baseline on any matrix.
3. **Given** a frontal matrix with dimension <= 256, **When** the APTP kernel is invoked, **Then** the kernel processes a single outer block (functionally equivalent to single-level) with no performance overhead from the two-level dispatch logic.

---

### User Story 2 - Factor/Apply/Update Decomposition for BLAS-3 (Priority: P1)

The dense APTP kernel decomposes each outer block into three phases: Factor (inner APTP on the diagonal sub-block), Apply (update panel rows below the diagonal using TRSM), and Update (rank-nb Schur complement update on the trailing submatrix using GEMM). This decomposition replaces the current per-column rank-1/rank-2 Schur updates with blocked BLAS-3 operations.

**Why this priority**: The Factor/Apply/Update decomposition is inseparable from the two-level algorithm — it IS the two-level algorithm. Without it, there is no performance benefit to blocking.

**Independent Test**: Can be tested by comparing the numerical output (L, D, permutation, delays) of the two-level kernel against the single-level kernel on identical input matrices of various sizes.

**Acceptance Scenarios**:

1. **Given** a frontal matrix of size 512 with 512 fully-summed columns, **When** factored with outer block size 256 and inner block size 32, **Then** the Factor phase produces L11/D11 satisfying reconstruction error < 1e-12 on the diagonal block, the Apply phase produces L21 consistent with the factored L11/D11, and the Update phase produces a Schur complement such that the overall factorization has reconstruction error < 1e-12.
2. **Given** the Factor/Apply/Update decomposition, **When** some columns in an outer block fail the stability check, **Then** the failed columns are correctly tracked, backup is restored for the affected regions, and the effective elimination count (nelim) is reduced to exclude failed columns before Apply and Update proceed.

---

### User Story 3 - Per-Block Backup and Restore (Priority: P2)

The two-level kernel uses a per-block backup strategy: before factoring each outer block, the relevant matrix entries are backed up. If the a posteriori stability check reveals that some columns in the block produced L entries exceeding the threshold bound, the backup is restored for the failed portion and the failed columns are handled (delayed or re-attempted). This replaces the current per-column backup used in the single-level kernel.

**Why this priority**: Per-block backup is necessary for the two-level algorithm to work correctly but is an implementation concern subordinate to the algorithmic decomposition itself.

**Independent Test**: Can be tested by constructing matrices that deliberately trigger pivot failures at various positions within an outer block and verifying that the correct columns are restored and delayed.

**Acceptance Scenarios**:

1. **Given** an outer block where columns 5-7 (out of 32) fail the stability check after Apply, **When** the Adjust phase runs, **Then** only the first 5 columns' factorization is retained, columns 5-7 are restored from backup, and the contribution of columns 0-4 is correctly applied to the trailing submatrix.
2. **Given** a matrix where no columns fail within any block, **When** factored with two-level APTP, **Then** backup memory is allocated but never restored, and results satisfy reconstruction error < 1e-12 (note: results will not be bitwise identical to single-level due to different BLAS-3 operation ordering).

---

### Edge Cases

- What happens when the frontal matrix dimension is exactly equal to the outer block size (e.g., 256)? The kernel should handle a single outer block correctly (degenerate case of the two-level loop — equivalent to single-level).
- What happens when the frontal matrix dimension is slightly larger than a multiple of the outer block size (e.g., 257)? The final outer block is smaller than nb and must be handled correctly.
- What happens when all columns in an outer block fail the stability check? The entire block's backup is restored and all columns are delayed; the trailing matrix update for that block is skipped.
- What happens when 2x2 pivots straddle an outer block boundary? This cannot happen — 2x2 pivots are chosen within the inner Factor phase on the diagonal sub-block, which operates entirely within the current outer block. The inner kernel's swap-delayed-to-end architecture ensures 2x2 partners are adjacent.
- What happens when the number of delayed columns from previous outer blocks accumulates? Delayed columns from earlier outer blocks become "failed entries" in the terminology of Algorithm 3.1 (Duff, Hogg & Lopez 2020) and receive UpdateNT/UpdateTN updates from subsequent blocks.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The APTP kernel MUST support a two-level blocking strategy with configurable outer block size (default: 256) and inner block size (default: 32).

- **FR-002**: For frontal matrices with dimension > outer block size, the kernel MUST use the two-level algorithm with multiple outer blocks. For dimensions <= outer block size, the kernel MUST process a single outer block (functionally equivalent to single-level behavior). No separate single-level code path is maintained.

- **FR-003**: Each outer block MUST be processed in three phases: Factor (inner APTP on diagonal sub-block), Apply (TRSM-based panel solve for rows below the diagonal), and Update (GEMM-based rank-nb Schur complement update on the trailing submatrix).

- **FR-004**: The Apply phase MUST perform an a posteriori threshold check on all L entries produced for the current outer block and reduce the effective elimination count (nelim) to the first column where any entry exceeds the stability bound 1/u.

- **FR-005**: The kernel MUST back up matrix entries per outer block before the Factor phase and restore failed columns from backup when the a posteriori check reduces nelim.

- **FR-006**: After all outer blocks are processed, failed/delayed columns MUST be permuted to the end of the fully-summed range, consistent with the existing single-level kernel's contract (swap-delayed-to-end).

- **FR-007**: The two-level kernel MUST produce the same output types as the single-level kernel (`AptpFactorResult` with `MixedDiagonal`, column permutation, delayed columns, and statistics) so that the multifrontal factorization (`AptpNumeric::factor`) requires no changes to its kernel call interface.

- **FR-008**: The two-level kernel MUST handle previously-delayed columns from earlier outer blocks by applying updates from each newly-factored block to the delayed portion (the UpdateNT/UpdateTN operations from Algorithm 3.1).

- **FR-009**: The kernel MUST handle 2x2 pivots chosen by the inner Factor phase correctly, including proper D storage in `MixedDiagonal` and correct application during the Apply phase (the Apply must account for both 1x1 and 2x2 pivot structures in the D solve).

- **FR-010**: The kernel MUST fall back gracefully when the remaining columns in an outer block are fewer than the configured outer block size (last-block handling).

- **FR-011**: The innermost ib x ib diagonal blocks within the inner Factor phase MUST use complete pivoting (Algorithm 4.1 from Duff, Hogg & Lopez 2020), which searches the entire remaining submatrix for the largest entry and uses 1x1 or 2x2 pivots accordingly. This replaces the current 1x1 threshold + 2x2 Bunch-Kaufman fallback at the leaf level of the recursion.

### Key Entities

- **Outer Block**: A contiguous group of nb columns processed together. The two-level loop iterates over outer blocks left-to-right.
- **Inner Sub-Block**: A contiguous group of ib columns within an outer block's diagonal. The inner Factor phase processes sub-blocks using single-level APTP, with the innermost ib x ib diagonal blocks factored using complete pivoting (Algorithm 4.1 from Duff, Hogg & Lopez 2020).
- **Factor Phase**: Inner APTP factorization of the nb x nb diagonal sub-block, producing L11, D11, and a local permutation.
- **Apply Phase**: Panel solve L21 = A21 * (L11 * D11)^{-T}, with a posteriori stability check. Also ApplyT for columns to the left of the current block (failed entries from prior blocks).
- **Update Phase**: Rank-nb Schur complement update A22 -= L21 * D11 * L21^T using GEMM. Also UpdateNT/UpdateTN for previously-delayed regions.
- **Backup**: Saved copy of matrix entries for a block, enabling restore on pivot failure.
- **nelim**: The effective number of successfully eliminated columns in an outer block, determined by the minimum across all Apply checks.

## Algorithm References

The following academic papers and reference implementations inform this feature:

- **Duff, Hogg & Lopez (2020)**, "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting", SIAM J. Sci. Comput. 42(4).
  - **Section 3**: Algorithm 3.1 — Full APTP with 2D block partitioning (Factor/ApplyN/ApplyT/Adjust/UpdateNN/UpdateNT/UpdateTN)
  - **Section 3, Figure 2**: Two-level blocking strategy (outer nb, inner ib)
  - **Section 4**: Complete pivoting within inner blocks, stability analysis (threshold u=0.25 equivalence)
  - **Section 5**: Multifrontal integration
  - **Markdown**: `/workspace/rivrs-linalg/references/ssids/duff2020.md`

- **Hogg, Scott & Sherlock (2016)**, implementation notes for SSIDS supernodal factorization
  - **Markdown**: `/workspace/rivrs-linalg/references/ssids/hogg2016.md`

- **Liu (1992)**, "The Multifrontal Method for Sparse Matrix Solution: Theory and Practice"
  - Elimination tree, level-set scheduling
  - **Markdown**: `/workspace/rivrs-linalg/references/ssids/liu1992.md`

- **SPRAL source** (BSD-3): `spral/src/ssids/cpu/kernels/ldlt_app.cxx`
  - CopyBackup/PoolBackup classes for backup strategy
  - Factor/ApplyN/ApplyT/UpdateNN decomposition
  - INNER_BLOCK_SIZE = 32, outer BLOCK_SIZE as template parameter
  - Located at `/opt/references/spral/src/ssids/cpu/kernels/ldlt_app.cxx`

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Two-level APTP produces reconstruction error (||P^T A P - L D L^T|| / ||A||) within 1e-12 on all hand-constructed test matrices and random matrices up to size 2048.

- **SC-002**: All 9 CI SuiteSparse matrices achieve backward error no worse than the Phase 7 single-level baseline when using MatchOrderMetis ordering. Target: improvement on currently-failing matrices (bratu3d, cvxqp3, stokes128, sparsine).

- **SC-003**: Factorization of dense frontal matrices with dimension 512+ completes faster with two-level APTP than with single-level, as measured by criterion benchmarks on isolated kernel calls.

- **SC-004**: Full SuiteSparse solver pipeline (67 matrices, MatchOrderMetis) shows no correctness regressions — every matrix that passed Phase 7's backward error threshold continues to pass.

- **SC-005**: The crossover point where two-level blocking outperforms single-level is identified through benchmarking and documented. The dispatch threshold is set at or near this crossover.

- **SC-006**: The two-level kernel uses BLAS-3 operations (faer's matmul/triangular_solve) for the Apply and Update phases, not per-column rank-1/rank-2 updates.

## Clarifications

### Session 2026-02-16

- Q: Should the innermost ib x ib diagonal blocks use complete pivoting (Algorithm 4.1) or the current 1x1 threshold + 2x2 BK fallback? → A: Complete pivoting (Algorithm 4.1) for innermost ib x ib diagonal blocks, matching the paper and SPRAL's `block_ldlt()`.

## Assumptions

- The outer block size default of 256 and inner block size default of 32 from SPRAL and the Duff, Hogg & Lopez (2020) paper are reasonable starting points. These may be tuned based on benchmarking but the algorithm should work correctly for any valid block sizes.
- The multifrontal factorization (`AptpNumeric::factor` in numeric.rs) will continue to pass the entire frontal matrix to the kernel. The two-level decomposition is internal to the kernel — the caller does not need to know whether single-level or two-level was used.
- faer's `faer::linalg::matmul::matmul` and `faer::linalg::triangular_solve` provide sufficient BLAS-3 performance for the Apply (TRSM) and Update (GEMM) phases. No external BLAS library is needed.
- Per-block backup (backing up the relevant portions of the matrix before each outer block's Factor phase) is sufficient. The more aggressive "backup entire matrix" strategy (SPRAL's APP_AGGRESSIVE) is not needed for correctness and is deferred unless benchmarking reveals excessive pivot failures requiring global rollback.
- The single-level kernel code will be replaced by the two-level implementation (not maintained as a separate code path) once two-level is validated. The single-level behavior is recovered when the frontal matrix is smaller than one outer block.
- Parallelism within the two-level decomposition (e.g., parallelizing Apply across block rows, or Update across block pairs) is deferred to Phase 8.2. Phase 8.1 implements the sequential two-level algorithm.
