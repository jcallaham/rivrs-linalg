# Feature Specification: Triangular Solve & Solver API

**Feature Branch**: `016-triangular-solve-api`
**Created**: 2026-02-15
**Status**: Draft
**Input**: User description: "Implement Phase 7: Triangular Solve & Solver API — per-supernode forward/backward solve, diagonal solve, SparseLDLT user-facing struct with analyze/factor/solve pipeline, scaling integration, backward error validation"

## User Scenarios & Testing *(mandatory)*

### User Story 1 — End-to-End Sparse Solve (Priority: P1)

A researcher or engineer has a sparse symmetric indefinite matrix and a right-hand side vector. They want to solve `Ax = b` through a single high-level call, obtaining an accurate solution validated by backward error.

**Why this priority**: This is the core deliverable — a working sparse symmetric indefinite solver. Without a correct end-to-end solve pipeline, the factorization from Phases 2-6 has no user-facing value.

**Independent Test**: Construct `b = A * x_exact` for test matrices (hand-constructed + SuiteSparse), solve via the full pipeline, verify `||Ax - b|| / (||A||_F ||x|| + ||b||) < 10^-10`.

**Acceptance Scenarios**:

1. **Given** a sparse symmetric positive definite matrix and a known exact solution, **When** the user calls the one-shot solve method, **Then** the backward error is below `10^-10`.
2. **Given** a sparse symmetric indefinite matrix from the SuiteSparse collection, **When** the user calls analyze → factor → solve, **Then** the backward error is below `10^-10` and inertia matches the reference (when available).
3. **Given** a matrix with 2x2 Bunch-Kaufman pivots and delayed columns (from Phase 5 test suite), **When** the user solves with the full pipeline, **Then** the solution is correct (backward error < `10^-10`), confirming that APTP pivot handling is end-to-end sound.

---

### User Story 2 — Three-Phase API with Reuse (Priority: P1)

A user solving multiple systems with the same sparsity pattern but different numeric values (e.g., interior point iteration) needs to analyze once, then factor and solve repeatedly. The symbolic analysis must be reusable, and refactoring must produce correct results.

**Why this priority**: Reusable symbolic analysis is the primary performance advantage of direct solvers in iterative contexts. This is co-equal with Story 1 because the API design directly affects all downstream usage.

**Independent Test**: Analyze a matrix, factor and solve, then refactor with different numeric values and solve again — both solutions must meet backward error tolerance.

**Acceptance Scenarios**:

1. **Given** a completed symbolic analysis, **When** the user calls `factor()` followed by `solve()`, **Then** the solution has backward error < `10^-10`.
2. **Given** a completed factorization, **When** the user calls `refactor()` with a matrix having the same sparsity pattern but different values, **Then** the new solve produces correct results.
3. **Given** a completed factorization, **When** the user solves with multiple different right-hand side vectors reusing the same factorization, **Then** all solutions are correct.

---

### User Story 3 — Per-Supernode Triangular Solve Correctness (Priority: P1)

The internal per-supernode forward/backward substitution must correctly gather, permute, solve, and scatter entries through the assembly tree. This is the numerically sensitive core that makes Story 1 and 2 possible.

**Why this priority**: Per-supernode index mapping (gather/scatter/local_perm) is explicitly identified in the development plan as "the most error-prone part of the solve." Incorrect index mapping silently produces wrong solutions. This must be unit-tested independently before integration.

**Independent Test**: For small hand-constructed matrices with known supernode structure, verify that each per-supernode forward solve step produces the analytically expected local vector.

**Acceptance Scenarios**:

1. **Given** a factored matrix with known per-supernode factors, **When** forward solve gathers entries from the global RHS using `col_indices`, **Then** the local vector contains exactly the expected entries in the correct order.
2. **Given** a supernode with APTP-reordered columns, **When** entries are gathered from the global RHS via `col_indices`, **Then** the local vector aligns with L11's column ordering (because `col_indices` encodes the APTP reordering from Phase 6's `extract_front_factors`).
3. **Given** a factored supernode with L21, **When** the scatter step updates non-fully-summed rows, **Then** `rhs[row_indices[i]] -= L21[i,:] * y_local` is applied correctly for all rows.
4. **Given** the complete forward → D → backward solve pipeline, **When** run on a 2-supernode hand-constructed matrix, **Then** the intermediate vectors at each stage match analytically computed values.

---

### User Story 4 — Scaling Integration (Priority: P2)

When the user requests MC64 matching/scaling during analysis, the scaling factors must be correctly applied during factorization and solve so that the final solution is in the original (unscaled) coordinate system.

**Why this priority**: MC64 scaling improves pivot quality for hard indefinite problems. However, the solver works without scaling (AMD/METIS ordering alone), making this secondary to the core solve pipeline.

**Independent Test**: For matrices where MC64 scaling is applied, verify that: (a) the factored matrix corresponds to the scaled system, and (b) the solve correctly unscales the solution. Compare backward error with and without scaling.

**Acceptance Scenarios**:

1. **Given** a matrix analyzed with MC64 ordering/scaling, **When** the user factors and solves, **Then** the backward error is computed against the original (unscaled) matrix and is below `10^-10`.
2. **Given** scaling factors `S`, **When** the scaling round-trip `S * (S^{-1} * x)` is applied, **Then** the result equals `x` to machine precision.
3. **Given** a hard indefinite matrix, **When** solved with MC64 scaling vs. without, **Then** the MC64 variant achieves equal or better backward error.

---

### User Story 5 — Workspace-Efficient Solve (Priority: P2)

The solve phase must use stack-based workspace allocation (faer's `MemStack`) to avoid heap allocation in the hot path. Users must be able to query workspace requirements upfront and reuse a single buffer across multiple solves.

**Why this priority**: Stack-based allocation is a faer convention and important for performance in repeated-solve workloads. However, correctness comes first — this is an efficiency concern.

**Independent Test**: Verify that `solve_scratch()` returns a valid `StackReq`, and that solve succeeds with a `MemStack` of exactly that size.

**Acceptance Scenarios**:

1. **Given** a factored matrix, **When** the user queries `solve_scratch(1)`, **Then** a valid `StackReq` is returned reflecting the maximum supernode dimensions.
2. **Given** a pre-allocated `MemBuffer` from the returned `StackReq`, **When** the user passes the corresponding `MemStack` to `solve_in_place()`, **Then** the solve succeeds without additional allocation.

---

### Edge Cases

- **Solve before factor**: Calling `solve()` on a `SparseLDLT` that has been analyzed but not factored must return a clear error (not panic).
- **Zero-dimension matrix**: A 0x0 matrix should be handled gracefully (trivial solve).
- **Rank-deficient matrix**: When `FactorizationStats::zero_pivots > 0`, the D-solve must handle zero pivots gracefully (zeroing the corresponding solution components, following SPRAL's convention) rather than panicking. Users detect rank deficiency via `stats().zero_pivots` and `inertia()`.
- **Single-column supernodes (simplicial fallback)**: The solve must work correctly when faer chooses simplicial decomposition (each column is a 1-column supernode).
- **Empty supernodes (ne = 0)**: If a supernode eliminates zero columns (all delayed), the solve must skip it gracefully.
- **RHS dimension mismatch**: Passing a RHS vector whose length doesn't match the matrix dimension must produce a clear error.
- **Sparse backward error for large matrices**: The current `validate::backward_error()` uses dense conversion (O(n^2) memory). Phase 7 must provide a sparse-aware backward error computation for SuiteSparse validation.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The solver MUST implement a per-supernode forward substitution that traverses supernodes in postorder, gathering entries from the global RHS via `col_indices` (which encodes the APTP column reordering from Phase 6), solving `L11 * y = rhs_local` via dense triangular solve, and scattering updates via `L21` to non-fully-summed rows.
- **FR-002**: The solver MUST implement a per-supernode diagonal solve that applies `d11.solve_in_place()` to each supernode's eliminated entries (gathered via `col_indices`), handling both 1x1 and 2x2 pivot blocks.
- **FR-003**: The solver MUST implement a per-supernode backward substitution that traverses supernodes in reverse postorder, applying `L21^T` updates and solving `L11^T * z = y_local`, with correct inverse-permutation and scatter.
- **FR-004**: The solver MUST apply the fill-reducing permutation (`AptpSymbolic::perm()`) to the RHS before forward solve and the inverse permutation to the solution after backward solve.
- **FR-005**: The solver MUST apply MC64 scaling factors (when present) as `S * rhs` before forward solve and `S * solution` after backward solve.
- **FR-006**: The solver MUST apply MC64 scaling to matrix entries during factorization as `scaled_value = s[i] * a[i][j] * s[j]` when scaling factors are present. Following SPRAL's approach, scaling factors are passed into the factorization and applied at the assembly (scatter) layer; the APTP pivot kernel itself remains scaling-unaware. No pre-scaled copy of the matrix is created.
- **FR-007**: The solver MUST expose a user-facing solver struct with the three-phase API: `analyze()` → `factor()` → `solve()`, where the symbolic analysis is reusable across factorizations with the same sparsity pattern.
- **FR-008**: The solver MUST expose both an allocating solve method (returning a new solution vector) and an in-place solve method (modifying the RHS).
- **FR-009**: The solver MUST expose a workspace-query method returning the stack requirement for solve, and the in-place solve method MUST accept a pre-allocated workspace stack.
- **FR-010**: The solver MUST expose a one-shot convenience method that performs analyze → factor → solve in a single call, allocating workspace internally.
- **FR-011**: The solver MUST expose a refactor method which reuses the symbolic analysis to factor a matrix with the same sparsity pattern but different numeric values.
- **FR-012**: The solver MUST expose inertia (eigenvalue sign counts) and factorization statistics from the most recent factorization.
- **FR-013**: The solver MUST return clear errors for: solve before factor, RHS dimension mismatch, and matrix dimension mismatch with symbolic analysis. For rank-deficient matrices (zero pivots), the solver MUST handle them gracefully by default (zeroing solution components at zero-pivot entries, following SPRAL's `action=true` convention) and expose rank deficiency through factorization statistics (`zero_pivots`) and `inertia()` rather than returning an error.
- **FR-014**: The validation module MUST provide a sparse-aware backward error computation that does not require O(n^2) dense conversion, enabling backward error validation on large SuiteSparse matrices.
- **FR-015**: The core per-supernode solve logic MUST be implemented as a separate internal function that is independent of scaling and permutation, with the user-facing solver wrapping it with the scaling/permutation layer.

### Key Entities

- **Solver struct**: The user-facing solver. Wraps symbolic analysis, optional numeric factors, and optional scaling factors. Provides the analyze/factor/solve API.
- **Core solve function**: Internal solve routine taking symbolic analysis, numeric factors, and a mutable RHS. Performs forward → D → backward solve through the supernode tree. Scaling- and permutation-unaware.
- **Analysis options**: Configuration for the analysis phase — ordering strategy (AMD, METIS, match-order-METIS, user-supplied), whether to request MC64 scaling.
- **Factorization options**: Configuration for factorization — APTP threshold, fallback strategy. Wraps existing factorization configuration.
- **Solver options**: Configuration for the one-shot solve — combines analysis and factorization options.

### Algorithm References

The following academic references and code sources inform the solve algorithm design:

1. **Duff, Hogg & Lopez (2020)** — "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting" (`/workspace/rivrs-linalg/references/ssids/duff2020.md`): Primary reference for the APTP factorization and solve algorithm. Section 3 describes how per-supernode factors integrate with the solve pipeline.

2. **Hogg, Ovtchinnikov & Scott (2016)** — "A sparse symmetric indefinite direct solver for GPU architectures" (`/workspace/rivrs-linalg/references/ssids/hogg2016.md`): Section 3 describes the SSIDS solve phase implementation including partial forward/backward solve, scaling integration, and performance considerations.

3. **Liu (1992)** — "The Multifrontal Method for Sparse Matrix Solution: Theory and Practice" (`/workspace/rivrs-linalg/references/ssids/liu1992.md`): Sections 4-5 on frontal matrix organization and postorder traversal, which determines the supernode visit order for forward/backward solve.

4. **Duff & Pralet (2005)** — "Strategies for Scaling and Pivoting for Symmetric Indefinite Problems" (`/workspace/rivrs-linalg/references/ssids/duff2005.md`): MC64SYM scaling formula and how scaling integrates with the solve pipeline (apply S before factor, apply/unapply S around solve).

5. **Duff & Reid (1984)** — "The Multifrontal Solution of Unsymmetric Sets of Linear Equations" (`/workspace/rivrs-linalg/references/ssids/duff1984.md`): Original multifrontal solve phase description.

6. **Duff & Reid (1999)** — "The Design of MA57" (`/workspace/rivrs-linalg/references/ssids/duff1999.md`): MA57's solve phase is a direct antecedent to SSIDS and provides additional reference for the triangular solve mechanics.

7. **SPRAL source code** (BSD-3): `spral/src/ssids/cpu/kernels/ldlt_app.cxx` (`ldlt_app_solve_fwd`, `ldlt_app_solve_bwd`) and `spral/src/ssids/cpu/subtree.hxx` (`solve()`) — per-supernode solve kernels and assembly-tree traversal reference. Documented in `/workspace/rivrs-linalg/sparse/dev/references/notes/SPRAL-CODE-REVIEW.md`.

8. **faer integration notes** (`/workspace/rivrs-linalg/sparse/dev/references/notes/FAER-INTEGRATION-NOTES.md`): Section 3.5 covers faer's dense triangular solve API; Section 5 outlines the solve phase flow. Note: faer's sparse triangular solves operate on CSC matrices and are NOT directly applicable — our solve uses dense TRSV within each supernode.

### Open Design Questions

The following areas from the development plan contain potential ambiguity or tension with the current codebase. These should be resolved during the planning/clarification phase:

1. **Scaling integration with factorization (resolved)**: SPRAL passes scaling factors through the call chain down to the assembly kernel (`add_a_block` in `assemble.hxx`), where each entry is scaled at scatter time: `node.lcol[k] = rscale * aval[src] * cscale`. No pre-scaled copy of A is created. Following SPRAL's approach, `AptpNumeric::factor()` will accept an optional scaling slice (`Option<&[f64]>`) and apply it during `scatter_original_entries`. The APTP pivot kernel itself remains scaling-unaware — only the assembly layer touches scaling factors. This avoids the O(nnz) matrix copy while keeping the coupling minimal and localized to a single function.

2. **Rank-deficient matrix handling in solve (resolved)**: SPRAL uses an `options.action` flag. When `action=true` (default), zero pivots are stored as `d[2*i] = 0.0` and the D-solve naturally computes `x[i] *= 0.0`, zeroing the corresponding solution component. The factorization returns a warning flag and records `num_zero` in statistics. When `action=false`, the factorization aborts with an error. Following SPRAL's approach: by default, the solver will allow rank-deficient solves — the D-solve will set `x[i] = 0.0` for zero-pivot entries (modifying `MixedDiagonal::solve_in_place()` to handle zero 1x1 pivots and zero-determinant 2x2 blocks gracefully instead of panicking). The factorization statistics already record `zero_pivots`, and `inertia()` reports zero eigenvalues. Users can inspect these to detect rank deficiency. An explicit "strict" mode (returning an error on singular matrices) may be added to factorization options if needed, but the default is SPRAL's permissive behavior.

3. **Multiple RHS (batch solve)**: The plan mentions "Multiple RHS handled efficiently" and `solve_scratch(rhs_ncols: usize)` suggests matrix-RHS support. However, the notional core solve signature uses a single-column RHS. Assumption: Phase 7 implements single-column solve; multi-column support is a performance optimization deferred to Phase 8.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Backward error `||Ax - b|| / (||A||_F ||x|| + ||b||) < 10^-10` on all hand-constructed test matrices (15 matrices).
- **SC-002**: Backward error `< 10^-10` on all SuiteSparse CI subset matrices (10 matrices) and `>95%` of the full SuiteSparse collection (67 matrices).
- **SC-003**: Median backward error across all test matrices is `< 10^-9`.
- **SC-004**: Inertia computed from the factorization matches the reference inertia for all test matrices where reference data is available.
- **SC-005**: Per-supernode index mapping operations (gather, local_perm application, scatter) produce analytically correct results on hand-constructed test cases.
- **SC-006**: Solve with MC64 scaling produces correct results on all test matrices that use MC64 ordering (backward error < `10^-10` against the original unscaled matrix).
- **SC-007**: Refactoring with a different numeric matrix (same sparsity) produces correct solutions for the new system.
- **SC-008**: All solve operations complete using only the pre-allocated workspace stack (no additional heap allocation in the solve hot path).
- **SC-009**: The three-phase API (analyze → factor → solve), one-shot API, and refactor API all produce correct, equivalent results.
- **SC-010**: Clear, actionable errors are returned for: solve-before-factor and dimension mismatch. Rank-deficient matrices produce a solution (with zeroed components at zero-pivot entries) and are detectable via `stats().zero_pivots > 0`.

### Assumptions

- **Single-column RHS for Phase 7**: Multi-column batch RHS is deferred to Phase 8. Phase 7 supports solving one column at a time; "multiple RHS" means calling solve repeatedly with different vectors reusing the same factorization.
- **Sequential solve**: Parallelism (tree-level or within-supernode) is deferred to Phase 8.2. Phase 7 is single-threaded.
- **Factorization allocates internally**: The solver's `factor()` does not take a workspace stack; frontal matrix allocation remains heap-based per supernode (arena allocation deferred to Phase 9.1).
- **Sparse backward error for validation**: The existing dense `backward_error()` in `validate.rs` is insufficient for large SuiteSparse matrices. Phase 7 adds a sparse-aware variant.
- **Analysis options cover ordering selection**: The analysis options struct must specify the ordering strategy (defaulting to METIS when available) and whether MC64 scaling is requested. This design will be refined during planning.
