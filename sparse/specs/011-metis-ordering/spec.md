# Feature Specification: METIS Nested Dissection Ordering

**Feature Branch**: `011-metis-ordering`
**Created**: 2026-02-12
**Status**: Draft
**Input**: User description: "implement Phase 4.1 in ssids-plan.md"

## Clarifications

### Session 2026-02-12

- Q: Should METIS be a required or optional (feature-flagged) dependency? → A: Required dependency. The `metis` Rust crate vendors METIS C source by default (compiled via `build.rs`), so no system library installation is needed — only a C compiler, which is already required. No feature flag necessary.

## Context

This feature implements Phase 4.1 of the SSIDS development plan. Phase 3 (AptpSymbolic analysis) is complete and accepts fill-reducing orderings as input via faer's `SymmetricOrdering` enum. However, the only ordering currently available is AMD (Approximate Minimum Degree), which produces catastrophically poor fill on many real-world matrices — 2-20x more nonzeros in L than METIS on matrices with geometric structure (FEM, quantum chemistry, optimization).

SPRAL's SSIDS uses METIS as its default ordering (`options%ordering = 1`), and all benchmark results in the primary SSIDS paper (Hogg et al. 2016) use METIS v4. Without METIS, the solver cannot match published performance baselines and is effectively limited to small or structurally simple matrices.

The existing `AptpSymbolic::analyze()` function already accepts `SymmetricOrdering::Custom(perm)`, so this feature requires no API changes to the symbolic analysis pipeline — it produces a fill-reducing permutation that feeds into the existing interface.

### Algorithm References

- **George (1973)** — "Nested Dissection of a Regular Finite Element Mesh", SIAM J. Numer. Anal. — foundational theory proving O(n^3) optimality of nested dissection for 2D meshes.
  Reference: `/workspace/rivrs-linalg/references/ssids/george1973.md`

- **Karypis & Kumar (1998)** — "A Fast and High Quality Multilevel Scheme for Partitioning Irregular Graphs", SIAM J. Sci. Comput. — the METIS algorithm: coarsening, initial partitioning, and uncoarsening with Kernighan-Lin refinement.

- **Hogg, Ovtchinnikov & Scott (2016)** — "A Sparse Symmetric Indefinite Direct Solver for GPU Architectures", ACM TOMS 42(1) — SSIDS paper using METIS v4 ordering; Table III provides benchmark nnz(L) values.
  Reference: `/workspace/rivrs-linalg/references/ssids/hogg2016.md`

- **Davis (2016)** — "Direct Methods for Sparse Linear Systems" — comprehensive survey of nested dissection and METIS in Section 8.6.
  Reference: `/workspace/rivrs-linalg/references/ssids/davis2016.md`

### Current State

- `AptpSymbolic::analyze(matrix, SymmetricOrdering::Custom(perm))` is implemented and tested (Phase 3)
- `perm_from_forward()` helper converts a forward permutation array to `Perm<usize>` (Phase 2)
- Full SuiteSparse test suite exists with `MAX_DIM_FOR_AMD = 30_000` guard that should be removable after METIS is integrated
- AMD is the only ordering currently available; no METIS dependency exists yet

## User Scenarios & Testing

### User Story 1 — Compute METIS Ordering for a Symmetric Sparse Matrix (Priority: P1)

A solver developer has a symmetric sparse matrix in CSC format and needs a fill-reducing ordering that produces significantly less fill than AMD, particularly for matrices with geometric structure (e.g., from finite element discretizations). They call a METIS ordering function, receive a valid permutation, and pass it to `AptpSymbolic::analyze()` via `SymmetricOrdering::Custom`.

**Why this priority**: This is the core deliverable — without it, the SSIDS solver cannot match published performance on real-world matrices. Every downstream phase (numeric factorization, solve) depends on ordering quality.

**Independent Test**: Can be fully tested by computing a METIS ordering on any SuiteSparse matrix and verifying the result is a valid permutation of the correct dimension.

**Acceptance Scenarios**:

1. **Given** a symmetric sparse matrix in CSC format, **When** the METIS ordering function is called, **Then** it returns a valid permutation where every index 0..n-1 appears exactly once.
2. **Given** a symmetric sparse matrix, **When** METIS ordering is computed and passed to `AptpSymbolic::analyze()` as `SymmetricOrdering::Custom`, **Then** the symbolic analysis completes successfully and reports predicted fill statistics.
3. **Given** the `sparsine` matrix (50K dimension), **When** METIS ordering is used instead of AMD, **Then** the predicted nnz(L) is substantially reduced (historically 10-20x less fill with METIS).

---

### User Story 2 — Reproduce Published Fill Predictions (Priority: P2)

A solver developer wants confidence that the METIS integration produces orderings of comparable quality to METIS v4 as used in the SSIDS reference paper. They run symbolic analysis with METIS ordering on the SuiteSparse benchmark matrices and compare predicted nnz(L) against the values reported in Hogg et al. (2016) Table III.

**Why this priority**: Validates that the METIS integration is correct and produces competitive orderings. Without this, we cannot claim equivalence with the published SSIDS results.

**Independent Test**: Can be tested by running symbolic analysis with METIS on a subset of Table III matrices and comparing predicted nnz(L) against paper-reported values.

**Acceptance Scenarios**:

1. **Given** a matrix from Hogg et al. (2016) Table III (e.g., `bmwcra_1`), **When** METIS ordering is computed and symbolic analysis is performed, **Then** predicted nnz(L) is within 20% of the paper-reported value.
2. **Given** the full set of SuiteSparse CI-subset matrices, **When** METIS ordering is used, **Then** METIS produces equal or less fill than AMD on at least 80% of matrices.

---

### User Story 3 — Analyze Large Matrices Without AMD Limitations (Priority: P3)

A solver developer wants to run symbolic analysis on all 67 SuiteSparse benchmark matrices without the current `MAX_DIM_FOR_AMD = 30_000` guard. With METIS ordering, large matrices that AMD handles poorly should produce reasonable fill predictions, allowing the full test suite to run end-to-end.

**Why this priority**: Unblocks full-collection testing and removes an artificial limitation. The current guard exists solely because AMD ordering quality degrades catastrophically on large matrices.

**Independent Test**: Can be tested by running the full SuiteSparse symbolic analysis test suite (currently `#[ignore]`d tests) without the dimension guard and verifying all matrices complete successfully.

**Acceptance Scenarios**:

1. **Given** the full SuiteSparse collection (67 matrices), **When** symbolic analysis is run with METIS ordering, **Then** all matrices complete successfully without the `MAX_DIM_FOR_AMD` guard.
2. **Given** the full SuiteSparse collection, **When** symbolic analysis is run with METIS ordering, **Then** the total wall-clock time is less than 2 minutes.

---

### Edge Cases

- What happens when the input matrix has dimension 0 or 1? The ordering function should handle trivially-sized matrices gracefully (identity permutation or equivalent).
- What happens when the input matrix is diagonal (no off-diagonal entries)? METIS should still produce a valid permutation, though the graph has no edges.
- What happens when the input matrix is dense (fully connected graph)? METIS should still produce a valid ordering, though nested dissection offers no fill advantage for dense matrices.
- What happens if the graph is disconnected (block-diagonal matrix)? METIS should handle disconnected components and produce a valid permutation for the full matrix.

## Requirements

### Functional Requirements

- **FR-001**: The system MUST provide a function that accepts a symmetric sparse matrix (CSC format) and returns a fill-reducing permutation computed via METIS nested dissection.
- **FR-002**: The returned permutation MUST be a valid permutation of {0, 1, ..., n-1} for an n-dimensional matrix.
- **FR-003**: The returned permutation MUST be directly usable with the existing `AptpSymbolic::analyze()` function via `SymmetricOrdering::Custom` — no changes to the symbolic analysis API.
- **FR-004**: The system MUST extract the adjacency structure from the symmetric sparse matrix for METIS consumption, using only the structural pattern (not numerical values).
- **FR-005**: The system MUST handle matrices of all dimensions representable in the test suite (up to ~150K, the largest SuiteSparse benchmark matrix), without artificial dimension guards.
- **FR-006**: The system MUST return an appropriate error for invalid inputs (e.g., non-square matrix, structural issues that prevent METIS from computing an ordering).
- **FR-007**: The system MUST handle edge cases gracefully: dimension-0 matrices, dimension-1 matrices, diagonal matrices (no off-diagonal structure), and disconnected graphs.

### Key Entities

- **METIS Permutation**: A fill-reducing reordering of matrix rows/columns, represented as faer's `Perm<usize>`. Computed from the graph structure of the symmetric matrix via multilevel nested dissection.
- **Adjacency Structure**: The graph representation of the symmetric sparse matrix extracted for METIS: each matrix row/column is a vertex, each off-diagonal nonzero is an edge. Self-loops (diagonal entries) are excluded.

## Success Criteria

### Measurable Outcomes

- **SC-001**: METIS ordering produces valid permutations on 100% of the 82 test matrices (15 hand-constructed + 10 CI-subset + 67 full SuiteSparse, where available).
- **SC-002**: Predicted nnz(L) using METIS ordering is within 20% of Hogg et al. (2016) Table III reported values for all overlapping benchmark matrices.
- **SC-003**: METIS ordering produces equal or less predicted fill than AMD on at least 80% of the SuiteSparse matrix collection.
- **SC-004**: Full SuiteSparse symbolic analysis (67 matrices) with METIS ordering completes in under 2 minutes total wall-clock time.
- **SC-005**: Integration requires no changes to the `AptpSymbolic::analyze()` function signature — METIS ordering is consumed via the existing `SymmetricOrdering::Custom` mechanism.
- **SC-006**: The `MAX_DIM_FOR_AMD` guard in the full SuiteSparse tests can be removed (or made ordering-conditional) after METIS is integrated.

## Assumptions

- METIS is a required (non-optional) dependency. The `metis` Rust crate (v0.2.x) vendors METIS 5.x C source and compiles it automatically via `build.rs` — no system library installation needed, only a C compiler.
- The `metis` crate is Apache-2.0 compatible (dual MIT/Apache-2.0 wrapper, Apache-2.0 METIS source).
- Differences in predicted nnz(L) between our METIS v5 and the paper's METIS v4 are expected to be small (within the 20% tolerance specified in SC-002). METIS v5 generally produces orderings of comparable or better quality than v4.
- The ordering function only needs to work with `f64` matrices (the only numeric type currently supported by the sparse solver).
- METIS options (e.g., number of separators, compression, etc.) will use sensible defaults initially. Configuration of METIS parameters is not in scope for this feature but could be added later.

## Dependencies

- **Phase 3 (AptpSymbolic)**: Complete. Provides the `SymmetricOrdering::Custom` integration point.
- **Phase 2 (APTP Data Structures)**: Complete. Provides `perm_from_forward()` for constructing `Perm<usize>` from forward arrays.
- **Existing test infrastructure**: Full SuiteSparse test suite and CI-subset tests are in place.
- **External dependency**: `metis` Rust crate (v0.2.x) — required dependency, vendors METIS C source (no system library needed).

## Out of Scope

- MC64 matching and scaling (Phase 4.2 — separate feature).
- Combined matching + METIS ordering pipeline (Phase 4.2 scope).
- METIS parameter tuning or configuration API (could be added in a future enhancement).
- Parallel ordering computation (METIS is already internally parallelized where beneficial).
- Alternative graph partitioning libraries (Scotch, etc.) — METIS only for now.
- Changes to the `AptpSymbolic::analyze()` API signature.
