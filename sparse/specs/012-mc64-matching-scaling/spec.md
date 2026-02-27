# Feature Specification: MC64 Matching & Scaling

**Feature Branch**: `012-mc64-matching-scaling`
**Created**: 2026-02-13
**Status**: Draft
**Input**: Phase 4.2 of SSIDS plan — weighted bipartite matching and symmetric scaling for indefinite matrix preprocessing

## Clarifications

### Session 2026-02-13

- Q: When the input matrix is structurally singular (no perfect matching exists), should MC64 return a partial matching or an error? → A: Return partial matching with Duff-Pralet scaling correction (domain best practice from Duff & Pralet 2005 Section 4.2.3 and SPRAL `match_order.f90`). Structural singularity is a warning, not an error.
- Q: What should MC64's permutation output represent — a matching, a symmetric ordering, or a combined matching+ordering? → A: Matching and scaling only. MC64's matching permutation is NOT a fill-reducing symmetric ordering. Fill-reducing ordering (METIS) operates independently on the original graph. Condensed ordering (SPRAL-style `mo_split()` composing matching + METIS via graph condensation) is deferred as a lightweight follow-up (~150-200 lines bookkeeping). MC64 is validated in isolation first, following SPRAL's own testing pattern where matching is tested independently in `tests/scaling.f90`.

## Context

MC64 matching and scaling is Phase 4.2 of the SSIDS development plan. It complements Phase 4.1 (METIS fill-reducing ordering, complete) by providing numerical preprocessing for indefinite matrices:

- **Phase 4.1 (METIS)**: Reduces fill-in in the sparsity structure (structural quality)
- **Phase 4.2 (MC64)**: Improves diagonal dominance for numerical stability (numerical quality)

Together, these form SPRAL's `ordering=2` strategy: matching followed by fill-reducing ordering. In SPRAL, the matching and ordering are bundled into a single `match_order_metis()` call that includes a graph condensation step (collapsing matched 2-cycle pairs into single nodes before running METIS). This phase implements matching and scaling independently of ordering; condensed ordering composition is a lightweight follow-up (see Assumptions).

### Algorithm References

- **Duff & Koster (2001)**: "On algorithms for permuting large entries to the diagonal of a sparse matrix" — core MC64 weighted matching algorithm (shortest augmenting paths with Dijkstra on reduced costs). `/workspace/rivrs-linalg/references/ssids/duff2001.md`
- **Duff & Koster (1999)**: "The design and use of algorithms for permuting large entries to the diagonal of sparse matrices" — bipartite matching foundations, transversal theory. `/workspace/rivrs-linalg/references/ssids/duff1999.md`
- **Duff & Pralet (2005)**: "Strategies for scaling and pivoting for sparse symmetric indefinite problems" — symmetric adaptation of MC64 (MC64SYM), scaling symmetrization, experimental results showing 2-100x speedups on difficult indefinite problems. `/workspace/rivrs-linalg/references/ssids/duff2005.md`
- **Duff & Lopez (2020)**: APTP algorithm and multifrontal factorization context. `/workspace/rivrs-linalg/references/ssids/duff2020.md`

### SPRAL Implementation References

- **`/opt/references/spral/src/scaling.f90`**: Hungarian matching algorithm (`hungarian_match`, `hungarian_scale_sym`). Core weighted matching implementation.
- **`/opt/references/spral/src/match_order.f90`**: Combined matching+ordering pipeline (`match_order_metis`). Contains `mo_match()` (matching wrapper), `mo_split()` (graph condensation — ~115 lines, purely bookkeeping: cycle detection → index remapping → condensed CSC construction → METIS on condensed graph → expansion back to original indices). Reference for future condensed ordering work.
- **`/opt/references/spral/tests/scaling.f90`**: Standalone matching tests — 4 Hungarian test subroutines that validate matching + scaling independently of ordering. Pattern for our test design.

### Codebase Status

- **Phase 3 (AptpSymbolic)**: Complete. Accepts `SymmetricOrdering::Custom(PermRef)` for external orderings.
- **Phase 4.1 (METIS)**: Complete. `metis_ordering()` in `src/aptp/ordering.rs` returns `Perm<usize>`.
- **Phase 2 (Data structures)**: Complete. `MixedDiagonal`, `Inertia`, `PivotType` available.
- **Error types**: `SparseError` with `InvalidInput`, `AnalysisFailure`, `StructurallySingular` variants.
- **faer has no matching or scaling functionality** — MC64 is genuinely new code, not a wrapper around faer.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Preprocess an Indefinite Matrix for Better Pivoting (Priority: P1)

A solver developer has a sparse symmetric indefinite matrix that produces many delayed pivots during APTP factorization. They want to preprocess the matrix by permuting large entries onto the diagonal and computing symmetric scaling factors, so that the subsequent factorization encounters fewer numerical difficulties.

**Why this priority**: This is the core purpose of MC64 in the APTP pipeline. Without matching and scaling, indefinite matrices with small or zero diagonal entries cause excessive pivot delays, increasing memory usage and factorization time. Duff & Pralet (2005) report 2-100x speedups on difficult problems.

**Independent Test**: Can be tested by applying MC64 to known indefinite test matrices and verifying that (a) a valid matching is produced, (b) positive scaling factors are computed, and (c) the scaled matrix has improved diagonal dominance. Following SPRAL's testing pattern (`tests/scaling.f90`), matching correctness is validated through the scaling properties it produces, not in isolation.

**Acceptance Scenarios**:

1. **Given** a sparse symmetric indefinite matrix with poor diagonal dominance, **When** MC64 matching is computed with the maximum-product objective, **Then** a matching is returned that identifies large-magnitude entry pairs, and symmetric scaling factors are returned such that the scaled matrix has unit diagonal entries and off-diagonal entries with absolute value at most 1.
2. **Given** the same matrix, **When** the diagonal dominance metric is measured before and after scaling, **Then** the scaled matrix shows improved diagonal dominance (ratio of diagonal magnitude to off-diagonal row sum).
3. **Given** the matching result, **When** the matching structure is examined, **Then** it decomposes into singletons (σ(i)=i, suitable for 1x1 pivots) and 2-cycles (σ(i)=j, σ(j)=i, suitable for 2x2 pivot blocks), with no longer cycles. This classification is useful metadata for downstream phases.

---

### User Story 2 — Use MC64 Scaling with Independent METIS Ordering (Priority: P2)

A solver developer wants to use MC64 scaling for numerical quality and METIS ordering for fill reduction — the two preprocessing strategies operating independently. MC64 produces scaling factors (consumed by numeric factorization in Phase 5), while METIS produces a fill-reducing ordering (consumed by symbolic analysis in Phase 3). The matching and ordering are independent: METIS orders the original graph, not a condensed graph.

**Why this priority**: This validates that MC64 and METIS compose cleanly as independent preprocessing steps. SPRAL bundles these via graph condensation (`mo_split()` in `match_order.f90`), where matched 2-cycle pairs are collapsed into single nodes before METIS ordering. That condensation is a quality optimization (~150-200 lines of bookkeeping) that keeps paired entries adjacent in the elimination order; it is deferred as a follow-up once MC64 is validated independently.

**Independent Test**: Can be tested by computing MC64 scaling on a test matrix, independently computing METIS ordering on the same matrix, feeding the METIS ordering into symbolic analysis, and verifying that the scaling factors are valid for the independently-ordered matrix.

**Acceptance Scenarios**:

1. **Given** a sparse symmetric indefinite matrix, **When** MC64 matching + scaling is computed and METIS ordering is computed independently on the original matrix, **Then** the METIS ordering produces valid fill estimates via `AptpSymbolic::analyze`, and the MC64 scaling factors remain valid (the scaling is structure-independent).
2. **Given** the MC64 matching result, **When** the 2-cycle pair count is examined, **Then** this count provides useful metadata for assessing whether future condensed ordering would improve quality (many 2-cycles → condensation likely beneficial).

---

### User Story 3 — Validate Matching on SuiteSparse Test Suite (Priority: P3)

A solver developer wants confidence that MC64 produces correct and useful matchings across a broad range of real-world sparse matrices from the SuiteSparse collection, including matrices with varying sizes, densities, and spectral properties.

**Why this priority**: Broad validation ensures the implementation is robust and not over-fitted to specific matrix structures. The existing SuiteSparse test infrastructure (67 matrices, CI subset of 10) can be reused from Phase 4.1.

**Independent Test**: Can be tested by running MC64 on all indefinite matrices in the SuiteSparse collection and verifying matching validity, scaling positivity, and diagonal dominance improvement.

**Acceptance Scenarios**:

1. **Given** the full SuiteSparse indefinite test matrix collection, **When** MC64 matching is computed on each matrix, **Then** a valid matching is produced for every matrix (matching count equals matrix dimension for structurally nonsingular matrices, and matching count equals structural rank for singular matrices).
2. **Given** each matching result, **When** scaling factors are verified, **Then** all scaling factors are positive and finite.
3. **Given** a structurally singular test matrix, **When** MC64 matching is computed, **Then** a partial matching is returned with `matched < n`, unmatched indices are placed at the end of the permutation, and the Duff-Pralet scaling correction is applied to unmatched rows.

---

### Edge Cases

- What happens when the matrix is structurally singular (no full matching exists)? Per domain best practice (Duff & Pralet 2005 Section 4.2.3; SPRAL `match_order.f90`), structural singularity is a **warning, not an error**. MC64 returns a partial matching (`matched < n`), extracts the structurally nonsingular submatrix (Property 4.2), re-runs the weighted matching on the nonsingular block, applies the Duff-Pralet scaling correction for unmatched indices (scale based on largest entry connected to matched set, or 1.0 if isolated), and places unmatched rows/columns at the end of the permutation. The `matched` count in the result distinguishes full from partial matchings.
- What happens when the matrix has zero entries on the diagonal? MC64 should still find a matching that permutes nonzero entries onto the diagonal if one exists.
- What happens when the matrix has entries spanning many orders of magnitude (badly scaled)? The logarithmic cost formulation (`c_ij = log|a_ij|`) should handle this naturally, but scaling factors should remain finite and well-conditioned.
- What happens when the matrix is already well-scaled with dominant diagonal? MC64 should produce an identity or near-identity permutation and scaling factors near unity, without degrading the matrix.
- What happens for a 1x1 or 2x2 matrix? Matching should be trivial and correct.
- What happens when all off-diagonal entries are zero (diagonal matrix)? Identity matching and unit scaling.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST compute a weighted bipartite matching that maximizes the product of diagonal entry magnitudes for a given sparse symmetric matrix — the "maximum product" objective from Duff & Koster (2001, Section 4).
- **FR-002**: System MUST compute symmetric scaling factors from the matching dual variables, such that the scaled permuted matrix has diagonal entries of absolute value 1 and all entries with absolute value at most 1 — the MC64SYM symmetrization from Duff & Pralet (2005, Property 4.1).
- **FR-003**: System MUST accept a sparse symmetric matrix with numeric values as input (not just the symbolic sparsity pattern), since the matching depends on entry magnitudes.
- **FR-004**: System MUST return the matching as a permutation vector (row i is matched to column σ(i)). For symmetric matrices, the matching decomposes into singletons (σ(i)=i) and 2-cycles (σ(i)=j, σ(j)=i) — no longer cycles occur. Note: this matching permutation is NOT a fill-reducing symmetric ordering and should NOT be used directly with `SymmetricOrdering::Custom`. Fill-reducing ordering is computed independently by METIS. A future condensed ordering utility will compose matching + METIS properly.
- **FR-005**: System MUST return scaling factors as a separate output, suitable for later consumption by the numeric factorization phase (Phase 5). Scaling does not affect sparsity structure and is not consumed by the symbolic phase.
- **FR-006**: System MUST report the number of entries successfully matched, to allow callers to distinguish full matchings from partial matchings.
- **FR-007**: System MUST validate that the input matrix is square and has a symmetric sparsity pattern (consistent with existing validation patterns in the codebase).
- **FR-008**: System MUST handle matrices where the existing diagonal contains zero entries, finding a matching that permutes nonzeros onto the diagonal when structurally possible.
- **FR-009**: System MUST preserve the pattern that faer types are used at the boundary — permutations as `Perm<usize>`, no custom permutation wrappers.
- **FR-010**: System MUST handle structurally singular matrices by returning a partial matching with the Duff-Pralet scaling correction (Duff & Pralet 2005, Section 4.2.3). Specifically: (a) return a partial matching with `matched < n`, (b) extract the structurally nonsingular submatrix via Property 4.2, (c) compute scaling factors for matched indices from the weighted matching on the nonsingular block, (d) compute scaling factors for unmatched indices based on their largest entry connected to the matched set (default to 1.0 if isolated), and (e) place unmatched rows/columns at the end of the permutation. This follows SPRAL's behavior where structural singularity produces a warning (`SSIDS_WARNING_ANAL_SINGULAR`), not an error.

### Key Entities

- **Matching Result**: The output of the MC64 algorithm, comprising a matching (as a permutation vector), symmetric scaling factors, and the count of matched entries. The matching classifies pivot structure (singletons vs 2-cycle pairs). The scaling feeds into numeric factorization (Phase 5). The matching does NOT serve as a fill-reducing ordering — that role belongs to METIS (Phase 4.1).
- **Matching Objective**: The optimization criterion for the matching. The primary objective is maximizing the product of diagonal entry magnitudes (equivalent to minimizing the sum of `-log|a_ij|` costs). A sum-maximization variant may also be useful.
- **Scaling Factors**: A vector of positive real values `s_i` such that the symmetrically scaled matrix `A_scaled[i,j] = s_i * A[i,j] * s_j` has improved diagonal dominance. For MC64SYM, these are derived from the dual variables of the assignment problem: `s_i = exp(-(u_i + v_i) / 2)` where `u, v` are the row and column duals from the minimum-cost assignment (Duff & Pralet 2005, Property 4.1).
- **Matching Cycle Structure**: For symmetric matrices, the matching decomposes into singletons (1-cycles, σ(i)=i) and 2-cycles (σ(i)=j, σ(j)=i). Singletons indicate 1x1 pivot candidates; 2-cycles indicate 2x2 pivot block candidates. This metadata is useful for Phase 5/6 (pivot pre-classification) and for future condensed ordering (2-cycle count indicates potential benefit of graph condensation).

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: MC64 produces a valid full matching (matched count equals matrix dimension) on all structurally nonsingular indefinite matrices in the SuiteSparse CI subset (10 matrices).
- **SC-002**: Scaling factors are positive and finite for all test matrices. The scaled permuted matrix has diagonal entries with absolute value 1 (within floating-point tolerance of 1e-12).
- **SC-003**: Diagonal dominance (measured as minimum ratio of `|d_ii| / sum_j |a_ij|` for j != i) improves or stays the same on at least 80% of indefinite test matrices after MC64 preprocessing.
- **SC-004**: The MC64 matching decomposes into only singletons and 2-cycles (no longer cycles) for all symmetric test matrices.
- **SC-005**: MC64 scaling and METIS ordering compose independently: METIS ordering on the original graph produces valid fill estimates, and the MC64 scaling factors remain valid regardless of the ordering applied.
- **SC-006**: MC64 completes on all CI subset matrices (largest ~100K dimension) without timeout or excessive memory usage.

## Testing Approach

Following SPRAL's testing pattern from `tests/scaling.f90`, matching correctness is validated **through the scaling properties it produces**, not by testing matching structure in isolation.

### SPRAL's Validation Criteria (numerical specifics)

SPRAL's Hungarian matching tests check these properties on the scaled matrix:

1. **Scaling constraint**: All entries of the scaled matrix satisfy `|s_i * a_ij * s_j| <= 1.0`
2. **Row maximum constraint**: For each row, `max_j |s_i * a_ij * s_j| >= 0.75` (SPRAL uses 1.0 - 0.25 tolerance for symmetric, looser for unsymmetric)
3. **Matching validity**: Every matched pair (i, σ(i)) corresponds to an actual nonzero entry in the matrix; no column is matched twice
4. **Matching cardinality**: Hungarian matching achieves `matched == n` for structurally nonsingular matrices (unlike auction matching, which only requires ≥ 90%)

These criteria validate the matching implicitly: if the scaling satisfies the constraints, the underlying matching must be correct.

### Our Validation Strategy

- **Hand-constructed matrices**: Exact scaling property verification (constraints 1-4 above)
- **SuiteSparse CI subset**: Same scaling property checks at scale; additionally verify matching cycle structure (singletons + 2-cycles only)
- **Full SuiteSparse**: `#[ignore]` tests for comprehensive validation before phase completion

## Assumptions

- **A-001**: The primary matching objective is maximum product (SPRAL default). Sum-maximization is a secondary objective that may be deferred if it adds significant scope.
- **A-002**: Symmetric scaling is computed using the MC64SYM approach from Duff & Pralet (2005): apply the unsymmetric MC64 algorithm to the symmetric matrix, then symmetrize the scaling factors via geometric mean of row and column duals.
- **A-003**: The scaling vector is stored and returned but not applied during this phase. The actual scaling application (scale before factorize, unscale after solve) is deferred to Phase 5 (numeric factorization).
- **A-004**: MC64 is implemented from academic references (Duff & Koster 2001, Duff & Pralet 2005) and SPRAL (BSD-3) patterns, maintaining clean room status for Apache-2.0 licensing.
- **A-005**: No external C/Fortran library is used for MC64 (unlike Phase 4.1 which used `metis-sys`). MC64 is implemented in pure Rust as the algorithm is well-documented and no permissive-licensed C implementation is readily available.
- **A-006**: The existing test matrix infrastructure (registry, loading, filtering) from Phase 0.4/1.1 is reused for MC64 testing.
- **A-007**: MC64 matching and fill-reducing ordering (METIS) operate independently in this phase. MC64 produces scaling + matching metadata; METIS produces ordering on the original graph. A condensed ordering utility (composing matching + METIS via graph condensation, as SPRAL does in `mo_split()`) is a lightweight follow-up (~150-200 lines of index bookkeeping, O(nnz + n log n)) that can be added without changing MC64's API. The condensation is deferred because: (a) MC64 is the hard algorithmic work and should be validated in isolation first; (b) condensation only benefits matrices with many 2-cycle pairs (augmented/saddle-point systems); (c) APTP's pivot delay mechanism handles arbitrary pivot patterns regardless of ordering quality. Reference for future condensation: SPRAL's `mo_split()` in `/opt/references/spral/src/match_order.f90` (lines 220-396).

## Dependencies

- **Phase 3 (AptpSymbolic)**: Complete. Provides `SymmetricOrdering::Custom` integration point (used by METIS, not directly by MC64).
- **Phase 4.1 (METIS ordering)**: Complete. Provides `metis_ordering()` for independent fill-reducing ordering.
- **Phase 2 (Data structures)**: Complete. Provides `perm_from_forward()` utility.
- **Phase 0.4/1.1 (Test infrastructure)**: Complete. Provides test matrix loading, validation, and filtering.
- **Phase 5 (Numeric factorization)**: NOT yet started. Scaling consumption API will be defined in Phase 5. This phase only produces and validates the scaling factors.
