# Research: MC64 Matching & Scaling

**Feature**: 012-mc64-matching-scaling
**Date**: 2026-02-13

## Decision 1: Core Algorithm Selection

**Decision**: Implement the Duff & Koster (2001) Algorithm MPD — Dijkstra-based shortest augmenting path on reduced costs with greedy initialization heuristic.

**Rationale**: This is the standard MC64 algorithm used in SPRAL (`scaling.f90::hungarian_match`), HSL (MC64), and MUMPS. It has O(n(tau + n) log n) complexity and is well-documented with pseudocode in the reference literature. The greedy initialization typically achieves ~80% matching cardinality before the main loop, reducing practical cost significantly.

**Alternatives considered**:
- Auction algorithm (Bertsekas): SPRAL implements this too (`auction_match`) but only guarantees 90% matching cardinality. Hungarian achieves 100% for nonsingular matrices. Auction is faster on dense problems but less suitable for our sparse use case.
- Kuhn-Munkres (Hungarian): O(n^3) complexity, too expensive for large sparse matrices.
- Hopcroft-Karp: Finds maximum cardinality matching but not maximum weight. MC64 needs weighted matching.

## Decision 2: Symmetric Scaling Approach

**Decision**: Use the MC64SYM approach from Duff & Pralet (2005): apply the unsymmetric MC64 algorithm to the symmetric matrix (treating it as a bipartite graph), then symmetrize the scaling factors via `D[i] = exp(-(u[i] + v[i]) / 2)` where u, v are the row and column dual variables.

**Rationale**: This is exactly what SPRAL does in `match_order.f90::mo_match()`. The symmetrized scaling satisfies Property 4.1: scaled matrix DAD has unit diagonal entries and all off-diagonal entries <= 1. The alternative (symmetric matching on the undirected graph) is more complex and doesn't improve quality for our use case.

**Alternatives considered**:
- Symmetric matching via graph doubling (Duff & Pralet 2005, Section 5.1): Creates a doubled graph and finds matching on the undirected representation. More theoretically elegant but significantly more complex. SPRAL uses the simpler MC64SYM approach.
- MC30 (symmetric scaling without matching): Minimizes sum of (log|a_ij|)^2 iteratively. Simpler but doesn't produce a matching for pivot classification.
- MC77 (norm equilibration): Iterative scaling to doubly stochastic form. Complementary but doesn't maximize diagonal dominance.

## Decision 3: Input Format

**Decision**: Accept `&SparseColMat<usize, f64>` (upper-triangular CSC, faer convention for symmetric matrices). Internally expand to full bipartite graph for the matching algorithm.

**Rationale**: All existing APIs in the codebase (`AptpSymbolic::analyze`, `metis_ordering`) accept the upper-triangular representation. MC64 needs the full bipartite graph internally, but the expansion is a standard O(nnz) operation that `ordering.rs` already implements for METIS (the `extract_adjacency` helper). Requiring callers to provide full symmetric storage would break the faer convention.

**Alternatives considered**:
- Accept full symmetric CSC: Would require callers to preprocess their matrices. Inconsistent with faer convention and existing APIs.
- Accept `SymbolicSparseColMatRef` + separate values: Over-complicated. MC64 needs values tightly coupled with structure for the cost matrix computation.

## Decision 4: Scaling Storage Domain

**Decision**: Store scaling factors in log domain internally (`log_scale[i] = -(u[i] + v[i]) / 2`) and return in linear domain (`scale[i] = exp(log_scale[i])`). Expose a method to access log-domain scaling for callers that need it.

**Rationale**: SPRAL stores scaling in log domain internally (`match_order.f90` lines 458-484) and only converts to linear domain at the boundary. Log domain avoids overflow/underflow for extremely badly-scaled matrices where entries span hundreds of orders of magnitude. The dual variables u, v are naturally in log domain (the cost matrix uses `log|a_ij|`). However, the downstream consumer (Phase 5 numeric factorization) will apply `s_i * a_ij * s_j`, which is most natural in linear domain.

**Alternatives considered**:
- Linear domain only: Risk of overflow/underflow for extreme scaling factors. Would need special handling.
- Log domain only: Would require downstream consumers to work in log domain. Unnatural for Phase 5.

## Decision 5: Priority Queue Implementation

**Decision**: Use a standard binary heap (Rust's `BinaryHeap` from std) for the Dijkstra frontier. Consider the Q1/Q2 split optimization from Duff & Koster (2001, Section 4) as a follow-up if profiling reveals heap operations as a bottleneck.

**Rationale**: The Q1/Q2 split (Q1 = array of nodes within tolerance of d_min, Q2 = binary heap for the rest) reduces heap operations by 30-50% empirically. However, Rust's `BinaryHeap` is already well-optimized, and the added complexity of maintaining two data structures is not justified until profiling shows it's needed. The greedy initialization already eliminates ~80% of augmentation iterations, so the Dijkstra inner loop runs fewer times than worst case.

**Alternatives considered**:
- Fibonacci heap: O(1) amortized decrease-key, but high constant factors and complex implementation. Not worth it for sparse matrices where n << tau.
- Q1/Q2 split from day one: Premature optimization. Binary heap is simpler and likely fast enough.

## Decision 6: Structural Singularity Handling

**Decision**: Follow Duff & Pralet (2005, Section 4.2.3) and SPRAL (`match_order.f90`). When the initial matching has cardinality < n: (1) extract the structurally nonsingular submatrix via Property 4.2, (2) re-run MC64 on the submatrix, (3) apply the Duff-Pralet scaling correction for unmatched indices.

**Rationale**: This is the established domain best practice. SPRAL reports structural singularity as `WARNING_SINGULAR` (not an error), allowing the solver to continue with partial information. The Duff-Pralet correction ensures unmatched rows get reasonable scaling factors based on their connections to the matched set.

**Alternatives considered**:
- Return error on structural singularity: Overly restrictive. Structurally singular matrices arise naturally in saddle-point systems (optimization, Stokes equations).
- Ignore unmatched rows (scale to 1.0): Simpler but loses the benefit of scaling for those rows.

## Decision 7: Module Structure

**Decision**: Single file `src/aptp/matching.rs` containing all MC64 functionality: public API, internal algorithm, cost matrix computation, matching augmentation, and symmetric scaling. Add integration tests in `tests/mc64_matching.rs`.

**Rationale**: Follows the pattern established by `ordering.rs` (545 lines for METIS). MC64 is estimated at 600-900 lines including internal helpers, which is manageable in a single file. The algorithm is a cohesive unit — splitting it across files would add complexity without improving readability. If the file grows beyond ~1000 lines, internal helpers can be extracted later.

**Alternatives considered**:
- Multiple files (mc64.rs, scaling.rs, dijkstra.rs): Over-modularized for the scope. The Dijkstra kernel is tightly coupled with the matching state.
- Separate crate: Overkill. MC64 is part of the APTP preprocessing pipeline.
