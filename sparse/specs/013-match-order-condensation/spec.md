# Feature Specification: Match-Order Condensation Pipeline

**Feature Branch**: `013-match-order-condensation`
**Created**: 2026-02-13
**Status**: Draft
**Input**: User description: "Implement SPRAL-style condensed matching-ordering pipeline combining MC64 matching with METIS ordering via cycle condensation"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Combined Matching-Ordering with Guaranteed Pair Adjacency (Priority: P1)

A solver developer requests a fill-reducing ordering for a sparse symmetric indefinite matrix that respects the MC64 matching structure. The system computes MC64 matching, splits the matching into singletons and 2-cycles, builds a condensed graph where each 2-cycle pair becomes a single super-node, runs METIS on the smaller condensed graph, and expands the result back to a full-size ordering where matched 2-cycle pairs are guaranteed to be adjacent in the elimination order.

**Why this priority**: This is the core value of the feature. Without guaranteed adjacency, the APTP kernel (Phase 5) cannot efficiently exploit 2-cycle matches as natural 2x2 pivot candidates. This is the primary reason SPRAL uses condensation rather than naive permutation composition.

**Independent Test**: Can be validated by computing a condensed ordering on any indefinite SuiteSparse matrix and verifying that every 2-cycle pair in the matching occupies consecutive positions in the output ordering.

**Acceptance Scenarios**:

1. **Given** a symmetric indefinite matrix with MC64 matching containing 2-cycles, **When** the condensed match-order pipeline is invoked, **Then** every 2-cycle pair (i,j) occupies consecutive positions in the output ordering.
2. **Given** a symmetric positive definite matrix (matching is all singletons), **When** the condensed pipeline is invoked, **Then** the output is equivalent to running METIS directly on the original graph.
3. **Given** a structurally nonsingular matrix, **When** the pipeline completes, **Then** the output ordering is a valid permutation of {0, ..., n-1} and can be passed to `AptpSymbolic::analyze()` via `SymmetricOrdering::Custom`.

---

### User Story 2 - Structurally Singular Matrix Handling (Priority: P2)

A solver developer invokes the combined pipeline on a structurally singular matrix (MC64 returns `matched < n`). The system correctly places unmatched indices at the end of the ordering, condenses only the matched subgraph, and produces a valid permutation where the unmatched rows/columns are positioned last in the elimination order.

**Why this priority**: Structurally singular matrices occur in real applications (constrained optimization, saddle-point systems). The pipeline must handle them gracefully without panics or incorrect orderings.

**Independent Test**: Can be validated by constructing or loading a structurally singular matrix, running the pipeline, and verifying that unmatched indices appear at positions [matched..n) in the output ordering.

**Acceptance Scenarios**:

1. **Given** a structurally singular matrix where MC64 matches k < n entries, **When** the condensed pipeline runs, **Then** the first k positions in the ordering correspond to matched indices and the last n-k positions correspond to unmatched indices.
2. **Given** a structurally singular matrix, **When** the pipeline completes, **Then** the output is still a valid permutation of {0, ..., n-1} and the scaling vector is well-defined for all indices (using Duff-Pralet correction for unmatched entries).

---

### User Story 3 - Reduced METIS Input Size for Performance (Priority: P3)

A solver developer processes a large indefinite matrix (dimension > 100K) and benefits from condensation reducing the METIS input by up to ~50%. The system builds the condensed graph efficiently and passes it to METIS, resulting in faster ordering computation compared to running METIS on the full graph.

**Why this priority**: Performance improvement is a secondary benefit — the primary motivation is correctness (pair adjacency). However, for large matrices the ~2x reduction in METIS input size translates to measurable speedup, which matters for repeated analyses of matrices with the same sparsity pattern.

**Independent Test**: Can be validated by timing METIS on the condensed graph vs the full graph for large SuiteSparse matrices and comparing ordering computation time.

**Acceptance Scenarios**:

1. **Given** a large indefinite matrix where MC64 produces many 2-cycles, **When** the condensed graph is built, **Then** the condensed dimension is strictly less than the original dimension (n_condensed < n).
2. **Given** the condensed graph, **When** METIS ordering time is measured, **Then** METIS completes faster on the condensed graph than on the original graph for the majority (>= 60%) of test matrices.

---

### Edge Cases

- What happens when the matrix is diagonal (all singletons, no edges in condensed graph)? The pipeline should return the identity permutation or a valid trivial ordering.
- What happens when dimension is 0 or 1? The pipeline should handle these trivially without calling MC64 or METIS.
- What happens when the matching contains only singletons (no 2-cycles)? The condensed graph has the same dimension as the original — condensation adds negligible overhead and degrades to the standard METIS-only path.
- What happens when every pair is a 2-cycle (fully matched, no singletons except possibly one for odd n)? The condensed graph is approximately n/2 — maximum compression benefit.
- What happens when MC64 matching contains longer cycles (length > 2) in the raw matching? The cycle-splitting step must break these into 2-cycles and singletons (as SPRAL does), since only 2-cycles correspond to natural 2x2 pivot candidates.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST decompose an MC64 matching permutation into singletons (sigma(i)=i), 2-cycles (sigma(i)=j, sigma(j)=i), and unmatched indices (no match partner).
- **FR-002**: System MUST handle longer cycles (length > 2) in the raw MC64 matching by splitting them into 2-cycles and singletons, following SPRAL's `mo_split` approach — pair consecutive cycle members, leaving an odd-one-out as a singleton if the cycle has odd length.
- **FR-003**: System MUST build a condensed graph where each singleton maps to one super-node and each 2-cycle maps to one super-node, reducing the graph dimension. The condensed graph must include all edges from both members of a 2-cycle pair into a single column.
- **FR-004**: System MUST exclude unmatched indices from the condensed graph. Only matched super-nodes participate in the METIS ordering.
- **FR-005**: System MUST invoke METIS on the condensed graph and produce a fill-reducing ordering of the super-nodes.
- **FR-006**: System MUST expand the condensed ordering back to original variable indices, placing matched 2-cycle pairs in consecutive positions and appending unmatched indices at the end.
- **FR-007**: System MUST return a valid permutation (dimension n) compatible with `SymmetricOrdering::Custom` for input to `AptpSymbolic::analyze()`.
- **FR-008**: System MUST also return the MC64 scaling vector alongside the ordering, so downstream numeric factorization (Phase 5) can apply scaling.
- **FR-009**: System MUST handle trivial cases (n <= 1, diagonal matrices) without invoking MC64 or METIS.
- **FR-010**: System MUST provide bidirectional mappings (original-to-condensed, condensed-to-original) as intermediate results to enable validation and debugging.

### Key Entities

- **CycleDecomposition**: The result of splitting an MC64 matching into singletons, 2-cycles, and unmatched indices. Includes the old-to-new and new-to-old super-node mappings and the number of matched super-nodes.
- **CondensedGraph**: The compressed sparse graph with one node per singleton/2-cycle. Stores only the matched portion — unmatched nodes are excluded.
- **MatchOrderResult**: The final output combining the fill-reducing ordering, the MC64 scaling vector, the match count, and diagnostic information (condensed dimension, cycle statistics).

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Every 2-cycle pair in the MC64 matching occupies consecutive positions in the output ordering for 100% of test matrices (both hand-constructed and SuiteSparse collection).
- **SC-002**: Unmatched indices (for structurally singular matrices) appear in the final positions of the ordering, after all matched indices.
- **SC-003**: The condensed ordering produces fill predictions (nnz(L)) no more than 10% worse than running unconstrained METIS directly on the full graph, on the SuiteSparse CI-subset. One-sided tolerance: `condensed_nnz <= unconstrained_nnz * 1.10`. (Condensation constrains METIS by fusing paired nodes, so fill regression is the expected direction.)
- **SC-004**: The condensed graph dimension is strictly less than n for any matrix where the MC64 matching contains at least one 2-cycle.
- **SC-005**: Total pipeline time (MC64 + condense + METIS on condensed + expand) does not exceed 1.5x the time of MC64 + METIS on the full graph, ensuring condensation overhead is bounded.
- **SC-006**: All existing MC64 and METIS tests continue to pass — the new pipeline is additive, not replacing existing components.
- **SC-007**: The output ordering produces valid symbolic analysis results via `AptpSymbolic::analyze()` on all SuiteSparse CI-subset matrices.

## Assumptions

- The MC64 matching (`mc64_matching()`) and METIS ordering (`metis_ordering()`) functions from Phases 4.1 and 4.2 are correct and available as building blocks.
- MC64 matchings for symmetric matrices consist primarily of singletons and 2-cycles. Longer cycles are rare but must be handled (split into pairs).
- The condensed graph can reuse the existing CSC construction pattern from `metis_ordering()`'s `extract_adjacency()` function.
- METIS can handle the condensed graph (which may have multi-edges after merging 2-cycle columns) — METIS is documented to handle general graphs.
- The scaling vector is passed through unchanged; condensation affects only the ordering, not the scaling.
- The combined match-order pipeline is an explicit option alongside separate `mc64_matching()` and `metis_ordering()` calls. It is not the default for indefinite matrices — users explicitly choose which approach to use. It may become the recommended default in later phases once empirical benchmarks confirm its benefit.

## Scope Boundaries

**In scope:**
- Cycle decomposition of MC64 matchings
- Condensed graph construction
- METIS ordering on condensed graph
- Expansion back to full-size ordering
- Integration with existing `AptpSymbolic::analyze()` API
- Comprehensive validation and performance tests

**Out of scope:**
- Changes to MC64 matching algorithm (Phase 4.2)
- Changes to METIS ordering algorithm (Phase 4.1)
- Changes to symbolic analysis (Phase 3)
- Alternative ordering algorithms (AMD, SCOTCH) — the condensation approach works with any ordering backend, but only METIS is in scope
- Numeric factorization integration (Phase 5) — this phase produces the ordering; Phase 5 consumes it
- Parallel condensation — sequential implementation is sufficient given condensation is O(n+nnz) and dominated by METIS cost

## Clarifications

### Session 2026-02-13

- Q: Should SC-003 fill quality tolerance be one-sided (allow regression) or bidirectional? → A: One-sided — condensed_nnz <= unconstrained_nnz * 1.10. Condensation constrains METIS's freedom, so fill regression is the expected direction.
- Q: Should match_order be the recommended default for indefinite matrices or an explicit option? → A: Explicit option alongside separate MC64/METIS. Users choose which approach to use. May become default after Phase 5 benchmarks.
