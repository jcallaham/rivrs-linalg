# Feature Specification: APTP Symbolic Analysis

**Feature Branch**: `010-aptp-symbolic`
**Created**: 2026-02-10
**Status**: Draft
**Input**: Phase 3 of ssids-plan.md — Build `AptpSymbolic` as a thin composition over faer's symbolic Cholesky pipeline, adding APTP-specific delayed-pivot buffer estimation.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Perform Symbolic Analysis with Default Ordering (Priority: P1)

A solver developer provides a sparse symmetric matrix and receives a symbolic analysis result that captures the elimination tree structure, fill-in prediction, and permutation — all necessary inputs for the numeric factorization phase. The default ordering (AMD) is used when no custom ordering is specified.

**Why this priority**: This is the foundational capability. Without symbolic analysis, no numeric factorization or solve can proceed. Every downstream phase (Phases 5–8) depends on this result.

**Independent Test**: Can be fully tested by loading any test matrix (hand-constructed or SuiteSparse), running symbolic analysis with AMD ordering, and verifying that the result reports valid statistics (dimension, nonzero count, column structure). No numeric factorization is needed.

**Acceptance Scenarios**:

1. **Given** any symmetric sparse matrix from the test suite (hand-constructed or SuiteSparse CI subset), **When** symbolic analysis is run with default (AMD) ordering, **Then** the result contains a valid permutation, a positive predicted nonzero count, and statistics matching the matrix dimension.
2. **Given** the same matrix analyzed twice with the same ordering, **When** results are compared, **Then** they are identical (deterministic output).
3. **Given** a valid symmetric sparse matrix, **When** symbolic analysis completes, **Then** the predicted nonzero count is at least as large as the number of structural nonzeros in the lower triangle of the input matrix.

---

### User Story 2 — Accept Custom Ordering for Alternative Fill-Reducing Strategies (Priority: P2)

A solver developer provides a pre-computed fill-reducing permutation (e.g., from METIS or MC64 in a future phase) and the symbolic analysis uses that ordering instead of computing its own. This enables Phase 4 (alternative orderings) to plug in seamlessly.

**Why this priority**: Extensibility for future ordering strategies. Phase 4 (MC64) will produce custom orderings that must be accepted here. Without this, the symbolic analysis would be locked to AMD only.

**Independent Test**: Can be tested by constructing an identity permutation (or any valid permutation), passing it as a custom ordering, and verifying that the symbolic analysis completes successfully and the permutation is reflected in the result.

**Acceptance Scenarios**:

1. **Given** a symmetric sparse matrix and a valid custom permutation, **When** symbolic analysis is run with that custom ordering, **Then** the analysis completes and the stored permutation matches the one provided.
2. **Given** two different valid permutations for the same matrix, **When** symbolic analysis is run with each, **Then** the predicted nonzero counts may differ (since different orderings produce different fill-in), but both results are valid.

---

### User Story 3 — Access Supernodal Structure for Multifrontal Factorization (Priority: P3)

A developer building the multifrontal numeric factorization (Phase 6) accesses the supernodal decomposition from the symbolic analysis result — supernode column ranges, the assembly tree, and row structure per supernode. These are needed to organize the blocked dense operations within each front.

**Why this priority**: The multifrontal factorization (Phase 6) requires supernodal structure to organize its frontal matrices. Without supernodal accessors, the symbolic analysis result would need to be restructured later.

**Independent Test**: Can be tested by running symbolic analysis on a matrix known to have multiple supernodes (e.g., a banded matrix or a SuiteSparse matrix), then verifying that supernodal structure is accessible and structurally valid (supernode column ranges cover all columns, assembly tree has correct parent-child relationships).

**Acceptance Scenarios**:

1. **Given** a symbolic analysis result for a matrix with non-trivial structure, **When** supernodal accessors are queried, **Then** the supernode column ranges partition all columns of the factored matrix without gaps or overlaps.
2. **Given** a symbolic analysis result, **When** the assembly tree is queried, **Then** it forms a valid tree (each supernode except the root has exactly one parent).

---

### User Story 4 — Estimate APTP Pivot Buffer Requirements (Priority: P3)

The symbolic analysis provides per-column estimates of extra workspace needed for delayed pivots during the APTP numeric factorization. These buffer estimates allow the numeric phase to pre-allocate sufficient memory, avoiding costly reallocation during factorization.

**Why this priority**: Delayed pivots are an APTP-specific concern that faer's symbolic Cholesky does not model. Without buffer estimates, the numeric factorization would either over-allocate (wasting memory) or need dynamic reallocation (hurting performance).

**Independent Test**: Can be tested by running symbolic analysis and verifying that buffer estimates are non-negative, proportional to the column structure, and consistently reproducible.

**Acceptance Scenarios**:

1. **Given** a symbolic analysis result, **When** pivot buffer estimates are queried, **Then** each estimate is non-negative and the total buffer size is proportional to the predicted fill-in.
2. **Given** two analyses of the same matrix with the same ordering, **When** buffer estimates are compared, **Then** they are identical.

---

### Edge Cases

- What happens when the input matrix has dimension 1x1? The analysis should succeed, producing trivial structure (one column, one supernode, identity permutation or no permutation).
- What happens when the input matrix is diagonal (no off-diagonal entries)? The analysis should succeed with predicted nonzero count equal to the dimension and trivially ordered elimination tree.
- What happens when the input matrix has a disconnected sparsity graph (block-diagonal structure)? The analysis should succeed; the elimination tree should reflect the independent blocks.
- What happens when a custom permutation has incorrect dimension (does not match matrix size)? The system should return a descriptive error, not panic or produce silently wrong results.
- What happens when the input matrix is empty (0x0)? The system should either return a trivial result or a descriptive error.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST perform symbolic analysis on any valid sparse symmetric matrix, producing an analysis result that contains an elimination tree, predicted factor structure, and fill-reducing permutation.
- **FR-002**: System MUST accept a fill-reducing ordering as an input parameter, supporting at minimum (a) automatic AMD ordering and (b) a user-supplied custom permutation.
- **FR-003**: System MUST predict the number of nonzeros in the Cholesky factor L, consistent with the prediction made by the underlying symbolic factorization for the same ordering.
- **FR-004**: System MUST expose the fill-reducing permutation from the symbolic analysis, if one was computed or supplied.
- **FR-005**: System MUST estimate per-column pivot buffer requirements for APTP delayed pivots, based on heuristic analysis of the symbolic factor structure.
- **FR-006**: System MUST produce deterministic results — the same matrix and ordering always yield identical analysis output.
- **FR-007**: System MUST expose supernodal structure from the symbolic analysis: supernode column ranges, assembly (supernodal) tree, and row structure per supernode.
- **FR-008**: System MUST report symbolic analysis statistics: matrix dimension, predicted nonzero count, and average column count in the factor.
- **FR-009**: System MUST return a descriptive error (not panic) when the input matrix is structurally invalid (e.g., dimension mismatch with a custom permutation, non-square matrix).
- **FR-010**: System MUST compose with the existing APTP data structures from Phase 2 (specifically, the buffer estimates should be usable by the numeric factorization phase to size `MixedDiagonal` and related structures).

### Key Entities

- **AptpSymbolic**: The central analysis result. Composes a symbolic Cholesky factorization with APTP-specific metadata (pivot buffer estimates). Produced by the `analyze` operation and consumed by the numeric factorization phase.
- **SymbolicStatistics**: A summary of the symbolic analysis for diagnostics and debugging. Fields: `dimension` (matrix dimension), `predicted_nnz` (predicted nonzero count in L), `average_col_count` (mean nonzeros per column), `is_supernodal` (whether faer chose supernodal over simplicial), `n_supernodes` (number of supernodes if supernodal, `None` if simplicial), `total_pivot_buffer` (sum of all pivot buffer estimates).
- **Pivot Buffer Estimates**: Per-column workspace estimates for delayed pivots. A heuristic quantity specific to APTP that is not part of standard Cholesky symbolic analysis. Based on column counts from the symbolic structure.

### Assumptions

- The symbolic phase for indefinite LDL^T is structurally identical to SPD Cholesky — pivoting is purely a numeric-phase concern. The same elimination tree and sparsity structure apply. (Reference: Hogg et al. 2016 Section 2.2; ssids-plan.md Phase 3 design decisions.)
- AMD ordering is the appropriate default for the initial implementation. Alternative orderings (METIS, MC64-derived) will be provided in Phase 4 via the custom ordering input.
- The 10% heuristic for pivot buffer estimation (from ssids-plan.md) is a reasonable starting point. The actual buffer needs will be validated empirically during the numeric factorization phase (Phase 5–6) and the heuristic may be refined.
- faer's symbolic Cholesky computes both simplicial and supernodal structure internally. Both should be accessible through the analysis result for use by downstream phases.

### Algorithm References

The following academic references inform the symbolic analysis approach:

| Reference | File | Relevance |
|-----------|------|-----------|
| Liu (1990), "The role of elimination trees in sparse factorization" | (cited in plan; foundational theory) | Elimination tree structure, postorder traversal |
| Gilbert, Moler, Schreiber (1992), "Sparse Matrices in MATLAB" | `/workspace/rivrs-linalg/references/ssids/gilbert1992.md` | Section 3.3.4: elimination tree algorithms; Section 3.3.2: supernodes |
| Gilbert, Ng, Peyton (1994), "An Efficient Algorithm to Compute Row and Column Counts" | `/workspace/rivrs-linalg/references/ssids/gilbert1994.md` | Predicting factor sparsity structure (row/column counts) without numeric computation |
| Liu (1992), "The Multifrontal Method for Sparse Matrix Solution" | `/workspace/rivrs-linalg/references/ssids/liu1992.md` | Assembly tree construction, supernode clustering, formal definitions |
| Hogg, Ovtchinnikov, Scott (2016), "A Sparse Symmetric Indefinite Direct Solver for GPU Architectures" | `/workspace/rivrs-linalg/references/ssids/hogg2016.md` | Section 2.2: analyze phase for SSIDS; supernode detection; delayed pivot impact on structure |
| Duff, Reid (1984), "The Multifrontal Solution of Unsymmetric Sets of Linear Equations" | `/workspace/rivrs-linalg/references/ssids/duff1984.md` | Assembly tree creation, front matrix structure |
| George (1973), "Nested Dissection of a Regular Finite Element Mesh" | `/workspace/rivrs-linalg/references/ssids/george1973.md` | Ordering impact on elimination structure |

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Symbolic analysis completes successfully on 100% of the test matrix suite (15 hand-constructed + 10 SuiteSparse CI subset) with AMD ordering.
- **SC-002**: Symbolic analysis produces deterministic results — running the same analysis twice on any test matrix yields bit-identical output.
- **SC-003**: Custom ordering input is accepted and propagated correctly — the stored permutation matches the supplied one, and the analysis completes without error.
- **SC-004**: Predicted nonzero count is consistent with the underlying Cholesky prediction for the same ordering and matrix.
- **SC-005**: Pivot buffer estimates are non-negative for every column and reproducible across runs.
- **SC-006**: Symbolic analysis time is less than 5% of numeric factorization time on representative matrices (validated when numeric factorization exists; recorded as baseline now).
- **SC-007**: Supernodal structure is accessible and structurally valid — supernode column ranges partition all columns without gaps or overlaps.
- **SC-008**: Invalid inputs (dimension mismatch, invalid permutation) produce descriptive errors, not panics.
