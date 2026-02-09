# Feature Specification: APTP Data Structures

**Feature Branch**: `009-aptp-data-structures`
**Created**: 2026-02-08
**Status**: Draft
**Input**: User description: "Implement Phase 2 in ssids-plan.md — APTP data structures for indefinite LDL^T factorization"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Store Mixed 1x1/2x2 Block Diagonal (Priority: P1)

A solver developer building the numeric factorization kernel (Phase 5) needs a data structure to represent the D factor in P^T A P = L D L^T, where D contains a mix of 1x1 scalar pivots and 2x2 symmetric pivot blocks. The developer must be able to construct D incrementally (column by column, as pivots are decided) and then solve systems D x = b using the stored blocks.

**Why this priority**: The mixed diagonal D is the central data structure that distinguishes indefinite LDL^T from positive-definite Cholesky. Every downstream phase (numeric factorization, triangular solve) depends on it. Without correct D storage and solve, no solver output can be validated.

**Independent Test**: Can be fully tested by constructing known mixed diagonals (hand-computed examples with 1x1 and 2x2 blocks), solving D x = b, and verifying the result against analytical solutions. Delivers a self-contained, testable numeric kernel.

**Acceptance Scenarios**:

1. **Given** an empty MixedDiagonal of dimension n, **When** the developer sets a sequence of 1x1 and 2x2 pivots covering all n columns, **Then** the structure reports correct pivot types for each column and stores the correct values.
2. **Given** a fully populated MixedDiagonal with mixed 1x1 and 2x2 blocks, **When** `solve_in_place` is called with a known right-hand side b, **Then** the solution x satisfies ||D x - b|| / ||b|| < 10^-14.
3. **Given** a MixedDiagonal with some columns marked as Delayed, **When** the developer queries the number of delayed pivots, **Then** it returns the correct count.

---

### User Story 2 - Track Pivot Decisions During Factorization (Priority: P1)

A solver developer needs to record and query pivot decisions made during APTP factorization. Each column is classified as a 1x1 pivot, the first or second column of a 2x2 Bunch-Kaufman pivot, or delayed (failed stability check). This classification drives the factorization kernel's control flow and is needed for correct triangular solve and inertia computation.

**Why this priority**: Pivot tracking is co-equal with MixedDiagonal — they are built together during factorization. The PivotType classification is the control signal that determines how each column is processed.

**Independent Test**: Can be tested by constructing PivotType sequences for known factorization patterns and verifying that inertia (positive/negative/zero eigenvalue counts) can be derived from the pivot classifications.

**Acceptance Scenarios**:

1. **Given** a sequence of pivot decisions, **When** the developer queries the pivot type for each column, **Then** the correct classification (1x1, 2x2 with partner index, or Delayed) is returned.
2. **Given** a MixedDiagonal with known 1x1 values and 2x2 blocks, **When** inertia is computed from the pivot information, **Then** the result matches the analytically known eigenvalue sign counts.
3. **Given** a 2x2 pivot at columns (i, i+1), **When** the developer queries both columns, **Then** column i reports TwoByTwo with partner i+1 and column i+1 reports TwoByTwo with partner i.

---

### User Story 3 - Construct faer Permutations from Ordering Output (Priority: P2)

A solver developer implementing ordering algorithms (Phase 4) needs to construct faer's `Perm<usize>` from the forward permutation array produced by ordering algorithms. Most ordering algorithms (AMD, METIS) produce only the forward mapping; faer's `Perm::new_checked` requires both forward and inverse arrays.

**Why this priority**: This is a standalone utility needed before Phase 4 (ordering) but simpler than the diagonal structures. It is a bridge between external ordering output and faer's permutation infrastructure.

**Independent Test**: Can be tested by constructing permutations from forward arrays, verifying round-trip properties (apply then apply inverse = identity), and confirming compatibility with faer's permutation operations.

**Acceptance Scenarios**:

1. **Given** a valid forward permutation array of length n, **When** `perm_from_forward` is called, **Then** it returns a `Perm<usize>` where applying the permutation and then its inverse yields the identity.
2. **Given** two forward permutation arrays, **When** both are converted to `Perm<usize>` and composed, **Then** the result matches the mathematically composed permutation.
3. **Given** an invalid forward array (duplicates or out-of-bounds), **When** `perm_from_forward` is called, **Then** it returns an appropriate error.

---

### Edge Cases

- What happens when the matrix dimension is 0? MixedDiagonal should handle empty matrices gracefully — construction succeeds, solve is a no-op, delayed count is 0.
- What happens when all pivots are 2x2? The dimension n must be even for this to be valid; MixedDiagonal should reject an attempt to set a 2x2 block starting at column n-1 when n is odd.
- What happens when a 2x2 block's determinant is zero or near-zero? `solve_in_place` debug-asserts non-singularity; in release mode, division-by-zero produces Inf/NaN naturally. The public solve API (Phase 7) will add `Result`-returning checks at a higher layer.
- What happens when a 1x1 pivot value is exactly zero? Same as above — debug-assert catches this in development; the factorization kernel (Phase 5) is responsible for preventing singular pivots from reaching MixedDiagonal.
- What happens when `solve_in_place` is called with delayed columns still present? Debug-assert that no delayed columns exist — all delays must be resolved by the factorization kernel before solve. Consistent with the singular-pivot policy.
- What happens when the forward permutation is the identity? `perm_from_forward` should produce an identity `Perm`.
- What happens when `perm_from_forward` receives an empty array? It should produce a valid dimension-0 `Perm`.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide a PivotType classification that distinguishes three states: 1x1 pivot, 2x2 pivot (with partner column identification), and delayed column.
- **FR-002**: System MUST provide a Block2x2 structure that stores three independent values (a, b, c) of a 2x2 symmetric block [a, b; b, c] along with the column index of the first column.
- **FR-003**: System MUST provide a MixedDiagonal structure that stores a sequence of 1x1 and 2x2 diagonal blocks, indexed by column, representing the D factor in LDL^T.
- **FR-004**: MixedDiagonal MUST support incremental construction — setting 1x1 pivots and 2x2 pivot blocks one at a time in column order.
- **FR-005**: MixedDiagonal MUST provide a `solve_in_place` operation that solves D x = b where D is the stored mixed block diagonal, modifying x in place. The operation debug-asserts two preconditions: (1) all pivots are non-singular, and (2) no delayed columns remain. Both invariants are the responsibility of upstream phases; the public solve API (Phase 7) adds `Result`-returning validation at a higher layer.
- **FR-006**: MixedDiagonal MUST report the count of delayed pivots.
- **FR-007**: MixedDiagonal MUST provide query access to the pivot type of any column.
- **FR-007a**: MixedDiagonal MUST provide a `compute_inertia` method that returns an `Inertia` (positive/negative/zero eigenvalue counts) derived from the stored 1x1 values and 2x2 block eigenvalue signs.
- **FR-008**: System MUST provide a `perm_from_forward` function that constructs a faer `Perm<usize>` from a forward permutation array by computing the inverse automatically.
- **FR-009**: `perm_from_forward` MUST validate its input (no duplicates, no out-of-bounds values) and return an error for invalid inputs.
- **FR-010**: MixedDiagonal's `solve_in_place` MUST correctly solve 2x2 symmetric systems using the analytical solution for 2x2 linear systems.
- **FR-011**: All new types MUST follow the project's transparent composition principle — using faer types (`Perm`, `SparseColMat`) at the boundary rather than introducing custom wrappers for concepts faer already models.
- **FR-012**: The existing `DBlock` enum in `io/reference.rs` MUST NOT be replaced or modified — it serves a different purpose (deserializing reference factorization JSON files). The new types are for the live factorization data path.
- **FR-013**: The `Inertia` struct MUST be relocated from `io/reference.rs` to the new APTP module, since it is a domain concept (eigenvalue sign classification) now used by both reference data loading and live solver operations. A `pub use` re-export MUST be added in `io/reference.rs` to preserve backward compatibility for existing imports.

### Key Entities

- **PivotType**: Classification of a single column's pivot decision during APTP factorization. Three variants: 1x1 scalar pivot, 2x2 Bunch-Kaufman pivot (identifying the partner column), and Delayed (column failed stability check and must be deferred to an ancestor node).
- **Block2x2**: Storage for a single 2x2 symmetric diagonal block [a, b; b, c]. Identified by the index of its first column. Stores only three values (exploiting symmetry: b = off-diagonal element shared by both positions).
- **MixedDiagonal**: The complete D factor in P^T A P = L D L^T, containing a mix of 1x1 and 2x2 blocks. Provides construction, query, solve, and inertia computation operations. Tracks delayed columns for downstream use by the numeric factorization phase.
- **Inertia**: Eigenvalue sign classification (positive, negative, zero counts). Relocated from `io/reference.rs` to the APTP module as the canonical domain type, re-exported from the original location for backward compatibility.
- **perm_from_forward**: A standalone utility function bridging ordering algorithm output (forward permutation arrays) to faer's `Perm<usize>` type by computing the inverse mapping automatically.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: MixedDiagonal correctly stores and retrieves any valid pattern of mixed 1x1 and 2x2 pivots for dimensions up to at least 10,000.
- **SC-002**: `solve_in_place` produces solutions with relative error ||D x - b|| / ||b|| < 10^-14 for all hand-constructed test cases with known analytical solutions.
- **SC-003**: PivotType correctly classifies all three pivot states, and inertia derived from pivot classifications matches known values for the 15 hand-constructed reference factorizations.
- **SC-004**: `perm_from_forward` produces valid faer `Perm<usize>` objects that satisfy the round-trip identity property (apply then invert = original vector) for permutations up to dimension 10,000.
- **SC-005**: All new public types and functions pass linting (`cargo clippy`), formatting (`cargo fmt --check`), and documentation (`cargo doc`) with zero warnings. All public items include rustdoc citing the academic reference(s) they implement (per Constitution principles II and IV).
- **SC-006**: The new module integrates into the existing codebase with minimal changes — the only modification to existing modules is adding a `pub use` re-export of `Inertia` in `io/reference.rs` after relocating the struct to the new APTP module. No behavioral changes to existing code.

## Assumptions

- The numeric factorization (Phase 5) will construct MixedDiagonal incrementally in column order during the elimination process. The API should support this construction pattern.
- 2x2 pivot partners are always adjacent columns (column i and i+1), consistent with the Bunch-Kaufman strategy described in the APTP literature.
- The `solve_in_place` operation on MixedDiagonal is not performance-critical at this phase — it will be called once per solve, not in an inner loop. Optimization can be deferred.
- `perm_from_forward` will initially be used in tests and Phase 4 ordering integration. It does not need to handle permutations beyond the scale of the test matrix collection (~1.6M dimension).
- Delayed columns in MixedDiagonal are tracked but not "resolved" in this phase — resolution happens during the numeric factorization (Phase 5) when delayed columns are passed to ancestor nodes in the elimination tree.
- The existing `DBlock` enum serves a separate purpose (JSON deserialization of reference factorizations) and will coexist with the new PivotType/MixedDiagonal types. No unification is needed.

## Clarifications

### Session 2026-02-08

- Q: What should `solve_in_place` do when it encounters a singular pivot (zero 1x1 or zero-determinant 2x2)? → A: Debug-assert non-singularity. The factorization kernel (Phase 5) is responsible for preventing singular pivots. The public solve API (Phase 7) adds `Result`-returning error handling at a higher layer. This keeps the internal solve path simple while catching invariant violations during development.
- Q: What should `solve_in_place` do when delayed columns are present? → A: Debug-assert no delayed columns exist. All delays must be resolved by the factorization kernel before solve is called. Same invariant-enforcement strategy as singular pivots.
- Q: Should MixedDiagonal provide inertia computation, and should Inertia be relocated? → A: Yes to both. MixedDiagonal provides `compute_inertia() -> Inertia` since it has all required information (1x1 signs, 2x2 eigenvalue signs). `Inertia` is relocated from `io/reference.rs` to the new APTP module as the canonical domain type, with a `pub use` re-export for backward compatibility.
