# Feature Specification: Robustness — Testing & Hardening

**Feature Branch**: `026-robustness-hardening`
**Created**: 2026-02-25
**Status**: Draft
**Input**: User description: "Phase 9.2: Robustness — Testing & Hardening. SPRAL test parity audit, SPRAL-style torture testing, property-based testing (proptest), fuzzing/adversarial inputs."

## User Scenarios & Testing *(mandatory)*

### User Story 1 — SPRAL Test Parity Audit (Priority: P1)

A solver maintainer audits the existing test suite in two directions: (1) against SPRAL's test coverage to identify gaps where SPRAL tests something we don't, and (2) across the rivrs-sparse suite itself to assess each test's independent value. The audit compares SPRAL's kernel-level test scenarios (pivot handling, delayed columns, singular blocks, mixed 1×1/2×2 pivots) with the current rivrs-sparse test suite, producing a mapping document that shows which SPRAL scenarios are covered, which map to N/A (architecture differences), and which represent genuine gaps.

Separately, the audit evaluates each existing rivrs-sparse test on its own merits: does it test a meaningful behavior? Tests that cover scenarios SPRAL doesn't test are valuable and should be kept — the goal is SPRAL parity as a *floor*, not a ceiling. Only tests that are genuinely redundant (duplicating another test's coverage) or were scaffolding for intermediate TDD steps (testing internal details that are no longer independently meaningful) should be removed.

**Why this priority**: The audit is foundational — it identifies what gaps exist before adding new tests. Without it, new torture tests may duplicate existing coverage or miss critical scenarios.

**Independent Test**: Can be validated by reviewing the mapping document and confirming (a) every SPRAL test scenario has a corresponding entry (covered, N/A with rationale, or gap with new test added) and (b) every removed rivrs-sparse test has a documented rationale (redundant with specific other test, or TDD scaffolding no longer meaningful).

**Acceptance Scenarios**:

1. **Given** the SPRAL kernel test suites (`spral/tests/ssids/kernels/ldlt_app.cxx`, `ldlt_tpp.cxx`), **When** the audit is performed, **Then** every distinct test scenario in SPRAL maps to either an existing rivrs-sparse test, a documented N/A (with rationale), or a new test filling the gap.
2. **Given** the current rivrs-sparse test suite (~524 tests), **When** the audit evaluates each test's independent value, **Then** tests that are genuinely redundant or TDD scaffolding are removed, tests that cover meaningful scenarios (even ones SPRAL doesn't test) are kept, and the remaining suite still provides equivalent or better coverage (all SuiteSparse matrices pass, all hand-constructed validations pass).
3. **Given** the completed audit document, **When** a reviewer reads it, **Then** they can verify both (a) SPRAL coverage completeness and (b) that every removal has a documented rationale beyond "not in SPRAL."

---

### User Story 2 — SPRAL-Style Torture Testing (Priority: P1)

A solver developer adds probabilistic stress tests that exercise the dense APTP kernel in isolation (matching SPRAL's approach of testing `aptp_factor_in_place` / `tpp_factor_as_primary` with dense matrices, not through the full multifrontal pipeline). The tests create adversarial numerical conditions: forced pivot delays, rank deficiency, and singular diagonal blocks. They run hundreds of random instances per configuration, catching subtle numerical issues that deterministic tests miss. End-to-end coverage through `SparseLDLT` is already provided by the SuiteSparse suite.

**Why this priority**: SPRAL's torture tests (ldlt_app.cxx lines 423-451, ldlt_tpp.cxx lines 343-369) caught issues that deterministic tests missed. This is the highest-value new test category.

**Independent Test**: Can be validated by running the torture test suite and confirming all instances complete without panics or assertion failures, with backward error within tolerance.

**Acceptance Scenarios**:

1. **Given** a random symmetric indefinite matrix, **When** `cause_delays()` multiplies n/8 random rows by 1000, **Then** the solver handles the forced pivot delays and produces backward error < 5e-11 (or cleanly reports numerical singularity).
2. **Given** a random symmetric matrix, **When** `make_singular()` creates rank deficiency by copying one column to another (with scaling), **Then** the solver produces a valid factorization with delayed columns and inertia reflecting the rank deficiency (n_zero > 0). Solve may produce large backward error on singular systems — that is acceptable; the invariant is correct factorization and inertia, not solve accuracy.
3. **Given** a random symmetric matrix, **When** `make_dblk_singular()` makes a specific diagonal block singular, **Then** the solver handles the singular block via 2×2 pivoting or delay, without panics.
4. **Given** 500 random instances per configuration (matrix sizes 32–256, ~70% with delays, ~20% singular, ~10% singular diagonal blocks), **When** the torture test suite runs, **Then** zero instances panic and backward error is within tolerance for all non-singular systems.

---

### User Story 3 — Property-Based Testing (Priority: P2)

A solver developer adds property-based tests using a property-based testing framework that generates random symmetric matrices of varying size, density, and definiteness, then verifies structural invariants of the factorization output.

**Why this priority**: Property-based testing complements torture tests by exploring a broader parameter space (size, density, definiteness) with automatically-shrunk counterexamples on failure.

**Independent Test**: Can be validated by running the property-based test suite and confirming all properties hold across thousands of generated instances.

**Acceptance Scenarios**:

1. **Given** a randomly generated symmetric positive-definite matrix (size 5–500, varying density), **When** the solver runs analyze/factor/solve, **Then** backward error < 5e-11.
2. **Given** a randomly generated symmetric indefinite matrix (size 5–500), **When** the solver runs, **Then** either backward error < 5e-11 or a clean error is returned (no panics).
3. **Given** a randomly generated symmetric matrix, **When** the solver produces inertia counts, **Then** `n_positive + n_negative + n_zero == matrix_dimension`.
4. **Given** a randomly generated symmetric matrix, **When** the solver produces a permutation, **Then** the permutation is valid (each index appears exactly once).

---

### User Story 4 — Adversarial & Edge-Case Input Testing (Priority: P2)

A solver developer adds tests for malformed, extreme, and degenerate inputs to ensure the solver never panics and always returns clean error values.

**Why this priority**: Production robustness requires graceful handling of all inputs, not just well-formed ones. This prevents downstream users from encountering panics.

**Independent Test**: Can be validated by running the adversarial test suite and confirming every case either succeeds or returns an appropriate error variant.

**Acceptance Scenarios**:

1. **Given** a 0×0 empty matrix, **When** the solver is invoked, **Then** it returns a clean error (not a panic).
2. **Given** a 1×1 matrix (including zero diagonal), **When** the solver is invoked, **Then** it produces a correct result or clean error.
3. **Given** a diagonal matrix, **When** the solver is invoked, **Then** the factorization is trivial and correct.
4. **Given** a matrix with extreme values (near `f64::MAX`, near `f64::MIN_POSITIVE`, exact zeros), **When** the solver is invoked, **Then** it either produces a valid result or returns a clean error.
5. **Given** a structurally non-symmetric sparsity pattern, **When** passed to the solver, **Then** the solver returns an appropriate error (not a panic or silent corruption).

---

### Edge Cases

- Matrix where all off-diagonal entries are zero (pure diagonal)
- Matrix with a single dense row/column (arrowhead pattern)
- Matrix where every pivot is delayed (worst-case for APTP)
- Matrix requiring all 2×2 pivots (no stable 1×1 pivots)
- Matrix with exact numerical cancellation during elimination
- Matrix dimension at power-of-2 boundaries (32, 64, 128, 256, 512) that interact with BLAS blocking
- Near-singular matrix where condition number approaches `1/eps`
- Matrix where MC64 matching fails to find a perfect matching
- Matrix with disconnected components in the sparsity graph

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The project MUST include a documented audit that (a) maps each distinct SPRAL kernel test scenario to a corresponding rivrs-sparse test, N/A rationale, or newly-added test, and (b) evaluates each existing rivrs-sparse test on its own merits — retaining tests that cover meaningful scenarios regardless of whether SPRAL tests them.
- **FR-002**: The test suite MUST include probabilistic perturbation helpers (`cause_delays`, `make_singular`, `make_dblk_singular`) that transform well-formed matrices into adversarial test cases following SPRAL's conceptual patterns.
- **FR-003**: The torture test suite MUST run at least 500 random instances per configuration (matching SPRAL's convention) across multiple matrix sizes.
- **FR-004**: The torture test distribution MUST include approximately 70% delay-inducing, 20% singular, and 10% singular-diagonal-block instances (matching SPRAL's probabilistic generation pattern).
- **FR-005**: The project MUST add property-based testing and include property tests for the solver that verify backward error, inertia consistency, and permutation validity.
- **FR-006**: The test suite MUST include adversarial input tests covering: empty matrices (0×0), trivial matrices (1×1), diagonal matrices, extreme floating-point values, and structurally invalid inputs.
- **FR-007**: The solver MUST NOT panic on any input covered by FR-006 — all failure modes MUST return errors with appropriate error variants.
- **FR-008**: Tests identified during the audit as genuinely redundant (duplicating another test's coverage) or as TDD scaffolding (testing intermediate internal details no longer independently meaningful) MUST be removed or consolidated. Tests covering meaningful scenarios not tested by SPRAL MUST be retained. The pruned suite MUST maintain equivalent or better coverage (all SuiteSparse matrices pass, all hand-constructed validations pass).
- **FR-009**: All torture test and property-based test modules MUST be gated behind the `test-util` feature flag, consistent with the existing test infrastructure.
- **FR-010**: The torture test suite MUST be runnable as `#[ignore]` tests (consistent with the SuiteSparse convention) since they are long-running.

### Key Entities

- **Perturbation Helper**: A function that transforms a well-formed symmetric matrix into an adversarial variant (delays, singularity, singular blocks). Belongs alongside existing random matrix generators.
- **Audit Mapping Document**: A document recording both the SPRAL-to-rivrs test coverage mapping and the per-test value assessment of the rivrs-sparse suite, living in the project's documentation directory as a reference artifact.
- **Torture Test Configuration**: Parameters controlling the probabilistic test generation — matrix sizes, perturbation probabilities, instance count — with SPRAL-matching defaults.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Audit covers 100% of distinct SPRAL test scenarios from `ldlt_app.cxx` and `ldlt_tpp.cxx` (each mapped to covered/N-A/new-test), and 100% of existing rivrs-sparse tests are evaluated (each marked as retained-with-rationale or removed-with-rationale).
- **SC-002**: Torture tests run 500+ random instances per configuration with zero panics across all configurations.
- **SC-003**: Property-based tests exercise at least 1000 generated matrices (across size range 5–500) with zero property violations.
- **SC-004**: Zero panics on any adversarial input (empty, trivial, extreme, malformed) — all return clean errors.
- **SC-005**: Existing SuiteSparse suite (65/65 pass) and hand-constructed validations remain passing after test suite pruning.
- **SC-006**: Test suite is leaner after audit — genuinely redundant and TDD-scaffolding tests removed, meaningful non-SPRAL tests retained, no coverage loss.

## Clarifications

### Session 2026-02-25

- Q: Torture tests target which API level — end-to-end (`SparseLDLT`), dense kernel, or both? → A: Dense kernel only (`aptp_factor_in_place` / `tpp_factor_as_primary`), matching SPRAL's `ldlt_app.cxx` / `ldlt_tpp.cxx` approach. End-to-end stress coverage is already provided by the SuiteSparse suite.
- Q: For rank-deficient matrices in torture tests, what is the expected outcome? → A: Factorization succeeds with delayed columns; inertia reflects the rank deficiency (n_zero > 0). Solve may produce large backward error on singular systems — that is acceptable.

## Assumptions

- SPRAL source code (`spral/tests/ssids/kernels/ldlt_app.cxx`, `ldlt_tpp.cxx`) is available under BSD-3 and can be freely consulted for test scenario analysis.
- Torture tests are long-running (minutes) and should be `#[ignore]` tests, run explicitly rather than in the default test suite.
- The existing `test-util` feature flag is the appropriate gate for new test infrastructure code.
- The perturbation helpers follow SPRAL's approach conceptually but are implemented from scratch in idiomatic Rust (clean room — we consult SPRAL's test design, not translate code).
- Matrix sizes for torture testing (32–256) are chosen to exercise the two-level APTP architecture (inner block size 32, outer block size 128/256) without excessive runtime.
- The existing random matrix generators (`generate_random_symmetric`, `generate_arrow`, `generate_tridiagonal` in `src/testing/generators.rs`) provide a foundation that the perturbation helpers build on.

## Algorithm References

The following academic papers and reference implementations inform this phase:

| Reference | Location | Relevance |
|-----------|----------|-----------|
| Duff, Hogg & Lopez 2020 | `/workspace/rivrs-linalg/references/ssids/duff2020.md` | Primary APTP algorithm reference; Section 5 discusses numerical experiments and test methodology |
| Hogg & Scott 2016 | `/workspace/rivrs-linalg/references/ssids/hogg2016.md` | SSIDS architecture; test matrix selection methodology |
| Duff & Pralet 2005 | `/workspace/rivrs-linalg/references/ssids/duff2005.md` | MC64 scaling effects on pivot stability; informs torture test design |
| Davis 2016 | `/workspace/rivrs-linalg/references/ssids/davis2016.md` | SuiteSparse collection methodology; test matrix classification |
| Schenk & Gärtner 2006 | `/workspace/rivrs-linalg/references/ssids/schenk2006.md` | Pivot strategy comparison; informs property-based test design for pivot behavior |
| Liu 1992 | `/workspace/rivrs-linalg/references/ssids/liu1992.md` | Multifrontal method theory; assembly correctness invariants |
| SPRAL ldlt_app.cxx tests | `/opt/references/spral/tests/ssids/kernels/ldlt_app.cxx` | Primary torture test reference (lines 78-134 perturbation helpers, lines 423-451 probabilistic generation, line 593 instance count) |
| SPRAL ldlt_tpp.cxx tests | `/opt/references/spral/tests/ssids/kernels/ldlt_tpp.cxx` | TPP-specific torture test patterns (lines 343-369) |
