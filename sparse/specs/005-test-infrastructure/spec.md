# Feature Specification: Core Test Infrastructure

**Feature Branch**: `005-test-infrastructure`
**Created**: 2026-02-06
**Status**: Draft
**Input**: User description: "implement phase 1.1 in docs/ssids-plan.md"

## Clarifications

### Session 2026-02-06

- Q: How should test infrastructure be made accessible to both unit tests and integration tests? → A: Cargo feature flag (`test-util`). Test infra lives in `src/`, gated behind `#[cfg(feature = "test-util")]`, enabled via `[dev-dependencies] rivrs-sparse = { path = ".", features = ["test-util"] }`.
- Q: Should the `NumericalValidator` replace or wrap the existing standalone `validate` functions? → A: Keep both. Existing `validate::reconstruction_error`, `validate::backward_error`, and `validate::check_inertia` remain public. `NumericalValidator` wraps them with configurable tolerances and structured `TestResult` reporting. No deprecation.
- Q: Should `test_roundtrip` aggregate individual phase results or run an independent end-to-end check? → A: Independent end-to-end. `test_roundtrip` runs the full analyze→factor→solve pipeline as one operation, then validates the final result (reconstruction error, backward error, inertia). It does not call `test_analyze`/`test_factor`/`test_solve` internally.

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Validate Solver Component Correctness (Priority: P1)

A developer implementing a new solver component (e.g., symbolic analysis, numeric factorization, triangular solve) uses the test harness to verify their implementation produces mathematically correct results. They write tests using the `SolverTest` trait, which provides standardized test methods for each solver phase (analyze, factor, solve, roundtrip). The harness automatically checks reconstruction error, backward error, and inertia against configurable tolerances.

**Why this priority**: This is the foundational value proposition — every subsequent development phase (2-11) depends on a reliable, standardized way to validate correctness. Without this, each component reinvents validation logic.

**Independent Test**: Can be fully tested by implementing a mock/trivial solver (identity permutation, no-op factor, direct solve) and running it through the harness. Delivers value as soon as Phase 2 begins.

**Acceptance Scenarios**:

1. **Given** a developer has implemented a solver component, **When** they create a `SolverTestCase` from any registry matrix and run `test_roundtrip`, **Then** the harness returns a structured `TestResult` with reconstruction error, backward error, and inertia correctness — all compared against configurable tolerances.
2. **Given** a solver produces incorrect results, **When** the harness runs validation, **Then** the failure message includes the specific metric that failed, the computed value, the threshold, and the matrix name — enabling rapid diagnosis.
3. **Given** a developer needs to test only the analysis phase, **When** they call `test_analyze` on a test case, **Then** only analysis-related checks run (valid permutation, correct elimination tree shape) without requiring a factorization implementation.

---

### User Story 2 - Numerical Validation with Configurable Tolerances (Priority: P1)

A developer configures the numerical validator with problem-specific tolerances. For double-precision well-conditioned problems, reconstruction error must be below 10^-12 and backward error below 10^-10 (per constitution). For ill-conditioned matrices, the developer relaxes tolerances appropriately. The validator provides clear pass/fail results with diagnostic detail. This consolidates and supersedes the existing `validate.rs` functions by wrapping them in a configurable, structured validation engine.

**Why this priority**: Configurable tolerances are essential because different matrix classes require different accuracy expectations. Hard-coding tolerances leads to either false passes (too loose) or false failures (too tight on ill-conditioned matrices). Consolidating the existing validate module into the new harness eliminates scattered threshold logic.

**Independent Test**: Can be tested immediately using existing hand-constructed matrices. The validator wraps existing reconstruction and backward error computations with configurable thresholds and structured result reporting.

**Acceptance Scenarios**:

1. **Given** default tolerance settings, **When** validating a hand-constructed matrix with its known factorization, **Then** reconstruction error check passes with error below 10^-12.
2. **Given** a deliberately perturbed factorization, **When** the validator checks reconstruction error, **Then** the result reports failure with the exact computed error and the threshold it exceeded.
3. **Given** custom tolerances set to 10^-6, **When** checking a moderately accurate result, **Then** the check passes even though it would fail with default 10^-12 tolerance.

---

### User Story 3 - Generate Random Test Matrices (Priority: P2)

A developer generates random sparse symmetric matrices with specified properties (size, density, definiteness) for property-based testing and stress testing. Generators produce matrices with known structural properties, enabling tests that verify solver invariants (symmetry preservation, valid permutations, correct inertia sign counts) across many random inputs.

**Why this priority**: Random matrix generation enables property-based testing that catches edge cases not covered by the fixed 82-matrix test collection. This becomes critical in Phases 2-8 where subtle bugs may only manifest on specific sparsity patterns.

**Independent Test**: Can be tested by generating random matrices and verifying they have the requested structural properties (symmetric, correct dimension, correct definiteness, nonzero count within bounds).

**Acceptance Scenarios**:

1. **Given** a request for a 100x100 sparse symmetric positive definite matrix with ~500 nonzeros, **When** the generator runs, **Then** the output matrix is symmetric, has the requested dimension, is positive definite, and has nonzero count within 20% of the target.
2. **Given** a request for an indefinite matrix, **When** the generator runs, **Then** the output matrix has both positive and negative eigenvalues (verified via inertia of a dense LDL^T factorization).
3. **Given** a request for a matrix with specific sparsity pattern (arrow, tridiagonal, banded), **When** the generator runs, **Then** the nonzero pattern matches the requested structure.

---

### User Story 4 - Unified Test Case Management (Priority: P2)

A developer loads test cases from the existing matrix registry into the standardized `SolverTestCase` format, which bundles the matrix, its properties, and optional reference results. This supersedes the current pattern of calling `registry::load_test_matrix()` directly in each test file and manually checking properties. Test cases can be iterated, filtered by category (hand-constructed, easy-indefinite, hard-indefinite, positive-definite), and used with any `SolverTest` implementation.

**Why this priority**: Standardized test case management reduces boilerplate in integration tests and ensures consistent validation across all solver components. It consolidates the existing registry loading patterns into a uniform interface.

**Independent Test**: Can be tested by loading test cases from the registry and verifying they produce valid `SolverTestCase` instances with correct properties and reference data.

**Acceptance Scenarios**:

1. **Given** the existing matrix registry, **When** a developer loads all hand-constructed test cases, **Then** 15 `SolverTestCase` instances are returned, each with matrix data, properties, and reference factorization.
2. **Given** a filter for "hard-indefinite" category, **When** loading test cases, **Then** only matrices classified as hard-indefinite are returned.
3. **Given** any loaded test case, **When** accessing its properties, **Then** matrix size, nonzero count, definiteness, and difficulty are available without additional computation.

---

### Edge Cases

- What happens when a test matrix file is missing (gitignored SuiteSparse matrix)? The test case loader returns `None` for that matrix and tests skip gracefully with a clear skip message.
- How does the validator handle a zero matrix (all entries zero)? Returns 0.0 reconstruction error (degenerate but correct).
- What happens when generating a random matrix with contradictory properties (e.g., 1x1 indefinite)? Returns an error indicating the properties are infeasible.
- How does the harness handle NaN/Inf in solver output? Reports a clear validation failure identifying the specific entry.
- What happens when the random matrix generator targets an unreachable density (e.g., 1000 nnz for a 5x5 matrix)? Clamps to maximum possible nonzeros and documents the clamping.
- How does the harness handle a solver that panics? Tests use standard Rust panic handling — panics surface as test failures with backtraces.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST provide a `SolverTestCase` type that bundles a sparse matrix, its structural/numerical properties, and optional reference results into a single test entity
- **FR-002**: System MUST provide a `SolverTest` trait with methods `test_analyze`, `test_factor`, `test_solve`, and `test_roundtrip` that return structured `TestResult` values. `test_roundtrip` runs the full analyze→factor→solve pipeline as one independent end-to-end operation and validates the final result; it does not call the individual phase methods internally
- **FR-003**: System MUST provide a `NumericalValidator` with configurable tolerances for reconstruction error, backward error, and inertia comparison — wrapping the existing standalone functions in `validate.rs` (which remain public and unchanged)
- **FR-004**: The `NumericalValidator` MUST support default tolerances matching the constitution standards: reconstruction error < 10^-12, backward error < 10^-10
- **FR-005**: Every validation failure MUST include diagnostic context: metric name, computed value, threshold, matrix name, and matrix dimensions
- **FR-006**: System MUST provide a function to load all registry matrices matching a filter (source, category, difficulty) as `SolverTestCase` instances, superseding the pattern of calling `load_test_matrix` individually in each test file
- **FR-007**: System MUST provide random matrix generators for sparse symmetric matrices with configurable size, approximate nonzero count, and definiteness (positive definite, indefinite)
- **FR-008**: Random generators MUST produce matrices with verifiable structural properties (symmetry, correct dimension, sparsity within tolerance of target)
- **FR-009**: System MUST provide pattern generators for common sparsity structures: arrow, tridiagonal, and banded matrices
- **FR-010**: System MUST integrate with Rust's standard test framework (`#[test]`, `assert!`) and Criterion benchmarking framework
- **FR-011**: Loading any individual test matrix from the registry into a `SolverTestCase` MUST complete in under 100ms
- **FR-012**: Existing integration tests (`hand_constructed.rs`, `suitesparse_ci.rs`) MUST be refactored to use the new test harness types, demonstrating the harness works end-to-end and eliminating duplicated validation logic

### Key Entities

- **SolverTestCase**: A complete test scenario bundling a sparse matrix, its structural/numerical properties (size, nnz, definiteness, difficulty, structure), and optional reference results (factorization, inertia, expected solver statistics)
- **TestResult**: Structured outcome of a solver test, including pass/fail status, numeric metrics (reconstruction error, backward error), inertia comparison, and diagnostic messages
- **NumericalValidator**: Configurable validation engine that checks factorization quality against tolerance thresholds, wrapping reconstruction error, backward error, and inertia checking into a unified interface
- **SolverTest trait**: Interface that any solver implementation must satisfy to be tested by the harness — methods for each solver phase (analyze, factor, solve, full roundtrip)

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: All 15 hand-constructed matrices load into `SolverTestCase` format and validate with reconstruction error below 10^-12 using the `NumericalValidator`
- **SC-002**: Validation failures produce diagnostic messages that include the matrix name, the failing metric, the computed value, and the threshold — enabling a developer to diagnose the issue without additional debugging
- **SC-003**: Random matrix generators produce valid symmetric matrices for sizes from 5 to 10,000 with configurable density, completing generation in under 1 second for matrices up to size 1,000
- **SC-004**: The `SolverTest` trait can be implemented for a trivial pass-through solver (identity permutation, no factorization) and all four test methods (`test_analyze`, `test_factor`, `test_solve`, `test_roundtrip`) execute without errors on hand-constructed test cases
- **SC-005**: Test case filtering by category returns the correct subset (e.g., 15 hand-constructed, 9 CI-subset) with zero incorrect inclusions or exclusions
- **SC-006**: The test infrastructure adds no more than 2 seconds to the total test suite execution time for the existing 24 matrices (15 hand-constructed + 9 CI-subset)
- **SC-007**: Existing integration tests are successfully refactored to use harness types, with identical pass/fail behavior to their current form

## Assumptions

- The existing `io::registry`, `io::reference`, `io::mtx`, and `validate` modules from Phase 0.4 provide foundational I/O and validation functions. This feature consolidates and supersedes their direct use in tests by wrapping them in a higher-level harness, but the underlying functions remain available.
- The `SolverTest` trait is designed for testing future solver implementations (Phases 2-8). During this phase, the trait is defined and tested with a mock/trivial solver. Real implementations will satisfy the trait in later phases.
- Random matrix generation uses `rand` and `rand_distr` (already in dev-dependencies) and does not require external libraries for sparse matrix generation.
- Positive definite matrix generation uses diagonal dominance as a sufficient condition (no eigenvalue computation needed).
- The test infrastructure is a library-internal testing tool. It lives in `src/` behind a Cargo feature flag (`test-util`), gated with `#[cfg(feature = "test-util")]`. Integration tests access it via `[dev-dependencies] rivrs-sparse = { path = ".", features = ["test-util"] }`. This ensures the test harness is available to both unit and integration tests while excluded from production builds.

## Dependencies

- **Phase 0.4** (complete): Provides `io::registry`, `io::reference`, `io::mtx`, `validate` modules, test data (82 matrices), and basic integration tests that this feature consolidates.
- **faer 0.22**: Sparse matrix types (`SparseColMat`), dense matrix operations (`Mat`, `Col`), and linear algebra primitives.
- **rand + rand_distr** (dev-dependencies): Random number generation for matrix generators.
- **criterion** (dev-dependency): Benchmark framework integration.

## Out of Scope

- Benchmarking framework redesign (Phase 1.2 in ssids-plan.md — separate feature)
- CI pipeline setup (Phase 1.3 in ssids-plan.md)
- Profiling and debug visualization tools (Phase 1.4 in ssids-plan.md)
- SPRAL comparison testing infrastructure (deferred per constitution v1.1.0)
- Formal property-based testing via `proptest` or `quickcheck` — random generators provide the foundation, but crate integration is deferred
- Python bindings or external language interfaces
