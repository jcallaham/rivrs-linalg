# Feature Specification: Repository Setup for Solver Development

**Feature Branch**: `004-repo-setup`
**Created**: 2026-02-06
**Status**: Draft
**Input**: User description: "implement phase 0.4 in dev/ssids-plan.md"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Load Test Matrices in Rust Tests (Priority: P1)

A solver developer writes an integration test that loads a hand-constructed `.mtx` matrix file and its companion `.json` factorization data from `test-data/`, converts them into faer sparse matrix types, and asserts properties (dimensions, symmetry, known values). This is the foundational capability — all subsequent solver development depends on being able to load test matrices in Rust code.

**Why this priority**: Without the ability to load test matrices into Rust types, no solver testing is possible. Every subsequent phase depends on this.

**Independent Test**: Can be fully tested by running `cargo test` and verifying that hand-constructed matrices load correctly with expected dimensions and values.

**Acceptance Scenarios**:

1. **Given** the test-data directory contains Matrix Market (`.mtx`) files, **When** a test calls the matrix loading function with a matrix name, **Then** the function returns a faer-compatible sparse symmetric matrix with correct dimensions and nonzero count.
2. **Given** a hand-constructed matrix has a companion `.json` file with known factorization data, **When** a test loads the reference data, **Then** the expected inertia, permutation, L factor entries, and D diagonal are accessible as typed Rust structures.
3. **Given** a matrix name that does not exist in test-data, **When** a test attempts to load it, **Then** a clear error message is returned indicating the matrix was not found.

---

### User Story 2 - CI Pipeline Validates Sparse Code (Priority: P2)

When a developer pushes code or opens a pull request, the CI pipeline automatically runs the sparse domain's tests, clippy lints, rustfmt checks, and documentation build — catching regressions before merge. This extends the existing CI (which only covers `control/`) to also cover `sparse/`.

**Why this priority**: Automated quality gates prevent regressions and enforce code standards as solver development begins in earnest.

**Independent Test**: Can be tested by pushing a commit to a PR branch and verifying that sparse CI jobs run and report pass/fail status.

**Acceptance Scenarios**:

1. **Given** a pull request targeting `main`, **When** the CI pipeline runs, **Then** `cargo test`, `cargo clippy`, `cargo fmt --check`, and `cargo doc` are all executed for the sparse domain.
2. **Given** a test failure in the sparse domain, **When** CI runs, **Then** the pipeline reports failure with clear output identifying the failing test.
3. **Given** the sparse domain's tests depend on hand-constructed test matrices in `test-data/hand-constructed/` and CI-subset SuiteSparse matrices in `test-data/suitesparse-ci/`, **When** CI runs, **Then** those files are available (checked into the repository) and tests that load them pass.

---

### User Story 3 - Numerical Validation Utilities (Priority: P3)

A solver developer writing a factorization test uses shared validation utilities to check reconstruction error (`||P^T A P - L D L^T|| / ||A||`), backward error (`||Ax - b|| / (||A|| ||x|| + ||b||)`), and inertia correctness. These utilities are reusable across all future solver phases.

**Why this priority**: Standardized numerical validation prevents inconsistent or incorrect test assertions across phases. Having these ready before Phase 1 means solver tests can be written immediately.

**Independent Test**: Can be tested by applying validation functions to hand-constructed matrices with known factorizations and confirming they report pass/fail correctly.

**Acceptance Scenarios**:

1. **Given** a hand-constructed matrix with a known LDL^T factorization, **When** the reconstruction error function is called with the correct factors, **Then** the relative error is below machine epsilon scale (< 10^-12).
2. **Given** an intentionally wrong factorization (e.g., perturbed L factor), **When** the reconstruction error function is called, **Then** it reports an error above the threshold.
3. **Given** a known solution to Ax = b, **When** the backward error function is called, **Then** it returns a value consistent with the expected residual.

---

### User Story 4 - Benchmark Scaffold (Priority: P4)

A developer can run `cargo bench` against placeholder benchmarks that load test matrices and time a no-op analysis pass, establishing the benchmarking infrastructure for future phases to fill in with actual solver operations.

**Why this priority**: Lower priority because actual benchmarking only becomes meaningful once solver components exist, but the scaffold ensures the infrastructure is ready.

**Independent Test**: Can be tested by running `cargo bench` and confirming it completes without error, producing timing output.

**Acceptance Scenarios**:

1. **Given** the benchmark suite is configured, **When** `cargo bench` is run, **Then** at least one benchmark executes and produces timing output.
2. **Given** a test matrix is loaded in a benchmark, **When** the benchmark runs, **Then** the matrix loading time is reported separately from the (placeholder) operation time.

---

### Edge Cases

- What happens when a `.mtx` file is malformed (e.g., missing header, wrong format)? The loader returns a descriptive error rather than panicking.
- What happens when a `.json` factorization file references dimensions inconsistent with its companion `.mtx`? Validation catches the mismatch at load time.
- What happens when SuiteSparse CI-subset matrices are not present (e.g., shallow clone)? Tests that require them fail with a clear message pointing to the data source, since these are expected to be in-repo for CI.
- What happens when the full SuiteSparse collection (gitignored, 67 matrices) is not extracted? Tests using those matrices are gated behind a runtime path check and skipped gracefully.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The project MUST provide a lightweight custom Matrix Market parser for loading symmetric coordinate-format `.mtx` files into faer-compatible sparse matrix types (no external parsing crate).
- **FR-002**: The project MUST provide a Rust module for loading companion `.json` factorization data (inertia, permutation, L factor, D diagonal) into typed Rust structures.
- **FR-003**: The project MUST provide a test matrix registry that maps matrix names to their file paths and metadata, using `test-data/metadata.json` as the source of truth.
- **FR-004**: The project MUST provide fully working numerical validation functions — reconstruction error, backward error, and inertia comparison — tested end-to-end against hand-constructed matrices with known factorizations.
- **FR-005**: The CI pipeline MUST run test, lint (clippy), format check (rustfmt), and documentation build for the sparse domain on every pull request.
- **FR-006**: The project MUST have at least one integration test that loads a hand-constructed matrix, verifies its properties, and runs a validation function.
- **FR-007**: The project MUST have a benchmark scaffold using criterion that loads a test matrix.
- **FR-008**: Tests that depend on the full SuiteSparse collection (67 gitignored matrices) MUST be skipped gracefully when the data is not present, rather than failing. Tests against hand-constructed and CI-subset matrices MUST always run (these are committed to the repository).

### Key Entities

- **TestMatrix**: A named test case consisting of a sparse matrix (from `.mtx`), its metadata (from `metadata.json`), and optional reference results (from `.json` factorization files).
- **ReferenceFactorization**: The known-correct LDL^T factorization of a hand-constructed matrix, including L factor, D diagonal, permutation, and inertia.
- **NumericalValidator**: A collection of validation functions with configurable tolerances for reconstruction error, backward error, and inertia checks.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: `cargo test` in the sparse directory passes, including at least one test that loads a hand-constructed matrix and validates its properties.
- **SC-002**: `cargo clippy --all-targets -- -D warnings` produces zero warnings in the sparse domain.
- **SC-003**: `cargo fmt --check` passes in the sparse domain.
- **SC-004**: `cargo doc --no-deps` builds without warnings in the sparse domain.
- **SC-005**: CI pipeline includes sparse domain jobs and they pass on the main branch.
- **SC-006**: `cargo bench` executes at least one benchmark without error.
- **SC-007**: All 15 hand-constructed matrices can be loaded and their properties verified in under 1 second total.
- **SC-008**: Reconstruction error validation returns < 10^-12 for all hand-constructed matrices with known-correct factorizations.

## Clarifications

### Session 2026-02-06

- Q: Should Matrix Market parsing use an external crate or a custom parser? → A: Lightweight custom parser scoped to symmetric coordinate format (no external crate).
- Q: Which test-data tiers should CI run tests against? → A: Both hand-constructed (15) and CI-subset SuiteSparse (10) matrices. Full collection (67) skipped gracefully when absent.
- Q: Should numerical validation utilities be fully implemented or stubbed in Phase 0.4? → A: Full working implementations (reconstruction error, backward error, inertia comparison) tested end-to-end against hand-constructed known factorizations.

## Assumptions

- The faer crate (v0.22) provides sufficient sparse matrix types (CSC format) for representing loaded matrices, plus the dense arithmetic needed for validation: `permute_self_adjoint` for P^T A P, `.to_dense()` for sparse-to-dense conversion, operator overloads for dense multiply/subtract, `.norm_l2()` for Frobenius norm, and `sparse_dense_matmul` for SpMV. Validation functions will leverage these rather than implementing matrix arithmetic from scratch. Dense conversion is appropriate for hand-constructed matrices (max 20x20); large-matrix validation will be revisited in later phases.
- Matrix Market parsing will use a lightweight custom parser scoped to symmetric coordinate format (the only format used by the project's test data). No external crate dependency for this.
- The existing `test-data/metadata.json` schema is stable and will not change during this phase.
- The existing `error.rs` module provides a foundation for loader error types; new error variants will extend `SparseError`.
- CI runs on GitHub Actions using the existing workflow file structure (extending `ci.yml`).

## Dependencies

- **Phase 0.2 (complete)**: Test matrices and metadata.json must be in place — they are.
- **faer 0.22**: Must be available as a dependency — it is (already in Cargo.toml).
- **criterion 0.5**: Must be available as a dev-dependency — it is (already in Cargo.toml).
