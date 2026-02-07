# Feature Specification: Continuous Integration Setup

**Feature Branch**: `007-ci-setup`
**Created**: 2026-02-07
**Status**: Draft
**Input**: User description: "Implement Phase 1.3 in ssids-plan.md — CI pipeline automation for sparse solver testing and validation"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Automated Test Validation on Every PR (Priority: P1)

A developer pushes a commit or opens a pull request against the main branch. The CI pipeline automatically runs the full sparse test suite — including unit tests, integration tests with hand-constructed matrices, and CI-subset SuiteSparse integration tests — across multiple Rust toolchain versions. The developer sees a clear pass/fail status without any manual intervention.

**Why this priority**: Without automated test validation, correctness regressions can silently enter the codebase. Given the project's "Correctness First" principle (Constitution I), this is the most critical CI capability.

**Independent Test**: Can be verified by opening a PR with a passing test suite and confirming that CI reports green status; then opening a PR with a deliberately broken test and confirming CI reports failure with clear diagnostic output.

**Acceptance Scenarios**:

1. **Given** a PR with all tests passing, **When** CI runs, **Then** the test job reports success for each configured toolchain version
2. **Given** a PR that introduces a test failure, **When** CI runs, **Then** the test job reports failure with output identifying the failing test and its assertion message
3. **Given** a PR that modifies only files outside the sparse directory, **When** CI runs, **Then** the sparse test jobs still run (ensuring no accidental breakage from monorepo-level changes)

---

### User Story 2 - Code Quality Enforcement (Priority: P1)

A developer submits a PR. The CI pipeline automatically checks formatting consistency and runs lint analysis, treating lint warnings as errors. This prevents gradual code quality degradation and enforces the project's Rust best practices (Constitution VII).

**Why this priority**: Equal priority with testing — formatting drift and lint warnings left unchecked compound into maintenance burden and can mask real issues.

**Independent Test**: Can be verified by submitting a PR with a formatting violation (e.g., inconsistent indentation) and confirming CI blocks it; then submitting a PR with a clippy warning and confirming CI blocks it.

**Acceptance Scenarios**:

1. **Given** a PR with correctly formatted code and no lint warnings, **When** CI runs, **Then** the lint job reports success
2. **Given** a PR with a formatting violation, **When** CI runs, **Then** the lint job reports failure identifying the file and nature of the violation
3. **Given** a PR with a clippy warning, **When** CI runs, **Then** the lint job reports failure identifying the warning

---

### User Story 3 - Documentation Build Verification (Priority: P2)

A developer adds or modifies public API documentation. The CI pipeline builds rustdoc and treats documentation warnings as errors, ensuring that documentation stays complete and correctly linked as the codebase evolves.

**Why this priority**: Documentation quality is important (Constitution IV, VII) but is lower impact than correctness and lint enforcement. A broken doc link won't cause incorrect solver output.

**Independent Test**: Can be verified by introducing a broken intra-doc link in a doc comment and confirming CI catches it.

**Acceptance Scenarios**:

1. **Given** a PR with well-formed documentation, **When** CI runs, **Then** the doc job reports success
2. **Given** a PR with a broken intra-doc link or missing required doc section, **When** CI runs, **Then** the doc job reports failure identifying the issue

---

### User Story 4 - Feature-Gated Test Coverage (Priority: P2)

The sparse crate has a `test-util` feature flag that gates testing and benchmarking infrastructure. CI validates that the crate compiles and tests pass both with and without this feature enabled, ensuring feature-gated code doesn't break the default build and that the test-util module itself stays functional.

**Why this priority**: Feature flag correctness prevents users who depend on the library (default features) from encountering compilation errors, and ensures the test-util infrastructure remains usable for development.

**Independent Test**: Can be verified by introducing a compilation error only visible under `--features test-util` and confirming CI catches it.

**Acceptance Scenarios**:

1. **Given** a PR where the crate compiles with default features, **When** CI runs, **Then** the default-features test job passes
2. **Given** a PR where the crate compiles with the `test-util` feature, **When** CI runs, **Then** the test-util test job passes
3. **Given** a PR that breaks compilation under `test-util` but not under default features, **When** CI runs, **Then** the test-util job fails while the default job passes, clearly identifying the feature-gated issue

---

### User Story 5 - Performance Regression Awareness (Priority: P3)

After the solver has measurable components (future phases), the CI pipeline provides a mechanism for running benchmarks and comparing against baselines. For now, the pipeline validates that the benchmark harness compiles and can execute, so that when solver code arrives it can immediately be performance-tracked.

**Why this priority**: No solver code exists yet, so benchmark results are not meaningful today. However, ensuring the benchmark infrastructure stays compilable prevents it from bit-rotting before it's needed.

**Independent Test**: Can be verified by confirming the benchmark binary compiles and the benchmark job completes successfully (even if all benchmarks are skipped due to mock solver).

**Acceptance Scenarios**:

1. **Given** a PR that modifies benchmark infrastructure, **When** CI runs, **Then** the benchmark compilation check passes
2. **Given** a PR that breaks the benchmark binary, **When** CI runs, **Then** the benchmark check fails with a clear compilation error

---

### Edge Cases

- What happens when the CI-subset SuiteSparse test matrices are not present in the checkout? The CI pipeline must ensure these matrices are available (they are committed to git, so a standard checkout suffices).
- What happens when a new Rust stable version introduces a breaking change? The pipeline runs against the project's MSRV (1.87) and stable, catching incompatibilities in either direction.
- What happens when a test is flaky (non-deterministic failure)? The pipeline should not retry tests automatically — flaky tests in numerical code indicate a real problem that must be investigated (Constitution I).

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The CI pipeline MUST run the sparse crate's full test suite (`cargo test --all-targets`) on every push to main and every pull request targeting main
- **FR-002**: The CI pipeline MUST test against at least two Rust toolchain versions: the project's MSRV (currently 1.87) and stable
- **FR-003**: The CI pipeline MUST check code formatting via `cargo fmt --check`
- **FR-004**: The CI pipeline MUST run clippy with warnings treated as errors (`-D warnings`) on all targets
- **FR-005**: The CI pipeline MUST build rustdoc with warnings treated as errors (`-D warnings`)
- **FR-006**: The CI pipeline MUST test the crate with the `test-util` feature enabled in addition to default features
- **FR-007**: The CI pipeline MUST verify that the benchmark binary compiles successfully
- **FR-008**: The CI pipeline MUST cache build artifacts between runs to keep execution time reasonable
- **FR-009**: The CI pipeline MUST produce clear, actionable output on failure — identifying the specific test, file, or lint rule that failed
- **FR-010**: The CI pipeline MUST run the sparse jobs independently from control-domain jobs, so that a failure in one domain does not block or obscure results from the other

### Key Entities

- **CI Pipeline**: The automated workflow configuration that triggers on repository events and orchestrates validation jobs
- **Test Suite**: The collection of unit tests, integration tests (hand-constructed and CI-subset SuiteSparse matrices), and feature-gated tests that validate correctness
- **Lint Suite**: The combination of formatting checks and static analysis that enforces code quality standards
- **MSRV (Minimum Supported Rust Version)**: The oldest Rust toolchain version the project guarantees compatibility with (currently 1.87)

## Assumptions

- The CI runs on Linux (ubuntu-latest) only. Multi-OS support (macOS, Windows) is deferred until the solver has platform-specific code paths that require cross-platform validation.
- The MSRV is 1.87, matching the Rust 2024 edition requirement already specified in Cargo.toml.
- SuiteSparse CI-subset matrices (~73MB) are committed directly to git and available in a standard checkout without additional download steps.
- SPRAL comparison tests are deferred to Phases 2-8 per the Phase 0.3 decision; the CI pipeline does not need to build or invoke SPRAL.
- Benchmark execution (not just compilation) is deferred until there is actual solver code producing meaningful performance data.
- The existing `.github/workflows/ci.yml` already contains sparse test, lint, and doc jobs. This feature enhances and extends that existing configuration rather than creating it from scratch.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Every pull request targeting main receives automated pass/fail feedback on test correctness within 10 minutes of being opened
- **SC-002**: The CI pipeline catches 100% of formatting violations and clippy warnings before merge
- **SC-003**: The CI pipeline validates the crate under both MSRV and stable Rust, ensuring no accidental use of unstable or too-new features
- **SC-004**: The CI pipeline validates both default-feature and test-util-feature builds, catching feature-gating errors
- **SC-005**: A developer can determine the cause of any CI failure from the CI output alone, without needing to reproduce locally (clear error messages with file paths and line numbers)
- **SC-006**: The CI pipeline completes all jobs for the sparse domain within 10 minutes on a clean cache, and within 5 minutes with a warm cache
