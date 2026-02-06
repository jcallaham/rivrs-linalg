# Feature Specification: SPRAL Golden Results Generation

**Feature Branch**: `003-spral-golden-results`
**Created**: 2026-02-06
**Status**: Deferred (see decision.md)
**Input**: User description: "implement phase 0.3 in docs/ssids-plan.md"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Generate Reference Results for All Test Matrices (Priority: P1)

A solver developer wants to run SPRAL on every test matrix collected in Phase 0.2 (82 matrices total: 15 hand-constructed + 67 SuiteSparse) and capture comprehensive reference outputs. These "golden results" will serve as the ground truth for validating the Rust SSIDS implementation across all future development phases.

**Why this priority**: Without golden results, there is no way to verify correctness of the Rust solver. Every subsequent phase depends on having known-good reference outputs to compare against. This is the single most critical deliverable.

**Independent Test**: Can be fully tested by running the SPRAL driver program on a single test matrix and verifying the output JSON contains all required fields (analysis, factorization, solve metrics). Delivers immediate value as a reference point.

**Acceptance Scenarios**:

1. **Given** the 82 test matrices from Phase 0.2 and a working SPRAL build, **When** the driver program is executed against each matrix, **Then** a structured result file is produced for each matrix containing analysis, factorization, and solve phase outputs.
2. **Given** a generated result file, **When** an inspector examines the file, **Then** it contains the matrix name, SPRAL options used, timing data, pivot statistics, inertia, memory usage, and forward/backward error measures.
3. **Given** any hand-constructed matrix with a known factorization, **When** SPRAL results are compared to the known factorization, **Then** SPRAL's inertia and backward error are consistent with the known answer.

---

### User Story 2 - Verify Reproducibility of Results (Priority: P2)

A solver developer wants assurance that the golden results are deterministic and reproducible. Running SPRAL twice on the same matrix with the same options must yield identical numerical results (bit-for-bit on analysis outputs, within rounding tolerance on numerical outputs).

**Why this priority**: Reproducibility is essential for the results to function as a reliable test oracle. If results vary between runs, they cannot serve as a stable reference.

**Independent Test**: Can be tested by running the driver program twice on a representative subset of matrices and comparing outputs. Delivers confidence in the test infrastructure.

**Acceptance Scenarios**:

1. **Given** a test matrix and fixed SPRAL options, **When** the driver is run twice, **Then** the structural outputs (elimination tree, supernodes, predicted fill) are bit-for-bit identical.
2. **Given** a test matrix and fixed SPRAL options, **When** the driver is run twice, **Then** the numerical outputs (residuals, inertia, pivot counts) agree within machine epsilon.

---

### User Story 3 - Produce Summary Report and Identify Challenging Cases (Priority: P3)

A solver developer wants a summary view across all test matrices showing timing statistics, failure modes, and which matrices are "hard" (high delayed pivots, near-singular, or requiring 2x2 pivots). This summary guides prioritization of which cases to focus on first during Rust implementation.

**Why this priority**: A summary provides strategic insight into the test suite — identifying easy wins for early validation and hard cases that will stress-test the APTP implementation. Lower priority because individual results are sufficient for basic correctness checking.

**Independent Test**: Can be tested by generating the summary from existing result files and verifying it contains aggregated statistics. Delivers a high-level map of the test landscape.

**Acceptance Scenarios**:

1. **Given** result files for all successfully-factored matrices, **When** a summary is generated, **Then** it contains per-matrix rows with size, nnz, factor time, solve time, backward error, delayed pivots, and inertia.
2. **Given** the summary, **When** a developer reviews it, **Then** matrices are categorized by difficulty (trivial, easy, moderate, hard) and failure modes (if any) are clearly noted.
3. **Given** any matrix that SPRAL fails to factor, **When** the failure is recorded, **Then** the result file captures the error code and diagnostic message from SPRAL.

---

### Edge Cases

- What happens when SPRAL fails to factorize a matrix (e.g., structural singularity)?
  The driver captures the error code, diagnostic message, and partial results (analysis phase data), and writes a result file marked as failed.
- What happens when a matrix is too large for available memory?
  The driver captures the out-of-memory condition and records it as a known limitation in the result file.
- What happens for the singular hand-constructed matrix (`singular-3.mtx`)?
  The result captures SPRAL's behavior on singular input (expected: detected singularity with zero eigenvalues reported in inertia).
- What happens when solving with a random right-hand side produces poor backward error?
  The result records the backward error as-is; the summary flags matrices where backward error exceeds 1e-10 as requiring investigation.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The system MUST build SPRAL from source as a Dockerfile build stage (CPU-only, meson, release mode) so the library is pre-installed in the development container image. Build dependencies (python3, meson, ninja, libopenblas-dev, libmetis-dev, libhwloc-dev) are installed in the Dockerfile. The build configuration is documented and reproducible via the Dockerfile itself.
- **FR-002**: The system MUST provide a driver program that accepts a Matrix Market file and SPRAL options as input, runs the full analyze-factorize-solve pipeline, and writes structured output.
- **FR-003**: For each test matrix, the driver MUST capture analysis phase outputs: predicted fill-in (nnz), predicted floating-point operations, number of partitions, maximum tree depth, and number of supernodes. (Note: full elimination tree parent pointers and supernode membership arrays are internal to SPRAL's `akeep` and not exposed via the C API; only aggregate statistics are available.)
- **FR-004**: For each test matrix, the driver MUST capture factorization phase outputs: actual fill-in (nnz), actual floating-point operations, number of delayed pivots, and inertia (positive, negative, zero eigenvalue counts).
- **FR-005**: For each test matrix, the driver MUST generate a random right-hand side vector, solve the system, and capture forward error and backward error.
- **FR-006**: For each test matrix, the driver MUST capture timing data for each phase (analyze, factorize, solve) and peak memory usage.
- **FR-007**: The driver MUST record the SPRAL version, build configuration, and options used, so results can be tied to a specific SPRAL configuration.
- **FR-008**: All results MUST be written in a structured format (one file per matrix) to a designated output directory.
- **FR-009**: The system MUST produce a summary file aggregating key metrics across all matrices (size, nnz, timing, errors, pivot statistics, difficulty category).
- **FR-010**: When SPRAL fails on a matrix, the driver MUST capture the failure mode (error code, diagnostic message) and still write a result file indicating the failure.
- **FR-011**: The system MUST support running on the CI subset (10 matrices) as a fast validation path, in addition to the full set.
- **FR-012**: Results for the hand-constructed matrices and CI subset MUST be stored in the git repository for use in automated testing.
- **FR-013**: Results for the full SuiteSparse set MUST be stored outside git (in the gitignored `test-data/` area or as a reference archive), consistent with the three-tier storage strategy from Phase 0.2.

### Key Entities

- **SPRAL Reference Result**: Per-matrix output containing analysis outputs, factorization outputs, solve outputs, timing, memory, and SPRAL configuration. Keyed by matrix name.
- **SPRAL Build Configuration**: The compiler, flags, SPRAL version/commit, and runtime options used to produce the results. Recorded once and referenced by all result files.
- **Summary Report**: Aggregated view of all results with per-matrix rows and difficulty categorization.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: SPRAL successfully factors at least 95% of the 82 test matrices (allowing for a small number of expected failures on singular or extreme matrices).
- **SC-002**: Complete result files are produced for every matrix attempted (including failed ones), with all required data fields populated.
- **SC-003**: Results are reproducible: running the driver twice on the same matrix produces identical structural outputs and numerical outputs that agree within 1e-14 relative tolerance.
- **SC-004**: A summary report is produced that categorizes all matrices by difficulty and highlights any failure modes.
- **SC-005**: The golden results for the CI subset (10 matrices) and hand-constructed matrices (15 matrices) are committed to the repository and accessible for automated testing.

## Assumptions

- SPRAL is built from source as a Dockerfile build stage, so it is pre-installed in the container image. The Dockerfile already has gfortran and gcc; additional deps (python3, meson, ninja, libopenblas-dev, libmetis-dev, libhwloc-dev) are added in the build stage.
- The SPRAL version used is the one available at `references/spral/` (BSD-3 licensed source).
- "Reproducibility" refers to deterministic results on the same hardware with the same build; cross-platform reproducibility is not required.
- The right-hand side vector for solve testing is generated deterministically as `b = A * ones(n)` (known exact solution `x = ones(n)`), matching SPRAL's own driver approach.
- Results for the full SuiteSparse collection require the archive to be extracted (same as Phase 0.2's extraction mechanism).
- The driver program is a C program that links against SPRAL via its C API; it is test infrastructure, not part of the Rust solver.

## Clarifications

### Session 2026-02-06

- Q: Should SPRAL be built inside the running container at runtime, as a Dockerfile build stage, as a separate Docker image, or manually? → A: Dockerfile build stage (Option B). SPRAL is built as a new multi-stage step in the existing Dockerfile, installing deps and compiling with meson. The library is pre-installed in the image, cached by Docker layers, and available immediately on container start.

## Dependencies

- **Phase 0.2 (Complete)**: Test matrix collection (82 matrices in three-tier storage).
- **SPRAL Source**: Available at `references/spral/` (BSD-3 licensed).
- **Build Toolchain**: Fortran compiler (gfortran), C compiler, and build system (meson) for building SPRAL.
