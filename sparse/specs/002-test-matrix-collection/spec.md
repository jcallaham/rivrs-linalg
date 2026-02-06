# Feature Specification: Test Matrix Collection Assembly

**Feature Branch**: `002-test-matrix-collection`
**Created**: 2026-02-05
**Status**: Draft
**Input**: User description: "Review Phase 0.2 in the ssids-plan.md doc and determine what needs to be done to compile a comprehensive test suite."

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Run Solver Tests Against Known Matrices (Priority: P1)

A developer working on the SSIDS solver needs to validate factorization correctness against matrices with known properties and known correct results. They load a hand-constructed matrix (e.g., a 10x10 arrow matrix), run the solver, and compare the output against pre-computed reference values. This catches regressions and verifies basic algorithmic correctness.

**Why this priority**: Without small, well-understood test matrices with known factorizations, there is no way to verify that any solver implementation is correct. This is the foundation of all testing.

**Independent Test**: Can be fully tested by loading any hand-constructed matrix from the collection, verifying it parses correctly, confirming its documented properties match, and comparing against known factorization results.

**Acceptance Scenarios**:

1. **Given** a hand-constructed 10x10 arrow matrix in Matrix Market format, **When** the matrix is loaded through the test infrastructure, **Then** it parses into the correct dimensions and nonzero pattern, and its metadata accurately describes its properties (symmetric, indefinite, arrow structure).
2. **Given** a hand-constructed tridiagonal matrix with a known LDL^T factorization, **When** the solver computes the factorization, **Then** the result matches the pre-computed reference to within documented tolerance.
3. **Given** a matrix designed to force delayed pivots, **When** the solver processes it, **Then** the expected number of delayed pivots occurs and the factorization remains numerically stable.

---

### User Story 2 - Validate Solver on Real-World Sparse Problems (Priority: P2)

A developer needs to verify the solver handles real-world sparse symmetric indefinite systems from established collections (SuiteSparse Matrix Collection). They select matrices from categorized difficulty levels (easy indefinite, hard indefinite) and run the solver, comparing residuals and inertia against reference results.

**Why this priority**: Real-world matrices expose failure modes that hand-constructed matrices cannot, including fill-in patterns, conditioning issues, and scaling challenges. This is essential for production-quality validation.

**Independent Test**: Can be fully tested by downloading a SuiteSparse matrix, loading it via the test infrastructure, verifying its metadata, and running a factorization with residual checks against documented reference values.

**Acceptance Scenarios**:

1. **Given** a SuiteSparse matrix classified as "easy indefinite," **When** it is loaded and factored, **Then** the backward error is below 1e-10 and the computed inertia matches the reference.
2. **Given** a SuiteSparse matrix classified as "hard indefinite," **When** it is loaded and factored, **Then** either the factorization succeeds with backward error below 1e-8, or a well-documented failure mode is reported with diagnostics.
3. **Given** a positive definite matrix from SuiteSparse, **When** the solver runs in its fast-path (positive definite) mode, **Then** the factorization completes without any pivot delays and with backward error below 1e-12.

---

### User Story 3 - Browse and Select Matrices by Properties (Priority: P2)

A developer diagnosing a specific solver behavior (e.g., excessive delayed pivots, poor conditioning) needs to find test matrices that exhibit that property. They query the metadata index by properties such as size range, definiteness, difficulty, expected delayed pivots, or source, and get a filtered list of suitable matrices.

**Why this priority**: As the solver matures, targeted testing against matrices with specific characteristics becomes essential for debugging and optimization. This enables efficient root-cause analysis.

**Independent Test**: Can be fully tested by querying the metadata index for matrices matching specific property criteria and verifying the returned list is non-empty and all results satisfy the query.

**Acceptance Scenarios**:

1. **Given** the complete metadata index, **When** a developer filters for symmetric indefinite matrices smaller than 1000x1000, **Then** the results include all hand-constructed indefinite matrices and any SuiteSparse matrices meeting the criteria.
2. **Given** the metadata index, **When** a developer filters for matrices with "high" expected delayed pivots, **Then** at least 3 matrices are returned, each with documented justification for the classification.

---

### User Story 4 - Reproduce Academic Paper Results (Priority: P3)

A developer implementing the APTP algorithm wants to validate against the specific test problems used in the original papers. They load matrices mentioned in key APTP publications and compare solver output against results reported in those papers.

**Why this priority**: Reproducing published results provides the strongest evidence that the algorithm is correctly implemented. However, this depends on first having basic correctness (P1) and real-world validation (P2).

**Independent Test**: Can be fully tested by loading each paper-referenced matrix, running the solver, and comparing key metrics (residual, inertia, pivot delay count) against values reported in the publications.

**Acceptance Scenarios**:

1. **Given** a matrix used in Hogg, Duff & Lopez (2020), **When** the solver processes it with comparable settings, **Then** the backward error and number of delayed pivots are within the same order of magnitude as reported in the paper.

---

### Edge Cases

- What happens when a Matrix Market file is malformed (truncated, wrong header, inconsistent dimensions)?
- How does the system handle matrices that are listed in metadata but whose files are missing from disk?
- What happens when a matrix is labeled symmetric but the file contains asymmetric entries?
- How does the system handle a zero-dimension or empty matrix?
- What happens when two matrices from different sources have the same name?

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The collection MUST include at least 15 hand-constructed matrices ranging from 5x5 to 20x20 with known exact factorizations, covering arrow, block diagonal, tridiagonal, and stress-test structures.
- **FR-002**: The collection MUST include at least 30 "easy indefinite" symmetric matrices sourced from the SuiteSparse Matrix Collection.
- **FR-003**: The collection MUST include at least 15 "hard indefinite" symmetric matrices sourced from SuiteSparse, representing challenging real-world problems.
- **FR-004**: The collection MUST include at least 20 positive definite matrices from SuiteSparse for fast-path validation.
- **FR-005**: All matrices MUST be stored in standard Matrix Market (.mtx) format.
- **FR-006**: Every matrix in the collection MUST have a corresponding metadata entry describing its name, source, dimensions, nonzero count, symmetry, definiteness, difficulty classification, and any known reference results.
- **FR-007**: The metadata index MUST be stored in a single machine-readable file (JSON) that can be loaded programmatically by test infrastructure.
- **FR-008**: Hand-constructed matrices MUST include pre-computed exact factorizations (L, D, permutation) documented alongside the matrix.
- **FR-009**: The collection MUST include matrices specifically mentioned in the core APTP papers (Hogg, Duff & Lopez 2020; Hogg, Ovtchinnikov & Scott 2016).
- **FR-010**: The collection MUST include at least 5 matrices known to cause failures with static pivoting strategies ("killer" cases).
- **FR-011**: The collection MUST span a size range from 5x5 to at least 100,000x100,000 problems.
- **FR-012**: The collection MUST include matrices extracted from or inspired by SPRAL's test suite.
- **FR-013**: Hand-constructed stress-test matrices MUST include cases designed for: maximum delayed pivots, worst-case fill-in, and ill-conditioning.
- **FR-014**: The test data directory MUST be organized into clearly separated subdirectories by source and category (hand-constructed, suitesparse/easy-indefinite, suitesparse/hard-indefinite, suitesparse/positive-definite, spral-tests, interior-point).

### Key Entities

- **Test Matrix**: A sparse matrix stored in Matrix Market format with associated metadata. Key attributes: name, source, dimensions (n), nonzero count (nnz), symmetry type, definiteness, difficulty classification.
- **Matrix Metadata**: Structured record for each matrix containing identification, structural properties, difficulty classification, and any available reference results (residuals, inertia, timing).
- **Metadata Index**: The complete collection of all matrix metadata entries, stored as a single queryable JSON file.
- **Reference Result**: Known-correct factorization output for a matrix, including inertia (positive/negative/zero eigenvalue counts), backward error, number of delayed pivots, and factorization statistics.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: The collection contains a minimum of 70 test matrices across all categories.
- **SC-002**: All matrices in the collection can be loaded and parsed without errors by standard Matrix Market readers.
- **SC-003**: Every matrix has a complete metadata entry with all required fields populated (no missing or placeholder values).
- **SC-004**: The collection covers a size range spanning at least 4 orders of magnitude (from single-digit to 100,000+ dimensions).
- **SC-005**: At least 5 "killer" matrices are identified and documented — matrices known to cause failures in solvers that use only static pivoting.
- **SC-006**: Hand-constructed matrices include verified exact factorizations that can be independently checked by hand or symbolic computation.
- **SC-007**: A developer unfamiliar with the collection can find a matrix matching specific property criteria (size, definiteness, difficulty) within 2 minutes using the metadata index.

## Assumptions

- The SuiteSparse Matrix Collection is publicly accessible and matrices can be downloaded in Matrix Market format without licensing restrictions.
- SPRAL test matrices can be extracted from the SPRAL repository (BSD-3 licensed) without licensing issues.
- Matrices referenced in academic papers are either available in SuiteSparse or can be reconstructed from the paper descriptions.
- The test data will be stored in the repository (or fetched by a script) and does not need a separate hosting solution for matrices under ~100MB total.
- "Easy" vs "hard" indefinite classification follows the conventions used in the SPRAL papers and SuiteSparse metadata (based on conditioning, fill-in ratio, and pivot behavior).

## Dependencies

- SuiteSparse Matrix Collection web API or download infrastructure.
- SPRAL repository (already available at `references/spral/`).
- Matrix Market format parser (Rust crate `matrixmarket-rs` or equivalent).
- Phase 0.1 reference library (complete) — provides paper references and algorithm documentation needed to identify paper-specific test matrices.

## Scope Boundaries

**In scope:**
- Collecting and organizing test matrices from the four sources defined in ssids-plan.md Phase 0.2
- Creating metadata for every matrix in the collection
- Building hand-constructed matrices with known exact factorizations
- Organizing the test-data directory structure
- Documenting "killer" cases and difficulty classifications

**Out of scope:**
- Running SPRAL to generate golden reference results (that is Phase 0.3)
- Building the Rust test harness or test runner infrastructure (that is Phase 0.4)
- Performance benchmarking
- Implementing any solver code
- Downloading matrices larger than the repository can reasonably host (>100MB individual files); such matrices should be documented in metadata with download instructions instead
