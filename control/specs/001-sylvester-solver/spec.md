# Feature Specification: Sylvester Equation Solver

**Feature Branch**: `001-sylvester-solver`
**Created**: 2026-01-27
**Status**: Draft
**Input**: User description: "Implement the Sylvester solver as outlined in docs/sylvester-plan.md. Start with a proper setup for a Rust project and then proceed with the rest of the plan"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Solve Standard Sylvester Equations (Priority: P1)

A control systems engineer needs to solve the continuous-time Sylvester equation AX + XB = C for matrix X, where A, B, and C are real matrices. This is the fundamental operation required for many control theory applications including observer design, model reduction, and controller synthesis.

**Why this priority**: This is the most basic and essential functionality. Without this, the library cannot solve any Sylvester equations. This forms the foundation for all other features.

**Independent Test**: Can be fully tested by providing known test matrices (A, B, C) from academic literature (e.g., Bartels-Stewart benchmark problems) and verifying the solution X satisfies the equation AX + XB = C within numerical tolerance.

**Acceptance Scenarios**:

1. **Given** three real matrices A (n×n), B (m×m), and C (n×m) where (A, -B) has no common eigenvalues, **When** the solver is invoked with these matrices, **Then** it returns a solution matrix X (n×m) such that ||AX + XB - C|| < ε for specified tolerance ε
2. **Given** matrices with small dimensions (e.g., 2×2, 3×3), **When** the solver is invoked, **Then** the solution is computed in under 10 milliseconds
3. **Given** well-conditioned matrices with known analytical solutions from textbooks (e.g., Golub & Van Loan examples), **When** the solver is invoked, **Then** the computed solution matches the analytical solution within 10⁻¹² relative error

---

### User Story 2 - Handle Discrete-Time Sylvester Equations (Priority: P2)

A control systems engineer needs to solve the discrete-time Sylvester equation AXB + X = C for matrix X, which arises in discrete-time control systems and digital signal processing applications.

**Why this priority**: Discrete-time systems are common in digital control implementations. This extends the solver's applicability to a major class of control problems without requiring a completely separate implementation.

**Independent Test**: Can be tested independently by providing discrete-time test cases (e.g., from SLICOT benchmark suite) and verifying the solution satisfies AXB + X = C within numerical tolerance.

**Acceptance Scenarios**:

1. **Given** three real matrices A (n×n), B (m×m), and C (n×m) for a discrete-time problem, **When** the discrete-time solver is invoked, **Then** it returns a solution matrix X such that ||AXB + X - C|| < ε
2. **Given** benchmark problems from control systems literature, **When** the discrete-time solver is invoked, **Then** the solution accuracy matches or exceeds SLICOT reference implementation

---

### User Story 3 - Detect and Report Singular Cases (Priority: P3)

A user attempts to solve a Sylvester equation where the uniqueness condition is violated (A and -B have common eigenvalues), and the system provides clear diagnostic information about why the equation cannot be solved uniquely.

**Why this priority**: Error detection and reporting is critical for numerical software, but users can work around it manually if needed. This prevents silent failures and helps users diagnose problem formulations.

**Independent Test**: Can be tested by constructing matrices with known common eigenvalues and verifying that the solver detects this condition and returns an informative error.

**Acceptance Scenarios**:

1. **Given** matrices A and B that share a common eigenvalue, **When** the solver is invoked, **Then** it returns an error indicating the uniqueness condition is violated and reports an estimate of the problematic eigenvalue
2. **Given** nearly-singular matrix pencils (A, -B) with small separation constants, **When** the solver is invoked, **Then** it provides a warning about potential numerical ill-conditioning with an estimated condition number

---

### User Story 4 - Solve Large-Scale Problems Efficiently (Priority: P4)

A researcher needs to solve Sylvester equations with matrices of dimension 500×500 or larger, requiring efficient memory usage and computational performance comparable to established libraries.

**Why this priority**: Performance is important for production use but not critical for initial validation. Users can start with smaller problems while performance is optimized.

**Independent Test**: Can be tested using large randomly generated matrices and comparing execution time and memory usage against SLICOT benchmarks.

**Acceptance Scenarios**:

1. **Given** matrices of dimension 500×500, **When** the solver is invoked, **Then** the solution is computed in under 5 seconds on standard hardware
2. **Given** matrices of dimension 1000×1000, **When** the solver is invoked, **Then** peak memory usage does not exceed 100 MB beyond input storage
3. **Given** the same test cases used to benchmark SLICOT, **When** the Rust solver is invoked, **Then** execution time is within 2× of SLICOT performance

---

### Edge Cases

- What happens when input matrices have NaN or Inf values?
- How does the solver handle matrices with extremely small or large condition numbers (10⁻¹⁵ to 10¹⁵)?
- What happens when A or B matrices are exactly singular?
- How does the solver perform when matrices have repeated eigenvalues?
- What happens when input matrices have inconsistent dimensions?
- How does the solver handle empty matrices or matrices with zero dimensions?
- What happens when the right-hand side matrix C has special structure (zero, identity, sparse)?

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST solve the continuous-time Sylvester equation AX + XB = C for real matrices A, B, and C
- **FR-002**: System MUST solve the discrete-time Sylvester equation AXB + X = C for real matrices A, B, and C
- **FR-003**: System MUST validate that input matrix dimensions are compatible (A is n×n, B is m×m, C is n×m)
- **FR-004**: System MUST detect when the uniqueness condition is violated (common eigenvalues between A and -B) and return an appropriate error
- **FR-005**: System MUST use numerically stable algorithms based on Schur decomposition (following Bartels-Stewart method from LAPACK documentation)
- **FR-006**: System MUST provide error estimates or residual norms for computed solutions
- **FR-007**: System MUST support both f32 and f64 floating-point precision through generic implementations
- **FR-008**: System MUST return informative errors for invalid inputs (NaN, Inf, dimension mismatches, null pointers)
- **FR-009**: System MUST document the numerical algorithm used and cite academic references (e.g., Bartels-Stewart algorithm, Golub & Van Loan section 7.6)
- **FR-010**: System MUST implement workspace allocation patterns compatible with faer's memory management
- **FR-011**: System MUST provide test cases comparing outputs against SLICOT benchmark results to validate correctness
- **FR-012**: System MUST validate solutions by computing residual norms ||AX + XB - C|| or ||AXB + X - C||

### Key Entities

- **Sylvester Solver**: Represents the computational routine that solves AX + XB = C or AXB + X = C, takes input matrices and configuration parameters, returns solution matrix and diagnostic information
- **Matrix A**: Square coefficient matrix (n×n) on the left side of the equation, must be real-valued, controls the dynamics of the first term
- **Matrix B**: Square coefficient matrix (m×m) on the right side of the equation, must be real-valued, controls the dynamics of the second term
- **Matrix C**: Rectangular right-hand side matrix (n×m), must be real-valued, represents the forcing term or target
- **Matrix X**: Solution matrix (n×m) to be computed, satisfies the Sylvester equation within numerical tolerance
- **Residual**: Computed quantity ||AX + XB - C|| or ||AXB + X - C|| that measures solution accuracy, provided as diagnostic output
- **Condition Estimate**: Numerical estimate of problem conditioning that helps users assess solution reliability, based on eigenvalue separation for (A, -B)

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Solver produces solutions for standard benchmark problems from academic literature (e.g., Bartels-Stewart test cases) with residual norms below 10⁻¹² in double precision
- **SC-002**: Solver handles matrix dimensions from 2×2 up to 1000×1000 without numerical failures
- **SC-003**: Computed solutions for test cases match SLICOT reference outputs within 10⁻¹⁰ relative difference
- **SC-004**: Solver execution time for 100×100 matrices is under 100 milliseconds on standard hardware
- **SC-005**: Solver correctly detects and reports error conditions (singular cases, dimension mismatches) in 100% of invalid test cases
- **SC-006**: All public functions have documentation citing academic references (papers, textbooks, or LAPACK documentation) used for the implementation
- **SC-007**: Test suite achieves 95% code coverage for the solver module
- **SC-008**: Solver passes all test cases when using both f32 and f64 precision

## Assumptions *(optional)*

- Users are familiar with basic linear algebra and control systems theory (understand what a Sylvester equation is)
- Input matrices are stored in standard dense format (not sparse), with column-major or row-major layout compatible with faer
- The primary use case is small to medium-sized problems (n, m < 1000) typical of control systems applications
- Users can provide their own matrices and are responsible for ensuring numerical well-conditioning when possible
- The Rust ecosystem has mature linear algebra libraries (faer) that provide the necessary primitive operations (QR, Schur decomposition, triangular solves)
- SLICOT benchmark data and test cases are available for validation (test data is non-copyrightable)
- Academic literature (Bartels-Stewart, Golub & Van Loan, Laub) provides sufficient algorithmic detail for clean room implementation
- Target platforms have IEEE 754 compliant floating-point arithmetic

## Dependencies *(optional)*

- **faer**: Required for high-performance linear algebra operations (Schur decomposition, QR factorization, triangular system solves)
- **ndarray**: Required for array structure and Python interoperability via PyO3
- **LAPACK source code**: BSD-licensed reference implementation (dtrsyl.f, strsyl.f) - safe to consult during implementation
- **LAPACK documentation**: Algorithmic reference for Bartels-Stewart method and numerical stability considerations
- **SLICOT test suite**: Required for validation and accuracy benchmarking (test data only, not source code)
- **Academic literature**:
  - Golub & Van Loan "Matrix Computations" Chapter 7.6.3 (Algorithm 7.6.2)
  - Bartels & Stewart (1972) "Solution of the Matrix Equation AX + XB = C" (CACM 15(9):820-826)
  - TOMS Algorithm 432 (for test case generation and validation only, not implementation reference)

## Out of Scope *(optional)*

- Sparse matrix support (initial implementation focuses on dense matrices)
- Generalized Sylvester equations (AXB + CXD = E) - deferred to future enhancement
- Complex-valued matrices (initial implementation is real-valued only)
- Parallel or distributed computation for very large problems
- Iterative solvers for large-scale problems (Krylov methods)
- Automatic preconditioning or scaling
- Sylvester observers and Lyapunov equations (separate features)
- GUI or visualization tools
- Automatic optimization of workspace sizes
- Support for structured matrices (Toeplitz, Hankel, etc.) beyond standard dense format
