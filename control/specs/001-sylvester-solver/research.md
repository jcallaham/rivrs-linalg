# Research: Sylvester Equation Solver Implementation

**Feature**: 001-sylvester-solver
**Date**: 2026-01-28
**Status**: Complete

## Overview

This document consolidates research findings for implementing Sylvester equation solvers (continuous-time and discrete-time) in Rust using the faer linear algebra library. All technical unknowns from the planning phase have been resolved.

---

## 1. Bartels-Stewart Algorithm

### Algorithm Description

The Bartels-Stewart algorithm solves the continuous-time Sylvester equation **AX + XB = C** by reducing it to a triangular system through Schur decomposition.

**Source**: Golub & Van Loan "Matrix Computations" (4th Ed), Section 7.6.3, Algorithm 7.6.2

### Three-Step Process

**Step 1: Compute Schur Decompositions**
```
A = U₁ T U₁ᵀ  (A is n×n, U₁ orthogonal, T upper quasi-triangular)
B = U₂ S U₂ᵀ  (B is m×m, U₂ orthogonal, S upper quasi-triangular)
```

**Step 2: Transform Right-Hand Side**
```
F = U₁ᵀ C U₂  (Transform C using Schur vectors)
```

**Step 3: Solve Triangular Sylvester System**
```
TY + YS = F  (T and S are quasi-triangular, solve via back-substitution)
```

**Step 4: Back-Transform Solution**
```
X = U₁ Y U₂ᵀ  (Recover original solution)
```

### Quasi-Triangular Structure

Real Schur form is "quasi-triangular":
- Upper triangular with 1×1 and 2×2 blocks on the diagonal
- 2×2 blocks represent complex conjugate eigenvalue pairs
- 2×2 blocks have structure: `[[a, b], [c, a]]` where `b·c < 0`

### Algorithm Adaptation

G&VL Algorithm 7.6.2 presents the equation as **FZ - ZG = C**. To adapt for **AX + XB = C**:
- Replace F with A, Z with X, G with -B
- Algorithm structure remains identical
- Triangular solve becomes: **TX - X(-B) = F** or **TX + XB = F**

### Reference Implementation

LAPACK dtrsyl.f implements the triangular solve (Step 3) with parameters:
- `TRANA, TRANB`: Handle transpose operations
- `ISGN`: Handle sign (±) in equation
- `SCALE`: Prevent overflow via automatic rescaling

---

## 2. faer Library Capabilities

### Matrix Types

**Core Types**:
```rust
Mat<T>           // Owned matrix (heap-allocated)
MatRef<'a, T>    // Immutable view (like &[T])
MatMut<'a, T>    // Mutable view (like &mut [T])
```

**Key Operations**:
- Create: `Mat::zeros(m, n)`, `Mat::identity(n, n)`, `Mat::from_fn(m, n, |i, j| ...)`
- Access: `matrix[(i, j)]` for element indexing
- Views: `.as_ref()`, `.as_mut()`, `.rb()`, `.rb_mut()` for reborrowing
- Layout: Column-major by default (Fortran-style)

### Schur Decomposition

**Module**: `faer::linalg::evd::schur`

**Key Functions**:
```rust
use faer::linalg::evd::schur::real_schur;

// Compute real Schur decomposition
let (schur_form, schur_vectors) = real_schur::compute(
    matrix.as_ref(),
    params,        // SchurParams for tuning
    parallelism,   // Par::Rayon or Par::None
    stack,         // Workspace
);
```

**Outputs**:
- `schur_form`: Quasi-triangular matrix (T or S in algorithm)
- `schur_vectors`: Orthogonal matrix (U₁ or U₂ in algorithm)

**Configuration**: `SchurParams` allows tuning QR algorithm parameters (shift counts, deflation windows)

### Triangular System Solvers

**Module**: `faer::linalg::triangular_solve`

**Available Functions**:
```rust
use faer::linalg::triangular_solve::*;

// Lower triangular: L*X = B
solve_lower_triangular_in_place_with_conj(
    triangular_lower: MatRef<T>,
    conj_lhs: Conj,
    rhs: MatMut<T>,  // Overwritten with solution
    par: Par
);

// Upper triangular: U*X = B
solve_upper_triangular_in_place_with_conj(
    triangular_upper: MatRef<T>,
    conj_lhs: Conj,
    rhs: MatMut<T>,  // Overwritten with solution
    par: Par
);
```

**Notes**:
- In-place solvers modify RHS directly
- Conjugation support via `Conj` enum
- Parallelism via `Par` parameter
- Equivalent to BLAS `trsm` operations

**Limitation**: These solve **matrix equations** (AX = B), not Sylvester equations. We'll need to implement custom triangular Sylvester solver.

### Matrix Multiplication

**Module**: `faer::linalg::matmul`

**Key Function**:
```rust
use faer::linalg::matmul::matmul_with_conj;

matmul_with_conj(
    dst: MatMut<T>,     // Output: C
    accum: Accum,       // Replace or Add
    lhs: MatRef<T>,     // A
    conj_lhs: Conj,
    rhs: MatRef<T>,     // B
    conj_rhs: Conj,
    alpha: T,           // Scale factor
    par: Par
);
```

**Accum Options**:
- `Accum::Replace`: Overwrites dst: `C = α·A·B`
- `Accum::Add`: Adds to dst: `C = C + α·A·B`

### Workspace Management

**Stack-Based Allocation** (preferred for performance):
```rust
use dyn_stack::{MemBuffer, MemStack, StackReq};

// 1. Query workspace requirements
fn my_algorithm_scratch<T>(m: usize, n: usize) -> StackReq {
    faer::linalg::temp_mat_scratch::<T>(m, n)
}

// 2. Allocate workspace
let scratch = my_algorithm_scratch::<f64>(100, 100);
let mut mem = MemBuffer::new(scratch);
let mut stack = MemStack::new(&mut mem);

// 3. Create temporary matrices
let (mut tmp, stack) = unsafe {
    faer::linalg::temp_mat_uninit(m, n, stack)
};
let mut tmp: MatMut<f64> = tmp.as_mat_mut();
```

**Pattern**: Functions compute workspace size, caller allocates, temporaries borrowed from stack.

### Generic Type System

**Traits**:
```rust
use faer::prelude::*;

// Real numbers (f32, f64)
trait RealField: ComplexField<Real = Self> + PartialOrd { ... }

// Complex and real numbers
trait ComplexField: Clone + PartialEq {
    type Real: RealField;
    const IS_REAL: bool;
}
```

**Usage**:
```rust
fn sylvester<T: RealField>(a: MatRef<T>, b: MatRef<T>, c: MatMut<T>) {
    // Works for f32 and f64 automatically
}
```

**Utilities**: `faer_traits::math_utils` provides `zero()`, `one()`, `from_f64()`, arithmetic ops

---

## 3. LAPACK dtrsyl Implementation Insights

### Function Parameters

**TRANA, TRANB**: Operation on A and B ('N' = no transpose, 'T' = transpose)
**ISGN**: Sign in equation (+1 for AX + XB = C, -1 for AX - XB = C)
**M, N**: Matrix dimensions (A is M×M, B is N×N, C is M×N)
**A, B**: Upper quasi-triangular Schur forms
**C**: Right-hand side (overwritten with solution)
**SCALE**: Output scale factor (0 < scale ≤ 1)
**INFO**: Return code (0 = success, 1 = near-singular, <0 = invalid argument)

### Overflow Prevention via SCALE

**Initialization**:
```fortran
SMLNUM = machine_safe_min * (M*N) / machine_precision
BIGNUM = 1 / SMLNUM
SMIN = MAX(SMLNUM, precision * ||A||_max, precision * ||B||_max)
```

**Local Scaling** (SCALOC):
- For each block solve, check if `|diagonal_element| < SMIN`
- If near-singular: perturb to SMIN, set INFO = 1
- If overflow risk: compute SCALOC < 1 to prevent overflow
- Rescale all computed blocks when SCALOC < 1

**Global Scaling**:
```fortran
SCALE = SCALE * SCALOC  // Accumulate all local scaling
```

**Result**: Returned solution is X/SCALE (user must apply final scaling)

### Four Solution Cases

Algorithm branches based on transpose flags (TRANA, TRANB):

1. **TRANA='N', TRANB='N'**: Solve AX + ISGN·XB = C (bottom-left to top-right ordering)
2. **TRANA='T', TRANB='N'**: Solve AᵀX + ISGN·XB = C (top-left to bottom-right)
3. **TRANA='T', TRANB='T'**: Solve AᵀX + ISGN·XBᵀ = C (top-right to bottom-left)
4. **TRANA='N', TRANB='T'**: Solve AX + ISGN·XBᵀ = C (bottom-right to top-left)

**Block Detection**: Scans off-diagonal to find 2×2 blocks (A(K+1,K) ≠ 0 indicates 2×2 block at K)

### Subsystem Solving

Four block combinations at each step:
1. **1×1 × 1×1**: Direct division with overflow check
2. **2×1 × 1×1**: Call DLALN2 (solves 2×1 systems with scaling)
3. **1×1 × 2×1**: Call DLALN2
4. **2×2 × 2×2**: Call DLASY2 (specialized 2×2 Sylvester solver)

### Workspace Requirements

**dtrsyl** (Level-2 BLAS):
- O(1) constant workspace (local arrays VEC, X, DUM)
- Suitable for embedded systems

**dtrsyl3** (Level-3 BLAS blocked):
- O(M+N) workspace (IWORK, SWORK arrays)
- 10-50× faster for large matrices (M, N > 100)
- Uses blocked DGEMM for off-diagonal updates

### Error Codes

**INFO = 0**: Success
**INFO = 1**: Near-singular (eigenvalues close), solution computed with perturbation
**INFO < 0**: Invalid argument (INFO = -k means k-th parameter invalid)

**Note**: SCALE < 1 is normal and does NOT cause INFO = 1

---

## 4. Discrete-Time Sylvester Equation

### Problem Statement

Solve **AXB + X = C** (or equivalently **X + AXB = C**)

### Uniqueness Condition

**Continuous-time** (AX + XB = C): A and -B have no common eigenvalues
**Discrete-time** (AXB + X = C): A and -B have no common **reciprocal** eigenvalues (λᵢ ≠ 1/μⱼ)

### G&VL Coverage

**Finding**: Golub & Van Loan Chapter 7 focuses on continuous-time equations. Discrete-time variant is not explicitly covered in Chapter 7.

### SLICOT Approach

**Separate Implementations**:
- **SB04MD**: Continuous-time (AX + XB = C)
- **SB04QD**: Discrete-time (X + AXB = C)
- **SB04PD**: Unified solver with DICO parameter ('C' or 'D')

**Algorithm** (from SB04QD documentation):
1. Transform A to upper Hessenberg form: H = UᵀAU
2. Transform Bᵀ to real Schur form: S = ZᵀBᵀZ
3. Transform RHS: F = UᵀCZ
4. Solve transformed system: Y + HYSᵀ = F (back-substitution)
5. Recover solution: X = UYZᵀ

**References**:
- Golub, Nash, Van Loan (1979) "A Hessenberg-Schur method for the problem AX + XB = C"
- Sima (1996) "Algorithms for Linear-Quadratic Optimization"

### SciPy and Other Libraries

**SciPy**: Only supports continuous-time (`solve_sylvester`), no discrete-time solver
**Julia MatrixEquations.jl**: Provides `sylvd()` for discrete-time using extended Bartels-Stewart
**Scilab**: Provides `sylv(A,B,C,'d')` for discrete-time
**Modelica/Wolfram**: Uses Hessenberg-Schur method

### Implementation Decision

**Recommendation**: Implement **separate discrete-time solver**, do NOT transform to continuous form.

**Rationale**:
1. ✅ Standard practice: SLICOT, Julia, Scilab all use separate solvers
2. ✅ Numerical stability: Avoids conditioning issues from transformation
3. ✅ Clean room compliance: Can reference academic sources directly
4. ✅ Performance: Direct solution more efficient
5. ❌ No standard transformation exists in literature
6. ❌ Transformation would require original research

**Algorithm Structure** (parallel to continuous-time):
```
1. Compute Schur decompositions: A = USUᵀ, B = VTVᵀ
2. Transform RHS: F = UᵀCV
3. Solve triangular discrete system: Y + SYT = F  (KEY DIFFERENCE)
4. Back-transform: X = UYVᵀ
```

**Back-Substitution Difference**:
- Continuous: TY + YS = F
- Discrete: Y + TYS = F (different recursion pattern)

### Condition Number Estimation

**Continuous**: sep(A,B) = min ||AX - XB||_F / ||X||_F
**Discrete**: sep(A,B) = min ||AXB + X||_F / ||X||_F

Need separate estimation for discrete-time warnings.

---

## 5. API Design Decisions

### Public API Functions

**Continuous-Time Solver**:
```rust
pub fn solve_continuous<T: RealField>(
    a: MatRef<T>,     // n×n square matrix
    b: MatRef<T>,     // m×m square matrix
    c: MatRef<T>,     // n×m right-hand side
) -> Result<SylvesterSolution<T>, SylvesterError>
```

**Discrete-Time Solver**:
```rust
pub fn solve_discrete<T: RealField>(
    a: MatRef<T>,
    b: MatRef<T>,
    c: MatRef<T>,
    sign: Sign,       // enum Sign { Plus, Minus } for ±X term
) -> Result<SylvesterSolution<T>, SylvesterError>
```

**Schur Form Solver** (advanced users with pre-computed Schur forms):
```rust
pub fn solve_continuous_schur<T: RealField>(
    schur_a: MatRef<T>,  // A in Schur form
    schur_b: MatRef<T>,  // B in Schur form
    u: MatRef<T>,        // Schur vectors for A
    v: MatRef<T>,        // Schur vectors for B
    c: MatRef<T>,
) -> Result<SylvesterSolution<T>, SylvesterError>
```

### Return Types

**Solution Structure**:
```rust
pub struct SylvesterSolution<T> {
    pub solution: Mat<T>,      // Computed solution X
    pub scale: T,              // Scale factor (0 < scale ≤ 1)
    pub residual_norm: T,      // ||AX + XB - C|| for validation
    pub near_singular: bool,   // Eigenvalues nearly overlapping
}
```

**Error Enumeration**:
```rust
pub enum SylvesterError {
    DimensionMismatch {
        expected: (usize, usize),
        got: (usize, usize)
    },
    NotSquare {
        matrix: char,
        dims: (usize, usize)
    },
    InvalidInput {
        reason: String  // NaN, Inf, etc.
    },
    CommonEigenvalues {
        separation: f64  // Estimated sep(A,B)
    },
    ConvergenceFailure {
        algorithm: String
    },
}
```

### Generic vs Concrete Types

**Decision**: Use generic implementations with `T: RealField` trait bound.

**Rationale**:
- Single implementation works for f32 and f64
- Follows faer's design patterns
- Enables future extension to fixed-point or custom types
- No code duplication

**Example**:
```rust
pub fn solve_continuous<T: RealField>(
    a: MatRef<T>,
    b: MatRef<T>,
    c: MatRef<T>
) -> Result<SylvesterSolution<T>, SylvesterError>
```

### In-Place vs Allocating

**Decision**: Provide allocating variants initially, add in-place variants later if needed.

**Rationale**:
- Allocating variants simpler to implement and use
- Matches LAPACK pattern (C matrix overwritten, but we return new Mat<T>)
- In-place optimization deferred to Priority 4 (performance)

**Future Enhancement**:
```rust
pub fn solve_continuous_inplace<T: RealField>(
    a: MatRef<T>,
    b: MatRef<T>,
    c: MatMut<T>,  // Overwritten with solution
) -> Result<SylvesterInfo<T>, SylvesterError>
```

### Workspace Management

**Decision**: Internal allocation for initial implementation.

**Rationale**:
- Simplifies API (users don't manage workspace)
- faer's stack allocation is internal anyway
- LAPACK-style workspace parameters complicate API

**Internal Pattern**:
```rust
fn solve_continuous_impl<T: RealField>(/* ... */) -> Result<...> {
    let scratch = compute_workspace_size(m, n);
    let mut mem = MemBuffer::new(scratch);
    let mut stack = MemStack::new(&mut mem);

    // Use stack for temporaries
    // ...
}
```

### Error Handling Strategy

**Panic vs Result**:
- Dimension mismatches → `Result::Err` (recoverable, user error)
- Non-square matrices → `Result::Err` (user error)
- NaN/Inf inputs → `Result::Err` (user error)
- Common eigenvalues → `Result::Err` with warning (ill-conditioned but detectable)
- Internal assertions (unreachable code paths) → `panic!` (library bug)

**Informative Errors**:
```rust
SylvesterError::CommonEigenvalues { separation: 1e-14 }
// Error message: "A and -B have nearly common eigenvalues (separation: 1e-14).
//                 Consider preconditioning or checking input conditioning."
```

### Discrete-Time Sign Parameter

**Decision**: Separate parameter for ±X term.

**Options Considered**:
1. Single function with `sign` parameter ✅ **CHOSEN**
2. Separate functions `solve_discrete_plus`, `solve_discrete_minus`
3. Implicit (always +X)

**Rationale**:
- LAPACK uses ISGN parameter (proven API pattern)
- Both forms common in control theory applications
- Single function reduces code duplication
- Clear parameter name avoids confusion

---

## 6. Test Case Strategy

### Unit Tests (Small Analytical Cases)

**2×2 Examples**:
```rust
// Known solution: A = [[1, 0], [0, -1]], B = [[2, 0], [0, 3]], C = eye(2)
// Solution: X = [[1/3, 0], [0, 1/2]]
```

**3×3 Examples**:
- Diagonal matrices (trivial to verify)
- Upper triangular with known eigenvalues

**Source**: Create by hand or use Golub & Van Loan examples if available

### Integration Tests (G&VL Textbook Examples)

**Location**: docs/gvl-ch7.md (scan for worked examples with numerical values)

**Approach**:
1. Extract example matrices from text
2. Verify our solution matches published solution
3. Check residual norm

### Validation Tests (SLICOT Benchmark Problems)

**Source**: slicot/examples/ directory

**Test Data Extraction**:
1. Identify relevant example files:
   - SB04MD examples (continuous-time)
   - SB04QD examples (discrete-time)
2. Extract test matrices (A, B, C) from example data files
3. Compare our solution to SLICOT's reference solution
4. Verify numerical accuracy within tolerance (e.g., 10⁻¹⁰ relative error)

**Important**: Use test DATA only, never read SLICOT source code.

### Edge Cases

**Singular Matrices**:
```rust
// A is singular (det(A) = 0)
let a = mat![[1.0, 2.0], [2.0, 4.0]];  // Rank-deficient
// Expect: SylvesterError::CommonEigenvalues or convergence failure
```

**Ill-Conditioned**:
```rust
// Near-singular matrix (condition number ~ 10^15)
let a = mat![[1.0, 1.0], [1.0, 1.0 + 1e-15]];
// Expect: Large residual or near-singular warning
```

**Dimension Mismatches**:
```rust
let a = Mat::zeros(3, 3);
let b = Mat::zeros(2, 2);
let c = Mat::zeros(3, 3);  // Wrong! Should be 3×2
// Expect: SylvesterError::DimensionMismatch
```

**NaN/Inf Inputs**:
```rust
let a = mat![[f64::NAN, 1.0], [0.0, 1.0]];
// Expect: SylvesterError::InvalidInput
```

**Empty Matrices**:
```rust
let a = Mat::zeros(0, 0);
// Expect: Valid (degenerate case) or early return
```

**Special Structure** (deferred testing):
- Zero RHS (C = 0)
- Identity matrices
- Sparse patterns (outside scope but note behavior)

### Coverage Target

**Goal**: 95% code coverage (per SC-007)

**Priority**:
1. All branching logic (transpose cases, block sizes)
2. Error paths (each error variant triggered)
3. Numerical edge cases (near-singular, overflow scenarios)
4. Both f32 and f64 precision

---

## 7. Academic References for Documentation

Every implemented algorithm must cite sources. Key references:

### Primary Algorithm Sources

1. **Bartels, R.H. and Stewart, G.W. (1972)**
   "Solution of the Matrix Equation AX + XB = C"
   *Communications of the ACM*, Vol. 15, No. 9, pp. 820-826
   DOI: 10.1145/361573.361582

2. **Golub, G.H. and Van Loan, C.F. (2013)**
   "Matrix Computations" (4th Edition)
   Johns Hopkins University Press
   Section 7.6.3: "Block Diagonalization" (Algorithm 7.6.2)

3. **Golub, G.H., Nash, S., and Van Loan, C.F. (1979)**
   "A Hessenberg-Schur method for the problem AX + XB = C"
   *IEEE Transactions on Automatic Control*, AC-24, pp. 909-913
   DOI: 10.1109/TAC.1979.1102170

4. **Sima, V. (1996)**
   "Algorithms for Linear-Quadratic Optimization"
   Marcel Dekker, Inc., New York
   ISBN: 0-8247-9612-4

### Numerical Stability and Conditioning

5. **Higham, N.J. (1993)**
   "Perturbation theory and backward error for AX - XB = C"
   *BIT Numerical Mathematics*, Vol. 33, pp. 124-136
   DOI: 10.1007/BF01990349

6. **Kågström, B. and Westin, L. (1989)**
   "Generalized Schur Methods with Condition Estimators for Solving the Generalized Sylvester Equation"
   *IEEE Transactions on Automatic Control*, AC-34, pp. 745-751

### Implementation References (BSD-Licensed)

7. **LAPACK dtrsyl.f**
   Reference implementation (BSD-3-Clause)
   Location: lapack/SRC/dtrsyl.f
   Contributors: Sven Hammarling (NAG), Anderson et al.

8. **LAPACK dtrsyl3.f**
   Blocked Level-3 BLAS variant (BSD-3-Clause)
   Location: lapack/SRC/dtrsyl3.f
   Based on: Quintana-Orti & Van De Geijn (2003), Schwarz & Kjelgaard Mikkelsen (2020)

### Cross-Reference to SLICOT

Document equivalent SLICOT routines for user migration:

- **Continuous-time**: rivrs-linalg `solve_continuous` ↔ SLICOT SB04MD/SB04ND
- **Discrete-time**: rivrs-linalg `solve_discrete` ↔ SLICOT SB04QD/SB04PD
- **Triangular solver**: rivrs-linalg `solve_triangular_sylvester` ↔ LAPACK DTRSYL

**Important**: Note that rivrs-linalg is independently implemented from academic sources, not derived from SLICOT code.

---

## 8. Implementation Roadmap

### Phase 1: Core Infrastructure (Priority 1)

**Tasks**:
1. Define error types (`error.rs`)
2. Implement input validation (dimension checks, NaN/Inf detection)
3. Set up test infrastructure and first unit tests

**Deliverables**:
- `src/error.rs` with `SylvesterError` enum
- `src/utils/validation.rs` with validation functions
- Basic test harness

### Phase 2: Triangular Sylvester Solver (Priority 1)

**Tasks**:
1. Implement 1×1 block solver (direct division with SCALE)
2. Implement 2×2 block solver (following LAPACK DLASY2 pattern)
3. Implement back-substitution loop (4 transpose cases)
4. Add overflow prevention (SCALE factor logic)

**Deliverables**:
- `src/sylvester/triangular.rs` with `solve_triangular_sylvester` function
- Unit tests for each block size
- Edge case tests (near-singular, overflow)

**Reference**: LAPACK dtrsyl.f (BSD-licensed)

### Phase 3: Continuous-Time High-Level Solver (Priority 1)

**Tasks**:
1. Integrate faer Schur decomposition
2. Implement Bartels-Stewart algorithm (Steps 1-4)
3. Compute residual norms for validation
4. Handle workspace allocation

**Deliverables**:
- `src/sylvester/continuous.rs` with `solve_continuous` function
- Integration tests with G&VL examples
- Validation tests with SLICOT data

**Reference**: G&VL Algorithm 7.6.2

### Phase 4: Discrete-Time Solver (Priority 2)

**Tasks**:
1. Implement discrete-time triangular solver (Y + SYT = F)
2. Implement high-level discrete solver
3. Add sign parameter (±X term)
4. Test against SLICOT SB04QD examples

**Deliverables**:
- `src/sylvester/discrete.rs` with `solve_discrete` function
- Discrete-time test suite

**Reference**: Golub, Nash, Van Loan (1979), SLICOT SB04QD documentation

### Phase 5: Condition Estimation (Priority 3)

**Tasks**:
1. Implement sep(A,B) estimation for continuous-time
2. Implement discrete-time separation estimation
3. Add warnings for near-singular cases
4. Document conditioning in error messages

**Deliverables**:
- `src/sylvester/condition.rs` with separation estimation
- Enhanced error reporting

**Reference**: G&VL Section 7.6.4, Higham (1993)

### Phase 6: Performance Optimization (Priority 4)

**Two-Phase Strategy**: After verifying correctness with unblocked algorithm (Phases 1-3), implement blocked Level-3 BLAS optimization for competitive performance.

**Tasks**:
1. Benchmark against SLICOT and SciPy (establish baseline)
2. Profile bottlenecks in unblocked implementation
3. Implement blocked algorithm (dtrsyl3 style, 10-50× faster for n > 500)
4. Optimize workspace usage for blocked variant
5. Parallelize via faer's Par parameter

**Deliverables**:
- Benchmark suite in `benches/`
- Performance report comparing to reference implementations
- Blocked algorithm maintaining numerical accuracy
- Performance competitive with SLICOT for all matrix sizes

### Phase 7: Documentation and Polish

**Tasks**:
1. Write comprehensive rustdoc comments
2. Add usage examples to documentation
3. Create quickstart guide
4. Document numerical properties and failure modes
5. Add migration guide from MATLAB/SLICOT

**Deliverables**:
- API documentation with examples
- User guide in docs/
- Academic attribution in all function comments

---

## 9. Open Questions / Future Research

### Resolved

All critical unknowns have been resolved:
- ✅ faer provides necessary linear algebra primitives
- ✅ LAPACK dtrsyl patterns understood and documented
- ✅ Discrete-time implementation approach decided
- ✅ API design finalized
- ✅ Test case strategy established

### Deferred to Future Work

1. **Complex-valued matrices**: Current scope is real-valued only (f64, f32). Complex support (c64, c32) deferred.

2. **Generalized Sylvester equations**: AXB + CXD = E form deferred to future enhancement.

3. **Sparse matrices**: Current implementation assumes dense matrices. Sparse solvers are out of scope.

4. **Iterative methods**: Direct method only. Krylov subspace methods for large-scale problems deferred.

5. **Level-3 BLAS blocked algorithm**: Two-phase approach: unblocked (dtrsyl style) for correctness verification, then blocked variant (dtrsyl3 style) for performance optimization (Priority 4).

6. **Python bindings**: Implement core Rust library first, add PyO3 bindings in separate feature.

---

## 10. Summary and Recommendations

### Key Decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Algorithm | Bartels-Stewart via Schur decomposition | Standard method, numerically stable, well-documented |
| Library | faer for linear algebra | Provides all necessary primitives, pure Rust, high performance |
| API Style | Generic over `T: RealField` | Single implementation for f32/f64, idiomatic Rust |
| Error Handling | `Result<T, E>` with informative errors | Recoverable errors, clear messages guide users |
| Workspace | Internal allocation | Simplifies API, can optimize later |
| Discrete-Time | Separate algorithm implementation | Standard practice, numerically stable, clean room compliant |
| Testing | Unit + Integration + Validation + Benchmarks | Comprehensive coverage, reference comparisons |

### Implementation Confidence

**High Confidence**:
- ✅ All necessary library features available in faer
- ✅ Algorithm well-documented in academic literature
- ✅ LAPACK provides BSD-licensed reference implementation
- ✅ SLICOT test cases available for validation
- ✅ No licensing conflicts with clean room approach

**Risks**:
- ⚠️ 2×2 block handling complexity (mitigate: extensive testing)
- ⚠️ Overflow prevention correctness (mitigate: follow LAPACK SCALE logic closely)
- ⚠️ Performance vs SLICOT (mitigate: defer optimization to P4, benchmark early)

### Next Steps

1. ✅ Research complete (this document)
2. ⏭️ Proceed to Phase 1: Data Model & API Contracts (create data-model.md, contracts/)
3. ⏭️ Generate implementation tasks (run `/speckit.tasks`)
4. ⏭️ Begin implementation following task dependencies

---

## References

### Academic Papers

- Bartels & Stewart (1972) - Original Sylvester solver algorithm
- Golub & Van Loan (2013) - Matrix Computations textbook, Algorithm 7.6.2
- Golub, Nash & Van Loan (1979) - Hessenberg-Schur method
- Higham (1993) - Perturbation theory and conditioning
- Sima (1996) - Discrete-time algorithms

### Implementation References

- LAPACK dtrsyl.f, dtrsyl3.f (BSD-3-Clause)
- SLICOT SB04MD, SB04QD documentation (HTML only, not source)
- faer-rs library source code (MIT/Apache-2.0)

### Validation Resources

- docs/gvl-ch7.md - Golub & Van Loan Chapter 7
- docs/432.f - TOMS Algorithm 432 (validation reference only)
- slicot/examples/ - Test case data (not source code)
- lapack/TESTING/ - LAPACK test suite (if needed)

---

**Research Status**: ✅ **COMPLETE**
**Blocking Issues**: None
**Ready for**: Phase 1 (Data Model & API Contracts)
