# Data Model: Sylvester Equation Solver

**Feature**: 001-sylvester-solver
**Date**: 2026-01-28

## Overview

This document defines the data structures, types, and validation rules for the Sylvester equation solver implementation.

---

## Core Types

### 1. Solution Type

```rust
/// Result of solving a Sylvester equation
pub struct SylvesterSolution<T> {
    /// Computed solution matrix X
    pub solution: Mat<T>,

    /// Scale factor applied during computation (0 < scale ≤ 1)
    /// Actual mathematical solution is solution/scale
    /// SCALE < 1 indicates overflow prevention was needed
    pub scale: T,

    /// Residual norm ||AX + XB - C|| or ||AXB + X - C||
    /// Measures solution accuracy
    pub residual_norm: T,

    /// True if A and B have nearly overlapping eigenvalues
    /// Indicates ill-conditioned problem
    pub near_singular: bool,
}
```

**Invariants**:
- `0 < scale <= 1` always
- `solution.nrows() == A.nrows()` (n)
- `solution.ncols() == B.ncols()` (m)
- `residual_norm >= 0`

**Usage**:
```rust
let result = solve_continuous(&a, &b, &c)?;
let x_true = &result.solution / result.scale;  // Apply scaling
println!("Residual: {}", result.residual_norm);
if result.near_singular {
    eprintln!("Warning: Problem is ill-conditioned");
}
```

---

### 2. Error Type

```rust
/// Errors that can occur during Sylvester equation solving
#[derive(Debug, Clone, PartialEq)]
pub enum SylvesterError {
    /// Matrix dimensions incompatible
    DimensionMismatch {
        expected: (usize, usize),
        got: (usize, usize),
        context: String,  // e.g., "Matrix C must be n×m"
    },

    /// Matrix is not square when it should be
    NotSquare {
        matrix: char,  // 'A' or 'B'
        dims: (usize, usize),
    },

    /// Input contains NaN, Inf, or other invalid values
    InvalidInput {
        reason: String,
    },

    /// A and -B have common or nearly common eigenvalues
    CommonEigenvalues {
        separation: f64,  // Estimated sep(A, B)
        threshold: f64,   // Threshold used
    },

    /// Algorithm failed to converge
    ConvergenceFailure {
        algorithm: String,  // e.g., "Schur decomposition"
        iterations: usize,
    },

    /// Not in quasi-triangular Schur form
    NotQuasiTriangular {
        matrix: char,
        location: (usize, usize),  // Where violation detected
    },
}
```

**Error Messages**:
```rust
impl Display for SylvesterError {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        match self {
            Self::DimensionMismatch { expected, got, context } => {
                write!(f, "Dimension mismatch: expected {:?}, got {:?}. {}",
                       expected, got, context)
            }
            Self::CommonEigenvalues { separation, threshold } => {
                write!(f, "A and -B have nearly common eigenvalues \
                          (separation: {:.2e} < threshold: {:.2e}). \
                          Problem is ill-conditioned. Consider \
                          preconditioning or regularization.",
                       separation, threshold)
            }
            // ... other variants
        }
    }
}
```

---

### 3. Configuration Types

```rust
/// Equation sign for discrete-time Sylvester equations
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sign {
    /// Solve AXB + X = C
    Plus,

    /// Solve AXB - X = C
    Minus,
}

/// Transpose operation on matrix
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Transpose {
    /// No transpose
    NoTrans,

    /// Transpose
    Trans,
}
```

---

### 4. Input Types

All inputs are **faer matrix types**:

```rust
use faer::Mat;          // Owned matrix
use faer::MatRef;       // Immutable view
use faer::MatMut;       // Mutable view
```

**Matrix A**:
- Type: `MatRef<T>` where `T: RealField`
- Constraints: Square (n×n), finite values
- Role: Left coefficient matrix

**Matrix B**:
- Type: `MatRef<T>` where `T: RealField`
- Constraints: Square (m×m), finite values
- Role: Right coefficient matrix

**Matrix C**:
- Type: `MatRef<T>` where `T: RealField`
- Constraints: Rectangular (n×m), finite values
- Role: Right-hand side

---

## Validation Rules

### Dimension Compatibility

**Continuous-Time (AX + XB = C)**:
```
A.nrows() == A.ncols() == n  (square)
B.nrows() == B.ncols() == m  (square)
C.nrows() == n
C.ncols() == m
```

**Discrete-Time (AXB + X = C)**:
```
A.nrows() == A.ncols() == n  (square)
B.nrows() == B.ncols() == m  (square)
C.nrows() == n
C.ncols() == m
```

**Validation Function**:
```rust
fn validate_dimensions<T>(
    a: MatRef<T>,
    b: MatRef<T>,
    c: MatRef<T>,
) -> Result<(), SylvesterError> {
    // Check A is square
    if a.nrows() != a.ncols() {
        return Err(SylvesterError::NotSquare {
            matrix: 'A',
            dims: (a.nrows(), a.ncols()),
        });
    }

    // Check B is square
    if b.nrows() != b.ncols() {
        return Err(SylvesterError::NotSquare {
            matrix: 'B',
            dims: (b.nrows(), b.ncols()),
        });
    }

    // Check C dimensions
    if c.nrows() != a.nrows() || c.ncols() != b.ncols() {
        return Err(SylvesterError::DimensionMismatch {
            expected: (a.nrows(), b.ncols()),
            got: (c.nrows(), c.ncols()),
            context: "Matrix C must be n×m where A is n×n and B is m×m".into(),
        });
    }

    Ok(())
}
```

### Finite Value Check

```rust
fn validate_finite<T: RealField>(
    matrix: MatRef<T>,
    name: char,
) -> Result<(), SylvesterError> {
    for i in 0..matrix.nrows() {
        for j in 0..matrix.ncols() {
            let val = matrix[(i, j)];
            if !val.is_finite() {
                return Err(SylvesterError::InvalidInput {
                    reason: format!(
                        "Matrix {} contains non-finite value at ({}, {}): {:?}",
                        name, i, j, val
                    ),
                });
            }
        }
    }
    Ok(())
}
```

### Quasi-Triangular Validation (for Schur forms)

```rust
fn validate_quasi_triangular<T: RealField>(
    matrix: MatRef<T>,
    name: char,
) -> Result<(), SylvesterError> {
    let n = matrix.nrows();
    let tol = T::from_f64(1e-14).unwrap();

    // Check strictly lower triangle is zero (except 2×2 blocks)
    let mut i = 0;
    while i < n {
        // Check if 2×2 block starts at i
        let is_2x2_block = if i < n - 1 {
            matrix[(i + 1, i)].abs() > tol
        } else {
            false
        };

        if is_2x2_block {
            // Check 2×2 block structure: [[a, b], [c, a]]
            let a11 = matrix[(i, i)];
            let a22 = matrix[(i + 1, i + 1)];
            let a12 = matrix[(i, i + 1)];
            let a21 = matrix[(i + 1, i)];

            if (a11 - a22).abs() > tol {
                return Err(SylvesterError::NotQuasiTriangular {
                    matrix: name,
                    location: (i, i),
                });
            }

            if (a12 * a21).abs() < tol || (a12 * a21) > T::zero() {
                return Err(SylvesterError::NotQuasiTriangular {
                    matrix: name,
                    location: (i, i),
                });
            }

            // Check below 2×2 block
            for k in i+2..n {
                if matrix[(k, i)].abs() > tol || matrix[(k, i+1)].abs() > tol {
                    return Err(SylvesterError::NotQuasiTriangular {
                        matrix: name,
                        location: (k, i),
                    });
                }
            }

            i += 2;  // Skip next row
        } else {
            // 1×1 block - check below diagonal
            for k in i+1..n {
                if matrix[(k, i)].abs() > tol {
                    return Err(SylvesterError::NotQuasiTriangular {
                        matrix: name,
                        location: (k, i),
                    });
                }
            }

            i += 1;
        }
    }

    Ok(())
}
```

---

## State Transitions

### Solver State Machine

```
┌─────────────┐
│   INPUT     │  User provides A, B, C
│ VALIDATION  │
└──────┬──────┘
       │ validate_dimensions, validate_finite
       │ ✓ Pass / ✗ Return SylvesterError
       ▼
┌─────────────┐
│   SCHUR     │  Compute A = U₁TU₁ᵀ, B = U₂SU₂ᵀ
│DECOMPOSITION│
└──────┬──────┘
       │ faer::linalg::evd::schur::real_schur
       │ ✓ Converged / ✗ Return ConvergenceFailure
       ▼
┌─────────────┐
│ TRANSFORM   │  F = U₁ᵀCU₂
│     RHS     │
└──────┬──────┘
       │ matmul operations
       ▼
┌─────────────┐
│  TRIANGULAR │  Solve TY + YS = F
│    SOLVE    │  (or Y + TYS = F for discrete)
└──────┬──────┘
       │ solve_triangular_sylvester
       │ Tracks SCALE, detects near-singular
       │ ✓ Solution found / ✗ CommonEigenvalues
       ▼
┌─────────────┐
│   BACK-     │  X = U₁YU₂ᵀ
│  TRANSFORM  │
└──────┬──────┘
       │ matmul operations
       ▼
┌─────────────┐
│  COMPUTE    │  residual = ||AX + XB - C||
│  RESIDUAL   │
└──────┬──────┘
       │ matmul, norm computation
       ▼
┌─────────────┐
│   RETURN    │  SylvesterSolution<T>
│   SOLUTION  │
└─────────────┘
```

---

## Memory Layout

### Stack Allocation Pattern

```rust
// Workspace requirements
fn sylvester_workspace<T: RealField>(n: usize, m: usize) -> StackReq {
    let schur_a = schur_workspace::<T>(n);
    let schur_b = schur_workspace::<T>(m);
    let tmp_nm = temp_mat_scratch::<T>(n, m);

    // Sum of all requirements
    schur_a.and(schur_b).and(tmp_nm)
}

// Algorithm with workspace
fn solve_continuous_impl<T: RealField>(
    a: MatRef<T>,
    b: MatRef<T>,
    c: MatRef<T>,
    stack: &mut MemStack,
) -> Result<SylvesterSolution<T>, SylvesterError> {
    let n = a.nrows();
    let m = b.nrows();

    // Allocate temporaries from stack
    let (mut f, stack) = temp_mat_zeroed::<T>(n, m, stack);
    let (mut y, stack) = temp_mat_zeroed::<T>(n, m, stack);

    // ... algorithm using f, y ...
}
```

### Column-Major Layout

faer uses **column-major** layout (Fortran-style):
- Element `(i, j)` at memory offset: `j * nrows + i`
- Columns are contiguous in memory
- Compatible with BLAS/LAPACK conventions

**Implication**: Matrix-matrix products, Schur decompositions are optimized for column access.

---

## Relationships Between Entities

```
┌────────────┐
│  Matrix A  │─────┐
│   (n×n)    │     │
└────────────┘     │
                   │ Coefficients
┌────────────┐     │ for equation
│  Matrix B  │─────┤
│   (m×m)    │     │
└────────────┘     │
                   ▼
┌────────────┐   ┌────────────────────┐
│  Matrix C  │──▶│ Sylvester Equation │
│   (n×m)    │   │   AX + XB = C      │
└────────────┘   │  or AXB + X = C    │
                 └──────────┬─────────┘
                            │
                            │ Solve via
                            │ Bartels-Stewart
                            ▼
                 ┌────────────────────┐
                 │ SylvesterSolution  │
                 │  - solution (X)    │
                 │  - scale           │
                 │  - residual_norm   │
                 │  - near_singular   │
                 └────────────────────┘
```

---

## Type Constraints

### Generic Bounds

```rust
pub fn solve_continuous<T>(
    a: MatRef<T>,
    b: MatRef<T>,
    c: MatRef<T>,
) -> Result<SylvesterSolution<T>, SylvesterError>
where
    T: RealField,  // Supports f32, f64
```

**RealField Requirements**:
- Implements field operations (+, -, *, /)
- Has ordering (PartialOrd)
- Has zero, one, epsilon
- Is finite-testable (is_finite, is_nan, is_infinite)
- Compatible with faer's linear algebra routines

**Supported Types**:
- ✅ `f64` (primary target)
- ✅ `f32` (secondary target)
- ❌ Complex types (c32, c64) - out of scope for initial implementation
- ❌ Fixed-point types - possible future extension

---

## Invariants and Contracts

### Pre-conditions (checked by validation)

1. **A is square**: `A.nrows() == A.ncols()`
2. **B is square**: `B.nrows() == B.ncols()`
3. **C dimensions match**: `C.nrows() == A.nrows() && C.ncols() == B.ncols()`
4. **All values finite**: No NaN, no Inf in A, B, C
5. **Eigenvalue separation**: λ(A) ∩ λ(-B) = ∅ (checked during solve, not validated upfront)

### Post-conditions (guaranteed by solver)

1. **Solution dimensions**: `X.nrows() == C.nrows() && X.ncols() == C.ncols()`
2. **Scale validity**: `0 < solution.scale <= 1`
3. **Residual computed**: `solution.residual_norm >= 0`
4. **Approximate solution**: `||AX + XB - C||_F ≈ solution.residual_norm`

### Invariants During Computation

1. **Schur forms preserve eigenvalues**: `λ(A) = λ(T)`, `λ(B) = λ(S)`
2. **Orthogonal transformations preserve norms**: `||U₁ᵀCU₂||_F = ||C||_F`
3. **SCALE monotonically decreases**: Each rescaling multiplies SCALE by factor ≤ 1
4. **Block structure preserved**: Schur form 2×2 blocks remain intact during triangular solve

---

## Example Usage

### Basic Solve

```rust
use rivrs_control::sylvester::solve_continuous;
use faer::prelude::*;

let a = mat![[1.0, 0.5], [0.0, -1.0]];
let b = mat![[2.0, 0.0], [0.0, 3.0]];
let c = mat![[1.0, 2.0], [3.0, 4.0]];

let result = solve_continuous(&a, &b, &c)?;
let x = &result.solution / result.scale;

println!("Solution:\n{}", x);
println!("Residual: {:.2e}", result.residual_norm);
```

### Error Handling

```rust
let result = solve_continuous(&a, &b, &c);

match result {
    Ok(solution) => {
        if solution.near_singular {
            eprintln!("Warning: Problem is ill-conditioned");
        }
        // Use solution
    }
    Err(SylvesterError::DimensionMismatch { expected, got, context }) => {
        eprintln!("Invalid dimensions: expected {:?}, got {:?}", expected, got);
    }
    Err(SylvesterError::CommonEigenvalues { separation, .. }) => {
        eprintln!("Eigenvalues too close (sep = {:.2e})", separation);
    }
    Err(e) => {
        eprintln!("Solver failed: {}", e);
    }
}
```

---

## Testing Checklist

### Type Safety Tests

- [ ] Generic implementation works with f64
- [ ] Generic implementation works with f32
- [ ] Compilation fails for non-RealField types
- [ ] Lifetime constraints prevent use-after-free

### Validation Tests

- [ ] DimensionMismatch caught for non-square A
- [ ] DimensionMismatch caught for non-square B
- [ ] DimensionMismatch caught for incompatible C
- [ ] InvalidInput caught for NaN in A, B, C
- [ ] InvalidInput caught for Inf in A, B, C

### Invariant Tests

- [ ] Solution dimensions always match expected
- [ ] Scale always in range (0, 1]
- [ ] Residual always non-negative
- [ ] Near-singular flag set when eigenvalues overlap

---

**Document Status**: ✅ Complete
**Next**: Generate API contracts (contracts/sylvester.yaml)
