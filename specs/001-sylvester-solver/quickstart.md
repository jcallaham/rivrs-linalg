# Quickstart Guide: Sylvester Equation Solver

**Feature**: 001-sylvester-solver
**Audience**: Users and developers getting started with the rivrs-linalg Sylvester solver

---

## What is a Sylvester Equation?

A **Sylvester equation** is a matrix equation of the form:

**Continuous-time**: AX + XB = C
**Discrete-time**: AXB + X = C

Where A, B, and C are known matrices, and X is the unknown matrix to solve for.

**Applications**:
- Observer design in control systems
- Model reduction and balancing
- Controller synthesis (LQR, H-infinity)
- Signal processing and filtering
- Numerical linear algebra algorithms

---

## Installation

Add rivrs-linalg to your `Cargo.toml`:

```toml
[dependencies]
rivrs_linalg = "0.1"  # Check for latest version
faer = "0.19"
```

---

## Basic Usage

### Continuous-Time Equation (AX + XB = C)

```rust
use rivrs_linalg::sylvester::solve_continuous;
use faer::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Define matrices
    let a = mat![
        [1.0, 0.5],
        [0.0, -1.0]
    ];

    let b = mat![
        [2.0, 0.0],
        [0.0, 3.0]
    ];

    let c = mat![
        [1.0, 2.0],
        [3.0, 4.0]
    ];

    // Solve AX + XB = C
    let result = solve_continuous(&a, &b, &c)?;

    // Apply scale factor (important!)
    let x = &result.solution / result.scale;

    println!("Solution X:\n{}", x);
    println!("Residual norm: {:.2e}", result.residual_norm);

    Ok(())
}
```

**Key Points**:
- Import `solve_continuous` function
- Create matrices using `faer::mat!` macro or `Mat::from_fn`
- Call solver with references (`&a`, `&b`, `&c`)
- **Always divide by `scale` to get the true solution**
- Check `residual_norm` to assess accuracy

---

### Discrete-Time Equation (AXB + X = C)

```rust
use rivrs_linalg::sylvester::{solve_discrete, Sign};
use faer::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Matrices for discrete-time system
    let a = mat![
        [0.5, 0.0],
        [0.0, 0.8]
    ];

    let b = mat![
        [0.6, 0.0],
        [0.0, 0.9]
    ];

    let c = mat![
        [1.0, 0.0],
        [0.0, 1.0]
    ];

    // Solve AXB + X = C
    let result = solve_discrete(&a, &b, &c, Sign::Plus)?;

    // Apply scale factor
    let x = &result.solution / result.scale;

    println!("Solution X:\n{}", x);

    Ok(())
}
```

**Discrete-Time Notes**:
- Use `Sign::Plus` for AXB + X = C
- Use `Sign::Minus` for AXB - X = C
- Common in digital control systems and discrete-time filters

---

## Error Handling

The solver returns `Result<SylvesterSolution<T>, SylvesterError>`. Always handle errors:

```rust
use rivrs_linalg::sylvester::{solve_continuous, SylvesterError};

match solve_continuous(&a, &b, &c) {
    Ok(result) => {
        // Check for warnings
        if result.near_singular {
            eprintln!("⚠ Warning: Problem is ill-conditioned");
            eprintln!("  Eigenvalues of A and -B are nearly overlapping");
        }

        let x = &result.solution / result.scale;
        println!("✓ Solution computed successfully");
        println!("  Residual: {:.2e}", result.residual_norm);
    }

    Err(SylvesterError::DimensionMismatch { expected, got, context }) => {
        eprintln!("✗ Dimension error:");
        eprintln!("  Expected: {:?}", expected);
        eprintln!("  Got: {:?}", got);
        eprintln!("  {}", context);
    }

    Err(SylvesterError::CommonEigenvalues { separation, threshold }) => {
        eprintln!("✗ Eigenvalues too close:");
        eprintln!("  Separation: {:.2e} < threshold: {:.2e}", separation, threshold);
        eprintln!("  Consider preconditioning or regularization");
    }

    Err(SylvesterError::InvalidInput { reason }) => {
        eprintln!("✗ Invalid input: {}", reason);
    }

    Err(e) => {
        eprintln!("✗ Solver failed: {}", e);
    }
}
```

---

## Understanding the Result

The `SylvesterSolution` struct contains:

```rust
pub struct SylvesterSolution<T> {
    pub solution: Mat<T>,      // Computed solution (may be scaled)
    pub scale: T,              // Scale factor (0 < scale ≤ 1)
    pub residual_norm: T,      // Accuracy measure
    pub near_singular: bool,   // Ill-conditioning warning
}
```

### Scale Factor

**Why does it exist?**
- Prevents numerical overflow during computation
- LAPACK and SLICOT use the same mechanism
- `scale` is typically 1.0 for well-conditioned problems
- `scale < 1` indicates overflow prevention was triggered

**How to use it:**
```rust
let result = solve_continuous(&a, &b, &c)?;

// CORRECT: Apply scaling
let x_true = &result.solution / result.scale;

// INCORRECT: Forgetting to scale
let x_wrong = result.solution;  // This is X * scale, not X!
```

### Residual Norm

The residual measures how well the solution satisfies the equation:

```rust
if result.residual_norm < 1e-10 {
    println!("✓ High accuracy solution");
} else if result.residual_norm < 1e-6 {
    println!("⚠ Moderate accuracy");
} else {
    println!("✗ Low accuracy - problem may be ill-conditioned");
}
```

---

## Creating Matrices

### Using mat! Macro (Small Matrices)

```rust
use faer::mat;

let a = mat![
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0]
];
```

### Using Mat::from_fn (Large Matrices)

```rust
use faer::Mat;

// Create 100×100 random matrix
let a = Mat::from_fn(100, 100, |i, j| {
    if i == j {
        1.0  // Diagonal
    } else {
        0.0  // Off-diagonal
    }
});
```

### Using Mat::zeros / Mat::identity

```rust
let zeros = Mat::zeros(10, 10);
let identity = Mat::identity(10, 10);
```

### From Existing Data

```rust
// From flat array (column-major)
let data: Vec<f64> = vec![1.0, 2.0, 3.0, 4.0];
let a = Mat::from_fn(2, 2, |i, j| data[j * 2 + i]);
```

---

## Validating Solutions

Always verify your solution for critical applications:

```rust
use rivrs_linalg::sylvester::{solve_continuous, compute_residual, EquationType};

let result = solve_continuous(&a, &b, &c)?;
let x = &result.solution / result.scale;

// Compute residual independently
let residual = compute_residual(&a, &b, &c, &x, EquationType::Continuous);

if residual < 1e-10 {
    println!("✓ Solution validated (residual: {:.2e})", residual);
} else {
    eprintln!("✗ Solution accuracy issue (residual: {:.2e})", residual);
}
```

---

## Common Patterns

### Solving Multiple Equations with Same A

When solving many equations AX₁ + X₁B₁ = C₁, AX₂ + X₂B₂ = C₂, ... with the same A matrix:

```rust
use rivrs_linalg::sylvester::{solve_continuous_schur};
use faer::linalg::evd::schur::real_schur;

// Pre-compute Schur decomposition of A (expensive, done once)
let (schur_a, u) = real_schur::compute(
    a.as_ref(),
    Default::default(),  // Use default parameters
    Par::None,           // No parallelism
    &mut stack,          // Workspace
)?;

// Solve multiple equations efficiently
for (b, c) in equations {
    let (schur_b, v) = real_schur::compute(b.as_ref(), ...)?;
    let result = solve_continuous_schur(&schur_a, &schur_b, &u, &v, &c)?;
    // ...
}
```

### Working with f32 (Single Precision)

```rust
use rivrs_linalg::sylvester::solve_continuous;
use faer::prelude::*;

let a: Mat<f32> = mat![
    [1.0f32, 0.5],
    [0.0, -1.0]
];

let b: Mat<f32> = mat![
    [2.0f32, 0.0],
    [0.0, 3.0]
];

let c: Mat<f32> = mat![
    [1.0f32, 2.0],
    [3.0, 4.0]
];

// Generic implementation works for f32
let result = solve_continuous(&a, &b, &c)?;
let x: Mat<f32> = &result.solution / result.scale;
```

---

## Troubleshooting

### Problem: DimensionMismatch Error

**Cause**: Matrix dimensions incompatible

**Solution**: Check that:
- A is square (n×n)
- B is square (m×m)
- C is rectangular (n×m)

```rust
println!("A: {}×{}", a.nrows(), a.ncols());
println!("B: {}×{}", b.nrows(), b.ncols());
println!("C: {}×{}", c.nrows(), c.ncols());
```

### Problem: CommonEigenvalues Error

**Cause**: A and -B have overlapping eigenvalues (ill-conditioned)

**Solutions**:
1. Check input conditioning: compute condition numbers
2. Scale matrices to similar magnitude
3. Add small regularization: A + εI instead of A
4. Use preconditioning techniques

### Problem: Large Residual Norm

**Cause**: Solution is inaccurate (poorly conditioned problem)

**Diagnosis**:
```rust
if result.near_singular {
    println!("Problem is ill-conditioned - expected");
}

// Check separation
let sep = /* compute eigenvalue separation */;
if sep < 1e-10 {
    println!("Eigenvalues nearly overlapping: sep = {:.2e}", sep);
}
```

**Solutions**:
- Use higher precision (f64 instead of f32)
- Improve problem conditioning
- Accept reduced accuracy or reformulate problem

### Problem: Residual Doesn't Match SylvesterSolution.residual_norm

**Cause**: Forgot to apply scale factor

**Solution**:
```rust
// WRONG
let x_wrong = result.solution;
let residual = compute_residual(&a, &b, &c, &x_wrong, ...);  // Won't match!

// CORRECT
let x_correct = &result.solution / result.scale;
let residual = compute_residual(&a, &b, &c, &x_correct, ...);  // Matches
```

---

## Performance Tips

### Matrix Size Impact

| Size (n×n) | Typical Time | Memory  |
|------------|--------------|---------|
| 10×10      | < 1ms        | < 1KB   |
| 100×100    | < 100ms      | < 1MB   |
| 500×500    | < 5s         | < 10MB  |
| 1000×1000  | < 30s        | < 100MB |

### Optimization Strategies

1. **Pre-compute Schur Forms**: For repeated solves with same matrices
2. **Use f32 for Speed**: If accuracy permits (8× faster, 4× less memory)
3. **Batch Processing**: Solve multiple independent equations in parallel
4. **Profile First**: Use `cargo-flamegraph` to find actual bottlenecks

---

## Next Steps

- **Read the API Documentation**: `cargo doc --open`
- **Explore Examples**: Check `examples/` directory for advanced usage
- **Run Benchmarks**: `cargo bench` to test performance on your hardware
- **Read Research Document**: `specs/001-sylvester-solver/research.md` for algorithm details

---

## Quick Reference

### Imports

```rust
use rivrs_linalg::sylvester::{
    solve_continuous,
    solve_discrete,
    solve_continuous_schur,
    compute_residual,
    Sign,
    EquationType,
    SylvesterSolution,
    SylvesterError,
};
use faer::prelude::*;
```

### Continuous-Time Template

```rust
let result = solve_continuous(&a, &b, &c)?;
let x = &result.solution / result.scale;
```

### Discrete-Time Template

```rust
let result = solve_discrete(&a, &b, &c, Sign::Plus)?;
let x = &result.solution / result.scale;
```

### Error Handling Template

```rust
match solve_continuous(&a, &b, &c) {
    Ok(result) => { /* use result */ }
    Err(e) => { /* handle error */ }
}
```

---

## Support

- **Documentation**: https://docs.rs/rivrs_linalg
- **Repository**: https://github.com/yourusername/rivrs_linalg
- **Issues**: https://github.com/yourusername/rivrs_linalg/issues
- **Academic References**: See `research.md` for algorithm citations

---

**Document Version**: 1.0.0
**Last Updated**: 2026-01-28
