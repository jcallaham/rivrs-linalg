# rivrs-linalg — Numerical Linear Algebra for Rivrs

A scientific computing library providing numerical linear algebra implementations
for the [Rivrs](https://github.com/jcallaham/rivrs) symbolic-numeric framework.
Currently focused on control systems algorithms similar in scope to
[SLICOT](http://www.slicot.org/), with plans to expand to sparse solvers and
other numerical methods.

## Current Features

### Sylvester Equation Solvers

| Equation | Function | Algorithm |
|---|---|---|
| Continuous: `AX + XB = C` | `solve_continuous` | Bartels-Stewart (Schur-Schur) |
| Discrete: `AXB + X = C` | `solve_discrete` | Modified Bartels-Stewart |

Both solvers include:

- Eigenvalue separation estimation for detecting near-singular problems
- Overflow prevention via scale factors
- Blocked Level-3 BLAS triangular solver for large matrices (continuous)
- Pre-computed Schur form API (`solve_continuous_schur`, `solve_discrete_schur`)

## Quick Start

Add to your `Cargo.toml`:

```toml
[dependencies]
rivrs-linalg = { git = "https://github.com/jcallaham/rivrs-linalg" }
faer = "0.22"
```

```rust
use rivrs_linalg::sylvester::solve_continuous;
use faer::mat;

let a = mat![[1.0, 0.0], [0.0, 2.0f64]];
let b = mat![[3.0, 0.0], [0.0, 4.0f64]];
let c = mat![[4.0, 5.0], [6.0, 12.0f64]];

let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
let x = &result.solution * (1.0 / result.scale);
// x ≈ [[1.0, 1.0], [2.0, 2.0]]
```

## Performance Comparison with SLICOT

Median wall-clock times (milliseconds) for solving Sylvester equations at
various matrix sizes N, measured on the same machine in a single run.

### Continuous: AX + XB = C

| N | SLICOT SB04MD (ms) | rivrs-linalg (ms) | Ratio |
|---:|---:|---:|---:|
| 10 | 0.031 | 0.025 | 0.8x |
| 50 | 2.59 | 1.18 | 0.5x |
| 100 | 11.2 | 7.00 | 0.6x |
| 200 | 139 | 50.5 | 0.4x |
| 500 | 1154 | 669 | 0.6x |

### Discrete: AXB + X = C

| N | SLICOT SB04QD (ms) | rivrs-linalg (ms) | Ratio |
|---:|---:|---:|---:|
| 10 | 0.041 | 0.040 | 1.0x |
| 50 | 1.47 | 2.21 | 1.5x |
| 100 | 11.1 | 18.0 | 1.6x |
| 200 | 83.3 | 191 | 2.3x |
| 500 | 1203 | 5311 | 4.4x |

Both implementations produce residuals well below 1e-10 at all sizes.

### Methodology and Caveats

This is an **apples-to-oranges comparison** — useful for confirming we are in
the right ballpark, not for drawing precise conclusions. Key differences:

- **SLICOT** uses the **Hessenberg-Schur** method (Golub, Nash & Van Loan 1979):
  reduces A to Hessenberg form and B to Schur form. Lower theoretical flop
  count but relies on Level-2 BLAS column solves.
- **rivrs-linalg** uses the **Bartels-Stewart** method (Bartels & Stewart 1972):
  reduces both A and B to Schur form. The continuous triangular solve uses a
  blocked Level-3 BLAS variant (Jonsson & Kagstrom 2002) for N > 64.
- **LAPACK version**: SLICOT is linked against **reference (Netlib) LAPACK/BLAS**,
  not an optimized implementation like OpenBLAS or MKL. rivrs-linalg uses **faer**
  for matrix operations and **nalgebra** for Schur decomposition, both of which
  include their own optimized kernels.
- **Discrete gap**: rivrs-linalg's discrete solver does not yet have a blocked
  triangular variant. The trilinear AXB term requires 6 matrix multiplications
  per block step (vs 2 for the continuous AX + XB case), making a blocked
  implementation more involved.
- **Random matrices** with diagonal shifts for conditioning (continuous) or
  scaling (discrete). Different RNGs between Fortran and Rust, but same
  statistical properties.

## Implementation Sources and Attribution

The Sylvester equation solvers are based on the following academic and open-source references:

### Primary Algorithm Sources
- **Bartels & Stewart (1972)**, "Solution of the Matrix Equation AX + XB = C",
  *Communications of the ACM* 15(9):820-826 — fundamental Schur-based algorithm
- **Golub & Van Loan (2013)**, *Matrix Computations* (4th Ed), Section 7.6.3 —
  algorithm description and numerical considerations
- **Jonsson & Kågström (2002)**, "Recursive blocked algorithms for solving
  triangular systems", *ACM TOMS* 28(4):416-435 — Level-3 BLAS blocked variant

### Reference Implementations
- **LAPACK** — implementation closely follows `dtrsyl`/`dtrsyl3`, including
  stability patterns, overflow prevention via scaling, and handling of quasi-triangular
  Schur form. Source code available at https://netlib.org/lapack/
- **SLICOT documentation and test suite** (SB04MD, SB04QD) — used only for understanding
  algorithm purpose, API, and expected numerical accuracy.

See [NOTICE](NOTICE) for complete licensing information.

## License

Licensed under Apache License, Version 2.0 ([LICENSE](LICENSE))
