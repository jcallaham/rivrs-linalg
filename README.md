# CSRRS — Control Systems Routines in Rust

A scientific computing library implementing control systems algorithms, similar
in scope to [SLICOT](http://www.slicot.org/) but with permissive
(MIT/Apache-2.0) licensing.

This is a **clean room implementation** — algorithms are implemented from
academic papers and textbooks, not from SLICOT source code, to avoid GPL
contamination.

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
csrrs = { git = "https://github.com/..." }  # update with actual URL
faer = "0.22"
```

```rust
use csrrs::sylvester::solve_continuous;
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

| N | SLICOT SB04MD (ms) | CSRRS (ms) | Ratio |
|---:|---:|---:|---:|
| 10 | 0.031 | 0.025 | 0.8x |
| 50 | 2.59 | 1.18 | 0.5x |
| 100 | 11.2 | 7.00 | 0.6x |
| 200 | 139 | 50.5 | 0.4x |
| 500 | 1154 | 669 | 0.6x |

### Discrete: AXB + X = C

| N | SLICOT SB04QD (ms) | CSRRS (ms) | Ratio |
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
- **CSRRS** uses the **Bartels-Stewart** method (Bartels & Stewart 1972):
  reduces both A and B to Schur form. The continuous triangular solve uses a
  blocked Level-3 BLAS variant (Jonsson & Kagstrom 2002) for N > 64.
- **LAPACK version**: SLICOT is linked against **reference (Netlib) LAPACK/BLAS**,
  not an optimized implementation like OpenBLAS or MKL. CSRRS uses **faer** for
  matrix operations and **nalgebra** for Schur decomposition, both of which
  include their own optimized kernels.
- **Discrete gap**: CSRRS's discrete solver does not yet have a blocked
  triangular variant. The trilinear AXB term requires 6 matrix multiplications
  per block step (vs 2 for the continuous AX + XB case), making a blocked
  implementation more involved.
- **Random matrices** with diagonal shifts for conditioning (continuous) or
  scaling (discrete). Different RNGs between Fortran and Rust, but same
  statistical properties.

## Academic References

- Bartels & Stewart (1972), "Solution of the Matrix Equation AX + XB = C",
  CACM 15(9):820-826
- Golub & Van Loan (2013), *Matrix Computations* (4th Ed), Chapter 7
- Jonsson & Kagstrom (2002), "Recursive blocked algorithms for solving
  triangular systems", ACM TOMS 28(4):416-435
- LAPACK `dtrsyl`/`dtrsyl3` (BSD-3-Clause) — consulted for numerical stability
  patterns

## License

Licensed under either of

- Apache License, Version 2.0 ([LICENSE-APACHE](LICENSE-APACHE))
- MIT License ([LICENSE-MIT](LICENSE-MIT))

at your option.
