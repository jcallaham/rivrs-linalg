//! rivrs-linalg — Numerical Linear Algebra for Rivrs
//!
//! A scientific computing library providing numerical linear algebra implementations
//! for the Rivrs symbolic-numeric framework. Currently focused on control systems
//! algorithms, with plans to expand to sparse solvers and other numerical methods.
//! Licensed under permissive MIT/Apache-2.0 terms.
//!
//! # Sylvester Equation Solvers
//!
//! The library currently provides solvers for the Sylvester matrix equation:
//!
//! - **Continuous-time**: `AX + XB = C` via [`sylvester::solve_continuous`]
//! - **Discrete-time**: `AXB + X = C` via [`sylvester::solve_discrete`]
//!
//! Both solvers use the Bartels-Stewart algorithm (Schur decomposition followed
//! by triangular back-substitution) with overflow prevention and condition
//! estimation for detecting near-singular problems.
//!
//! # Example
//!
//! ```
//! use rivrs_linalg::sylvester::solve_continuous;
//! use faer::mat;
//!
//! let a = mat![[1.0, 0.0], [0.0, 2.0f64]];
//! let b = mat![[3.0, 0.0], [0.0, 4.0f64]];
//! let c = mat![[4.0, 5.0], [6.0, 12.0f64]];
//!
//! let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
//! let x = &result.solution * (1.0 / result.scale);
//! // x ≈ [[1.0, 1.0], [2.0, 2.0]]
//! ```
//!
//! # Clean Room Implementation
//!
//! This library is independently implemented from academic sources to maintain
//! permissive licensing. The primary references are:
//!
//! - Bartels & Stewart (1972), "Solution of the Matrix Equation AX + XB = C",
//!   CACM 15(9):820-826
//! - Golub & Van Loan (2013), "Matrix Computations" (4th Ed), Section 7.6.3
//! - LAPACK dtrsyl/dtrsyl3 (BSD-3-Clause) for numerical stability patterns
//! - Jonsson & Kågström (2002), "Recursive blocked algorithms for solving
//!   triangular systems", ACM TOMS 28(4):416-435
//!
//! See [NOTICE](../NOTICE) for complete attribution and licensing information.

pub mod error;
pub mod sylvester;
