//! Sylvester equation solvers.
//!
//! This module provides solvers for continuous-time (AX + XB = C) and
//! discrete-time (AXB + X = C) Sylvester equations using the Bartels-Stewart
//! algorithm based on Schur decomposition.
//!
//! # References
//!
//! - Bartels & Stewart (1972), "Solution of the Matrix Equation AX + XB = C",
//!   CACM 15(9):820-826
//! - Golub & Van Loan (2013), "Matrix Computations" (4th Ed), Section 7.6.3

pub mod types;
pub mod validation;

pub use types::{EquationType, Sign, SylvesterSolution, Transpose};
