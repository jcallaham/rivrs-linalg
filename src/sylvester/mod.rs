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

pub mod continuous;
pub mod discrete;
pub mod triangular;
pub mod triangular_discrete;
pub mod types;
pub mod utils;
pub mod validation;

pub use continuous::{solve_continuous, solve_continuous_schur};
pub use discrete::{solve_discrete, solve_discrete_schur};
pub use types::{EquationType, Sign, SylvesterSolution, Transpose};
pub use utils::compute_residual;
