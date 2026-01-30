use faer::prelude::*;

/// Result of solving a Sylvester equation.
///
/// The solution may be scaled to prevent overflow during computation.
/// The actual mathematical solution is `solution / scale`.
///
/// # References
///
/// - Bartels & Stewart (1972), "Solution of the Matrix Equation AX + XB = C",
///   CACM 15(9):820-826
/// - LAPACK dtrsyl uses the same SCALE factor convention
pub struct SylvesterSolution {
    /// Computed solution matrix X (may be scaled).
    ///
    /// Divide by `scale` to obtain the true mathematical solution.
    pub solution: Mat<f64>,

    /// Scale factor applied during computation (0 < scale <= 1).
    ///
    /// When `scale < 1`, overflow prevention was needed during the triangular
    /// solve phase. The actual solution is `solution / scale`.
    pub scale: f64,

    /// Residual norm measuring solution accuracy.
    ///
    /// For continuous-time: `||AX + XB - C||_F`
    /// For discrete-time: `||AXB + X - C||_F`
    pub residual_norm: f64,

    /// True if A and B have nearly overlapping eigenvalues.
    ///
    /// This indicates an ill-conditioned problem where the solution
    /// may be inaccurate.
    pub near_singular: bool,
}

/// Sign of the X term in discrete-time Sylvester equations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sign {
    /// Solve AXB + X = C
    Plus,
    /// Solve AXB - X = C
    Minus,
}

/// Transpose operation on a matrix in the triangular solver.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Transpose {
    /// No transpose: use A as-is
    NoTrans,
    /// Transpose: use A^T
    Trans,
}

/// Type of Sylvester equation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquationType {
    /// Continuous-time: AX + XB = C
    Continuous,
    /// Discrete-time: AXB + X = C
    Discrete,
}
