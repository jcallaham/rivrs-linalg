use std::fmt;

/// Errors that can occur during Sylvester equation solving.
///
/// These cover input validation failures, numerical issues during computation,
/// and structural requirements for advanced APIs.
#[derive(Debug, Clone, PartialEq)]
pub enum SylvesterError {
    /// Matrix dimensions are incompatible for the equation.
    DimensionMismatch {
        expected: (usize, usize),
        got: (usize, usize),
        context: String,
    },

    /// A matrix that should be square is not.
    NotSquare { matrix: char, dims: (usize, usize) },

    /// Input contains NaN, Inf, or other invalid floating-point values.
    InvalidInput { reason: String },

    /// A and -B have common or nearly common eigenvalues, making the
    /// equation ill-conditioned or singular.
    CommonEigenvalues { separation: f64, threshold: f64 },

    /// An iterative algorithm (e.g. Schur decomposition) failed to converge.
    ConvergenceFailure { algorithm: String },

    /// A matrix is not in quasi-triangular (real Schur) form when it should be.
    NotQuasiTriangular {
        matrix: char,
        location: (usize, usize),
    },
}

impl fmt::Display for SylvesterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::DimensionMismatch {
                expected,
                got,
                context,
            } => {
                write!(
                    f,
                    "Dimension mismatch: expected {:?}, got {:?}. {}",
                    expected, got, context
                )
            }
            Self::NotSquare { matrix, dims } => {
                write!(
                    f,
                    "Matrix {} must be square, but has dimensions {}x{}",
                    matrix, dims.0, dims.1
                )
            }
            Self::InvalidInput { reason } => {
                write!(f, "Invalid input: {}", reason)
            }
            Self::CommonEigenvalues {
                separation,
                threshold,
            } => {
                write!(
                    f,
                    "A and -B have nearly common eigenvalues \
                     (separation: {:.2e} < threshold: {:.2e}). \
                     Problem is ill-conditioned. Consider preconditioning or regularization.",
                    separation, threshold
                )
            }
            Self::ConvergenceFailure { algorithm } => {
                write!(f, "Algorithm failed to converge: {}", algorithm)
            }
            Self::NotQuasiTriangular { matrix, location } => {
                write!(
                    f,
                    "Matrix {} is not in quasi-triangular (Schur) form at position ({}, {})",
                    matrix, location.0, location.1
                )
            }
        }
    }
}

impl std::error::Error for SylvesterError {}
