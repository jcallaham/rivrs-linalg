use std::fmt;

/// Errors that can occur during sparse solver operations.
#[derive(Debug, Clone, PartialEq)]
pub enum SparseError {
    /// Matrix dimensions are incompatible.
    DimensionMismatch {
        expected: (usize, usize),
        got: (usize, usize),
        context: String,
    },

    /// A matrix that should be square is not.
    NotSquare { dims: (usize, usize) },

    /// Input contains NaN, Inf, or other invalid floating-point values.
    InvalidInput { reason: String },

    /// The matrix is structurally singular (symbolic analysis detected zero diagonal).
    StructurallySingular { column: usize },

    /// Numeric factorization encountered a zero or near-zero pivot.
    NumericalSingularity { pivot_index: usize, value: f64 },

    /// An ordering or analysis algorithm failed.
    AnalysisFailure { reason: String },

    /// An IO operation failed (file not found, permission denied, etc.).
    IoError { source: String, path: String },

    /// A file could not be parsed (malformed Matrix Market, invalid JSON, etc.).
    ParseError {
        reason: String,
        path: String,
        line: Option<usize>,
    },
}

impl fmt::Display for SparseError {
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
            Self::NotSquare { dims } => {
                write!(
                    f,
                    "Matrix must be square, but has dimensions {}x{}",
                    dims.0, dims.1
                )
            }
            Self::InvalidInput { reason } => {
                write!(f, "Invalid input: {}", reason)
            }
            Self::StructurallySingular { column } => {
                write!(f, "Matrix is structurally singular at column {}", column)
            }
            Self::NumericalSingularity { pivot_index, value } => {
                write!(
                    f,
                    "Numerical singularity at pivot index {} (value: {:.2e})",
                    pivot_index, value
                )
            }
            Self::AnalysisFailure { reason } => {
                write!(f, "Symbolic analysis failed: {}", reason)
            }
            Self::IoError { source, path } => {
                write!(f, "IO error for '{}': {}", path, source)
            }
            Self::ParseError { reason, path, line } => {
                if let Some(line) = line {
                    write!(f, "Parse error in '{}' at line {}: {}", path, line, reason)
                } else {
                    write!(f, "Parse error in '{}': {}", path, reason)
                }
            }
        }
    }
}

impl std::error::Error for SparseError {}

impl From<std::io::Error> for SparseError {
    fn from(e: std::io::Error) -> Self {
        Self::IoError {
            source: e.to_string(),
            path: String::new(),
        }
    }
}

impl From<serde_json::Error> for SparseError {
    fn from(e: serde_json::Error) -> Self {
        Self::ParseError {
            reason: e.to_string(),
            path: String::new(),
            line: None,
        }
    }
}
