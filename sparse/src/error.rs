//! Error types for the sparse solver pipeline.
//!
//! Provides [`SparseError`], the unified error enum for all phases of the
//! sparse solver (analysis, factorization, solve, and I/O). Implements
//! [`std::error::Error`] and [`Display`](std::fmt::Display) for integration
//! with standard Rust error handling.
//!
//! `PartialEq` is intentionally not derived because the `NumericalSingularity`
//! variant contains `f64`, where `NaN != NaN` would cause subtle comparison
//! bugs. Use `matches!()` for pattern-matching assertions in tests.

use std::fmt;

/// Errors that can occur during sparse solver operations.
///
/// Note: `PartialEq` is intentionally not derived because the `NumericalSingularity`
/// variant contains `f64`, where `NaN != NaN` would cause subtle comparison bugs.
/// Use `matches!()` for pattern-matching assertions in tests.
#[derive(Debug, Clone)]
pub enum SparseError {
    /// Matrix dimensions are incompatible.
    DimensionMismatch {
        /// Expected dimensions (rows, cols).
        expected: (usize, usize),
        /// Actual dimensions (rows, cols).
        got: (usize, usize),
        /// Description of where the mismatch occurred.
        context: String,
    },

    /// A matrix that should be square is not.
    NotSquare {
        /// Actual dimensions (rows, cols).
        dims: (usize, usize),
    },

    /// Input contains NaN, Inf, or other invalid floating-point values.
    InvalidInput {
        /// Description of the invalid input.
        reason: String,
    },

    /// The matrix is structurally singular (symbolic analysis detected zero diagonal).
    StructurallySingular {
        /// Zero-diagonal column index.
        column: usize,
    },

    /// Numeric factorization encountered a zero or near-zero pivot.
    NumericalSingularity {
        /// Index of the singular pivot.
        pivot_index: usize,
        /// The near-zero pivot value.
        value: f64,
    },

    /// An ordering or analysis algorithm failed.
    AnalysisFailure {
        /// Description of the failure.
        reason: String,
    },

    /// An IO operation failed (file not found, permission denied, etc.).
    IoError {
        /// Error description.
        source: String,
        /// File path that caused the error.
        path: String,
    },

    /// A file could not be parsed (malformed Matrix Market, invalid JSON, etc.).
    ParseError {
        /// Description of the parse failure.
        reason: String,
        /// File path that could not be parsed.
        path: String,
        /// Line number where parsing failed, if available.
        line: Option<usize>,
    },

    /// A named matrix was not found in the registry.
    MatrixNotFound {
        /// Name that was looked up.
        name: String,
    },

    /// Solve was called before factor() completed.
    SolveBeforeFactor {
        /// Description of the premature solve attempt.
        context: String,
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
            Self::MatrixNotFound { name } => {
                write!(f, "Matrix '{}' not found in registry", name)
            }
            Self::SolveBeforeFactor { context } => {
                write!(f, "Solve called before factor(): {}", context)
            }
        }
    }
}

impl std::error::Error for SparseError {}

// Note: blanket `From<std::io::Error>` and `From<serde_json::Error>` impls were
// intentionally removed. They produced errors with empty `path` fields, losing
// context. All call sites should construct `SparseError::IoError` or
// `SparseError::ParseError` explicitly with the relevant file path.
