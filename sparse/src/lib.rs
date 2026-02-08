//! rivrs-sparse — Sparse Linear Algebra for Rivrs
//!
//! A scientific computing library providing sparse linear algebra solvers
//! for the Rivrs symbolic-numeric framework.
//!
//! # Sparse Symmetric Indefinite Solver (SSIDS)
//!
//! The library will provide a sparse symmetric indefinite direct solver
//! based on the A Posteriori Threshold Pivoting (APTP) algorithm:
//!
//! - **Symbolic analysis**: Ordering, elimination tree, symbolic factorization
//! - **Numeric factorization**: LDL^T with APTP pivoting
//! - **Triangular solve**: Forward/backward substitution
//!
//! # Clean Room Implementation
//!
//! This library is independently implemented from academic sources to maintain
//! permissive licensing. The primary references are:
//!
//! - Hogg, Duff, & Lopez (2020), "A New Sparse LDL^T Solver Using
//!   A Posteriori Threshold Pivoting", SIAM J. Sci. Comput.
//! - Duff & Reid (1983), "The multifrontal solution of indefinite sparse
//!   symmetric linear equations"
//! - Liu (1992), "The Multifrontal Method for Sparse Matrix Solution:
//!   Theory and Practice", SIAM Review
//! - SPRAL (BSD-3-Clause) for reference implementation patterns
//!
//! See [NOTICE](../NOTICE) for complete attribution and licensing information.

use std::fmt;

use serde::{Deserialize, Serialize};

pub mod error;
pub mod io;
pub mod validate;

#[cfg(feature = "test-util")]
pub mod benchmarking;
#[cfg(feature = "test-util")]
pub mod debug;
#[cfg(feature = "test-util")]
pub mod profiling;
#[cfg(feature = "test-util")]
pub mod testing;

/// Solver phases corresponding to the three-phase API: analyze → factorize → solve.
///
/// Variants are ordered to match the natural pipeline order. `Roundtrip` represents
/// the full end-to-end pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum SolverPhase {
    Analyze,
    Factor,
    Solve,
    Roundtrip,
}

impl fmt::Display for SolverPhase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Analyze => write!(f, "analyze"),
            Self::Factor => write!(f, "factor"),
            Self::Solve => write!(f, "solve"),
            Self::Roundtrip => write!(f, "roundtrip"),
        }
    }
}

impl SolverPhase {
    /// All individual component phases (excludes Roundtrip).
    pub const fn components() -> &'static [Self] {
        &[Self::Analyze, Self::Factor, Self::Solve]
    }

    /// All phases including Roundtrip.
    pub const fn all() -> &'static [Self] {
        &[Self::Analyze, Self::Factor, Self::Solve, Self::Roundtrip]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_display() {
        assert_eq!(SolverPhase::Analyze.to_string(), "analyze");
        assert_eq!(SolverPhase::Factor.to_string(), "factor");
        assert_eq!(SolverPhase::Solve.to_string(), "solve");
        assert_eq!(SolverPhase::Roundtrip.to_string(), "roundtrip");
    }

    #[test]
    fn phase_serde_roundtrip() {
        for phase in SolverPhase::all() {
            let json = serde_json::to_string(phase).unwrap();
            let back: SolverPhase = serde_json::from_str(&json).unwrap();
            assert_eq!(*phase, back);
        }
    }

    #[test]
    fn phase_ordering_matches_pipeline() {
        assert!(SolverPhase::Analyze < SolverPhase::Factor);
        assert!(SolverPhase::Factor < SolverPhase::Solve);
        assert!(SolverPhase::Solve < SolverPhase::Roundtrip);
    }
}
