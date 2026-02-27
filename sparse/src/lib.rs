#![warn(missing_docs)]
//! rivrs-sparse — Sparse Linear Algebra for Rivrs
//!
//! A scientific computing library providing sparse linear algebra solvers
//! for the Rivrs symbolic-numeric framework.
//!
//! # Sparse Symmetric Indefinite Solver (SSIDS)
//!
//! The library provides a sparse symmetric indefinite direct solver
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
//!
//! # Quick Start
//!
//! ```
//! use faer::sparse::{SparseColMat, Triplet};
//! use faer::Col;
//! use rivrs_sparse::symmetric::{SparseLDLT, SolverOptions};
//!
//! // Symmetric indefinite 2x2: A = [[2, 1], [1, -1]]
//! let triplets = vec![
//!     Triplet::new(0, 0, 2.0),
//!     Triplet::new(1, 0, 1.0), Triplet::new(0, 1, 1.0),
//!     Triplet::new(1, 1, -1.0),
//! ];
//! let a = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
//! let b = Col::from_fn(2, |i| [3.0, 0.0][i]);
//!
//! let x = SparseLDLT::solve_full(&a, &b, &SolverOptions::default()).unwrap();
//! assert!((x[0] - 1.0).abs() < 1e-12);
//! assert!((x[1] - 1.0).abs() < 1e-12);
//! ```
//!
//! # Three-Phase API
//!
//! For repeated solves with the same sparsity pattern, use the three-phase API
//! to reuse the symbolic analysis:
//!
//! ```
//! use faer::sparse::{SparseColMat, Triplet};
//! use faer::{Col, Par};
//! use faer::dyn_stack::{MemBuffer, MemStack};
//! use rivrs_sparse::symmetric::{SparseLDLT, AnalyzeOptions, FactorOptions};
//!
//! // Build a symmetric matrix (full CSC — both triangles stored)
//! let triplets = vec![
//!     Triplet::new(0, 0, 4.0),
//!     Triplet::new(1, 0, 1.0), Triplet::new(0, 1, 1.0),
//!     Triplet::new(1, 1, 3.0),
//! ];
//! let a = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
//!
//! // Phase 1: Symbolic analysis (reusable for same sparsity pattern)
//! let mut solver = SparseLDLT::analyze_with_matrix(&a, &AnalyzeOptions::default())?;
//!
//! // Phase 2: Numeric factorization
//! solver.factor(&a, &FactorOptions::default())?;
//!
//! // Phase 3: Triangular solve
//! let b = Col::from_fn(2, |i| [5.0, 4.0][i]);
//! let req = solver.solve_scratch(1);
//! let mut mem = MemBuffer::new(req);
//! let x = solver.solve(&b, &mut MemStack::new(&mut mem), Par::Seq)?;
//! assert!((x[0] - 1.0).abs() < 1e-12);
//! assert!((x[1] - 1.0).abs() < 1e-12);
//! # Ok::<(), rivrs_sparse::error::SparseError>(())
//! ```
//!
//! # Ordering Strategies
//!
//! | Strategy | Best for | Notes |
//! |----------|----------|-------|
//! | [`MatchOrderMetis`](symmetric::OrderingStrategy::MatchOrderMetis) | Hard-indefinite (KKT, saddle-point) | Default. MC64 matching + METIS |
//! | [`Metis`](symmetric::OrderingStrategy::Metis) | Easy-indefinite (FEM, thermal) | Pure METIS nested dissection |
//! | [`Amd`](symmetric::OrderingStrategy::Amd) | Small matrices, unit tests | faer built-in AMD |
//!
//! # Feature Flags
//!
//! | Feature | Purpose |
//! |---------|---------|
//! | `diagnostic` | Per-supernode timing instrumentation. Zero overhead when disabled. |
//! | `test-util` | Test infrastructure: random matrix generators, property-based testing. |
//!
//! # Error Handling
//!
//! All fallible operations return [`Result<T, SparseError>`](error::SparseError).
//! Error variants are organized by solver phase:
//!
//! - **Analysis**: [`NotSquare`](error::SparseError::NotSquare),
//!   [`AnalysisFailure`](error::SparseError::AnalysisFailure)
//! - **Factorization**: [`NumericalSingularity`](error::SparseError::NumericalSingularity),
//!   [`StructurallySingular`](error::SparseError::StructurallySingular)
//! - **Solve**: [`SolveBeforeFactor`](error::SparseError::SolveBeforeFactor),
//!   [`DimensionMismatch`](error::SparseError::DimensionMismatch)
//! - **I/O**: [`IoError`](error::SparseError::IoError),
//!   [`ParseError`](error::SparseError::ParseError)
//!
//! # Matrix Storage Convention
//!
//! All input matrices must be stored as **full symmetric CSC** (both upper and
//! lower triangles). The [`io::mtx`] reader automatically mirrors entries from
//! Matrix Market files.
//!
//! # Parallelism
//!
//! Both factorization and solve support shared-memory parallelism via
//! `Par::Seq` (sequential) or `Par::rayon(n)` (parallel with `n` threads):
//!
//! ```no_run
//! # use faer::Par;
//! # use rivrs_sparse::symmetric::FactorOptions;
//! let opts = FactorOptions { par: Par::rayon(4), ..Default::default() };
//! ```
//!
//! Tree-level parallelism (independent subtrees via rayon) and intra-node
//! parallelism (TRSM/GEMM via faer `Par`) are both controlled by this setting.
//!
//! # Performance
//!
//! Factorization dominates total solve time for most matrices. On the
//! 65-matrix SuiteSparse benchmark suite, rivrs-sparse is competitive with
//! SPRAL (median 0.94x sequential, 0.89x at 8 threads). Tuning
//! [`FactorOptions::threshold`](symmetric::FactorOptions::threshold) (default
//! 0.01) trades off stability vs fill-in for specific problem classes.
//!
//! # Examples
//!
//! See the [`examples/`](https://github.com/jcallaham/rivrs-linalg/tree/main/sparse/examples)
//! directory for complete programs:
//!
//! - `basic_usage.rs` — Self-contained hello world
//! - `multiple_rhs.rs` — Solve for multiple right-hand sides, reusing the factorization
//! - `refactorization.rs` — Refactorize with different values on the same sparsity pattern
//! - `solve_timing.rs` — End-to-end solve timing on SuiteSparse matrices
//! - `profile_matrix.rs` — Per-supernode profiling with Chrome Trace export
//! - `parallel_scaling.rs` — Parallel speedup measurement

use std::fmt;

use serde::{Deserialize, Serialize};

pub mod error;
pub mod io;
pub mod ordering;
pub mod symmetric;
pub mod validate;

#[cfg(any(feature = "test-util", feature = "diagnostic"))]
pub mod benchmarking;
#[cfg(feature = "test-util")]
pub mod debug;
#[cfg(any(feature = "test-util", feature = "diagnostic"))]
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
    /// Symbolic analysis: ordering, elimination tree, supernode structure.
    Analyze,
    /// Numeric factorization: LDL^T with APTP pivoting.
    Factor,
    /// Triangular solve: forward/backward substitution.
    Solve,
    /// Full end-to-end pipeline (analyze + factor + solve).
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
