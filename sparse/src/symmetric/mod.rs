//! Sparse symmetric indefinite solver.
//!
//! This module provides a multifrontal LDL^T factorization with A Posteriori
//! Threshold Pivoting (APTP) for symmetric indefinite sparse matrices:
//!
//! - [`SparseLDLT`] — user-facing solver with three-phase API (analyze → factor → solve)
//! - [`AptpSymbolic`] — symbolic analysis result (reusable across factorizations)
//! - [`AptpNumeric`] — numeric factorization result
//! - [`PivotType`] — classification of column pivot decisions (1x1, 2x2, delayed)
//! - [`Block2x2`] — storage for 2x2 symmetric diagonal blocks
//! - [`MixedDiagonal`] — the D factor with mixed 1x1/2x2 blocks, solve, and inertia
//! - [`Inertia`] — eigenvalue sign classification (positive/negative/zero counts)
//!
//! Ordering and preprocessing utilities are in the sibling [`ordering`](crate::ordering) module.
//!
//! # References
//!
//! - Hogg, Duff & Lopez (2020), "A New Sparse LDL^T Solver Using A Posteriori
//!   Threshold Pivoting", SIAM J. Sci. Comput. 42(4)
//! - Bunch & Kaufman (1977), "Some Stable Methods for Calculating Inertia and
//!   Solving Symmetric Linear Systems", Math. Comp.

pub(crate) mod amalgamation;
pub mod diagonal;
pub mod factor;
pub mod inertia;
pub mod numeric;
pub mod pivot;
pub mod solve;
pub mod solver;
pub mod symbolic;

pub use diagonal::{MixedDiagonal, PivotEntry, PivotIter};
pub use inertia::Inertia;
pub use pivot::{Block2x2, PivotType};
pub use symbolic::{AptpSymbolic, SymbolicStatistics};

pub use factor::{
    AptpFactorResult, AptpFactorization, AptpFallback, AptpOptions, AptpPivotRecord,
    AptpStatistics, FailedPivotMethod, aptp_factor, aptp_factor_in_place,
};

pub use numeric::{AptpNumeric, FactorizationStats, FrontFactors, PerSupernodeStats};
pub use solver::{AnalyzeOptions, FactorOptions, OrderingStrategy, SolverOptions, SparseLDLT};
