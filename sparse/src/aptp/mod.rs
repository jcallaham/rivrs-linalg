//! APTP (A Posteriori Threshold Pivoting) data structures.
//!
//! This module provides the data structures unique to indefinite LDL^T
//! factorization with a posteriori threshold pivoting:
//!
//! - [`PivotType`] — classification of column pivot decisions (1x1, 2x2, delayed)
//! - [`Block2x2`] — storage for 2x2 symmetric diagonal blocks
//! - [`MixedDiagonal`] — the D factor with mixed 1x1/2x2 blocks, solve, and inertia
//! - [`Inertia`] — eigenvalue sign classification (positive/negative/zero counts)
//! - [`perm_from_forward`] — construct faer `Perm<usize>` from a forward permutation array
//!
//! # References
//!
//! - Hogg, Duff & Lopez (2020), "A New Sparse LDL^T Solver Using A Posteriori
//!   Threshold Pivoting", SIAM J. Sci. Comput. 42(4)
//! - Bunch & Kaufman (1977), "Some Stable Methods for Calculating Inertia and
//!   Solving Symmetric Linear Systems", Math. Comp.

pub mod diagonal;
pub mod inertia;
pub mod matching;
pub mod ordering;
pub mod perm;
pub mod pivot;
pub mod symbolic;

pub use diagonal::MixedDiagonal;
pub use inertia::Inertia;
pub use matching::{Mc64Job, Mc64Result, mc64_matching};
pub use ordering::metis_ordering;
pub use perm::perm_from_forward;
pub use pivot::{Block2x2, PivotType};
pub use symbolic::{AptpSymbolic, SymbolicStatistics};
