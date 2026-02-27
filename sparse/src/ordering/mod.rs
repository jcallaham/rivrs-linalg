//! Fill-reducing orderings and matrix preprocessing.
//!
//! This module provides reusable ordering and scaling algorithms for sparse matrices:
//!
//! - [`metis_ordering()`] — METIS nested dissection ordering
//! - [`mc64_matching()`] — MC64 weighted bipartite matching and symmetric scaling
//! - [`match_order_metis()`] — Combined MC64 matching + METIS ordering pipeline
//! - [`perm_from_forward()`] — Construct faer `Perm<usize>` from a forward permutation array
//!
//! These algorithms are general-purpose and independent of any specific solver.
//! They can be composed with any sparse factorization or iterative solver.
//!
//! # Algorithm References
//!
//! - Karypis & Kumar (1998), "A Fast and High Quality Multilevel Scheme for
//!   Partitioning Irregular Graphs", SIAM J. Sci. Comput. 20(1)
//! - Duff & Koster (2001), "On Algorithms for Permuting Large Entries to the
//!   Diagonal of a Sparse Matrix", SIAM J. Matrix Anal. Appl. 22(4)
//! - Duff & Pralet (2005), "Strategies for Scaling and Pivoting for Sparse
//!   Symmetric Indefinite Problems", RAL Technical Report

/// MC64 weighted bipartite matching and symmetric scaling.
pub mod matching;
mod metis;
mod perm;

pub use matching::{Mc64Job, Mc64Result, mc64_matching};
pub use metis::{MatchOrderResult, match_order_metis, metis_ordering};
pub use perm::perm_from_forward;
