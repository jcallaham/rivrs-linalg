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

pub mod error;
