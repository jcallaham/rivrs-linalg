//! rivrs-linalg — Numerical Linear Algebra for Rivrs
//!
//! Re-exports domain crates as feature-gated modules.
//! Enable features to select which modules to compile:
//!
//! - `sparse` (default) — Sparse symmetric indefinite solver (SSIDS)
//!
//! # Example
//!
//! ```ignore
//! use rivrs_linalg::sparse::symmetric::SparseLDLT;
//! ```

#[cfg(feature = "sparse")]
pub use rivrs_sparse as sparse;
