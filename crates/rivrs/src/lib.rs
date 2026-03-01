//! rivrs — Symbolic and Numeric Computing for Rust
//!
//! Currently only the `linalg` module is released publicly:
//!
//! - `linalg` (default) — Numerical linear algebra ([`rivrs-linalg`])
//!
//! Future modules: `ode`, `optimize`, `signal`, `mesh`, etc.

#[cfg(feature = "linalg")]
pub use rivrs_linalg as linalg;
