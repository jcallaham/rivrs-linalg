//! rivrs — Scientific Computing for Rust
//!
//! Umbrella crate re-exporting domain libraries.
//!
//! - `linalg` (default) — Numerical linear algebra ([`rivrs-linalg`])
//!
//! Future modules: `ode`, `optimize`, `signal`, `mesh`.

#[cfg(feature = "linalg")]
pub use rivrs_linalg as linalg;
