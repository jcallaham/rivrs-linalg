//! CSRRS - Control Systems Routines in Rust
//!
//! A scientific computing library implementing control systems algorithms,
//! starting with Sylvester equation solvers.
//!
//! # Clean Room Implementation
//!
//! This library is independently implemented from academic sources
//! (Bartels & Stewart 1972, Golub & Van Loan 2013, LAPACK BSD-licensed code)
//! to maintain MIT/Apache-2.0 licensing. SLICOT source code was never consulted.

pub mod error;
pub mod sylvester;
