//! Benchmarking utilities.
//!
//! Provides RSS memory tracking for solver profiling. Gated behind the
//! `test-util` or `diagnostic` Cargo feature flag (except `rss`, which is
//! always available when the module is compiled).

pub mod rss;
pub use rss::{read_current_rss_kb, read_peak_rss_kb};
