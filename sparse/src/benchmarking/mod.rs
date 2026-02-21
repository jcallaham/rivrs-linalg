//! Benchmarking infrastructure for SSIDS solver components.
//!
//! Provides a [`Benchmarkable`] trait for solver integration with Criterion,
//! data types for recording and comparing benchmark results, and utilities
//! for baseline management and regression detection.
//!
//! Solver phases are represented by [`crate::SolverPhase`], shared with the
//! testing module.
//!
//! Gated behind the `test-util` Cargo feature flag. Not included in production builds.

// RSS module is always available (no dependencies on other feature-gated modules)
pub mod rss;
pub use rss::{read_current_rss_kb, read_peak_rss_kb};

// Remaining benchmarking modules depend on the testing module (test-util)
#[cfg(feature = "test-util")]
pub mod baseline;
#[cfg(feature = "test-util")]
pub mod config;
#[cfg(feature = "test-util")]
pub mod report;
#[cfg(feature = "test-util")]
pub mod results;
#[cfg(feature = "test-util")]
pub mod traits;

#[cfg(feature = "test-util")]
pub use baseline::{Baseline, Improvement, Regression, RegressionReport, detect_regressions};
#[cfg(feature = "test-util")]
pub use config::BenchmarkConfig;
#[cfg(feature = "test-util")]
pub use report::{export_csv, export_json, generate_markdown_table};
#[cfg(feature = "test-util")]
pub use results::{BenchmarkResult, BenchmarkSuiteResult, SkippedBenchmark, collect_results};
#[cfg(feature = "test-util")]
pub use traits::{Benchmarkable, MockBenchmarkable};
