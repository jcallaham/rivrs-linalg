//! Benchmarking infrastructure for SSIDS solver components.
//!
//! Provides a [`Benchmarkable`] trait for solver integration with Criterion,
//! data types for recording and comparing benchmark results, and utilities
//! for baseline management and regression detection.
//!
//! Gated behind the `test-util` Cargo feature flag. Not included in production builds.

pub mod baseline;
pub mod config;
pub mod report;
pub mod results;
pub mod rss;
pub mod traits;

pub use baseline::{Baseline, Improvement, Regression, RegressionReport, detect_regressions};
pub use config::{BenchmarkConfig, BenchmarkPhase};
pub use report::{export_csv, export_json, generate_markdown_table};
pub use results::{BenchmarkResult, BenchmarkSuiteResult, SkippedBenchmark, collect_results};
pub use rss::read_peak_rss_kb;
pub use traits::{Benchmarkable, MockBenchmarkable};
