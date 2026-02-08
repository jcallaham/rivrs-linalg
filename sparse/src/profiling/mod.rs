//! Profiling infrastructure for SSIDS solver components.
//!
//! Provides hierarchical timing instrumentation with thread-safe recording,
//! Chrome Trace format export, and memory (RSS) tracking.
//!
//! Gated behind the `test-util` Cargo feature flag. Not included in production builds.
//!
//! # Examples
//!
//! ```
//! use rivrs_sparse::profiling::ProfileSession;
//!
//! let session = ProfileSession::new();
//! {
//!     let _guard = session.enter_section("analyze");
//!     // ... do work ...
//! }
//! let finished = session.finish();
//! let report = finished.summary_report();
//! assert!(report.contains("analyze"));
//! ```

pub mod memory;
mod report;
mod section;
mod session;

pub use report::format_chrome_trace;
pub use section::{ProfileEvent, ProfileSection};
pub use session::{FinishedSession, ProfileSession, SectionGuard, flush_thread_events};
