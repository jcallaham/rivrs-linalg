//! Memory tracking via RSS snapshots.
//!
//! Records point-in-time RSS measurements at labeled checkpoints during
//! solver execution, computes deltas between consecutive snapshots,
//! and generates formatted reports.

use std::fmt::Write;
use std::time::Instant;

use crate::benchmarking::rss::{read_current_rss_kb, read_peak_rss_kb};

/// A point-in-time RSS measurement annotated with a label.
#[derive(Debug, Clone)]
pub struct MemorySnapshot {
    /// Developer-provided section name.
    pub label: String,
    /// When the snapshot was taken.
    pub timestamp: Instant,
    /// Current RSS in KB (None if unavailable).
    pub rss_kb: Option<u64>,
    /// Peak RSS (VmHWM) in KB (None if unavailable).
    pub peak_rss_kb: Option<u64>,
}

/// Computed difference between two consecutive snapshots.
#[derive(Debug, Clone)]
pub struct MemoryDelta {
    /// Label of the earlier snapshot.
    pub from_label: String,
    /// Label of the later snapshot.
    pub to_label: String,
    /// Change in current RSS (can be negative).
    pub rss_delta_kb: Option<i64>,
}

/// Aggregated view of a sequence of memory snapshots.
#[derive(Debug, Clone)]
pub struct MemoryReport {
    /// Ordered sequence of snapshots.
    pub snapshots: Vec<MemorySnapshot>,
    /// Maximum peak_rss_kb across all snapshots.
    pub peak_rss_kb: Option<u64>,
    /// Computed differences between consecutive snapshots.
    pub deltas: Vec<MemoryDelta>,
}

/// Tracks RSS at labeled points during solver execution.
pub struct MemoryTracker {
    snapshots: Vec<MemorySnapshot>,
}

impl MemoryTracker {
    /// Create a new, empty memory tracker.
    pub fn new() -> Self {
        Self {
            snapshots: Vec::new(),
        }
    }

    /// Record a snapshot of current and peak RSS with the given label.
    pub fn snapshot(&mut self, label: &str) {
        self.snapshots.push(MemorySnapshot {
            label: label.to_string(),
            timestamp: Instant::now(),
            rss_kb: read_current_rss_kb(),
            peak_rss_kb: read_peak_rss_kb(),
        });
    }

    /// Compute deltas and produce a memory report.
    pub fn report(&self) -> MemoryReport {
        let deltas = self.compute_deltas();
        let peak_rss_kb = self.snapshots.iter().filter_map(|s| s.peak_rss_kb).max();

        MemoryReport {
            snapshots: self.snapshots.clone(),
            peak_rss_kb,
            deltas,
        }
    }

    fn compute_deltas(&self) -> Vec<MemoryDelta> {
        self.snapshots
            .windows(2)
            .map(|pair| {
                let rss_delta_kb = match (pair[0].rss_kb, pair[1].rss_kb) {
                    (Some(a), Some(b)) => Some(b as i64 - a as i64),
                    _ => None,
                };
                MemoryDelta {
                    from_label: pair[0].label.clone(),
                    to_label: pair[1].label.clone(),
                    rss_delta_kb,
                }
            })
            .collect()
    }
}

impl Default for MemoryTracker {
    fn default() -> Self {
        Self::new()
    }
}

impl MemoryReport {
    /// Get the peak RSS across all snapshots.
    pub fn peak_rss_kb(&self) -> Option<u64> {
        self.peak_rss_kb
    }

    /// Get the ordered snapshots.
    pub fn snapshots(&self) -> &[MemorySnapshot] {
        &self.snapshots
    }

    /// Get the computed deltas.
    pub fn deltas(&self) -> &[MemoryDelta] {
        &self.deltas
    }

    /// Format as a text table for terminal display.
    pub fn display_report(&self) -> String {
        let mut out = String::new();

        writeln!(out, "Memory Report").unwrap();
        writeln!(out, "{:─<56}", "").unwrap();
        writeln!(
            out,
            "  {:<22} {:>10} {:>10} {:>10}",
            "Label", "RSS (KB)", "Peak (KB)", "Delta (KB)"
        )
        .unwrap();

        for (i, snapshot) in self.snapshots.iter().enumerate() {
            let rss = format_opt_u64(snapshot.rss_kb);
            let peak = format_opt_u64(snapshot.peak_rss_kb);
            let delta = if i == 0 {
                "\u{2014}".to_string() // em dash
            } else {
                match &self.deltas.get(i - 1) {
                    Some(d) => format_opt_delta(d.rss_delta_kb),
                    None => "N/A".to_string(),
                }
            };

            writeln!(
                out,
                "  {:<22} {:>10} {:>10} {:>10}",
                snapshot.label, rss, peak, delta
            )
            .unwrap();
        }

        writeln!(out, "{:─<56}", "").unwrap();

        match self.peak_rss_kb {
            Some(peak) => {
                let mb = peak as f64 / 1024.0;
                writeln!(out, "  Peak RSS: {peak} KB ({mb:.1} MB)").unwrap();
            }
            None => {
                writeln!(out, "  Peak RSS: N/A").unwrap();
            }
        }

        out
    }
}

fn format_opt_u64(val: Option<u64>) -> String {
    match val {
        Some(v) => format_with_commas(v),
        None => "N/A".to_string(),
    }
}

fn format_opt_delta(val: Option<i64>) -> String {
    match val {
        Some(v) if v >= 0 => format!("+{}", format_with_commas(v as u64)),
        Some(v) => format!("-{}", format_with_commas((-v) as u64)),
        None => "N/A".to_string(),
    }
}

fn format_with_commas(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::new();
    for (i, c) in s.chars().rev().enumerate() {
        if i > 0 && i % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }
    result.chars().rev().collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn memory_snapshot_construction() {
        let snap = MemorySnapshot {
            label: "test".to_string(),
            timestamp: Instant::now(),
            rss_kb: Some(1024),
            peak_rss_kb: Some(2048),
        };
        assert_eq!(snap.label, "test");
        assert_eq!(snap.rss_kb, Some(1024));
        assert_eq!(snap.peak_rss_kb, Some(2048));
    }

    #[test]
    fn memory_snapshot_none_values() {
        let snap = MemorySnapshot {
            label: "no_rss".to_string(),
            timestamp: Instant::now(),
            rss_kb: None,
            peak_rss_kb: None,
        };
        assert_eq!(snap.rss_kb, None);
        assert_eq!(snap.peak_rss_kb, None);
    }

    #[test]
    fn memory_delta_construction() {
        let delta = MemoryDelta {
            from_label: "start".to_string(),
            to_label: "end".to_string(),
            rss_delta_kb: Some(1000),
        };
        assert_eq!(delta.from_label, "start");
        assert_eq!(delta.to_label, "end");
        assert_eq!(delta.rss_delta_kb, Some(1000));
    }

    #[test]
    fn memory_delta_none() {
        let delta = MemoryDelta {
            from_label: "a".to_string(),
            to_label: "b".to_string(),
            rss_delta_kb: None,
        };
        assert_eq!(delta.rss_delta_kb, None);
    }

    #[test]
    fn new_tracker_starts_empty() {
        let tracker = MemoryTracker::new();
        let report = tracker.report();
        assert!(report.snapshots().is_empty());
        assert!(report.deltas().is_empty());
        assert!(report.peak_rss_kb().is_none());
    }

    #[test]
    fn snapshot_adds_entries() {
        let mut tracker = MemoryTracker::new();
        tracker.snapshot("first");
        tracker.snapshot("second");
        let report = tracker.report();
        assert_eq!(report.snapshots().len(), 2);
        assert_eq!(report.snapshots()[0].label, "first");
        assert_eq!(report.snapshots()[1].label, "second");
    }

    #[test]
    fn snapshots_in_correct_order() {
        let mut tracker = MemoryTracker::new();
        tracker.snapshot("a");
        tracker.snapshot("b");
        tracker.snapshot("c");
        let report = tracker.report();
        let labels: Vec<&str> = report
            .snapshots()
            .iter()
            .map(|s| s.label.as_str())
            .collect();
        assert_eq!(labels, vec!["a", "b", "c"]);
    }

    #[test]
    fn report_computes_deltas_between_consecutive() {
        let mut tracker = MemoryTracker::new();
        tracker.snapshot("start");
        // Allocate some memory to increase RSS
        let _data: Vec<u8> = vec![0; 1024 * 1024]; // 1MB
        tracker.snapshot("after_alloc");
        let report = tracker.report();

        assert_eq!(report.deltas().len(), 1);
        assert_eq!(report.deltas()[0].from_label, "start");
        assert_eq!(report.deltas()[0].to_label, "after_alloc");
    }

    #[test]
    #[cfg(target_os = "linux")]
    fn peak_rss_returns_maximum() {
        let mut tracker = MemoryTracker::new();
        tracker.snapshot("s1");
        tracker.snapshot("s2");
        let report = tracker.report();
        assert!(report.peak_rss_kb().is_some());
        let peak = report.peak_rss_kb().unwrap();
        for snap in report.snapshots() {
            if let Some(p) = snap.peak_rss_kb {
                assert!(peak >= p);
            }
        }
    }

    #[test]
    fn empty_tracker_produces_empty_report() {
        let tracker = MemoryTracker::new();
        let report = tracker.report();
        assert!(report.snapshots().is_empty());
        assert!(report.deltas().is_empty());
        assert!(report.peak_rss_kb().is_none());
    }

    #[test]
    fn display_report_contains_header() {
        let mut tracker = MemoryTracker::new();
        tracker.snapshot("test");
        let report = tracker.report();
        let display = report.display_report();
        assert!(display.contains("Memory Report"));
        assert!(display.contains("Label"));
        assert!(display.contains("RSS (KB)"));
        assert!(display.contains("Peak (KB)"));
        assert!(display.contains("Delta (KB)"));
    }

    #[test]
    fn display_report_contains_snapshot_labels() {
        let mut tracker = MemoryTracker::new();
        tracker.snapshot("baseline");
        tracker.snapshot("after_load");
        let report = tracker.report();
        let display = report.display_report();
        assert!(display.contains("baseline"));
        assert!(display.contains("after_load"));
    }

    #[test]
    fn display_report_shows_peak_footer() {
        let mut tracker = MemoryTracker::new();
        tracker.snapshot("test");
        let report = tracker.report();
        let display = report.display_report();
        assert!(display.contains("Peak RSS:"));
    }

    #[test]
    fn display_report_handles_none_gracefully() {
        // Construct a report with None values manually
        let report = MemoryReport {
            snapshots: vec![MemorySnapshot {
                label: "no_data".to_string(),
                timestamp: Instant::now(),
                rss_kb: None,
                peak_rss_kb: None,
            }],
            peak_rss_kb: None,
            deltas: Vec::new(),
        };
        let display = report.display_report();
        assert!(display.contains("N/A"));
    }

    #[test]
    fn format_with_commas_works() {
        assert_eq!(format_with_commas(0), "0");
        assert_eq!(format_with_commas(999), "999");
        assert_eq!(format_with_commas(1000), "1,000");
        assert_eq!(format_with_commas(1000000), "1,000,000");
        assert_eq!(format_with_commas(12345), "12,345");
    }

    #[test]
    fn format_opt_delta_positive() {
        assert_eq!(format_opt_delta(Some(1000)), "+1,000");
    }

    #[test]
    fn format_opt_delta_negative() {
        assert_eq!(format_opt_delta(Some(-500)), "-500");
    }

    #[test]
    fn format_opt_delta_none() {
        assert_eq!(format_opt_delta(None), "N/A");
    }
}
