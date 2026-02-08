//! Summary report generation and Chrome Trace export.

use std::fmt::Write;
use std::time::Duration;

use serde_json::{Value, json};

use super::section::{ProfileEvent, ProfileSection};
use super::session::FinishedSession;

impl FinishedSession {
    /// Generate a hierarchical text summary report.
    pub fn summary_report(&self) -> String {
        let total_ms = self.total_duration.as_secs_f64() * 1000.0;
        let mut out = String::new();

        writeln!(
            out,
            "Profile Summary (total: {})",
            format_duration(self.total_duration)
        )
        .unwrap();
        writeln!(out, "{:─<70}", "").unwrap();
        writeln!(
            out,
            "  {:<20} {:>10} {:>8} {:>10} {:>10} {:>10} {:>8}",
            "Section", "Total", "Calls", "Mean", "Min", "Max", "Parent%"
        )
        .unwrap();

        for section in &self.sections {
            write_section(&mut out, section, 0, self.total_duration);
        }

        let _ = total_ms; // used indirectly via self.total_duration
        out
    }

    /// Export timing data as Chrome Trace Event format JSON.
    ///
    /// Individual invocations are preserved (not aggregated) in the trace export.
    pub fn export_chrome_trace(&self) -> String {
        format_chrome_trace(&self.events)
    }
}

/// Write a section and its children recursively to the output string.
fn write_section(
    out: &mut String,
    section: &ProfileSection,
    indent: usize,
    _parent_duration: Duration,
) {
    let prefix = "  ".repeat(indent + 1);
    writeln!(
        out,
        "{prefix}{:<width$} {:>10} {:>8} {:>10} {:>10} {:>10} {:>7.1}%",
        section.name,
        format_duration(section.total_duration),
        section.call_count,
        format_duration(section.mean_duration),
        format_duration(section.min_duration),
        format_duration(section.max_duration),
        section.parent_pct,
        width = 20 - indent * 2,
    )
    .unwrap();

    for child in &section.children {
        write_section(out, child, indent + 1, section.total_duration);
    }
}

/// Format a duration for display in the report.
fn format_duration(d: Duration) -> String {
    let us = d.as_secs_f64() * 1_000_000.0;
    if us < 1000.0 {
        format!("{us:.1} µs")
    } else if us < 1_000_000.0 {
        format!("{:.1} ms", us / 1000.0)
    } else {
        format!("{:.2} s", us / 1_000_000.0)
    }
}

/// Format events as Chrome Trace Event JSON.
///
/// Produces a JSON string with a `traceEvents` array of Complete Events (type "X").
pub fn format_chrome_trace(events: &[ProfileEvent]) -> String {
    let pid = std::process::id();

    let trace_events: Vec<Value> = events
        .iter()
        .map(|e| {
            json!({
                "name": e.name,
                "cat": "profiling",
                "ph": "X",
                "ts": e.start_us,
                "dur": e.duration.as_secs_f64() * 1_000_000.0,
                "pid": pid,
                "tid": e.thread_idx,
                "args": {}
            })
        })
        .collect();

    let trace = json!({ "traceEvents": trace_events });
    serde_json::to_string_pretty(&trace).unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use std::time::Duration;

    use super::*;
    use crate::profiling::ProfileSession;

    #[test]
    fn summary_report_contains_header() {
        let session = ProfileSession::new();
        {
            let _guard = session.enter_section("analyze");
            std::thread::sleep(Duration::from_millis(1));
        }
        let finished = session.finish();
        let report = finished.summary_report();

        assert!(report.contains("Profile Summary"), "should have header");
        assert!(report.contains("Section"), "should have column header");
        assert!(report.contains("Total"), "should have Total column");
        assert!(report.contains("Calls"), "should have Calls column");
        assert!(report.contains("Mean"), "should have Mean column");
        assert!(report.contains("Min"), "should have Min column");
        assert!(report.contains("Max"), "should have Max column");
        assert!(report.contains("Parent%"), "should have Parent% column");
    }

    #[test]
    fn summary_report_shows_section_names() {
        let session = ProfileSession::new();
        {
            let _guard = session.enter_section("factor");
            {
                let _inner = session.enter_section("dense_kernel");
                std::thread::sleep(Duration::from_millis(1));
            }
        }
        let finished = session.finish();
        let report = finished.summary_report();

        assert!(report.contains("factor"), "should contain outer section");
        assert!(
            report.contains("dense_kernel"),
            "should contain inner section"
        );
    }

    #[test]
    fn summary_report_hierarchical_indentation() {
        let session = ProfileSession::new();
        {
            let _guard = session.enter_section("outer");
            {
                let _inner = session.enter_section("inner");
                std::thread::sleep(Duration::from_millis(1));
            }
        }
        let finished = session.finish();
        let report = finished.summary_report();

        // Inner section should be more indented than outer
        let lines: Vec<&str> = report.lines().collect();
        let outer_line = lines.iter().find(|l| l.contains("outer")).unwrap();
        let inner_line = lines.iter().find(|l| l.contains("inner")).unwrap();

        let outer_indent = outer_line.len() - outer_line.trim_start().len();
        let inner_indent = inner_line.len() - inner_line.trim_start().len();
        assert!(
            inner_indent > outer_indent,
            "inner should be more indented than outer"
        );
    }

    #[test]
    fn chrome_trace_valid_json() {
        let session = ProfileSession::new();
        {
            let _guard = session.enter_section("test");
            std::thread::sleep(Duration::from_millis(1));
        }
        let finished = session.finish();
        let json_str = finished.export_chrome_trace();

        let parsed: Value = serde_json::from_str(&json_str).expect("should be valid JSON");
        assert!(
            parsed.get("traceEvents").is_some(),
            "should have traceEvents"
        );
    }

    #[test]
    fn chrome_trace_complete_events() {
        let session = ProfileSession::new();
        {
            let _guard = session.enter_section("analyze");
            std::thread::sleep(Duration::from_millis(1));
        }
        let finished = session.finish();
        let json_str = finished.export_chrome_trace();
        let parsed: Value = serde_json::from_str(&json_str).unwrap();
        let events = parsed["traceEvents"].as_array().unwrap();

        assert!(!events.is_empty());
        let event = &events[0];
        assert_eq!(event["name"], "analyze");
        assert_eq!(event["ph"], "X");
        assert!(event["ts"].as_f64().is_some(), "should have ts");
        assert!(event["dur"].as_f64().is_some(), "should have dur");
        assert!(event["pid"].as_u64().is_some(), "should have pid");
        assert!(event["tid"].as_u64().is_some(), "should have tid");
    }

    #[test]
    fn chrome_trace_nested_events_overlap() {
        let session = ProfileSession::new();
        {
            let _outer = session.enter_section("outer");
            {
                let _inner = session.enter_section("inner");
                std::thread::sleep(Duration::from_millis(1));
            }
        }
        let finished = session.finish();
        let json_str = finished.export_chrome_trace();
        let parsed: Value = serde_json::from_str(&json_str).unwrap();
        let events = parsed["traceEvents"].as_array().unwrap();

        let outer = events.iter().find(|e| e["name"] == "outer").unwrap();
        let inner = events.iter().find(|e| e["name"] == "inner").unwrap();

        let outer_start = outer["ts"].as_f64().unwrap();
        let outer_dur = outer["dur"].as_f64().unwrap();
        let inner_start = inner["ts"].as_f64().unwrap();

        // Inner should start after outer and end before outer ends
        assert!(inner_start >= outer_start, "inner should start after outer");
        assert!(
            inner_start < outer_start + outer_dur,
            "inner should start before outer ends"
        );
    }

    #[test]
    fn format_duration_microseconds() {
        assert_eq!(format_duration(Duration::from_micros(50)), "50.0 µs");
    }

    #[test]
    fn format_duration_milliseconds() {
        assert_eq!(format_duration(Duration::from_millis(12)), "12.0 ms");
    }

    #[test]
    fn format_duration_seconds() {
        assert_eq!(format_duration(Duration::from_secs(2)), "2.00 s");
    }

    #[test]
    fn chrome_trace_empty_session() {
        let session = ProfileSession::new();
        let finished = session.finish();
        let json_str = finished.export_chrome_trace();
        let parsed: Value = serde_json::from_str(&json_str).unwrap();
        let events = parsed["traceEvents"].as_array().unwrap();
        assert!(events.is_empty());
    }
}
