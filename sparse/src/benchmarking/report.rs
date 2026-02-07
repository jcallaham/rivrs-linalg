//! Export and reporting utilities for benchmark results.

use std::fmt::Write as FmtWrite;
use std::io::Write;
use std::path::Path;

use super::baseline::RegressionReport;
use super::results::BenchmarkSuiteResult;

/// Escape a string for CSV output by quoting if it contains commas, quotes, or newlines.
fn csv_escape(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Export benchmark results to a CSV file.
pub fn export_csv(suite: &BenchmarkSuiteResult, path: &Path) -> Result<(), std::io::Error> {
    let mut file = std::fs::File::create(path)?;
    writeln!(
        file,
        "matrix_name,phase,mean_ns,std_dev_ns,median_ns,iterations,matrix_size,matrix_nnz,throughput_nnz_per_sec"
    )?;
    for r in &suite.results {
        writeln!(
            file,
            "{},{},{},{},{},{},{},{},{}",
            csv_escape(&r.matrix_name),
            r.phase,
            r.mean_ns,
            r.std_dev_ns,
            r.median_ns,
            r.iterations.map_or("".to_string(), |v| v.to_string()),
            r.matrix_size.map_or("".to_string(), |v| v.to_string()),
            r.matrix_nnz,
            r.throughput_nnz_per_sec
                .map_or("".to_string(), |v| format!("{}", v))
        )?;
    }
    Ok(())
}

/// Export benchmark results to a JSON file.
pub fn export_json(suite: &BenchmarkSuiteResult, path: &Path) -> Result<(), std::io::Error> {
    let json = serde_json::to_string_pretty(suite).map_err(std::io::Error::other)?;
    std::fs::write(path, json)?;
    Ok(())
}

/// Generate a Markdown table from benchmark suite results.
pub fn generate_markdown_table(suite: &BenchmarkSuiteResult) -> String {
    let mut out = String::new();
    writeln!(
        out,
        "| Matrix | Phase | Mean (ms) | Std Dev (ms) | Throughput (nnz/s) |"
    )
    .unwrap();
    writeln!(
        out,
        "|--------|-------|---------:|------------:|------------------:|"
    )
    .unwrap();
    for r in &suite.results {
        let tp = r
            .throughput_nnz_per_sec
            .map_or("-".to_string(), |v| format!("{:.2e}", v));
        writeln!(
            out,
            "| {} | {} | {:.3} | {:.3} | {} |",
            r.matrix_name,
            r.phase,
            r.mean_ns / 1_000_000.0,
            r.std_dev_ns / 1_000_000.0,
            tp
        )
        .unwrap();
    }
    out
}

/// Generate a Markdown table from a regression report.
pub fn generate_regression_markdown_table(report: &RegressionReport) -> String {
    let mut out = String::new();
    writeln!(
        out,
        "| Matrix | Phase | Baseline (ms) | Current (ms) | Change (%) | Status |"
    )
    .unwrap();
    writeln!(
        out,
        "|--------|-------|-------------:|------------:|-----------:|--------|"
    )
    .unwrap();
    for r in &report.regressions {
        writeln!(
            out,
            "| {} | {} | {:.3} | {:.3} | +{:.1} | REGRESSION |",
            r.matrix_name,
            r.phase,
            r.baseline_mean_ns / 1_000_000.0,
            r.current_mean_ns / 1_000_000.0,
            r.change_pct
        )
        .unwrap();
    }
    for i in &report.improvements {
        writeln!(
            out,
            "| {} | {} | {:.3} | {:.3} | {:.1} | IMPROVED |",
            i.matrix_name,
            i.phase,
            i.baseline_mean_ns / 1_000_000.0,
            i.current_mean_ns / 1_000_000.0,
            i.change_pct
        )
        .unwrap();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SolverPhase;
    use crate::benchmarking::baseline::{Improvement, Regression};
    use crate::benchmarking::results::{BenchmarkResult, BenchmarkSuiteResult};

    fn sample_suite() -> BenchmarkSuiteResult {
        BenchmarkSuiteResult {
            results: vec![
                BenchmarkResult {
                    matrix_name: "matrix-a".to_string(),
                    phase: SolverPhase::Factor,
                    mean_ns: 1_234_567.0,
                    std_dev_ns: 12_345.0,
                    median_ns: 1_200_000.0,
                    iterations: Some(100),
                    throughput_nnz_per_sec: Some(1.5e6),
                    matrix_size: Some(100),
                    matrix_nnz: 500,
                },
                BenchmarkResult {
                    matrix_name: "matrix-b".to_string(),
                    phase: SolverPhase::Analyze,
                    mean_ns: 500_000.0,
                    std_dev_ns: 5_000.0,
                    median_ns: 490_000.0,
                    iterations: Some(200),
                    throughput_nnz_per_sec: None,
                    matrix_size: Some(50),
                    matrix_nnz: 200,
                },
            ],
            peak_rss_kb: Some(65536),
            skipped: vec![],
            timestamp: "2026-02-07T12:00:00Z".to_string(),
        }
    }

    #[test]
    fn csv_export_has_correct_columns() {
        let suite = sample_suite();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("results.csv");
        export_csv(&suite, &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert_eq!(lines.len(), 3); // header + 2 data rows
        assert!(lines[0].contains("matrix_name"));
        assert!(lines[0].contains("throughput_nnz_per_sec"));
        // Verify comma-separated column count
        assert_eq!(lines[0].split(',').count(), 9);
        assert_eq!(lines[1].split(',').count(), 9);
    }

    #[test]
    fn csv_escapes_commas_in_names() {
        let suite = BenchmarkSuiteResult {
            results: vec![BenchmarkResult {
                matrix_name: "matrix,with,commas".to_string(),
                phase: SolverPhase::Factor,
                mean_ns: 1_000.0,
                std_dev_ns: 10.0,
                median_ns: 1_000.0,
                iterations: None,
                throughput_nnz_per_sec: None,
                matrix_size: None,
                matrix_nnz: 100,
            }],
            peak_rss_kb: None,
            skipped: vec![],
            timestamp: "test".to_string(),
        };
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("results.csv");
        export_csv(&suite, &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let lines: Vec<&str> = content.lines().collect();
        assert!(lines[1].starts_with("\"matrix,with,commas\""));
    }

    #[test]
    fn json_export_roundtrip() {
        let suite = sample_suite();
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("results.json");
        export_json(&suite, &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        let back: BenchmarkSuiteResult = serde_json::from_str(&content).unwrap();
        assert_eq!(back.results.len(), 2);
        assert_eq!(back.peak_rss_kb, Some(65536));
    }

    #[test]
    fn markdown_table_structure() {
        let suite = sample_suite();
        let table = generate_markdown_table(&suite);
        let lines: Vec<&str> = table.lines().collect();
        assert!(lines.len() >= 4); // header + separator + 2 data rows
        assert!(lines[0].contains("Matrix"));
        assert!(lines[0].contains("Phase"));
        assert!(lines[0].contains("Mean (ms)"));
        // Data rows
        assert!(lines[2].contains("matrix-a"));
        assert!(lines[3].contains("matrix-b"));
    }

    #[test]
    fn regression_markdown_table() {
        let report = RegressionReport {
            baseline_name: "main".to_string(),
            regressions: vec![Regression {
                matrix_name: "matrix-a".to_string(),
                phase: SolverPhase::Factor,
                baseline_mean_ns: 1_000_000.0,
                current_mean_ns: 1_200_000.0,
                change_pct: 20.0,
            }],
            improvements: vec![Improvement {
                matrix_name: "matrix-b".to_string(),
                phase: SolverPhase::Analyze,
                baseline_mean_ns: 500_000.0,
                current_mean_ns: 400_000.0,
                change_pct: -20.0,
            }],
            unchanged: vec![],
            threshold_pct: 5.0,
        };
        let table = generate_regression_markdown_table(&report);
        assert!(table.contains("REGRESSION"));
        assert!(table.contains("IMPROVED"));
        assert!(table.contains("matrix-a"));
        assert!(table.contains("matrix-b"));
    }
}
