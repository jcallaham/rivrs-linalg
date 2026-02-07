//! Benchmark result data structures and Criterion output collection.

use std::fmt;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::SolverPhase;

/// A single measurement outcome from one (matrix, phase) pair.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult {
    pub matrix_name: String,
    pub phase: SolverPhase,
    pub mean_ns: f64,
    pub std_dev_ns: f64,
    pub median_ns: f64,
    pub iterations: Option<u64>,
    pub throughput_nnz_per_sec: Option<f64>,
    pub matrix_size: Option<usize>,
    pub matrix_nnz: usize,
}

impl fmt::Display for BenchmarkResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mean_ms = self.mean_ns / 1_000_000.0;
        let stddev_ms = self.std_dev_ns / 1_000_000.0;
        write!(
            f,
            "{:<20} {:<10} {:>10.3} ms (+/- {:.3} ms)  nnz={}",
            self.matrix_name, self.phase, mean_ms, stddev_ms, self.matrix_nnz
        )?;
        if let Some(n) = self.matrix_size {
            write!(f, "  n={}", n)?;
        }
        if let Some(tp) = self.throughput_nnz_per_sec {
            write!(f, "  {:.2e} nnz/s", tp)?;
        }
        Ok(())
    }
}

/// Aggregated results from a complete benchmark run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkSuiteResult {
    pub results: Vec<BenchmarkResult>,
    pub peak_rss_kb: Option<u64>,
    pub skipped: Vec<SkippedBenchmark>,
    pub timestamp: String,
}

impl fmt::Display for BenchmarkSuiteResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Benchmark Suite Results ({})", self.timestamp)?;
        writeln!(f, "{:-<80}", "")?;
        for result in &self.results {
            writeln!(f, "  {}", result)?;
        }
        if !self.skipped.is_empty() {
            writeln!(f, "\nSkipped:")?;
            for skip in &self.skipped {
                writeln!(f, "  {}", skip)?;
            }
        }
        if let Some(rss) = self.peak_rss_kb {
            writeln!(f, "\nPeak RSS: {} KB ({:.1} MB)", rss, rss as f64 / 1024.0)?;
        }
        Ok(())
    }
}

/// Records a benchmark that was not executed.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SkippedBenchmark {
    pub matrix_name: String,
    pub phase: SolverPhase,
    pub reason: String,
}

impl fmt::Display for SkippedBenchmark {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} / {}: {}", self.matrix_name, self.phase, self.reason)
    }
}

/// Criterion's estimates.json structure (subset of fields we need).
#[derive(Deserialize)]
struct CriterionEstimates {
    mean: CriterionEstimate,
    median: CriterionEstimate,
    std_dev: CriterionEstimate,
}

#[derive(Deserialize)]
struct CriterionEstimate {
    point_estimate: f64,
}

/// Criterion's benchmark.json structure.
#[derive(Deserialize)]
struct CriterionBenchmark {
    group_id: String,
    value_str: Option<String>,
    throughput: Option<CriterionThroughput>,
}

#[derive(Deserialize)]
enum CriterionThroughput {
    Elements(u64),
    // Variant required for serde deserialization of Criterion's throughput enum.
    #[allow(dead_code)]
    Bytes(u64),
}

/// Collect results from Criterion's output directory.
///
/// Walks `target/criterion/` looking for `new/estimates.json` and
/// `new/benchmark.json` files, parsing them into `BenchmarkResult` entries.
pub fn collect_results(
    criterion_dir: &Path,
    peak_rss_kb: Option<u64>,
) -> Result<BenchmarkSuiteResult, std::io::Error> {
    let mut results = Vec::new();

    if !criterion_dir.exists() {
        return Ok(BenchmarkSuiteResult {
            results,
            peak_rss_kb,
            skipped: vec![],
            timestamp: now_epoch_string(),
        });
    }

    // Walk top-level directories (group directories like ssids_analyze)
    for group_entry in std::fs::read_dir(criterion_dir)? {
        let group_entry = group_entry?;
        let group_path = group_entry.path();
        if !group_path.is_dir() || group_entry.file_name() == "report" {
            continue;
        }

        // Walk matrix directories within each group
        for matrix_entry in std::fs::read_dir(&group_path)? {
            let matrix_entry = matrix_entry?;
            let matrix_path = matrix_entry.path();
            if !matrix_path.is_dir() {
                continue;
            }

            let estimates_path = matrix_path.join("new").join("estimates.json");
            let benchmark_path = matrix_path.join("new").join("benchmark.json");

            if !estimates_path.exists() || !benchmark_path.exists() {
                continue;
            }

            let estimates: CriterionEstimates = match read_json(&estimates_path) {
                Ok(e) => e,
                Err(e) => {
                    eprintln!(
                        "WARNING: failed to parse {}: {}",
                        estimates_path.display(),
                        e
                    );
                    continue;
                }
            };
            let benchmark: CriterionBenchmark = match read_json(&benchmark_path) {
                Ok(b) => b,
                Err(e) => {
                    eprintln!(
                        "WARNING: failed to parse {}: {}",
                        benchmark_path.display(),
                        e
                    );
                    continue;
                }
            };

            let phase = match parse_phase_from_group_id(&benchmark.group_id) {
                Some(p) => p,
                None => continue,
            };

            let matrix_name = benchmark
                .value_str
                .unwrap_or_else(|| matrix_entry.file_name().to_string_lossy().to_string());

            let nnz = match &benchmark.throughput {
                Some(CriterionThroughput::Elements(n)) => *n as usize,
                _ => 0,
            };

            let throughput_nnz_per_sec = if nnz > 0 && estimates.mean.point_estimate > 0.0 {
                Some(nnz as f64 / (estimates.mean.point_estimate * 1e-9))
            } else {
                None
            };

            results.push(BenchmarkResult {
                matrix_name,
                phase,
                mean_ns: estimates.mean.point_estimate,
                std_dev_ns: estimates.std_dev.point_estimate,
                median_ns: estimates.median.point_estimate,
                iterations: None,
                throughput_nnz_per_sec,
                matrix_size: None,
                matrix_nnz: nnz,
            });
        }
    }

    // Sort by pipeline order (Analyze < Factor < Solve < Roundtrip), then by matrix name.
    results.sort_by(|a, b| {
        a.phase
            .cmp(&b.phase)
            .then(a.matrix_name.cmp(&b.matrix_name))
    });

    Ok(BenchmarkSuiteResult {
        results,
        peak_rss_kb,
        skipped: vec![],
        timestamp: now_epoch_string(),
    })
}

fn parse_phase_from_group_id(group_id: &str) -> Option<SolverPhase> {
    // group_id is like "ssids/analyze" or "ssids/roundtrip"
    let phase_str = group_id.rsplit('/').next()?;
    match phase_str {
        "analyze" => Some(SolverPhase::Analyze),
        "factor" => Some(SolverPhase::Factor),
        "solve" => Some(SolverPhase::Solve),
        "roundtrip" => Some(SolverPhase::Roundtrip),
        _ => None,
    }
}

fn read_json<T: serde::de::DeserializeOwned>(path: &Path) -> Result<T, std::io::Error> {
    let content = std::fs::read_to_string(path)?;
    serde_json::from_str(&content).map_err(std::io::Error::other)
}

fn now_epoch_string() -> String {
    let duration = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    format!("{}s-since-epoch", duration.as_secs())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_result() -> BenchmarkResult {
        BenchmarkResult {
            matrix_name: "test-matrix".to_string(),
            phase: SolverPhase::Factor,
            mean_ns: 1_234_567.0,
            std_dev_ns: 12_345.0,
            median_ns: 1_200_000.0,
            iterations: Some(100),
            throughput_nnz_per_sec: Some(1.5e6),
            matrix_size: Some(100),
            matrix_nnz: 500,
        }
    }

    fn sample_suite() -> BenchmarkSuiteResult {
        BenchmarkSuiteResult {
            results: vec![sample_result()],
            peak_rss_kb: Some(65536),
            skipped: vec![SkippedBenchmark {
                matrix_name: "missing-matrix".to_string(),
                phase: SolverPhase::Analyze,
                reason: "matrix file not found".to_string(),
            }],
            timestamp: "2026-02-07T12:00:00Z".to_string(),
        }
    }

    #[test]
    fn benchmark_result_serde_roundtrip() {
        let result = sample_result();
        let json = serde_json::to_string(&result).unwrap();
        let back: BenchmarkResult = serde_json::from_str(&json).unwrap();
        assert_eq!(back.matrix_name, "test-matrix");
        assert_eq!(back.phase, SolverPhase::Factor);
        assert!((back.mean_ns - 1_234_567.0).abs() < f64::EPSILON);
    }

    #[test]
    fn suite_result_serde_roundtrip() {
        let suite = sample_suite();
        let json = serde_json::to_string_pretty(&suite).unwrap();
        let back: BenchmarkSuiteResult = serde_json::from_str(&json).unwrap();
        assert_eq!(back.results.len(), 1);
        assert_eq!(back.skipped.len(), 1);
        assert_eq!(back.peak_rss_kb, Some(65536));
    }

    #[test]
    fn benchmark_result_display() {
        let result = sample_result();
        let s = format!("{}", result);
        assert!(s.contains("test-matrix"));
        assert!(s.contains("factor"));
        assert!(s.contains("ms"));
    }

    #[test]
    fn collect_results_from_empty_dir() {
        let dir = tempfile::tempdir().unwrap();
        let suite = collect_results(dir.path(), Some(1234)).unwrap();
        assert!(suite.results.is_empty());
        assert_eq!(suite.peak_rss_kb, Some(1234));
    }

    #[test]
    fn collect_results_from_nonexistent_dir() {
        let suite = collect_results(std::path::Path::new("/nonexistent"), None).unwrap();
        assert!(suite.results.is_empty());
    }

    #[test]
    fn skipped_display() {
        let skip = SkippedBenchmark {
            matrix_name: "big-matrix".to_string(),
            phase: SolverPhase::Solve,
            reason: "timeout exceeded".to_string(),
        };
        let s = format!("{}", skip);
        assert!(s.contains("big-matrix"));
        assert!(s.contains("solve"));
        assert!(s.contains("timeout exceeded"));
    }
}
