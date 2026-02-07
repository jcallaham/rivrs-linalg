//! Baseline management and regression detection.

use std::fmt;
use std::io::Write;
use std::path::Path;

use serde::{Deserialize, Serialize};

use super::config::BenchmarkPhase;
use super::results::BenchmarkSuiteResult;

/// A saved suite result used for regression detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Baseline {
    pub name: String,
    pub suite_result: BenchmarkSuiteResult,
    pub created: String,
}

/// A detected performance regression.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Regression {
    pub matrix_name: String,
    pub phase: BenchmarkPhase,
    pub baseline_mean_ns: f64,
    pub current_mean_ns: f64,
    pub change_pct: f64,
}

impl fmt::Display for Regression {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "REGRESSION {}/{}: {:.3}ms -> {:.3}ms (+{:.1}%)",
            self.matrix_name,
            self.phase,
            self.baseline_mean_ns / 1_000_000.0,
            self.current_mean_ns / 1_000_000.0,
            self.change_pct
        )
    }
}

/// A detected performance improvement.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Improvement {
    pub matrix_name: String,
    pub phase: BenchmarkPhase,
    pub baseline_mean_ns: f64,
    pub current_mean_ns: f64,
    pub change_pct: f64,
}

impl fmt::Display for Improvement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "IMPROVED {}/{}: {:.3}ms -> {:.3}ms ({:.1}%)",
            self.matrix_name,
            self.phase,
            self.baseline_mean_ns / 1_000_000.0,
            self.current_mean_ns / 1_000_000.0,
            self.change_pct
        )
    }
}

/// Comparison of current results against a baseline.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegressionReport {
    pub baseline_name: String,
    pub regressions: Vec<Regression>,
    pub improvements: Vec<Improvement>,
    pub unchanged: Vec<String>,
    pub threshold_pct: f64,
}

impl fmt::Display for RegressionReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "Regression Report (vs '{}', threshold: {:.1}%)",
            self.baseline_name, self.threshold_pct
        )?;
        writeln!(f, "{:-<60}", "")?;
        if self.regressions.is_empty() {
            writeln!(f, "No regressions detected.")?;
        } else {
            writeln!(f, "Regressions ({}):", self.regressions.len())?;
            for r in &self.regressions {
                writeln!(f, "  {}", r)?;
            }
        }
        if !self.improvements.is_empty() {
            writeln!(f, "Improvements ({}):", self.improvements.len())?;
            for i in &self.improvements {
                writeln!(f, "  {}", i)?;
            }
        }
        writeln!(f, "Unchanged: {}", self.unchanged.len())?;
        Ok(())
    }
}

/// Save a benchmark suite result as a named baseline.
pub fn save_baseline(
    suite_result: &BenchmarkSuiteResult,
    name: &str,
    output_dir: &Path,
) -> Result<(), std::io::Error> {
    std::fs::create_dir_all(output_dir)?;
    let path = output_dir.join(format!("{}.json", name));
    let baseline = Baseline {
        name: name.to_string(),
        suite_result: suite_result.clone(),
        created: suite_result.timestamp.clone(),
    };
    let json = serde_json::to_string_pretty(&baseline).map_err(std::io::Error::other)?;
    let mut file = std::fs::File::create(&path)?;
    file.write_all(json.as_bytes())?;
    Ok(())
}

/// Load a previously saved baseline by name.
pub fn load_baseline(name: &str, dir: &Path) -> Result<Baseline, std::io::Error> {
    let path = dir.join(format!("{}.json", name));
    let content = std::fs::read_to_string(&path)?;
    let baseline: Baseline = serde_json::from_str(&content)
        .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
    Ok(baseline)
}

/// Compare current results against a baseline and detect regressions.
///
/// Matches results by (matrix_name, phase) pairs. A regression is flagged
/// when the current mean exceeds the baseline mean by more than `threshold_pct`.
/// An improvement is flagged when the current mean is less than the baseline
/// mean by more than `threshold_pct`.
pub fn detect_regressions(
    current: &BenchmarkSuiteResult,
    baseline: &Baseline,
    threshold_pct: f64,
) -> RegressionReport {
    let mut regressions = Vec::new();
    let mut improvements = Vec::new();
    let mut unchanged = Vec::new();

    for current_result in &current.results {
        let key = (&current_result.matrix_name, &current_result.phase);

        if let Some(baseline_result) = baseline
            .suite_result
            .results
            .iter()
            .find(|b| b.matrix_name == *key.0 && b.phase == *key.1)
        {
            if baseline_result.mean_ns == 0.0 {
                unchanged.push(format!("{}/{}", key.0, key.1));
                continue;
            }

            let change_pct = ((current_result.mean_ns - baseline_result.mean_ns)
                / baseline_result.mean_ns)
                * 100.0;

            if change_pct > threshold_pct {
                regressions.push(Regression {
                    matrix_name: current_result.matrix_name.clone(),
                    phase: current_result.phase,
                    baseline_mean_ns: baseline_result.mean_ns,
                    current_mean_ns: current_result.mean_ns,
                    change_pct,
                });
            } else if change_pct < -threshold_pct {
                improvements.push(Improvement {
                    matrix_name: current_result.matrix_name.clone(),
                    phase: current_result.phase,
                    baseline_mean_ns: baseline_result.mean_ns,
                    current_mean_ns: current_result.mean_ns,
                    change_pct,
                });
            } else {
                unchanged.push(format!("{}/{}", key.0, key.1));
            }
        }
    }

    RegressionReport {
        baseline_name: baseline.name.clone(),
        regressions,
        improvements,
        unchanged,
        threshold_pct,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::benchmarking::results::BenchmarkResult;

    fn make_result(matrix: &str, phase: BenchmarkPhase, mean_ns: f64) -> BenchmarkResult {
        BenchmarkResult {
            matrix_name: matrix.to_string(),
            phase,
            mean_ns,
            std_dev_ns: mean_ns * 0.01,
            median_ns: mean_ns,
            iterations: 100,
            throughput_nnz_per_sec: None,
            matrix_size: 100,
            matrix_nnz: 500,
        }
    }

    fn make_suite(results: Vec<BenchmarkResult>) -> BenchmarkSuiteResult {
        BenchmarkSuiteResult {
            results,
            peak_rss_kb: None,
            skipped: vec![],
            timestamp: "2026-02-07T12:00:00Z".to_string(),
        }
    }

    #[test]
    fn save_and_load_baseline_roundtrip() {
        let suite = make_suite(vec![make_result(
            "matrix-a",
            BenchmarkPhase::Factor,
            1_000_000.0,
        )]);
        let dir = tempfile::tempdir().unwrap();
        save_baseline(&suite, "test-baseline", dir.path()).unwrap();
        let loaded = load_baseline("test-baseline", dir.path()).unwrap();
        assert_eq!(loaded.name, "test-baseline");
        assert_eq!(loaded.suite_result.results.len(), 1);
        assert_eq!(loaded.suite_result.results[0].matrix_name, "matrix-a");
    }

    #[test]
    fn detect_regression() {
        let baseline_suite = make_suite(vec![make_result(
            "matrix-a",
            BenchmarkPhase::Factor,
            1_000_000.0,
        )]);
        let baseline = Baseline {
            name: "main".to_string(),
            suite_result: baseline_suite,
            created: "2026-02-07T12:00:00Z".to_string(),
        };

        // 10% slower → regression at 5% threshold
        let current = make_suite(vec![make_result(
            "matrix-a",
            BenchmarkPhase::Factor,
            1_100_000.0,
        )]);

        let report = detect_regressions(&current, &baseline, 5.0);
        assert_eq!(report.regressions.len(), 1);
        assert!(report.improvements.is_empty());
        assert!(report.regressions[0].change_pct > 9.0);
    }

    #[test]
    fn detect_improvement() {
        let baseline_suite = make_suite(vec![make_result(
            "matrix-a",
            BenchmarkPhase::Factor,
            1_000_000.0,
        )]);
        let baseline = Baseline {
            name: "main".to_string(),
            suite_result: baseline_suite,
            created: "2026-02-07T12:00:00Z".to_string(),
        };

        // 20% faster → improvement
        let current = make_suite(vec![make_result(
            "matrix-a",
            BenchmarkPhase::Factor,
            800_000.0,
        )]);

        let report = detect_regressions(&current, &baseline, 5.0);
        assert!(report.regressions.is_empty());
        assert_eq!(report.improvements.len(), 1);
        assert!(report.improvements[0].change_pct < -15.0);
    }

    #[test]
    fn detect_unchanged() {
        let baseline_suite = make_suite(vec![make_result(
            "matrix-a",
            BenchmarkPhase::Factor,
            1_000_000.0,
        )]);
        let baseline = Baseline {
            name: "main".to_string(),
            suite_result: baseline_suite,
            created: "2026-02-07T12:00:00Z".to_string(),
        };

        // 2% change → within 5% threshold → unchanged
        let current = make_suite(vec![make_result(
            "matrix-a",
            BenchmarkPhase::Factor,
            1_020_000.0,
        )]);

        let report = detect_regressions(&current, &baseline, 5.0);
        assert!(report.regressions.is_empty());
        assert!(report.improvements.is_empty());
        assert_eq!(report.unchanged.len(), 1);
    }
}
