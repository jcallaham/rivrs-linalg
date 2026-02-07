//! SolverTest trait, TestResult, MetricResult types.

use std::fmt;

use faer::Col;

use crate::SolverPhase;
use crate::validate;

use super::cases::SolverTestCase;
use super::validator::NumericalValidator;

/// Outcome of a single numerical metric check.
#[derive(Debug, Clone)]
pub struct MetricResult {
    pub name: String,
    pub value: f64,
    pub threshold: f64,
    pub passed: bool,
}

impl fmt::Display for MetricResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let status = if self.passed { "PASS" } else { "FAIL" };
        write!(
            f,
            "{}: {:.2e} (threshold: {:.2e}) — {}",
            self.name, self.value, self.threshold, status
        )
    }
}

/// Complete outcome of a solver test.
#[derive(Debug, Clone)]
pub struct TestResult {
    pub passed: bool,
    pub phase: SolverPhase,
    pub matrix_name: String,
    pub metrics: Vec<MetricResult>,
    pub diagnostics: Vec<String>,
}

impl fmt::Display for TestResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let status = if self.passed { "PASS" } else { "FAIL" };
        write!(f, "[{:?}] {} — {}", self.phase, self.matrix_name, status)?;
        for metric in &self.metrics {
            write!(f, "\n  {}", metric)?;
        }
        for diag in &self.diagnostics {
            write!(f, "\n  {}", diag)?;
        }
        Ok(())
    }
}

/// Interface for testing solver implementations across all phases.
///
/// # Current limitations
///
/// Each method is independent and stateless — it receives a `SolverTestCase`
/// and returns a `TestResult`. This works for the `MockSolver` (which uses
/// reference data directly) but will need rethinking for the real solver,
/// where `test_factor` needs symbolic analysis output and `test_solve` needs
/// the factorization. Options for Phase 1+:
///
/// - Add associated types for intermediate state (analysis result, factorization)
/// - Split into separate traits per phase
/// - Use a builder/session pattern that accumulates state
///
/// The `test_roundtrip` method partially addresses this by combining all phases.
pub trait SolverTest {
    fn test_analyze(&self, case: &SolverTestCase) -> TestResult;
    fn test_factor(&self, case: &SolverTestCase) -> TestResult;
    fn test_solve(&self, case: &SolverTestCase) -> TestResult;
    fn test_roundtrip(&self, case: &SolverTestCase) -> TestResult;
}

/// Mock solver for testing the harness itself.
///
/// Uses the reference factorization directly (no actual computation).
/// Useful for verifying that the test infrastructure works correctly.
pub struct MockSolver {
    validator: NumericalValidator,
}

impl MockSolver {
    pub fn new() -> Self {
        Self {
            validator: NumericalValidator::new(),
        }
    }
}

impl Default for MockSolver {
    fn default() -> Self {
        Self::new()
    }
}

impl SolverTest for MockSolver {
    fn test_analyze(&self, case: &SolverTestCase) -> TestResult {
        let mut metrics = Vec::new();
        let mut diagnostics = Vec::new();
        let n = case.matrix.nrows();

        diagnostics.push(format!("matrix: {} ({}x{})", case.name, n, n));

        if let Some(ref reference) = case.reference {
            // Validate permutation using shared utility
            let valid = validate::validate_permutation(&reference.permutation, n).is_ok();
            metrics.push(MetricResult {
                name: "permutation_valid".to_string(),
                value: if valid { 0.0 } else { 1.0 },
                threshold: 0.5,
                passed: valid,
            });
        } else {
            diagnostics.push("no reference — using identity permutation".to_string());
            metrics.push(MetricResult {
                name: "permutation_valid".to_string(),
                value: 0.0,
                threshold: 0.5,
                passed: true,
            });
        }

        let passed = metrics.iter().all(|m| m.passed);
        TestResult {
            passed,
            phase: SolverPhase::Analyze,
            matrix_name: case.name.clone(),
            metrics,
            diagnostics,
        }
    }

    fn test_factor(&self, case: &SolverTestCase) -> TestResult {
        let mut metrics = Vec::new();
        let mut diagnostics = Vec::new();

        diagnostics.push(format!(
            "matrix: {} ({}x{})",
            case.name,
            case.matrix.nrows(),
            case.matrix.ncols()
        ));

        if let Some(ref reference) = case.reference {
            let recon = self.validator.check_reconstruction(&case.matrix, reference);
            let inertia = self
                .validator
                .check_inertia(&reference.inertia, &reference.inertia);
            metrics.push(recon);
            metrics.push(inertia);
        } else {
            diagnostics.push("no reference factorization — cannot validate factor".to_string());
        }

        let passed = metrics.iter().all(|m| m.passed);
        TestResult {
            passed,
            phase: SolverPhase::Factor,
            matrix_name: case.name.clone(),
            metrics,
            diagnostics,
        }
    }

    fn test_solve(&self, case: &SolverTestCase) -> TestResult {
        let mut metrics = Vec::new();
        let mut diagnostics = Vec::new();
        let n = case.matrix.nrows();

        diagnostics.push(format!("matrix: {} ({}x{})", case.name, n, n));

        // Generate a known x and compute b = A*x, then check backward error
        let x_exact = Col::<f64>::from_fn(n, |i| (i + 1) as f64);
        let b = validate::dense_matvec(&case.matrix, &x_exact);

        let berr = self
            .validator
            .check_backward_error(&case.matrix, &x_exact, &b);
        metrics.push(berr);

        let passed = metrics.iter().all(|m| m.passed);
        TestResult {
            passed,
            phase: SolverPhase::Solve,
            matrix_name: case.name.clone(),
            metrics,
            diagnostics,
        }
    }

    fn test_roundtrip(&self, case: &SolverTestCase) -> TestResult {
        let mut metrics = Vec::new();
        let mut diagnostics = Vec::new();
        let n = case.matrix.nrows();

        diagnostics.push(format!("matrix: {} ({}x{})", case.name, n, n));

        // Reconstruction check
        if let Some(ref reference) = case.reference {
            let recon = self.validator.check_reconstruction(&case.matrix, reference);
            metrics.push(recon);

            let inertia = self
                .validator
                .check_inertia(&reference.inertia, &reference.inertia);
            metrics.push(inertia);
        }

        // Backward error check with known solution
        let x_exact = Col::<f64>::from_fn(n, |i| (i + 1) as f64);
        let b = validate::dense_matvec(&case.matrix, &x_exact);
        let berr = self
            .validator
            .check_backward_error(&case.matrix, &x_exact, &b);
        metrics.push(berr);

        let passed = metrics.iter().all(|m| m.passed);
        TestResult {
            passed,
            phase: SolverPhase::Roundtrip,
            matrix_name: case.name.clone(),
            metrics,
            diagnostics,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::cases::{TestCaseFilter, load_test_cases};

    #[test]
    fn metric_result_display_pass() {
        let m = MetricResult {
            name: "reconstruction_error".to_string(),
            value: 1.5e-14,
            threshold: 1e-12,
            passed: true,
        };
        let s = format!("{}", m);
        assert!(s.contains("reconstruction_error"));
        assert!(s.contains("1.50e-14"));
        assert!(s.contains("1.00e-12"));
        assert!(s.contains("PASS"));
    }

    #[test]
    fn metric_result_display_fail() {
        let m = MetricResult {
            name: "backward_error".to_string(),
            value: 2.3e-8,
            threshold: 1e-10,
            passed: false,
        };
        let s = format!("{}", m);
        assert!(s.contains("backward_error"));
        assert!(s.contains("FAIL"));
    }

    #[test]
    fn test_result_passed_reflects_metrics() {
        let result = TestResult {
            passed: false,
            phase: SolverPhase::Factor,
            matrix_name: "test-matrix".to_string(),
            metrics: vec![
                MetricResult {
                    name: "recon".to_string(),
                    value: 1e-14,
                    threshold: 1e-12,
                    passed: true,
                },
                MetricResult {
                    name: "inertia".to_string(),
                    value: 1.0,
                    threshold: 0.5,
                    passed: false,
                },
            ],
            diagnostics: vec!["inertia mismatch".to_string()],
        };
        assert!(!result.passed);
        let s = format!("{}", result);
        assert!(s.contains("FAIL"));
        assert!(s.contains("test-matrix"));
    }

    #[test]
    fn mock_solver_test_roundtrip_passes() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed()).expect("failed to load");
        let arrow = cases.iter().find(|c| c.name == "arrow-5-pd").unwrap();

        let solver = MockSolver::new();
        let result = solver.test_roundtrip(arrow);
        assert!(result.passed, "roundtrip should pass: {}", result);
    }

    #[test]
    fn mock_solver_test_analyze_passes() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed()).expect("failed to load");
        let arrow = cases.iter().find(|c| c.name == "arrow-5-pd").unwrap();

        let solver = MockSolver::new();
        let result = solver.test_analyze(arrow);
        assert!(result.passed, "analyze should pass: {}", result);
    }

    #[test]
    fn mock_solver_test_factor_passes() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed()).expect("failed to load");
        let arrow = cases.iter().find(|c| c.name == "arrow-5-pd").unwrap();

        let solver = MockSolver::new();
        let result = solver.test_factor(arrow);
        assert!(result.passed, "factor should pass: {}", result);
    }

    #[test]
    fn mock_solver_test_solve_passes() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed()).expect("failed to load");
        let arrow = cases.iter().find(|c| c.name == "arrow-5-pd").unwrap();

        let solver = MockSolver::new();
        let result = solver.test_solve(arrow);
        assert!(result.passed, "solve should pass: {}", result);
    }

    #[test]
    fn solver_test_all_hand_constructed() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed()).expect("failed to load");
        assert_eq!(cases.len(), 15);

        let solver = MockSolver::new();
        for case in &cases {
            let result = solver.test_roundtrip(case);
            assert!(result.passed, "{} failed: {}", case.name, result);
        }
    }
}
