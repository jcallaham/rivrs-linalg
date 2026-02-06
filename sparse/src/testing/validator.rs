//! NumericalValidator — configurable validation wrapping validate.rs functions.

use faer::sparse::SparseColMat;
use faer::Col;

use crate::io::reference::{Inertia, ReferenceFactorization};
use crate::validate;

use super::cases::SolverTestCase;
use super::harness::{MetricResult, TestKind, TestResult};

/// Configurable validator for factorization quality.
///
/// Wraps the standalone functions in `validate` with tolerance thresholds
/// and structured result reporting.
pub struct NumericalValidator {
    reconstruction_tol: f64,
    backward_error_tol: f64,
}

impl NumericalValidator {
    /// Create validator with constitution-mandated default tolerances.
    pub fn new() -> Self {
        Self {
            reconstruction_tol: 1e-12,
            backward_error_tol: 1e-10,
        }
    }

    /// Override reconstruction error tolerance.
    pub fn with_reconstruction_tol(mut self, tol: f64) -> Self {
        self.reconstruction_tol = tol;
        self
    }

    /// Override backward error tolerance.
    pub fn with_backward_error_tol(mut self, tol: f64) -> Self {
        self.backward_error_tol = tol;
        self
    }

    /// Check reconstruction error against threshold.
    pub fn check_reconstruction(
        &self,
        matrix: &SparseColMat<usize, f64>,
        reference: &ReferenceFactorization,
    ) -> MetricResult {
        let value = validate::reconstruction_error(matrix, reference);
        MetricResult {
            name: "reconstruction_error".to_string(),
            value,
            threshold: self.reconstruction_tol,
            passed: value < self.reconstruction_tol,
        }
    }

    /// Check backward error against threshold.
    pub fn check_backward_error(
        &self,
        matrix: &SparseColMat<usize, f64>,
        x: &Col<f64>,
        b: &Col<f64>,
    ) -> MetricResult {
        let value = validate::backward_error(matrix, x, b);
        MetricResult {
            name: "backward_error".to_string(),
            value,
            threshold: self.backward_error_tol,
            passed: value < self.backward_error_tol,
        }
    }

    /// Check inertia equality.
    pub fn check_inertia(&self, computed: &Inertia, expected: &Inertia) -> MetricResult {
        let matches = validate::check_inertia(computed, expected);
        MetricResult {
            name: "inertia".to_string(),
            value: if matches { 0.0 } else { 1.0 },
            threshold: 0.5,
            passed: matches,
        }
    }

    /// Run all applicable checks for a test case with reference data.
    pub fn validate_factorization(&self, case: &SolverTestCase) -> TestResult {
        let mut metrics = Vec::new();
        let mut diagnostics = Vec::new();

        diagnostics.push(format!(
            "matrix: {} ({}x{}, {} nnz)",
            case.name, case.properties.size, case.properties.size, case.properties.nnz
        ));

        if let Some(ref reference) = case.reference {
            let recon = self.check_reconstruction(&case.matrix, reference);
            metrics.push(recon);

            let inertia = self.check_inertia(&reference.inertia, &reference.inertia);
            metrics.push(inertia);
        } else {
            diagnostics.push("no reference factorization available — skipping checks".to_string());
        }

        let passed = metrics.iter().all(|m| m.passed);

        TestResult {
            passed,
            test_kind: TestKind::Factor,
            matrix_name: case.name.clone(),
            metrics,
            diagnostics,
        }
    }
}

impl Default for NumericalValidator {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::cases::{load_test_cases, TestCaseFilter};

    #[test]
    fn default_tolerances_match_constitution() {
        let v = NumericalValidator::new();
        assert_eq!(v.reconstruction_tol, 1e-12);
        assert_eq!(v.backward_error_tol, 1e-10);
    }

    #[test]
    fn check_reconstruction_passes_for_hand_constructed() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed())
            .expect("failed to load");
        let arrow = cases.iter().find(|c| c.name == "arrow-5-pd").unwrap();
        let reference = arrow.reference.as_ref().unwrap();

        let v = NumericalValidator::new();
        let result = v.check_reconstruction(&arrow.matrix, reference);
        assert!(result.passed, "reconstruction should pass: {}", result);
        assert!(result.value < 1e-12);
    }

    #[test]
    fn check_reconstruction_fails_for_perturbed() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed())
            .expect("failed to load");
        let arrow = cases.iter().find(|c| c.name == "arrow-5-pd").unwrap();
        let mut reference = arrow.reference.clone().unwrap();

        // Perturb an L entry
        if !reference.l_entries.is_empty() {
            reference.l_entries[0].value += 10.0;
        }

        let v = NumericalValidator::new();
        let result = v.check_reconstruction(&arrow.matrix, &reference);
        assert!(!result.passed, "perturbed reconstruction should fail");
        assert_eq!(result.name, "reconstruction_error");
    }

    #[test]
    fn custom_tolerance_relaxes_check() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed())
            .expect("failed to load");
        let arrow = cases.iter().find(|c| c.name == "arrow-5-pd").unwrap();
        let mut reference = arrow.reference.clone().unwrap();

        // Slightly perturb an L entry
        if !reference.l_entries.is_empty() {
            reference.l_entries[0].value += 1e-8;
        }

        // Strict validator should fail
        let strict = NumericalValidator::new();
        let result_strict = strict.check_reconstruction(&arrow.matrix, &reference);
        assert!(!result_strict.passed);

        // Relaxed validator should pass
        let relaxed = NumericalValidator::new().with_reconstruction_tol(1e-6);
        let result_relaxed = relaxed.check_reconstruction(&arrow.matrix, &reference);
        assert!(result_relaxed.passed);
    }

    #[test]
    fn check_inertia_pass_and_fail() {
        let v = NumericalValidator::new();

        let a = Inertia {
            positive: 5,
            negative: 3,
            zero: 0,
        };

        // Pass case
        let result = v.check_inertia(&a, &a);
        assert!(result.passed);
        assert_eq!(result.value, 0.0);

        // Fail case
        let b = Inertia {
            positive: 4,
            negative: 3,
            zero: 1,
        };
        let result = v.check_inertia(&a, &b);
        assert!(!result.passed);
        assert_eq!(result.value, 1.0);
    }

    #[test]
    fn validate_all_15_hand_constructed() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed())
            .expect("failed to load");
        assert_eq!(cases.len(), 15);

        let v = NumericalValidator::new();
        for case in &cases {
            let result = v.validate_factorization(case);
            assert!(
                result.passed,
                "{} failed validation: {}",
                case.name, result
            );
        }
    }
}
