//! The `Benchmarkable` trait and `MockBenchmarkable` for testing.

use std::any::Any;

use faer::sparse::SparseColMat;

use crate::SolverPhase;

/// Interface for benchmarking solver phases.
///
/// Each method returns `None` if the phase is not yet implemented,
/// causing the benchmark harness to skip it with a diagnostic message.
/// Opaque `Box<dyn Any>` return types prevent the compiler from optimizing
/// away computation and allow passing intermediate state between phases.
pub trait Benchmarkable {
    /// Check whether the given phase is supported without running it.
    ///
    /// Default implementation returns `true` for all phases. Override to
    /// indicate which phases are available without incurring the cost of
    /// a full execution.
    fn supports_phase(&self, phase: SolverPhase) -> bool {
        let _ = phase;
        true
    }

    /// Run symbolic analysis on the matrix.
    fn bench_analyze(&self, matrix: &SparseColMat<usize, f64>) -> Option<Box<dyn Any>>;

    /// Run numeric factorization, optionally using a prior analysis result.
    fn bench_factor(
        &self,
        matrix: &SparseColMat<usize, f64>,
        analysis: Option<&dyn Any>,
    ) -> Option<Box<dyn Any>>;

    /// Run triangular solve using a factorization result and right-hand side.
    fn bench_solve(
        &self,
        matrix: &SparseColMat<usize, f64>,
        factorization: Option<&dyn Any>,
        rhs: &[f64],
    ) -> Option<Vec<f64>>;

    /// Run the full pipeline: analyze -> factor -> solve.
    ///
    /// Default implementation chains the three individual phases.
    fn bench_roundtrip(&self, matrix: &SparseColMat<usize, f64>, rhs: &[f64]) -> Option<Vec<f64>> {
        let analysis = self.bench_analyze(matrix)?;
        let factorization = self.bench_factor(matrix, Some(analysis.as_ref()))?;
        self.bench_solve(matrix, Some(factorization.as_ref()), rhs)
    }
}

/// Mock solver for testing the benchmark harness.
///
/// Returns synthetic results for all phases without performing real computation.
/// Configurable to return `None` for specific phases to test skip behavior.
pub struct MockBenchmarkable {
    pub analyze_enabled: bool,
    pub factor_enabled: bool,
    pub solve_enabled: bool,
}

impl MockBenchmarkable {
    /// Create a mock with all phases enabled.
    pub fn new() -> Self {
        Self {
            analyze_enabled: true,
            factor_enabled: true,
            solve_enabled: true,
        }
    }

    /// Create a mock with specific phases enabled or disabled for testing skip behavior.
    pub fn with_phases_enabled(analyze: bool, factor: bool, solve: bool) -> Self {
        Self {
            analyze_enabled: analyze,
            factor_enabled: factor,
            solve_enabled: solve,
        }
    }
}

impl Default for MockBenchmarkable {
    fn default() -> Self {
        Self::new()
    }
}

impl Benchmarkable for MockBenchmarkable {
    fn supports_phase(&self, phase: SolverPhase) -> bool {
        match phase {
            SolverPhase::Analyze => self.analyze_enabled,
            SolverPhase::Factor => self.factor_enabled,
            SolverPhase::Solve => self.solve_enabled,
            SolverPhase::Roundtrip => {
                self.analyze_enabled && self.factor_enabled && self.solve_enabled
            }
        }
    }

    fn bench_analyze(&self, matrix: &SparseColMat<usize, f64>) -> Option<Box<dyn Any>> {
        if !self.analyze_enabled {
            return None;
        }
        // Return a trivial identity permutation as "analysis result"
        let n = matrix.nrows();
        let perm: Vec<usize> = (0..n).collect();
        Some(Box::new(perm))
    }

    fn bench_factor(
        &self,
        matrix: &SparseColMat<usize, f64>,
        _analysis: Option<&dyn Any>,
    ) -> Option<Box<dyn Any>> {
        if !self.factor_enabled {
            return None;
        }
        // Return a trivial "factorization" (just the matrix dimension)
        let n = matrix.nrows();
        Some(Box::new(n))
    }

    fn bench_solve(
        &self,
        _matrix: &SparseColMat<usize, f64>,
        _factorization: Option<&dyn Any>,
        rhs: &[f64],
    ) -> Option<Vec<f64>> {
        if !self.solve_enabled {
            return None;
        }
        // Return rhs as "solution" (identity solve)
        Some(rhs.to_vec())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::testing::{TestCaseFilter, load_test_cases};

    #[test]
    fn mock_all_phases_enabled() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed()).unwrap();
        let case = &cases[0];
        let mock = MockBenchmarkable::new();

        let analysis = mock.bench_analyze(&case.matrix);
        assert!(analysis.is_some(), "analyze should return Some");

        let factor = mock.bench_factor(&case.matrix, analysis.as_deref());
        assert!(factor.is_some(), "factor should return Some");

        let n = case.matrix.nrows();
        let rhs: Vec<f64> = vec![1.0; n];
        let solve = mock.bench_solve(&case.matrix, factor.as_deref(), &rhs);
        assert!(solve.is_some(), "solve should return Some");
    }

    #[test]
    fn mock_disabled_phase_returns_none() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed()).unwrap();
        let case = &cases[0];
        let mock = MockBenchmarkable::with_phases_enabled(true, false, true);

        assert!(mock.bench_analyze(&case.matrix).is_some());
        assert!(mock.bench_factor(&case.matrix, None).is_none());
    }

    #[test]
    fn mock_roundtrip_chains_phases() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed()).unwrap();
        let case = &cases[0];
        let mock = MockBenchmarkable::new();
        let n = case.matrix.nrows();
        let rhs: Vec<f64> = vec![1.0; n];

        let result = mock.bench_roundtrip(&case.matrix, &rhs);
        assert!(result.is_some(), "roundtrip should return Some");
        assert_eq!(result.unwrap().len(), n);
    }

    #[test]
    fn mock_roundtrip_returns_none_when_phase_disabled() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed()).unwrap();
        let case = &cases[0];
        let mock = MockBenchmarkable::with_phases_enabled(true, false, true);
        let n = case.matrix.nrows();
        let rhs: Vec<f64> = vec![1.0; n];

        let result = mock.bench_roundtrip(&case.matrix, &rhs);
        assert!(
            result.is_none(),
            "roundtrip should return None if factor disabled"
        );
    }

    #[test]
    fn supports_phase_reflects_config() {
        let mock = MockBenchmarkable::with_phases_enabled(true, false, true);
        assert!(mock.supports_phase(SolverPhase::Analyze));
        assert!(!mock.supports_phase(SolverPhase::Factor));
        assert!(mock.supports_phase(SolverPhase::Solve));
        assert!(!mock.supports_phase(SolverPhase::Roundtrip));
    }
}
