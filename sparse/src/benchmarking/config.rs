//! Benchmark configuration.

use std::time::Duration;

use crate::SolverPhase;
use crate::testing::TestCaseFilter;

/// Configuration for a benchmark run.
#[derive(Debug, Clone)]
pub struct BenchmarkConfig {
    pub filter: TestCaseFilter,
    pub phases: Vec<SolverPhase>,
    pub sample_size: Option<usize>,
    pub measurement_time: Option<Duration>,
    pub warm_up_time: Option<Duration>,
}

impl BenchmarkConfig {
    /// Create a config that benchmarks hand-constructed matrices across all component phases.
    pub fn default_components() -> Self {
        Self {
            filter: TestCaseFilter::hand_constructed(),
            phases: SolverPhase::components().to_vec(),
            sample_size: None,
            measurement_time: None,
            warm_up_time: None,
        }
    }

    pub fn with_filter(mut self, filter: TestCaseFilter) -> Self {
        self.filter = filter;
        self
    }

    pub fn with_phases(mut self, phases: Vec<SolverPhase>) -> Self {
        self.phases = phases;
        self
    }

    pub fn with_sample_size(mut self, n: usize) -> Self {
        self.sample_size = Some(n);
        self
    }

    pub fn with_measurement_time(mut self, d: Duration) -> Self {
        self.measurement_time = Some(d);
        self
    }

    pub fn with_warm_up_time(mut self, d: Duration) -> Self {
        self.warm_up_time = Some(d);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn config_builder() {
        let config = BenchmarkConfig::default_components()
            .with_sample_size(50)
            .with_measurement_time(Duration::from_secs(10))
            .with_warm_up_time(Duration::from_secs(1));

        assert_eq!(config.sample_size, Some(50));
        assert_eq!(config.measurement_time, Some(Duration::from_secs(10)));
        assert_eq!(config.warm_up_time, Some(Duration::from_secs(1)));
        assert_eq!(config.phases.len(), 3);
    }

    #[test]
    fn components_excludes_roundtrip() {
        let components = SolverPhase::components();
        assert!(!components.contains(&SolverPhase::Roundtrip));
        assert_eq!(components.len(), 3);
    }

    #[test]
    fn all_includes_roundtrip() {
        let all = SolverPhase::all();
        assert!(all.contains(&SolverPhase::Roundtrip));
        assert_eq!(all.len(), 4);
    }
}
