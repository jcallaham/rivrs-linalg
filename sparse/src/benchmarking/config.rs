//! Benchmark configuration and phase definitions.

use std::fmt;
use std::time::Duration;

use serde::{Deserialize, Serialize};

use crate::testing::TestCaseFilter;

/// Measurable solver operations.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BenchmarkPhase {
    Analyze,
    Factor,
    Solve,
    Roundtrip,
}

impl fmt::Display for BenchmarkPhase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Analyze => write!(f, "analyze"),
            Self::Factor => write!(f, "factor"),
            Self::Solve => write!(f, "solve"),
            Self::Roundtrip => write!(f, "roundtrip"),
        }
    }
}

impl BenchmarkPhase {
    /// All individual component phases (excludes Roundtrip).
    pub fn components() -> Vec<Self> {
        vec![Self::Analyze, Self::Factor, Self::Solve]
    }

    /// All phases including Roundtrip.
    pub fn all() -> Vec<Self> {
        vec![Self::Analyze, Self::Factor, Self::Solve, Self::Roundtrip]
    }
}

/// Configuration for a benchmark run.
#[derive(Debug, Clone)]
pub struct BenchmarkConfig {
    pub filter: TestCaseFilter,
    pub phases: Vec<BenchmarkPhase>,
    pub sample_size: Option<usize>,
    pub measurement_time: Option<Duration>,
    pub warm_up_time: Option<Duration>,
    pub timeout_per_matrix: Option<Duration>,
}

impl BenchmarkConfig {
    /// Create a config that benchmarks hand-constructed matrices across all component phases.
    pub fn default_components() -> Self {
        Self {
            filter: TestCaseFilter::hand_constructed(),
            phases: BenchmarkPhase::components(),
            sample_size: None,
            measurement_time: None,
            warm_up_time: None,
            timeout_per_matrix: None,
        }
    }

    pub fn with_filter(mut self, filter: TestCaseFilter) -> Self {
        self.filter = filter;
        self
    }

    pub fn with_phases(mut self, phases: Vec<BenchmarkPhase>) -> Self {
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

    pub fn with_timeout_per_matrix(mut self, d: Duration) -> Self {
        self.timeout_per_matrix = Some(d);
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn phase_display() {
        assert_eq!(BenchmarkPhase::Analyze.to_string(), "analyze");
        assert_eq!(BenchmarkPhase::Factor.to_string(), "factor");
        assert_eq!(BenchmarkPhase::Solve.to_string(), "solve");
        assert_eq!(BenchmarkPhase::Roundtrip.to_string(), "roundtrip");
    }

    #[test]
    fn phase_serde_roundtrip() {
        for phase in BenchmarkPhase::all() {
            let json = serde_json::to_string(&phase).unwrap();
            let back: BenchmarkPhase = serde_json::from_str(&json).unwrap();
            assert_eq!(phase, back);
        }
    }

    #[test]
    fn config_builder() {
        let config = BenchmarkConfig::default_components()
            .with_sample_size(50)
            .with_measurement_time(Duration::from_secs(10))
            .with_warm_up_time(Duration::from_secs(1))
            .with_timeout_per_matrix(Duration::from_secs(60));

        assert_eq!(config.sample_size, Some(50));
        assert_eq!(config.measurement_time, Some(Duration::from_secs(10)));
        assert_eq!(config.warm_up_time, Some(Duration::from_secs(1)));
        assert_eq!(config.timeout_per_matrix, Some(Duration::from_secs(60)));
        assert_eq!(config.phases.len(), 3);
    }

    #[test]
    fn components_excludes_roundtrip() {
        let components = BenchmarkPhase::components();
        assert!(!components.contains(&BenchmarkPhase::Roundtrip));
        assert_eq!(components.len(), 3);
    }

    #[test]
    fn all_includes_roundtrip() {
        let all = BenchmarkPhase::all();
        assert!(all.contains(&BenchmarkPhase::Roundtrip));
        assert_eq!(all.len(), 4);
    }
}
