//! SolverTestCase, TestMatrixProperties, TestCaseFilter, and loading functions.

use faer::sparse::SparseColMat;

use crate::error::SparseError;
use crate::io::reference::ReferenceFactorization;
use crate::io::registry;

/// Structural and numerical properties of a test matrix.
///
/// Flattened view of [`registry::MatrixMetadata`] + [`registry::MatrixProperties`]
/// for convenient access in test code.
#[derive(Debug, Clone)]
pub struct TestMatrixProperties {
    pub size: usize,
    pub nnz: usize,
    pub symmetric: bool,
    pub positive_definite: bool,
    pub indefinite: bool,
    pub difficulty: String,
    pub structure: Option<String>,
    pub source: String,
    pub category: String,
}

impl From<&registry::MatrixMetadata> for TestMatrixProperties {
    fn from(meta: &registry::MatrixMetadata) -> Self {
        Self {
            size: meta.size,
            nnz: meta.nnz,
            symmetric: meta.properties.symmetric,
            positive_definite: meta.properties.positive_definite,
            indefinite: meta.properties.indefinite,
            difficulty: meta.properties.difficulty.clone(),
            structure: meta.properties.structure.clone(),
            source: meta.source.clone(),
            category: meta.category.clone(),
        }
    }
}

/// A complete test scenario for solver validation.
#[derive(Debug, Clone)]
pub struct SolverTestCase {
    pub name: String,
    pub matrix: SparseColMat<usize, f64>,
    pub properties: TestMatrixProperties,
    pub reference: Option<ReferenceFactorization>,
}

/// Criteria for loading subsets of test matrices from the registry.
#[derive(Debug, Clone)]
pub struct TestCaseFilter {
    pub source: Option<String>,
    pub category: Option<String>,
    pub difficulty: Option<String>,
    pub ci_only: bool,
    pub require_reference: bool,
}

impl TestCaseFilter {
    /// No filtering — return all matrices.
    pub fn all() -> Self {
        Self {
            source: None,
            category: None,
            difficulty: None,
            ci_only: false,
            require_reference: false,
        }
    }

    /// Only hand-constructed matrices (with reference factorizations).
    pub fn hand_constructed() -> Self {
        Self {
            source: Some("hand-constructed".to_string()),
            category: None,
            difficulty: None,
            ci_only: false,
            require_reference: false,
        }
    }

    /// Only CI-subset SuiteSparse matrices.
    pub fn ci_subset() -> Self {
        Self {
            source: None,
            category: None,
            difficulty: None,
            ci_only: true,
            require_reference: false,
        }
    }

    /// Filter by source.
    pub fn with_source(mut self, source: &str) -> Self {
        self.source = Some(source.to_string());
        self
    }

    /// Filter by category.
    pub fn with_category(mut self, category: &str) -> Self {
        self.category = Some(category.to_string());
        self
    }

    /// Filter by difficulty.
    pub fn with_difficulty(mut self, difficulty: &str) -> Self {
        self.difficulty = Some(difficulty.to_string());
        self
    }

    /// Only CI-subset matrices.
    pub fn ci_only(mut self) -> Self {
        self.ci_only = true;
        self
    }

    /// Only matrices with reference factorizations.
    pub fn require_reference(mut self) -> Self {
        self.require_reference = true;
        self
    }

    fn matches(&self, meta: &registry::MatrixMetadata) -> bool {
        if let Some(ref src) = self.source {
            if meta.source != *src {
                return false;
            }
        }
        if let Some(ref cat) = self.category {
            if meta.category != *cat {
                return false;
            }
        }
        if let Some(ref diff) = self.difficulty {
            if meta.properties.difficulty != *diff {
                return false;
            }
        }
        if self.ci_only && !meta.ci_subset {
            // Hand-constructed matrices are always "in CI" (committed to git)
            if meta.source != "hand-constructed" {
                return false;
            }
        }
        if self.require_reference && meta.factorization_path.is_none() {
            return false;
        }
        true
    }
}

/// Load test cases matching filter criteria.
///
/// Loads the registry once and reuses it for all matching entries, avoiding
/// redundant metadata.json parsing. Matrices whose .mtx file is missing on
/// disk are silently skipped.
pub fn load_test_cases(filter: &TestCaseFilter) -> Result<Vec<SolverTestCase>, SparseError> {
    let all_meta = registry::load_registry()?;
    let mut cases = Vec::new();

    for meta in &all_meta {
        if !filter.matches(meta) {
            continue;
        }

        // Load directly from entry to avoid re-parsing metadata.json per matrix
        match registry::load_test_matrix_from_entry(meta)? {
            Some(test_matrix) => {
                let properties = TestMatrixProperties::from(meta);
                let case = SolverTestCase {
                    name: meta.name.clone(),
                    matrix: test_matrix.matrix,
                    properties,
                    reference: test_matrix.reference,
                };

                // Post-filter: require_reference checks actual loaded data
                if filter.require_reference && case.reference.is_none() {
                    continue;
                }

                cases.push(case);
            }
            None => {
                // .mtx file not on disk — silently skip
                continue;
            }
        }
    }

    Ok(cases)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn load_hand_constructed_returns_15() {
        let cases = load_test_cases(&TestCaseFilter::hand_constructed())
            .expect("failed to load test cases");
        assert_eq!(cases.len(), 15, "expected 15 hand-constructed cases");
        for case in &cases {
            assert!(
                case.reference.is_some(),
                "{} should have a reference factorization",
                case.name
            );
            assert_eq!(
                case.properties.source, "hand-constructed",
                "{} should have source 'hand-constructed'",
                case.name
            );
        }
    }

    #[test]
    fn load_ci_subset_returns_10() {
        let cases =
            load_test_cases(&TestCaseFilter::ci_subset()).expect("failed to load test cases");
        // CI subset = 10 suitesparse matrices + 15 hand-constructed (always available)
        // But the filter says ci_only, so only ci_subset=true suitesparse + hand-constructed
        let suitesparse_cases: Vec<_> = cases
            .iter()
            .filter(|c| c.properties.source == "suitesparse")
            .collect();
        assert_eq!(
            suitesparse_cases.len(),
            10,
            "expected 10 CI-subset suitesparse cases"
        );
    }

    #[test]
    fn filter_by_category() {
        let cases = load_test_cases(&TestCaseFilter::all().with_category("hand-constructed"))
            .expect("failed to load");
        for case in &cases {
            assert_eq!(case.properties.category, "hand-constructed");
        }
        assert!(!cases.is_empty());

        // Filter by hard-indefinite
        let hard = load_test_cases(&TestCaseFilter::all().with_category("hard-indefinite"))
            .expect("failed to load");
        for case in &hard {
            assert_eq!(case.properties.category, "hard-indefinite");
            assert!(
                case.properties.indefinite,
                "{} should be indefinite",
                case.name
            );
        }
    }

    #[test]
    fn require_reference_filters_correctly() {
        let cases =
            load_test_cases(&TestCaseFilter::all().require_reference()).expect("failed to load");
        for case in &cases {
            assert!(
                case.reference.is_some(),
                "{} should have a reference",
                case.name
            );
        }
    }

    #[test]
    #[ignore] // Wall-clock assertion may flake in constrained CI environments
    fn load_performance_under_100ms() {
        let start = std::time::Instant::now();
        let cases =
            load_test_cases(&TestCaseFilter::hand_constructed().with_category("hand-constructed"))
                .expect("failed to load");
        let elapsed = start.elapsed();
        let per_matrix = elapsed / cases.len() as u32;
        assert!(
            per_matrix.as_millis() < 100,
            "per-matrix load time: {:?} (expected < 100ms)",
            per_matrix
        );
    }
}
