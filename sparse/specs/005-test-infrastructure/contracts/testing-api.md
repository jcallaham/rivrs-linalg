# API Contract: Testing Module

**Branch**: `005-test-infrastructure` | **Date**: 2026-02-06
**Module**: `src/testing.rs` (gated behind `#[cfg(feature = "test-util")]`)

## Module Structure

```
src/
├── testing.rs          # Module root (re-exports)
├── testing/
│   ├── mod.rs          # Module declarations
│   ├── harness.rs      # SolverTest trait, TestResult, TestKind
│   ├── validator.rs    # NumericalValidator, MetricResult
│   ├── cases.rs        # SolverTestCase, TestMatrixProperties, TestCaseFilter, loading
│   └── generators.rs   # Random matrix generators, pattern generators
```

Or as a single flat module if the total LOC is small enough:
```
src/
├── testing.rs          # All test infrastructure types and functions
```

## SolverTest Trait

```rust
/// Interface for testing solver implementations across all phases.
///
/// Each method validates a specific solver phase against a test case.
/// Implementors provide a solver and the harness validates its output.
pub trait SolverTest {
    /// Validate symbolic analysis output.
    ///
    /// Checks: valid permutation (bijection on 0..n), elimination tree
    /// structure consistency. Does not require factorization.
    fn test_analyze(&self, case: &SolverTestCase) -> TestResult;

    /// Validate numeric factorization output.
    ///
    /// Checks: reconstruction error (||P^T A P - LDL^T|| / ||A||),
    /// inertia correctness (if reference available).
    fn test_factor(&self, case: &SolverTestCase) -> TestResult;

    /// Validate triangular solve output.
    ///
    /// Checks: backward error (||Ax - b|| / (||A|| ||x|| + ||b||))
    /// for a generated right-hand side.
    fn test_solve(&self, case: &SolverTestCase) -> TestResult;

    /// Validate full analyze → factor → solve pipeline.
    ///
    /// Runs the complete pipeline as one operation, then checks
    /// reconstruction error, backward error, and inertia.
    /// Independent of test_analyze/test_factor/test_solve.
    fn test_roundtrip(&self, case: &SolverTestCase) -> TestResult;
}
```

## NumericalValidator

```rust
/// Configurable validator for factorization quality.
///
/// Wraps the standalone functions in `validate` with tolerance thresholds
/// and structured result reporting.
pub struct NumericalValidator {
    reconstruction_tol: f64,  // default: 1e-12
    backward_error_tol: f64,  // default: 1e-10
}

impl NumericalValidator {
    /// Create validator with constitution-mandated default tolerances.
    pub fn new() -> Self;

    /// Override reconstruction error tolerance.
    pub fn with_reconstruction_tol(self, tol: f64) -> Self;

    /// Override backward error tolerance.
    pub fn with_backward_error_tol(self, tol: f64) -> Self;

    /// Check reconstruction error against threshold.
    pub fn check_reconstruction(
        &self,
        matrix: &SparseColMat<usize, f64>,
        reference: &ReferenceFactorization,
    ) -> MetricResult;

    /// Check backward error against threshold.
    pub fn check_backward_error(
        &self,
        matrix: &SparseColMat<usize, f64>,
        x: &Col<f64>,
        b: &Col<f64>,
    ) -> MetricResult;

    /// Check inertia equality.
    pub fn check_inertia(
        &self,
        computed: &Inertia,
        expected: &Inertia,
    ) -> MetricResult;

    /// Run all applicable checks for a test case with reference data.
    pub fn validate_factorization(
        &self,
        case: &SolverTestCase,
    ) -> TestResult;
}

impl Default for NumericalValidator { /* same as new() */ }
```

## Result Types

```rust
/// Outcome of a single numerical metric check.
pub struct MetricResult {
    pub name: String,       // e.g., "reconstruction_error"
    pub value: f64,         // computed metric
    pub threshold: f64,     // pass/fail boundary
    pub passed: bool,       // value < threshold (or exact match for inertia)
}

/// Which solver phase was tested.
pub enum TestKind {
    Analyze,
    Factor,
    Solve,
    Roundtrip,
}

/// Complete outcome of a solver test.
pub struct TestResult {
    pub passed: bool,              // all metrics passed
    pub test_kind: TestKind,
    pub matrix_name: String,
    pub metrics: Vec<MetricResult>,
    pub diagnostics: Vec<String>,  // human-readable messages
}
```

## Test Case Types

```rust
/// A complete test scenario for solver validation.
pub struct SolverTestCase {
    pub name: String,
    pub matrix: SparseColMat<usize, f64>,
    pub properties: TestMatrixProperties,
    pub reference: Option<ReferenceFactorization>,
}

/// Structural and numerical properties of a test matrix.
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

/// Criteria for loading subsets of test matrices from the registry.
pub struct TestCaseFilter {
    pub source: Option<String>,
    pub category: Option<String>,
    pub difficulty: Option<String>,
    pub ci_only: bool,
    pub require_reference: bool,
}

impl TestCaseFilter {
    pub fn all() -> Self;
    pub fn hand_constructed() -> Self;
    pub fn ci_subset() -> Self;
    pub fn with_source(self, source: &str) -> Self;
    pub fn with_category(self, category: &str) -> Self;
    pub fn with_difficulty(self, difficulty: &str) -> Self;
    pub fn ci_only(self) -> Self;
    pub fn require_reference(self) -> Self;
}

/// Load test cases matching filter criteria.
///
/// Returns Ok(cases) where cases may be empty if no matrices match.
/// Matrices whose .mtx file is missing on disk are silently skipped
/// (returns None for that entry).
pub fn load_test_cases(filter: &TestCaseFilter) -> Result<Vec<SolverTestCase>, SparseError>;
```

## Matrix Generators

```rust
/// Configuration for random sparse symmetric matrix generation.
pub struct RandomMatrixConfig {
    pub size: usize,
    pub target_nnz: usize,          // approximate; actual may differ by ±20%
    pub positive_definite: bool,     // if true, guarantees PD via diagonal dominance
}

/// Generate a random sparse symmetric matrix.
///
/// If positive_definite is true, uses diagonal dominance to guarantee PD.
/// If false, generates an indefinite matrix with mixed diagonal signs.
pub fn generate_random_symmetric(
    config: &RandomMatrixConfig,
    rng: &mut impl Rng,
) -> Result<SparseColMat<usize, f64>, SparseError>;

/// Generate a sparse symmetric arrow matrix of given size.
///
/// Arrow pattern: dense first row/column, diagonal remainder.
/// If positive_definite is true, diagonal dominance is applied.
pub fn generate_arrow(
    size: usize,
    positive_definite: bool,
    rng: &mut impl Rng,
) -> SparseColMat<usize, f64>;

/// Generate a sparse symmetric tridiagonal matrix of given size.
pub fn generate_tridiagonal(
    size: usize,
    positive_definite: bool,
    rng: &mut impl Rng,
) -> SparseColMat<usize, f64>;

/// Generate a sparse symmetric banded matrix of given size and bandwidth.
pub fn generate_banded(
    size: usize,
    bandwidth: usize,
    positive_definite: bool,
    rng: &mut impl Rng,
) -> SparseColMat<usize, f64>;
```
