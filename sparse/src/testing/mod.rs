//! Reusable test infrastructure for validating SSIDS solver components.
//!
//! This module provides four subsystems:
//!
//! - **[`cases`]**: Unified test case loading from the matrix registry with filtering.
//! - **[`harness`]**: `SolverTest` trait and result types for per-phase validation.
//! - **[`validator`]**: `NumericalValidator` with configurable tolerances.
//! - **[`generators`]**: Random and structured sparse symmetric matrix generators.
//!
//! Gated behind the `test-util` Cargo feature flag. Not included in production builds.
//!
//! # Quick Start
//!
//! ```rust,no_run
//! use rivrs_sparse::testing::{load_test_cases, NumericalValidator, TestCaseFilter};
//!
//! let cases = load_test_cases(&TestCaseFilter::hand_constructed()).unwrap();
//! let validator = NumericalValidator::new();
//! for case in &cases {
//!     let result = validator.validate_factorization(case);
//!     assert!(result.passed);
//! }
//! ```

pub mod cases;
pub mod generators;
pub mod harness;
pub mod mc64_validation;
pub mod validator;

pub use cases::{SolverTestCase, TestCaseFilter, TestMatrixProperties, load_test_cases};
pub use generators::{
    RandomMatrixConfig, generate_arrow, generate_banded, generate_random_symmetric,
    generate_tridiagonal,
};
pub use harness::{MetricResult, MockSolver, SolverTest, TestResult};
pub use mc64_validation::verify_spral_scaling_properties;
pub use validator::NumericalValidator;
