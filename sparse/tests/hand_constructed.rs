//! Integration tests for hand-constructed test matrices.
//!
//! Uses the test infrastructure harness to load all 15 hand-constructed
//! matrices and validate reconstruction error and inertia.

use rivrs_sparse::testing::{load_test_cases, NumericalValidator, TestCaseFilter};

#[test]
fn load_all_hand_constructed_matrices() {
    let start = std::time::Instant::now();

    let cases = load_test_cases(&TestCaseFilter::hand_constructed())
        .expect("failed to load hand-constructed test cases");

    assert_eq!(
        cases.len(),
        15,
        "expected 15 hand-constructed matrices"
    );

    let validator = NumericalValidator::new();

    for case in &cases {
        // Verify dimensions match properties
        assert_eq!(
            case.matrix.nrows(),
            case.properties.size,
            "dimension mismatch for '{}'",
            case.name
        );
        assert_eq!(
            case.matrix.ncols(),
            case.properties.size,
            "dimension mismatch for '{}'",
            case.name
        );

        // Verify reference factorization is available
        assert!(
            case.reference.is_some(),
            "hand-constructed matrix '{}' should have a reference factorization",
            case.name
        );

        let reference = case.reference.as_ref().unwrap();
        assert_eq!(
            reference.inertia.dimension(),
            case.properties.size,
            "inertia dimension mismatch for '{}'",
            case.name
        );

        // Validate reconstruction error < 10^-12
        let result = validator.validate_factorization(case);
        assert!(
            result.passed,
            "'{}' failed validation: {}",
            case.name, result
        );
    }

    let elapsed = start.elapsed();
    assert!(
        elapsed.as_secs_f64() < 1.0,
        "loading all 15 matrices took {:.3}s, expected < 1s",
        elapsed.as_secs_f64()
    );
}
