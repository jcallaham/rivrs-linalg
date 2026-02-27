//! Integration tests for hand-constructed test matrices.
//!
//! Uses the test infrastructure harness to load all 15 hand-constructed
//! matrices and validate reconstruction error and inertia.

use rivrs_sparse::debug::SparsityDisplay;
use rivrs_sparse::testing::{NumericalValidator, TestCaseFilter, load_test_cases};

#[test]
fn load_all_hand_constructed_matrices() {
    let start = std::time::Instant::now();

    let cases = load_test_cases(&TestCaseFilter::hand_constructed())
        .expect("failed to load hand-constructed test cases");

    assert_eq!(cases.len(), 15, "expected 15 hand-constructed matrices");

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

#[test]
fn sparsity_display_all_hand_constructed() {
    let cases = load_test_cases(&TestCaseFilter::hand_constructed())
        .expect("failed to load hand-constructed test cases");

    for case in &cases {
        let display = SparsityDisplay::from_sparse(&case.matrix);
        let rendered = display.render();

        // Every matrix should produce a non-empty visualization
        assert!(
            !rendered.is_empty(),
            "SparsityDisplay for '{}' should produce output",
            case.name
        );

        // Header should contain matrix dimensions
        let n = case.properties.size;
        let dim_str = format!("{n}x{n}");
        assert!(
            rendered.contains(&dim_str),
            "SparsityDisplay for '{}' should show dimensions '{dim_str}' in header, got:\n{rendered}",
            case.name
        );

        // Should contain nnz count
        assert!(
            rendered.contains("nnz="),
            "SparsityDisplay for '{}' should show nnz count",
            case.name
        );

        // Grid should have at least one line after the header
        let grid_lines: Vec<&str> = rendered.lines().skip(2).collect();
        assert!(
            !grid_lines.is_empty() || n == 0,
            "SparsityDisplay for '{}' should have grid lines",
            case.name
        );
    }
}
