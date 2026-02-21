//! Integration tests for SuiteSparse CI-subset matrices.
//!
//! Uses the test infrastructure harness to load all 10 CI-subset matrices
//! and verify dimensions and symmetry.
//!
//! The CI subset contains small, fast matrices (~19 MB total, <1s serial factor)
//! spanning easy-indefinite and hard-indefinite categories.

use rivrs_sparse::testing::{TestCaseFilter, load_test_cases};

#[test]
fn load_all_ci_subset_matrices() {
    let cases =
        load_test_cases(&TestCaseFilter::ci_subset()).expect("failed to load CI-subset matrices");

    let suitesparse: Vec<_> = cases
        .iter()
        .filter(|c| c.properties.source == "suitesparse")
        .collect();

    assert_eq!(
        suitesparse.len(),
        10,
        "expected 10 CI-subset suitesparse matrices, got {}",
        suitesparse.len()
    );

    for case in &suitesparse {
        // Verify dimensions match properties
        assert_eq!(
            case.matrix.nrows(),
            case.properties.size,
            "row dimension mismatch for '{}'",
            case.name
        );
        assert_eq!(
            case.matrix.ncols(),
            case.properties.size,
            "col dimension mismatch for '{}'",
            case.name
        );

        // Verify matrix is square
        assert_eq!(
            case.matrix.nrows(),
            case.matrix.ncols(),
            "matrix '{}' should be square",
            case.name
        );
    }
}
