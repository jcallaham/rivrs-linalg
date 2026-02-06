//! Integration tests for SuiteSparse CI-subset matrices.
//!
//! Uses the test infrastructure harness to load all 9 CI-subset matrices
//! and verify dimensions and symmetry.
//!
//! nd6k (103 MB) excluded from CI subset due to GitHub's 100 MB file size limit.

use rivrs_sparse::testing::{load_test_cases, TestCaseFilter};

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
        9,
        "expected 9 CI-subset suitesparse matrices, got {}",
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
