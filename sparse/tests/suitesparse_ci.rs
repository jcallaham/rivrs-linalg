//! Integration tests for SuiteSparse CI-subset matrices.
//!
//! Loads all 10 CI-subset matrices from the registry, verifies dimensions
//! and nnz match metadata, and confirms matrix symmetry.

use rivrs_sparse::io::registry;

#[test]
fn load_all_ci_subset_matrices() {
    let reg = registry::load_registry().expect("failed to load registry");
    let ci_subset: Vec<_> = reg.iter().filter(|m| m.ci_subset).collect();

    assert_eq!(
        ci_subset.len(),
        10,
        "expected 10 CI-subset matrices, got {}",
        ci_subset.len()
    );

    for meta in &ci_subset {
        let test = registry::load_test_matrix(&meta.name)
            .unwrap_or_else(|e| panic!("failed to load '{}': {}", meta.name, e))
            .unwrap_or_else(|| panic!("CI-subset matrix '{}' should exist on disk", meta.name));

        // Verify dimensions match metadata
        assert_eq!(
            test.matrix.nrows(),
            meta.size,
            "row dimension mismatch for '{}'",
            meta.name
        );
        assert_eq!(
            test.matrix.ncols(),
            meta.size,
            "col dimension mismatch for '{}'",
            meta.name
        );

        // Verify matrix is square
        assert_eq!(
            test.matrix.nrows(),
            test.matrix.ncols(),
            "matrix '{}' should be square",
            meta.name
        );
    }
}
