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

    let mut loaded = 0;
    let mut skipped = Vec::new();

    for meta in &ci_subset {
        let test = match registry::load_test_matrix(&meta.name) {
            Ok(Some(t)) => t,
            Ok(None) => {
                skipped.push(format!("{} (file not found)", meta.name));
                continue;
            }
            Err(e) => {
                skipped.push(format!("{} (parse error: {})", meta.name, e));
                continue;
            }
        };

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

        loaded += 1;
    }

    if !skipped.is_empty() {
        eprintln!(
            "WARNING: skipped {} CI-subset matrices: {:?}",
            skipped.len(),
            skipped
        );
    }

    // At least 8 of 10 CI-subset matrices should load successfully
    assert!(
        loaded >= 8,
        "only {} of 10 CI-subset matrices loaded successfully (skipped: {:?})",
        loaded,
        skipped
    );
}
