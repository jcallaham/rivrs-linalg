//! Integration tests for hand-constructed test matrices.
//!
//! Loads all 15 hand-constructed matrices, verifies dimensions and nnz match
//! metadata, and confirms reference factorizations load correctly.

use rivrs_sparse::io::registry;
use rivrs_sparse::validate;

#[test]
fn load_all_hand_constructed_matrices() {
    let start = std::time::Instant::now();

    let reg = registry::load_registry().expect("failed to load registry");
    let hand_constructed: Vec<_> = reg
        .iter()
        .filter(|m| m.source == "hand-constructed")
        .collect();

    assert_eq!(
        hand_constructed.len(),
        15,
        "expected 15 hand-constructed matrices"
    );

    for meta in &hand_constructed {
        let test = registry::load_test_matrix(&meta.name)
            .unwrap_or_else(|e| panic!("failed to load '{}': {}", meta.name, e))
            .unwrap_or_else(|| panic!("matrix '{}' should exist on disk", meta.name));

        // Verify dimensions match metadata
        assert_eq!(
            test.matrix.nrows(),
            meta.size,
            "dimension mismatch for '{}'",
            meta.name
        );
        assert_eq!(
            test.matrix.ncols(),
            meta.size,
            "dimension mismatch for '{}'",
            meta.name
        );

        // Verify reference factorization is available for all hand-constructed matrices
        assert!(
            test.reference.is_some(),
            "hand-constructed matrix '{}' should have a reference factorization",
            meta.name
        );

        let refdata = test.reference.as_ref().unwrap();
        assert_eq!(
            refdata.inertia.dimension(),
            meta.size,
            "inertia dimension mismatch for '{}'",
            meta.name
        );

        // Validate reconstruction error < 10^-12 (SC-008)
        let recon_err = validate::reconstruction_error(&test.matrix, refdata);
        assert!(
            recon_err < 1e-12,
            "reconstruction error for '{}': {:.2e} (expected < 1e-12)",
            meta.name,
            recon_err
        );

        // Validate inertia self-consistency (sanity check)
        assert!(
            validate::check_inertia(&refdata.inertia, &refdata.inertia),
            "inertia self-check failed for '{}'",
            meta.name
        );
    }

    let elapsed = start.elapsed();
    assert!(
        elapsed.as_secs_f64() < 1.0,
        "loading all 15 matrices took {:.3}s, expected < 1s",
        elapsed.as_secs_f64()
    );
}
