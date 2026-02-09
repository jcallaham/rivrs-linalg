//! Integration tests for APTP data structures.
//!
//! Tests that APTP types correctly interoperate with the existing test matrix
//! infrastructure and produce results consistent with reference factorizations.

use rivrs_sparse::aptp::{Block2x2, MixedDiagonal, perm_from_forward};
use rivrs_sparse::io::reference::DBlock;
use rivrs_sparse::io::registry;

/// Helper: convert reference `DBlock` entries into `MixedDiagonal` calls.
///
/// The reference factorization format (`DBlock`) and the live factorization
/// format (`MixedDiagonal`) are intentionally separate types (FR-012). This
/// test helper bridges the two for integration testing.
fn build_mixed_diagonal_from_dblocks(d_blocks: &[DBlock], n: usize) -> MixedDiagonal {
    let mut diag = MixedDiagonal::new(n);
    let mut col = 0;
    for block in d_blocks {
        match block {
            DBlock::OneByOne { value } => {
                diag.set_1x1(col, *value);
                col += 1;
            }
            DBlock::TwoByTwo { values } => {
                diag.set_2x2(Block2x2 {
                    first_col: col,
                    a: values[0][0],
                    b: values[0][1],
                    c: values[1][1],
                });
                col += 2;
            }
        }
    }
    diag
}

// ---- T018: Integration test — compute_inertia against 15 hand-constructed matrices ----

#[test]
fn inertia_matches_all_15_hand_constructed_references() {
    let entries = registry::load_registry().expect("failed to load registry");
    let hand_constructed: Vec<_> = entries
        .iter()
        .filter(|e| e.category == "hand-constructed")
        .collect();

    assert!(
        hand_constructed.len() >= 15,
        "expected at least 15 hand-constructed matrices, found {}",
        hand_constructed.len()
    );

    for entry in &hand_constructed {
        let test = registry::load_test_matrix(&entry.name)
            .expect("registry error")
            .expect("matrix should exist");
        let reference = test
            .reference
            .as_ref()
            .unwrap_or_else(|| panic!("no reference for {}", entry.name));

        let n = reference.permutation.len();
        let diag = build_mixed_diagonal_from_dblocks(&reference.d_blocks, n);
        let computed = diag.compute_inertia();

        assert_eq!(
            computed, reference.inertia,
            "inertia mismatch for '{}': computed {:?} vs reference {:?}",
            entry.name, computed, reference.inertia
        );
    }
}

// ---- T022: Integration test — perm_from_forward against hand-constructed references ----

#[test]
fn perm_from_forward_matches_hand_constructed_references() {
    let entries = registry::load_registry().expect("failed to load registry");
    let hand_constructed: Vec<_> = entries
        .iter()
        .filter(|e| e.category == "hand-constructed")
        .collect();

    assert!(
        hand_constructed.len() >= 15,
        "expected at least 15 hand-constructed matrices, found {}",
        hand_constructed.len()
    );

    for entry in &hand_constructed {
        let test = registry::load_test_matrix(&entry.name)
            .expect("registry error")
            .expect("matrix should exist");
        let reference = test
            .reference
            .as_ref()
            .unwrap_or_else(|| panic!("no reference for {}", entry.name));

        let perm = perm_from_forward(reference.permutation.clone())
            .unwrap_or_else(|e| panic!("perm_from_forward failed for '{}': {}", entry.name, e));
        let (fwd_arr, inv_arr) = perm.as_ref().arrays();

        // Verify forward array matches
        assert_eq!(
            fwd_arr,
            reference.permutation.as_slice(),
            "forward mismatch for '{}'",
            entry.name
        );

        // Verify forward/inverse relationship
        for i in 0..fwd_arr.len() {
            assert_eq!(
                inv_arr[fwd_arr[i]], i,
                "inv[fwd[{}]] != {} for '{}'",
                i, i, entry.name
            );
            assert_eq!(
                fwd_arr[inv_arr[i]], i,
                "fwd[inv[{}]] != {} for '{}'",
                i, i, entry.name
            );
        }
    }
}
