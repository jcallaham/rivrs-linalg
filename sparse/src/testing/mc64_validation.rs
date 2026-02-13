//! Shared validation helpers for MC64 matching & scaling tests.
//!
//! Used by both unit tests in `src/aptp/matching.rs` and integration tests
//! in `tests/mc64_matching.rs` to avoid duplicated verification logic.

use faer::sparse::SparseColMat;

use crate::aptp::Mc64Result;
use crate::aptp::matching::count_cycles;

/// Verify SPRAL-style scaling properties for an MC64 matching result.
///
/// Checks:
/// 1. Scaling bound: `|s_i * a_ij * s_j| <= 1.0 + tol` for all entries
/// 2. Row max quality: `max_j |s_i * a_ij * s_j| >= 1.0 - tol` for nonzero rows
/// 3. Matching is a valid permutation
/// 4. Scaling factors are positive and finite
/// 5. Cycle decomposition is consistent
///
/// Panics with a descriptive message if any check fails.
pub fn verify_spral_scaling_properties(
    name: &str,
    matrix: &SparseColMat<usize, f64>,
    result: &Mc64Result,
) {
    let n = matrix.nrows();
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let (fwd, _) = result.matching.as_ref().arrays();

    // Property 1: |s_i * a_ij * s_j| <= 1.0 for all stored entries
    let mut row_max = vec![0.0_f64; n];

    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for k in start..end {
            let i = row_indices[k];
            let scaled = (result.scaling[i] * values[k] * result.scaling[j]).abs();
            assert!(
                scaled <= 1.0 + 1e-10,
                "{}: |s[{}]*a[{},{}]*s[{}]| = {:.6e} > 1.0",
                name,
                i,
                i,
                j,
                j,
                scaled
            );
            if scaled > row_max[i] {
                row_max[i] = scaled;
            }
            if i != j && scaled > row_max[j] {
                row_max[j] = scaled;
            }
        }
    }

    // Property 2: row max >= 1.0 - epsilon for nonzero rows
    for (i, &rm) in row_max.iter().enumerate() {
        if rm > 0.0 {
            assert!(
                rm >= 1.0 - 1e-12,
                "{}: row_max[{}] = {:.6e} < 1.0 - 1e-12 (SPRAL expects ~1.0)",
                name,
                i,
                rm
            );
        }
    }

    // Property 3: matching is valid permutation
    let mut seen = vec![false; n];
    for i in 0..n {
        assert!(
            fwd[i] < n,
            "{}: matching[{}] = {} out of range",
            name,
            i,
            fwd[i]
        );
        assert!(
            !seen[fwd[i]],
            "{}: duplicate in matching at {}",
            name, fwd[i]
        );
        seen[fwd[i]] = true;
    }

    // Property 4: scaling factors positive and finite
    for (i, &s) in result.scaling.iter().enumerate() {
        assert!(s > 0.0, "{}: scaling[{}] = {} not positive", name, i, s);
        assert!(s.is_finite(), "{}: scaling[{}] = {} not finite", name, i, s);
    }

    // Property 5: cycle decomposition is consistent
    let (singletons, two_cycles, longer_cycles) = count_cycles(fwd);
    assert!(
        singletons + 2 * two_cycles <= n,
        "{}: cycle counts exceed dimension: {} singletons + {} 2-cycles > {}",
        name,
        singletons,
        two_cycles,
        n
    );
    if longer_cycles > 0 {
        eprintln!(
            "  {}: {} longer cycles (expected for asymmetric cost graph)",
            name, longer_cycles
        );
    }
}
