//! Shared validation helpers for MC64 matching & scaling tests.
//!
//! Used by both unit tests in `src/aptp/matching.rs` and integration tests
//! in `tests/mc64_matching.rs` to avoid duplicated verification logic.

use faer::sparse::SparseColMat;

use crate::aptp::Mc64Result;
use crate::aptp::matching::count_cycles;

/// Categorized scaling violation report.
///
/// Scaling violations (`|s_i * a_ij * s_j| > 1`) are classified by whether
/// the row and column are matched in the Hungarian algorithm. This matters
/// because:
///
/// - **Matched-matched**: Dual feasibility guarantees `|s_i * a_ij * s_j| <= 1`.
///   Violations indicate a bug in the Hungarian algorithm.
/// - **Matched-unmatched**: Duff-Pralet correction guarantees bounded products
///   with matched neighbors. Violations indicate a bug in Duff-Pralet.
/// - **Unmatched-unmatched**: No coupling between independent Duff-Pralet
///   corrections. Violations are an inherent mathematical limitation, not a bug.
///   SPRAL has the same behavior.
pub struct ScalingViolationReport {
    /// Max `|s_i * a_ij * s_j| - 1` for edges where both i and j are matched.
    /// Should be < 1e-10 (dual feasibility guarantee).
    pub max_matched_matched: f64,
    /// Max `|s_i * a_ij * s_j| - 1` for edges where exactly one of i, j is matched.
    /// Should be < 1e-10 (Duff-Pralet guarantee).
    pub max_matched_unmatched: f64,
    /// Max `|s_i * a_ij * s_j| - 1` for edges where neither i nor j is matched.
    /// May exceed 0 (inherent limitation).
    pub max_unmatched_unmatched: f64,
    /// Number of unmatched-unmatched edges with `|s_i * a_ij * s_j| > 1 + 1e-10`.
    pub unmatched_unmatched_violation_count: usize,
    /// Total number of unmatched-unmatched edges.
    pub unmatched_unmatched_edge_count: usize,
}

/// Classify scaling violations by edge type (matched/unmatched endpoints).
///
/// Iterates all stored entries in `matrix`, computes `|s_i * a_ij * s_j|`,
/// and categorizes violations by whether the endpoints are matched according
/// to `result.is_matched`.
pub fn classify_scaling_violations(
    matrix: &SparseColMat<usize, f64>,
    result: &Mc64Result,
) -> ScalingViolationReport {
    let symbolic = matrix.symbolic();
    let values = matrix.val();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let n = matrix.nrows();

    let mut report = ScalingViolationReport {
        max_matched_matched: 0.0,
        max_matched_unmatched: 0.0,
        max_unmatched_unmatched: 0.0,
        unmatched_unmatched_violation_count: 0,
        unmatched_unmatched_edge_count: 0,
    };

    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for k in start..end {
            let i = row_indices[k];
            let scaled = (result.scaling[i] * values[k] * result.scaling[j]).abs();
            let violation = scaled - 1.0;

            let i_matched = result.is_matched[i];
            let j_matched = result.is_matched[j];

            match (i_matched, j_matched) {
                (true, true) => {
                    if violation > report.max_matched_matched {
                        report.max_matched_matched = violation;
                    }
                }
                (true, false) | (false, true) => {
                    if violation > report.max_matched_unmatched {
                        report.max_matched_unmatched = violation;
                    }
                }
                (false, false) => {
                    report.unmatched_unmatched_edge_count += 1;
                    if violation > report.max_unmatched_unmatched {
                        report.max_unmatched_unmatched = violation;
                    }
                    if violation > 1e-10 {
                        report.unmatched_unmatched_violation_count += 1;
                    }
                }
            }
        }
    }

    report
}

/// Verify SPRAL-style scaling properties for an MC64 matching result.
///
/// Checks:
/// 1. Scaling bound (categorized): matched-matched and matched-unmatched edges
///    must satisfy `|s_i * a_ij * s_j| <= 1.0 + tol`. Unmatched-unmatched
///    violations are reported but not asserted (inherent limitation).
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

    // Property 1: Categorized scaling bound
    let report = classify_scaling_violations(matrix, result);

    assert!(
        report.max_matched_matched <= 1e-10,
        "{}: matched-matched violation = {:.6e} > 1e-10 (dual feasibility broken)",
        name,
        report.max_matched_matched,
    );
    assert!(
        report.max_matched_unmatched <= 1e-10,
        "{}: matched-unmatched violation = {:.6e} > 1e-10 (Duff-Pralet broken)",
        name,
        report.max_matched_unmatched,
    );
    if report.max_unmatched_unmatched > 1e-10 {
        eprintln!(
            "  {}: unmatched-unmatched max violation = {:.2e} ({}/{} edges) — inherent limitation",
            name,
            report.max_unmatched_unmatched,
            report.unmatched_unmatched_violation_count,
            report.unmatched_unmatched_edge_count,
        );
    }

    // Property 2: row max >= 1.0 - epsilon for nonzero rows
    let mut row_max = vec![0.0_f64; n];
    for j in 0..n {
        let start = col_ptrs[j];
        let end = col_ptrs[j + 1];
        for k in start..end {
            let i = row_indices[k];
            let scaled = (result.scaling[i] * values[k] * result.scaling[j]).abs();
            if scaled > row_max[i] {
                row_max[i] = scaled;
            }
            if i != j && scaled > row_max[j] {
                row_max[j] = scaled;
            }
        }
    }

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
