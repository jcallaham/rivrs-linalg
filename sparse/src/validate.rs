//! Numerical validation utilities for solver correctness testing.
//!
//! Provides functions to verify LDL^T factorization correctness via:
//! - **Reconstruction error**: `||P^T A P - L D L^T||_F / ||A||_F`
//! - **Backward error**: `||Ax - b|| / (||A||_F ||x|| + ||b||)`
//! - **Inertia comparison**: field-wise equality of eigenvalue sign counts

use faer::linalg::matmul::matmul;
use faer::sparse::SparseColMat;
use faer::{Accum, Col, Mat, Par};

use crate::io::reference::{DBlock, Inertia, ReferenceFactorization};

/// Compute the relative reconstruction error of a factorization.
///
/// Returns `||P^T A P - L D L^T||_F / ||A||_F` where L, D, P come from
/// the reference factorization.
///
/// Returns 0.0 if `||A||_F` is zero.
pub fn reconstruction_error(
    a: &SparseColMat<usize, f64>,
    reference: &ReferenceFactorization,
) -> f64 {
    let n = a.nrows();
    let a_dense = a.to_dense();
    let a_norm = a_dense.norm_l2();
    if a_norm == 0.0 {
        return 0.0;
    }

    // Build P^T A P by reindexing: (PAP)_{i,j} = A_{p[i], p[j]}
    let perm = &reference.permutation;
    let pap = Mat::<f64>::from_fn(n, n, |i, j| a_dense[(perm[i], perm[j])]);

    // Build L: n×n identity + strict lower triangle entries
    let mut l_dense = Mat::<f64>::identity(n, n);
    for entry in &reference.l_entries {
        l_dense[(entry.row, entry.col)] = entry.value;
    }

    // Build D: n×n block diagonal
    let mut d_dense = Mat::<f64>::zeros(n, n);
    let mut row_offset = 0;
    for block in &reference.d_blocks {
        match block {
            DBlock::OneByOne { value } => {
                d_dense[(row_offset, row_offset)] = *value;
                row_offset += 1;
            }
            DBlock::TwoByTwo { values } => {
                d_dense[(row_offset, row_offset)] = values[0][0];
                d_dense[(row_offset, row_offset + 1)] = values[0][1];
                d_dense[(row_offset + 1, row_offset)] = values[1][0];
                d_dense[(row_offset + 1, row_offset + 1)] = values[1][1];
                row_offset += 2;
            }
        }
    }

    // Compute L * D
    let mut ld = Mat::<f64>::zeros(n, n);
    matmul(
        ld.as_mut(),
        Accum::Replace,
        l_dense.as_ref(),
        d_dense.as_ref(),
        1.0,
        Par::Seq,
    );

    // Compute L * D * L^T
    let mut ldlt = Mat::<f64>::zeros(n, n);
    matmul(
        ldlt.as_mut(),
        Accum::Replace,
        ld.as_ref(),
        l_dense.as_ref().transpose(),
        1.0,
        Par::Seq,
    );

    // Compute ||P^T A P - L D L^T||_F / ||A||_F
    let diff = &pap - &ldlt;
    diff.norm_l2() / a_norm
}

/// Compute the scaled backward error for a linear solve.
///
/// Returns `||Ax - b|| / (||A||_F * ||x|| + ||b||)`.
///
/// Returns 0.0 if denominator is zero.
pub fn backward_error(a: &SparseColMat<usize, f64>, x: &Col<f64>, b: &Col<f64>) -> f64 {
    let n = a.nrows();
    let a_dense = a.to_dense();

    // Compute r = A*x - b using dense matmul
    let x_mat = x.as_mat();
    let mut ax = Mat::<f64>::zeros(n, 1);
    matmul(
        ax.as_mut(),
        Accum::Replace,
        a_dense.as_ref(),
        x_mat,
        1.0,
        Par::Seq,
    );

    // r = Ax - b
    let r = Col::<f64>::from_fn(n, |i| ax[(i, 0)] - b[i]);

    let r_norm = r.norm_l2();
    let a_norm = a_dense.norm_l2();
    let x_norm = x.norm_l2();
    let b_norm = b.norm_l2();

    let denom = a_norm * x_norm + b_norm;
    if denom == 0.0 {
        return 0.0;
    }

    r_norm / denom
}

/// Check whether two inertia values are identical.
pub fn check_inertia(computed: &Inertia, expected: &Inertia) -> bool {
    computed.positive == expected.positive
        && computed.negative == expected.negative
        && computed.zero == expected.zero
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::registry;

    // Helper: compute b = A * x using dense conversion
    fn dense_matvec(a: &SparseColMat<usize, f64>, x: &Col<f64>) -> Col<f64> {
        let a_dense = a.to_dense();
        let n = a.nrows();
        let x_mat = x.as_mat();
        let mut result = Mat::<f64>::zeros(n, 1);
        matmul(
            result.as_mut(),
            Accum::Replace,
            a_dense.as_ref(),
            x_mat,
            1.0,
            Par::Seq,
        );
        Col::from_fn(n, |i| result[(i, 0)])
    }

    // T021: reconstruction_error tests

    #[test]
    fn reconstruction_error_arrow_5_pd() {
        let test = registry::load_test_matrix("arrow-5-pd")
            .expect("registry error")
            .expect("matrix should exist");
        let refdata = test.reference.expect("reference should exist");
        let err = reconstruction_error(&test.matrix, &refdata);
        assert!(
            err < 1e-12,
            "reconstruction error for arrow-5-pd: {:.2e} (expected < 1e-12)",
            err
        );
    }

    #[test]
    fn reconstruction_error_stress_delayed_pivots() {
        let test = registry::load_test_matrix("stress-delayed-pivots")
            .expect("registry error")
            .expect("matrix should exist");
        let refdata = test.reference.expect("reference should exist");
        let err = reconstruction_error(&test.matrix, &refdata);
        assert!(
            err < 1e-12,
            "reconstruction error for stress-delayed-pivots: {:.2e} (expected < 1e-12)",
            err
        );
    }

    #[test]
    fn reconstruction_error_with_perturbed_l() {
        let test = registry::load_test_matrix("arrow-5-pd")
            .expect("registry error")
            .expect("matrix should exist");
        let mut refdata = test.reference.expect("reference should exist");
        // Perturb an L entry
        if !refdata.l_entries.is_empty() {
            refdata.l_entries[0].value += 10.0;
        }
        let err = reconstruction_error(&test.matrix, &refdata);
        assert!(
            err > 0.01,
            "perturbed reconstruction error should be large: {:.2e}",
            err
        );
    }

    // T022: backward_error tests

    #[test]
    fn backward_error_exact_solution() {
        let test = registry::load_test_matrix("arrow-5-pd")
            .expect("registry error")
            .expect("matrix should exist");
        let a = &test.matrix;
        let n = a.nrows();

        // Create known x_exact
        let x_exact = Col::<f64>::from_fn(n, |i| (i + 1) as f64);

        // Compute b = A * x_exact
        let b = dense_matvec(a, &x_exact);

        let err = backward_error(a, &x_exact, &b);
        assert!(
            err < 1e-14,
            "backward error for exact solution: {:.2e} (expected < 1e-14)",
            err
        );
    }

    #[test]
    fn backward_error_perturbed_solution() {
        let test = registry::load_test_matrix("arrow-5-pd")
            .expect("registry error")
            .expect("matrix should exist");
        let a = &test.matrix;
        let n = a.nrows();

        // Create known x_exact and compute b
        let x_exact = Col::<f64>::from_fn(n, |i| (i + 1) as f64);
        let b = dense_matvec(a, &x_exact);

        // Perturb x
        let x_perturbed = Col::<f64>::from_fn(n, |i| (i + 1) as f64 + 0.1);
        let err = backward_error(a, &x_perturbed, &b);
        assert!(
            err > 1e-4,
            "backward error for perturbed solution should be measurable: {:.2e}",
            err
        );
    }

    // T023: check_inertia tests

    #[test]
    fn check_inertia_matching() {
        let a = Inertia {
            positive: 5,
            negative: 3,
            zero: 0,
        };
        let b = Inertia {
            positive: 5,
            negative: 3,
            zero: 0,
        };
        assert!(check_inertia(&a, &b));
    }

    #[test]
    fn check_inertia_mismatched() {
        let a = Inertia {
            positive: 5,
            negative: 3,
            zero: 0,
        };
        let b = Inertia {
            positive: 4,
            negative: 3,
            zero: 1,
        };
        assert!(!check_inertia(&a, &b));
    }
}
