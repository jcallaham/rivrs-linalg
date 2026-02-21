//! Numerical validation utilities for solver correctness testing.
//!
//! Provides functions to verify LDL^T factorization correctness via:
//! - **Reconstruction error**: `||P^T A P - L D L^T||_F / ||A||_F`
//! - **Backward error**: `||Ax - b|| / (||A||_F ||x|| + ||b||)`
//! - **Inertia comparison**: field-wise equality of eigenvalue sign counts
//! - **Permutation validation**: checks that a permutation is a valid bijection

use faer::linalg::matmul::matmul;
use faer::sparse::SparseColMat;
use faer::{Accum, Col, Mat, Par};

use crate::error::SparseError;
use crate::io::reference::{DBlock, Inertia, ReferenceFactorization};

/// Validate that `perm` is a valid permutation of `0..n`.
///
/// Returns `Ok(())` if `perm` has exactly `n` elements and is a bijection on `0..n`.
///
/// # Errors
///
/// - `perm.len() != n`
/// - Any element is out of bounds (`>= n`)
/// - Any element appears more than once
pub fn validate_permutation(perm: &[usize], n: usize) -> Result<(), SparseError> {
    if perm.len() != n {
        return Err(SparseError::InvalidInput {
            reason: format!(
                "permutation length ({}) does not match expected size ({})",
                perm.len(),
                n
            ),
        });
    }
    let mut seen = vec![false; n];
    for (i, &p) in perm.iter().enumerate() {
        if p >= n {
            return Err(SparseError::InvalidInput {
                reason: format!("permutation[{}] = {} is out of bounds (n = {})", i, p, n),
            });
        }
        if seen[p] {
            return Err(SparseError::InvalidInput {
                reason: format!("permutation has duplicate index {} at position {}", p, i),
            });
        }
        seen[p] = true;
    }
    Ok(())
}

/// Compute dense matrix-vector product `A * x`.
///
/// Converts `A` to dense and uses BLAS-3 matmul. Intended for validation
/// of small-to-medium matrices in tests, not for production solves.
///
/// # Future optimization
///
/// For large matrices, use `faer::sparse::linalg::matmul::sparse_dense_matmul`
/// to avoid the O(n^2) dense conversion.
pub fn dense_matvec(a: &SparseColMat<usize, f64>, x: &Col<f64>) -> Col<f64> {
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

/// Compute the relative reconstruction error of a factorization.
///
/// Returns `||P^T A P - L D L^T||_F / ||A||_F` where L, D, P come from
/// the reference factorization. The Frobenius norm is computed via faer's
/// `norm_l2()` method on `Mat`, which returns `||M||_F` (not the spectral 2-norm).
///
/// Returns 0.0 if `||A||_F` is zero.
pub fn reconstruction_error(
    a: &SparseColMat<usize, f64>,
    reference: &ReferenceFactorization,
) -> f64 {
    let n = a.nrows();
    let a_dense = a.to_dense();
    // faer's Mat::norm_l2() computes the Frobenius norm: sqrt(sum of squared entries)
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
/// Returns 0.0 if the denominator is zero (which only occurs when both `A` and `b`
/// are zero, since the denominator is a sum of non-negative IEEE 754 values).
///
/// # Future optimization
///
/// Currently converts `A` to dense for the matrix-vector product. For large
/// matrices, use `faer::sparse::linalg::matmul::sparse_dense_matmul` to avoid
/// the O(n^2) dense conversion.
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

    // The denominator is a sum of non-negative IEEE 754 values, so it is
    // exactly 0.0 only when all terms are exactly 0.0 (i.e., A=0 and b=0).
    let denom = a_norm * x_norm + b_norm;
    if denom == 0.0 {
        return 0.0;
    }

    r_norm / denom
}

/// Compute backward error using sparse matrix-vector multiply.
///
/// Returns `||Ax - b|| / (||A||_F * ||x|| + ||b||)`.
/// Uses direct sparse iteration to avoid O(n^2) dense conversion.
///
/// **Matrix storage**: Expects a full symmetric CSC matrix (both triangles
/// stored), which is what our `.mtx` reader and `SparseColMat::try_new_from_triplets`
/// produce. Each off-diagonal entry (i,j) appears in both column i and column j.
///
/// The Frobenius norm and matrix-vector product are computed directly from
/// the stored entries. This is O(nnz) work and O(n) extra memory.
pub fn sparse_backward_error(a: &SparseColMat<usize, f64>, x: &Col<f64>, b: &Col<f64>) -> f64 {
    let n = a.nrows();

    // Compute A*x directly from the full symmetric CSC.
    // Since both triangles are stored, a regular matvec is correct.
    let symbolic = a.symbolic();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let values = a.val();

    let mut ax = vec![0.0f64; n];
    for j in 0..n {
        for idx in col_ptrs[j]..col_ptrs[j + 1] {
            let i = row_indices[idx];
            let v = values[idx];
            ax[i] += v * x[j];
        }
    }

    // Compute residual r = Ax - b
    let mut r_norm_sq = 0.0;
    for k in 0..n {
        let r = ax[k] - b[k];
        r_norm_sq += r * r;
    }
    let r_norm = r_norm_sq.sqrt();
    let x_norm = x.norm_l2();
    let b_norm = b.norm_l2();

    // Compute ||A||_F from stored values directly.
    // Full symmetric storage: each off-diagonal entry appears twice,
    // so sum(v^2) = sum(diag^2) + 2*sum(off_diag^2) = ||A||_F^2.
    let mut a_norm_sq = 0.0;
    for j in 0..n {
        for &v in &values[col_ptrs[j]..col_ptrs[j + 1]] {
            a_norm_sq += v * v;
        }
    }
    let a_norm = a_norm_sq.sqrt();

    let denom = a_norm * x_norm + b_norm;
    if denom == 0.0 {
        return 0.0;
    }

    r_norm / denom
}

/// Check whether two inertia values are identical and dimensionally consistent.
///
/// Returns `true` only if all three counts match AND both inertias have the
/// same total dimension (positive + negative + zero).
pub fn check_inertia(computed: &Inertia, expected: &Inertia) -> bool {
    computed.dimension() == expected.dimension()
        && computed.positive == expected.positive
        && computed.negative == expected.negative
        && computed.zero == expected.zero
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::registry;

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

    // T024: validate_permutation tests

    #[test]
    fn validate_permutation_valid() {
        assert!(validate_permutation(&[2, 0, 1], 3).is_ok());
        assert!(validate_permutation(&[0], 1).is_ok());
        assert!(validate_permutation(&[], 0).is_ok());
    }

    #[test]
    fn validate_permutation_wrong_length() {
        let result = validate_permutation(&[0, 1], 3);
        assert!(result.is_err());
    }

    #[test]
    fn validate_permutation_out_of_bounds() {
        let result = validate_permutation(&[0, 5, 2], 3);
        assert!(result.is_err());
    }

    #[test]
    fn validate_permutation_duplicate() {
        let result = validate_permutation(&[0, 1, 1], 3);
        assert!(result.is_err());
    }

    // T025: dense_matvec test

    #[test]
    fn dense_matvec_matches_manual() {
        let test = registry::load_test_matrix("arrow-5-pd")
            .expect("registry error")
            .expect("matrix should exist");
        let n = test.matrix.nrows();
        let x = Col::<f64>::from_fn(n, |i| (i + 1) as f64);
        let b = dense_matvec(&test.matrix, &x);
        assert_eq!(b.nrows(), n);

        // Verify backward error is tiny
        let err = backward_error(&test.matrix, &x, &b);
        assert!(err < 1e-14);
    }
}
