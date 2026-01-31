//! Shared test utilities for integration tests.

#![allow(dead_code)]

use faer::prelude::*;

/// Build a matrix with known eigenvalue structure that guarantees
/// 2x2 quasi-triangular blocks in its real Schur form.
/// Uses 2x2 rotation blocks on the diagonal to produce complex eigenvalue pairs.
pub fn build_quasi_triangular_test_matrix(n: usize) -> Mat<f64> {
    let mut m = Mat::zeros(n, n);
    let mut i = 0;
    let mut block_idx = 0;
    while i < n {
        if i + 1 < n && block_idx % 3 == 1 {
            // 2x2 rotation block: eigenvalues = center ± freq*i
            let center = (block_idx + 2) as f64;
            let freq = 0.5;
            m[(i, i)] = center;
            m[(i, i + 1)] = freq;
            m[(i + 1, i)] = -freq;
            m[(i + 1, i + 1)] = center;
            i += 2;
        } else {
            // 1x1 block with real eigenvalue
            m[(i, i)] = (block_idx + 1) as f64;
            i += 1;
        }
        block_idx += 1;
    }
    // Add some upper triangular fill to make it non-trivial
    for j in 0..n {
        for ii in 0..j.min(3) {
            m[(ii, j)] += 0.1 * ((ii + j) as f64);
        }
    }
    m
}

/// Compute Schur decomposition using nalgebra (same method as the solver internally).
/// Returns (schur_form T, orthogonal matrix U) where A = U*T*U^T.
pub fn compute_schur(a: MatRef<'_, f64>) -> (Mat<f64>, Mat<f64>) {
    let n = a.nrows();
    let na_matrix = nalgebra::DMatrix::from_fn(n, n, |i, j| a[(i, j)]);
    let schur = nalgebra::Schur::new(na_matrix);
    let (na_u, na_t) = schur.unpack();
    let mut t = Mat::zeros(n, n);
    let mut u = Mat::zeros(n, n);
    for j in 0..n {
        for i in 0..n {
            t[(i, j)] = na_t[(i, j)];
            u[(i, j)] = na_u[(i, j)];
        }
    }
    (t, u)
}

/// Check if a matrix's Schur form contains 2x2 blocks (complex eigenvalues).
pub fn schur_has_2x2_blocks(a: MatRef<'_, f64>) -> bool {
    let (schur, _) = compute_schur(a);
    (0..schur.nrows().saturating_sub(1)).any(|i| schur[(i + 1, i)].abs() > 1e-10)
}
