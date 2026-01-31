//! Utility functions for Sylvester equation solvers.

use faer::prelude::*;
use faer::Accum;

use super::types::EquationType;
use crate::error::SylvesterError;

/// Computes the residual norm for a Sylvester equation solution.
///
/// For continuous-time: `||AX + XB - C||_F`
/// For discrete-time: `||AXB + X - C||_F`
///
/// The Frobenius norm is used as it is both easy to compute and
/// provides a good measure of backward error.
pub fn compute_residual(
    a: MatRef<'_, f64>,
    b: MatRef<'_, f64>,
    c: MatRef<'_, f64>,
    x: MatRef<'_, f64>,
    equation_type: EquationType,
) -> f64 {
    let n = a.nrows();
    let m = b.nrows();

    match equation_type {
        EquationType::Continuous => {
            // R = AX + XB - C
            let mut r = Mat::zeros(n, m);
            // R = A*X
            faer::linalg::matmul::matmul(r.as_mut(), Accum::Replace, a, x, 1.0f64, Par::Seq);
            // R += X*B
            faer::linalg::matmul::matmul(r.as_mut(), Accum::Add, x, b, 1.0f64, Par::Seq);
            // R -= C
            for j in 0..m {
                for i in 0..n {
                    r[(i, j)] -= c[(i, j)];
                }
            }
            frobenius_norm(r.as_ref())
        }
        EquationType::Discrete => {
            // R = AXB + X - C
            let mut axb = Mat::zeros(n, m);
            // tmp = A*X
            let mut tmp = Mat::zeros(n, m);
            faer::linalg::matmul::matmul(tmp.as_mut(), Accum::Replace, a, x, 1.0f64, Par::Seq);
            // axb = tmp*B = A*X*B
            faer::linalg::matmul::matmul(
                axb.as_mut(),
                Accum::Replace,
                tmp.as_ref(),
                b,
                1.0f64,
                Par::Seq,
            );
            // R = AXB + X - C
            let mut r = Mat::zeros(n, m);
            for j in 0..m {
                for i in 0..n {
                    r[(i, j)] = axb[(i, j)] + x[(i, j)] - c[(i, j)];
                }
            }
            frobenius_norm(r.as_ref())
        }
    }
}

/// Computes the Frobenius norm of a matrix: sqrt(sum of squared elements).
pub fn frobenius_norm(m: MatRef<'_, f64>) -> f64 {
    let mut sum = 0.0f64;
    for j in 0..m.ncols() {
        for i in 0..m.nrows() {
            sum += m[(i, j)] * m[(i, j)];
        }
    }
    sum.sqrt()
}

/// Computes the real Schur decomposition of a matrix.
///
/// Returns `(schur_form, schur_vectors)` where:
/// - `schur_form` is the quasi-triangular Schur form T
/// - `schur_vectors` is the orthogonal matrix U such that A = U * T * U^T
///
/// Uses nalgebra's Schur decomposition internally, converting between
/// faer and nalgebra matrix types.
pub(crate) fn compute_real_schur(
    a: MatRef<'_, f64>,
) -> Result<(Mat<f64>, Mat<f64>), SylvesterError> {
    let n = a.nrows();

    if n == 0 {
        return Ok((Mat::zeros(0, 0), Mat::zeros(0, 0)));
    }

    if n == 1 {
        let mut t = Mat::zeros(1, 1);
        t[(0, 0)] = a[(0, 0)];
        let mut u = Mat::zeros(1, 1);
        u[(0, 0)] = 1.0;
        return Ok((t, u));
    }

    // Convert faer matrix to nalgebra DMatrix
    let na_matrix = nalgebra::DMatrix::from_fn(n, n, |i, j| a[(i, j)]);

    // Compute real Schur decomposition using nalgebra
    let schur = nalgebra::Schur::new(na_matrix);

    // Extract results: unpack() returns (Q, T) where A = Q * T * Q^T
    let (na_u, na_t) = schur.unpack();

    // Convert nalgebra matrices back to faer
    let mut t = Mat::zeros(n, n);
    let mut u = Mat::zeros(n, n);

    for j in 0..n {
        for i in 0..n {
            t[(i, j)] = na_t[(i, j)];
            u[(i, j)] = na_u[(i, j)];
        }
    }

    // Verify convergence by checking that U is approximately orthogonal
    // U^T * U should be close to identity
    let mut utu = Mat::zeros(n, n);
    faer::linalg::matmul::matmul(
        utu.as_mut(),
        Accum::Replace,
        u.as_ref().transpose(),
        u.as_ref(),
        1.0f64,
        Par::Seq,
    );

    let mut max_err = 0.0f64;
    for j in 0..n {
        for i in 0..n {
            let expected = if i == j { 1.0 } else { 0.0 };
            max_err = max_err.max((utu[(i, j)] - expected).abs());
        }
    }

    if max_err > 1e-10 {
        return Err(SylvesterError::ConvergenceFailure {
            algorithm: format!(
                "Schur decomposition: orthogonality check failed (error: {:.2e})",
                max_err
            ),
        });
    }

    Ok((t, u))
}

/// Transforms the right-hand side for the Bartels-Stewart algorithm.
///
/// Computes `F = U^T * C * V` where U and V are orthogonal Schur vectors.
pub(crate) fn transform_rhs(
    u: MatRef<'_, f64>,
    c: MatRef<'_, f64>,
    v: MatRef<'_, f64>,
) -> Mat<f64> {
    let n = u.nrows();
    let m = v.nrows();

    // tmp = U^T * C
    let mut tmp = Mat::zeros(n, m);
    faer::linalg::matmul::matmul(
        tmp.as_mut(),
        Accum::Replace,
        u.transpose(),
        c,
        1.0f64,
        Par::Seq,
    );
    // F = tmp * V
    let mut f = Mat::zeros(n, m);
    faer::linalg::matmul::matmul(
        f.as_mut(),
        Accum::Replace,
        tmp.as_ref(),
        v,
        1.0f64,
        Par::Seq,
    );
    f
}

/// Back-transforms the solution from Schur coordinates.
///
/// Computes `X = U * Y * V^T` where U and V are orthogonal Schur vectors.
pub(crate) fn back_transform(
    u: MatRef<'_, f64>,
    y: MatRef<'_, f64>,
    v: MatRef<'_, f64>,
) -> Mat<f64> {
    let n = u.nrows();
    let m = v.nrows();

    // tmp = U * Y
    let mut tmp = Mat::zeros(n, m);
    faer::linalg::matmul::matmul(tmp.as_mut(), Accum::Replace, u, y, 1.0f64, Par::Seq);
    // X = tmp * V^T
    let mut x = Mat::zeros(n, m);
    faer::linalg::matmul::matmul(
        x.as_mut(),
        Accum::Replace,
        tmp.as_ref(),
        v.transpose(),
        1.0f64,
        Par::Seq,
    );
    x
}

/// Reconstructs a matrix from its Schur decomposition.
///
/// Computes `A = U * T * U^T` where T is the Schur form and U is orthogonal.
pub(crate) fn reconstruct_from_schur(u: MatRef<'_, f64>, t: MatRef<'_, f64>) -> Mat<f64> {
    let n = u.nrows();

    // tmp = U * T
    let mut tmp = Mat::zeros(n, n);
    faer::linalg::matmul::matmul(tmp.as_mut(), Accum::Replace, u, t, 1.0f64, Par::Seq);
    // A = tmp * U^T
    let mut a = Mat::zeros(n, n);
    faer::linalg::matmul::matmul(
        a.as_mut(),
        Accum::Replace,
        tmp.as_ref(),
        u.transpose(),
        1.0f64,
        Par::Seq,
    );
    a
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_frobenius_norm_identity() {
        let m = Mat::<f64>::identity(3, 3);
        assert!((frobenius_norm(m.as_ref()) - 3.0f64.sqrt()).abs() < 1e-14);
    }

    #[test]
    fn test_frobenius_norm_single_element() {
        let m = mat![[5.0f64]];
        assert!((frobenius_norm(m.as_ref()) - 5.0).abs() < 1e-14);
    }

    #[test]
    fn test_frobenius_norm_zero() {
        let m = Mat::<f64>::zeros(3, 4);
        assert_eq!(frobenius_norm(m.as_ref()), 0.0);
    }

    #[test]
    fn test_frobenius_norm_rectangular() {
        // [[1, 2], [3, 4], [5, 6]] -> sqrt(1+4+9+16+25+36) = sqrt(91)
        let m = mat![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0f64]];
        assert!((frobenius_norm(m.as_ref()) - 91.0f64.sqrt()).abs() < 1e-14);
    }

    #[test]
    fn test_residual_continuous_exact_solution() {
        // A = [[1, 0], [0, 2]], B = [[3, 0], [0, 4]], X = [[1/4, 0], [0, 1/6]]
        // AX + XB = [[1/4, 0], [0, 2/6]] + [[3/4, 0], [0, 4/6]] = [[1, 0], [0, 1]] = C
        let a = mat![[1.0, 0.0], [0.0, 2.0f64]];
        let b = mat![[3.0, 0.0], [0.0, 4.0f64]];
        let c = mat![[1.0, 0.0], [0.0, 1.0f64]];
        let x = mat![[0.25, 0.0], [0.0, 1.0 / 6.0f64]];
        let r = compute_residual(
            a.as_ref(),
            b.as_ref(),
            c.as_ref(),
            x.as_ref(),
            EquationType::Continuous,
        );
        assert!(
            r < 1e-14,
            "Residual for exact solution should be ~0, got {:.2e}",
            r
        );
    }

    #[test]
    fn test_residual_continuous_wrong_solution() {
        // Using X = I (wrong) for AX + XB = C where A=I, B=I, C=[[1,0],[0,1]]
        // AX + XB = I + I = 2I, residual = ||2I - I|| = ||I|| = sqrt(2)
        let a = Mat::<f64>::identity(2, 2);
        let b = Mat::<f64>::identity(2, 2);
        let c = Mat::<f64>::identity(2, 2);
        let x = Mat::<f64>::identity(2, 2); // Wrong: correct is 0.5*I
        let r = compute_residual(
            a.as_ref(),
            b.as_ref(),
            c.as_ref(),
            x.as_ref(),
            EquationType::Continuous,
        );
        assert!(
            (r - 2.0f64.sqrt()).abs() < 1e-14,
            "Expected sqrt(2), got {}",
            r
        );
    }

    #[test]
    fn test_residual_discrete_exact_solution() {
        // A = [[2]], B = [[3]], X = [[1]], C = AXB + X = 2*1*3 + 1 = 7
        let a = mat![[2.0f64]];
        let b = mat![[3.0f64]];
        let c = mat![[7.0f64]];
        let x = mat![[1.0f64]];
        let r = compute_residual(
            a.as_ref(),
            b.as_ref(),
            c.as_ref(),
            x.as_ref(),
            EquationType::Discrete,
        );
        assert!(
            r < 1e-14,
            "Residual for exact solution should be ~0, got {:.2e}",
            r
        );
    }

    #[test]
    fn test_residual_discrete_wrong_solution() {
        // A = [[1]], B = [[1]], C = [[3]], X = [[1]] (wrong: correct is 3/2)
        // AXB + X = 1*1*1 + 1 = 2, residual = |2 - 3| = 1
        let a = mat![[1.0f64]];
        let b = mat![[1.0f64]];
        let c = mat![[3.0f64]];
        let x = mat![[1.0f64]];
        let r = compute_residual(
            a.as_ref(),
            b.as_ref(),
            c.as_ref(),
            x.as_ref(),
            EquationType::Discrete,
        );
        assert!((r - 1.0).abs() < 1e-14, "Expected 1.0, got {}", r);
    }

    #[test]
    fn test_residual_continuous_rectangular() {
        // A (3x3), B (2x2), X (3x2), C (3x2)
        let a = mat![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0f64]];
        let b = mat![[1.0, 0.0], [0.0, 1.0f64]];
        // AX + XB = IX + XI = 2X, so if C = 2X, residual = 0
        let x = mat![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0f64]];
        let c = &x * 2.0;
        let r = compute_residual(
            a.as_ref(),
            b.as_ref(),
            c.as_ref(),
            x.as_ref(),
            EquationType::Continuous,
        );
        assert!(r < 1e-14, "Residual should be ~0, got {:.2e}", r);
    }
}
