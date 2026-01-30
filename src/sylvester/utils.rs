//! Utility functions for Sylvester equation solvers.

use faer::prelude::*;
use faer::Accum;

use super::types::EquationType;

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
            faer::linalg::matmul::matmul(
                r.as_mut(),
                Accum::Replace,
                a,
                x,
                1.0f64,
                Par::Seq,
            );
            // R += X*B
            faer::linalg::matmul::matmul(
                r.as_mut(),
                Accum::Add,
                x,
                b,
                1.0f64,
                Par::Seq,
            );
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
            faer::linalg::matmul::matmul(
                tmp.as_mut(),
                Accum::Replace,
                a,
                x,
                1.0f64,
                Par::Seq,
            );
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
