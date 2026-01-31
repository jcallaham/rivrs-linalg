//! Triangular discrete-time Sylvester equation solver.
//!
//! Solves `X + sgn * A*X*B = scale*C` where A and B are quasi-triangular
//! (in real Schur form), and C is overwritten with the solution X.
//!
//! The discrete-time triangular Sylvester equation has a different
//! back-substitution structure than the continuous case because X appears
//! both standalone and in a product A*X*B.
//!
//! # Algorithm
//!
//! For `X + A*X*B = C`:
//! Process columns of B from left to right, rows of A from bottom to top.
//! At each (k, l) block:
//!   X(k,l) + A(k,k)*X(k,l)*B(l,l) = C(k,l) - sum_updates
//! which rearranges to:
//!   (I + A(k,k) ⊗ B(l,l)^T) vec(X(k,l)) = vec(C(k,l) - updates)
//!
//! # References
//!
//! - Golub, Nash & Van Loan (1979), "A Hessenberg-Schur method for AX + XB = C",
//!   IEEE Trans. Auto. Contr. AC-24:909-913
//! - Sima (1996), "Algorithms for Linear-Quadratic Optimization"

use faer::Accum;
use faer::prelude::*;

use super::triangular::{detect_blocks, matrix_max_abs};

/// Solves the triangular discrete-time Sylvester equation in-place.
///
/// Solves `X + sgn * A*X*B = scale*C` where A (m x m) and B (n x n) are
/// upper quasi-triangular (real Schur form). C (m x n) is overwritten
/// with the solution X. Returns `(scale, near_singular)`.
///
/// # References
///
/// - Sima (1996), "Algorithms for Linear-Quadratic Optimization"
pub fn solve_triangular_sylvester_discrete(
    a: MatRef<'_, f64>,
    b: MatRef<'_, f64>,
    c: MatMut<'_, f64>,
    sgn: f64,
) -> (f64, bool) {
    let m = a.nrows();
    let n = b.nrows();
    let mut c = c;
    let mut scale = 1.0f64;
    let mut near_singular = false;

    if m == 0 || n == 0 {
        return (scale, near_singular);
    }

    let eps = f64::EPSILON;
    let smlnum = f64::MIN_POSITIVE * (m * n) as f64 / eps;
    let a_max = matrix_max_abs(a);
    let b_max = matrix_max_abs(b);
    let smin = smlnum.max(eps * a_max).max(eps * b_max);

    let a_blocks = detect_blocks(a, m);
    let b_blocks = detect_blocks(b, n);

    // Solve X + sgn * A*X*B = scale*C
    // Process columns left to right, rows bottom to top
    for &(l1, l2) in &b_blocks {
        for k_idx in (0..a_blocks.len()).rev() {
            let (k1, k2) = a_blocks[k_idx];

            // Update C(k,l) by subtracting contributions from already-solved blocks
            //
            // For discrete case: X + sgn*A*X*B = C
            // => X(k,l) + sgn*A(k,k)*X(k,l)*B(l,l) = C(k,l)
            //    - sgn * sum_{i>k} A(k,i)*X(i,l)*B(l,l)     (A blocks below K in same B column)
            //    - sgn * sum_{j<l} A(k,k)*X(k,j)*B(j,l)     (B blocks left of L in same A row)
            //    - sgn * sum_{i>k,j<l} A(k,i)*X(i,j)*B(j,l) (cross terms)

            // Contribution from A blocks below K (for current B column block)
            if k2 + 1 < m {
                let a_sub = a.submatrix(k1, k2 + 1, k2 - k1 + 1, m - k2 - 1);
                let x_rows = m - k2 - 1;
                let x_cols = l2 - l1 + 1;

                // Compute A(k,i>k) * X(i>k, l) * B(l,l)
                let mut ax_tmp = Mat::zeros(k2 - k1 + 1, x_cols);
                {
                    let mut x_tmp = Mat::zeros(x_rows, x_cols);
                    for jj in 0..x_cols {
                        for ii in 0..x_rows {
                            x_tmp[(ii, jj)] = c[(k2 + 1 + ii, l1 + jj)];
                        }
                    }
                    // ax_tmp = A(k,i>k) * X(i>k, l)
                    faer::linalg::matmul::matmul(
                        ax_tmp.as_mut(),
                        Accum::Replace,
                        a_sub,
                        x_tmp.as_ref(),
                        1.0f64,
                        Par::Seq,
                    );
                }
                // Now multiply by B(l,l): ax_tmp * B(l,l)
                let b_ll = b.submatrix(l1, l1, l2 - l1 + 1, l2 - l1 + 1);
                let mut axb_tmp = Mat::zeros(k2 - k1 + 1, x_cols);
                faer::linalg::matmul::matmul(
                    axb_tmp.as_mut(),
                    Accum::Replace,
                    ax_tmp.as_ref(),
                    b_ll,
                    1.0f64,
                    Par::Seq,
                );

                // C(k,l) -= sgn * A(k,i>k) * X(i>k,l) * B(l,l)
                for jj in 0..x_cols {
                    for ii in 0..(k2 - k1 + 1) {
                        c[(k1 + ii, l1 + jj)] -= sgn * axb_tmp[(ii, jj)];
                    }
                }
            }

            // Contribution from B blocks left of L (for current A row block)
            if l1 > 0 {
                let a_kk = a.submatrix(k1, k1, k2 - k1 + 1, k2 - k1 + 1);
                let b_sub = b.submatrix(0, l1, l1, l2 - l1 + 1);

                // Compute A(k,k) * X(k, j<l) * B(j<l, l)
                let mut x_tmp = Mat::zeros(k2 - k1 + 1, l1);
                for jj in 0..l1 {
                    for ii in 0..(k2 - k1 + 1) {
                        x_tmp[(ii, jj)] = c[(k1 + ii, jj)];
                    }
                }
                // ax = A(k,k) * X(k, j<l)
                let mut ax = Mat::zeros(k2 - k1 + 1, l1);
                faer::linalg::matmul::matmul(
                    ax.as_mut(),
                    Accum::Replace,
                    a_kk,
                    x_tmp.as_ref(),
                    1.0f64,
                    Par::Seq,
                );
                // axb = ax * B(j<l, l)
                let mut axb = Mat::zeros(k2 - k1 + 1, l2 - l1 + 1);
                faer::linalg::matmul::matmul(
                    axb.as_mut(),
                    Accum::Replace,
                    ax.as_ref(),
                    b_sub,
                    1.0f64,
                    Par::Seq,
                );

                for jj in 0..(l2 - l1 + 1) {
                    for ii in 0..(k2 - k1 + 1) {
                        c[(k1 + ii, l1 + jj)] -= sgn * axb[(ii, jj)];
                    }
                }
            }

            // Cross terms: A blocks below K AND B blocks left of L
            if k2 + 1 < m && l1 > 0 {
                let a_sub = a.submatrix(k1, k2 + 1, k2 - k1 + 1, m - k2 - 1);
                let b_sub = b.submatrix(0, l1, l1, l2 - l1 + 1);

                // X(i>k, j<l) - already fully solved
                let x_rows = m - k2 - 1;
                let mut x_tmp = Mat::zeros(x_rows, l1);
                for jj in 0..l1 {
                    for ii in 0..x_rows {
                        x_tmp[(ii, jj)] = c[(k2 + 1 + ii, jj)];
                    }
                }
                // ax = A(k,i>k) * X(i>k, j<l)
                let mut ax = Mat::zeros(k2 - k1 + 1, l1);
                faer::linalg::matmul::matmul(
                    ax.as_mut(),
                    Accum::Replace,
                    a_sub,
                    x_tmp.as_ref(),
                    1.0f64,
                    Par::Seq,
                );
                // axb = ax * B(j<l, l)
                let mut axb = Mat::zeros(k2 - k1 + 1, l2 - l1 + 1);
                faer::linalg::matmul::matmul(
                    axb.as_mut(),
                    Accum::Replace,
                    ax.as_ref(),
                    b_sub,
                    1.0f64,
                    Par::Seq,
                );

                for jj in 0..(l2 - l1 + 1) {
                    for ii in 0..(k2 - k1 + 1) {
                        c[(k1 + ii, l1 + jj)] -= sgn * axb[(ii, jj)];
                    }
                }
            }

            // Now solve the small block system:
            // X(k,l) + sgn * A(k,k) * X(k,l) * B(l,l) = C(k,l)
            let k_size = k2 - k1 + 1;
            let l_size = l2 - l1 + 1;

            let (scaloc, ns) = solve_discrete_small_block(
                a.submatrix(k1, k1, k_size, k_size),
                b.submatrix(l1, l1, l_size, l_size),
                c.rb_mut().submatrix_mut(k1, l1, k_size, l_size),
                sgn,
                smin,
            );

            if ns {
                near_singular = true;
            }

            if scaloc != 1.0 {
                for j in 0..n {
                    for i in 0..m {
                        c[(i, j)] *= scaloc;
                    }
                }
                scale *= scaloc;
            }
        }
    }

    (scale, near_singular)
}

/// Solves the small discrete block Sylvester equation:
///   X + sgn * A_kk * X * B_ll = C_kl
///
/// This is equivalent to: (I + sgn * B_ll^T ⊗ A_kk) vec(X) = vec(C_kl)
#[allow(clippy::needless_range_loop)]
fn solve_discrete_small_block(
    a_kk: MatRef<'_, f64>,
    b_ll: MatRef<'_, f64>,
    mut c_kl: MatMut<'_, f64>,
    sgn: f64,
    smin: f64,
) -> (f64, bool) {
    let k_size = a_kk.nrows();
    let l_size = b_ll.nrows();

    match (k_size, l_size) {
        (1, 1) => {
            // x + sgn * a11 * x * b11 = c  =>  x * (1 + sgn*a11*b11) = c
            let denom = 1.0 + sgn * a_kk[(0, 0)] * b_ll[(0, 0)];
            let mut near_singular = false;
            let eff_denom = if denom.abs() <= smin {
                near_singular = true;
                smin
            } else {
                denom
            };
            c_kl[(0, 0)] /= eff_denom;
            (1.0, near_singular)
        }
        _ => {
            // General case: vectorize the equation
            // (I + sgn * B^T ⊗ A) vec(X) = vec(C)
            let total_size = k_size * l_size;
            assert!(total_size <= 4, "Block size too large");

            // Build coefficient matrix: I + sgn * (B^T ⊗ A)
            let mut mat = [[0.0f64; 4]; 4];
            for i in 0..total_size {
                mat[i][i] = 1.0; // Identity part
            }
            // Add sgn * (B^T ⊗ A)
            // (B^T ⊗ A)_{(j1*k_size+i1), (j2*k_size+i2)} = B[j2][j1] * A[i1][i2]
            for j1 in 0..l_size {
                for i1 in 0..k_size {
                    for j2 in 0..l_size {
                        for i2 in 0..k_size {
                            let row = j1 * k_size + i1;
                            let col = j2 * k_size + i2;
                            mat[row][col] += sgn * b_ll[(j2, j1)] * a_kk[(i1, i2)];
                        }
                    }
                }
            }

            // RHS: vec(C) in column-major order
            let mut rhs = [0.0f64; 4];
            for j in 0..l_size {
                for i in 0..k_size {
                    rhs[j * k_size + i] = c_kl[(i, j)];
                }
            }

            // Solve using Gaussian elimination with partial pivoting
            let mut near_singular = false;
            let mut ipiv = [0usize, 1, 2, 3];

            for col in 0..total_size {
                // Find pivot
                let mut max_val = 0.0f64;
                let mut max_row = col;
                for row in col..total_size {
                    let val = mat[ipiv[row]][col].abs();
                    if val > max_val {
                        max_val = val;
                        max_row = row;
                    }
                }
                ipiv.swap(col, max_row);

                let pivot = mat[ipiv[col]][col];
                if pivot.abs() <= smin {
                    near_singular = true;
                    mat[ipiv[col]][col] = smin;
                    let pivot = smin;
                    for row in (col + 1)..total_size {
                        let factor = mat[ipiv[row]][col] / pivot;
                        for k in (col + 1)..total_size {
                            mat[ipiv[row]][k] -= factor * mat[ipiv[col]][k];
                        }
                        rhs[ipiv[row]] -= factor * rhs[ipiv[col]];
                    }
                } else {
                    for row in (col + 1)..total_size {
                        let factor = mat[ipiv[row]][col] / pivot;
                        for k in (col + 1)..total_size {
                            mat[ipiv[row]][k] -= factor * mat[ipiv[col]][k];
                        }
                        rhs[ipiv[row]] -= factor * rhs[ipiv[col]];
                    }
                }
            }

            // Back-substitution
            let mut x = [0.0f64; 4];
            for i in (0..total_size).rev() {
                let mut sum = rhs[ipiv[i]];
                for j in (i + 1)..total_size {
                    sum -= mat[ipiv[i]][j] * x[j];
                }
                let diag = mat[ipiv[i]][i];
                if diag.abs() <= smin {
                    near_singular = true;
                    x[i] = sum / smin;
                } else {
                    x[i] = sum / diag;
                }
            }

            // Extract solution back to C (column-major)
            for j in 0..l_size {
                for i in 0..k_size {
                    c_kl[(i, j)] = x[j * k_size + i];
                }
            }

            (1.0, near_singular)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_discrete_1x1() {
        // X + A*X*B = C where A=[2], B=[3], C=[7]
        // x + 2*x*3 = 7 => x*(1+6) = 7 => x = 1
        let a = mat![[2.0f64]];
        let b = mat![[3.0f64]];
        let mut c = mat![[7.0f64]];
        let (scale, _ns) =
            solve_triangular_sylvester_discrete(a.as_ref(), b.as_ref(), c.as_mut(), 1.0);
        assert!((c[(0, 0)] / scale - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_discrete_2x2_diagonal() {
        // Diagonal A and B: each element independent
        // x_ij + a_ii * x_ij * b_jj = c_ij  =>  x_ij = c_ij / (1 + a_ii*b_jj)
        let a = mat![[1.0, 0.0], [0.0, 2.0f64]];
        let b = mat![[3.0, 0.0], [0.0, 4.0f64]];
        let c_orig = mat![[4.0, 5.0], [6.0, 9.0f64]];
        let mut c = c_orig.clone();
        let (scale, _ns) =
            solve_triangular_sylvester_discrete(a.as_ref(), b.as_ref(), c.as_mut(), 1.0);
        let s = scale;
        // x11 = 4/(1+1*3) = 1
        // x12 = 5/(1+1*4) = 1
        // x21 = 6/(1+2*3) = 6/7
        // x22 = 9/(1+2*4) = 1
        assert!((c[(0, 0)] / s - 1.0).abs() < 1e-12);
        assert!((c[(0, 1)] / s - 1.0).abs() < 1e-12);
        assert!((c[(1, 0)] / s - 6.0 / 7.0).abs() < 1e-12);
        assert!((c[(1, 1)] / s - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_discrete_upper_triangular() {
        let a = mat![[1.0, 0.5], [0.0, 2.0f64]];
        let b = mat![[3.0, 1.0], [0.0, 4.0f64]];
        let c_orig = mat![[10.0, 11.0], [12.0, 13.0f64]];
        let mut c = c_orig.clone();
        let (scale, _ns) =
            solve_triangular_sylvester_discrete(a.as_ref(), b.as_ref(), c.as_mut(), 1.0);

        // Verify: X + A*X*B should equal scale*C_orig
        let x = &c;
        let ax = &a * x;
        let axb = &ax * &b;
        let mut max_err = 0.0f64;
        for j in 0..2 {
            for i in 0..2 {
                let residual = (x[(i, j)] + axb[(i, j)] - scale * c_orig[(i, j)]).abs();
                max_err = max_err.max(residual);
            }
        }
        assert!(max_err < 1e-10, "Residual too large: {}", max_err);
    }

    #[test]
    fn test_discrete_3x3() {
        let a = mat![[0.5, 0.1, 0.0], [0.0, 0.8, 0.2], [0.0, 0.0, 0.3f64]];
        let b = mat![[0.6, 0.1, 0.0], [0.0, 0.9, 0.05], [0.0, 0.0, 0.4f64]];
        let c_orig = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0f64]];
        let mut c = c_orig.clone();
        let (scale, _ns) =
            solve_triangular_sylvester_discrete(a.as_ref(), b.as_ref(), c.as_mut(), 1.0);

        // Verify residual
        let x = &c;
        let axb = &(&a * x) * &b;
        let mut max_err = 0.0f64;
        for j in 0..3 {
            for i in 0..3 {
                let residual = (x[(i, j)] + axb[(i, j)] - scale * c_orig[(i, j)]).abs();
                max_err = max_err.max(residual);
            }
        }
        assert!(max_err < 1e-10, "Residual too large: {}", max_err);
    }

    #[test]
    fn test_discrete_negative_sign() {
        // Solve X - A*X*B = C (sgn = -1)
        let a = mat![[0.5, 0.0], [0.0, 0.8f64]];
        let b = mat![[0.6, 0.0], [0.0, 0.9f64]];
        let c_orig = mat![[1.0, 2.0], [3.0, 4.0f64]];
        let mut c = c_orig.clone();
        let (scale, _ns) =
            solve_triangular_sylvester_discrete(a.as_ref(), b.as_ref(), c.as_mut(), -1.0);

        // Verify: X - A*X*B = scale*C_orig
        let x = &c;
        let axb = &(&a * x) * &b;
        let mut max_err = 0.0f64;
        for j in 0..2 {
            for i in 0..2 {
                let residual = (x[(i, j)] - axb[(i, j)] - scale * c_orig[(i, j)]).abs();
                max_err = max_err.max(residual);
            }
        }
        assert!(max_err < 1e-10, "Residual too large: {}", max_err);
    }
}
