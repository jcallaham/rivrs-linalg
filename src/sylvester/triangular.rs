//! Triangular Sylvester equation solver.
//!
//! Solves `op(A)*X + ISGN*X*op(B) = scale*C` where A and B are quasi-triangular
//! (in real Schur form), and C is overwritten with the solution X.
//!
//! This is the core back-substitution algorithm used by the Bartels-Stewart method.
//!
//! # Algorithm
//!
//! Processes the equation block-by-block, handling 1x1 and 2x2 diagonal blocks
//! that arise from the real Schur decomposition. The ordering of processing
//! depends on the transpose flags.
//!
//! # References
//!
//! - Bartels & Stewart (1972), "Solution of the Matrix Equation AX + XB = C",
//!   CACM 15(9):820-826
//! - Golub & Van Loan (2013), "Matrix Computations" (4th Ed), Algorithm 7.6.2
//! - LAPACK dtrsyl (BSD-3-Clause) for parameter conventions and overflow prevention

use faer::prelude::*;
use faer::Accum;

/// Solves the triangular Sylvester equation in-place.
///
/// Solves `A*X + sgn*X*B = scale*C` where A (m x m) and B (n x n) are
/// upper quasi-triangular (real Schur form). C (m x n) is overwritten
/// with the solution X. Returns `(scale, near_singular)`.
///
/// The `sgn` parameter is +1.0 or -1.0.
///
/// # Algorithm
///
/// For the case A*X + X*B = C (sgn = +1, NoTrans/NoTrans):
/// Process columns of B from left to right, rows of A from bottom to top.
/// At each (k, l) block:
///   A(k,k)*X(k,l) + X(k,l)*B(l,l) = C(k,l) - sum_updates
///
/// # References
///
/// - Golub & Van Loan (2013), Algorithm 7.6.2
/// - LAPACK dtrsyl for overflow prevention via SCALE factor
pub fn solve_triangular_sylvester(
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

    // Overflow prevention constants (following LAPACK dtrsyl pattern)
    let eps = f64::EPSILON;
    let smlnum = f64::MIN_POSITIVE * (m * n) as f64 / eps;

    // Compute SMIN = max(smlnum, eps * max_element(A), eps * max_element(B))
    let a_max = matrix_max_abs(a);
    let b_max = matrix_max_abs(b);
    let smin = smlnum.max(eps * a_max).max(eps * b_max);

    // Detect block structure of A (row blocks)
    let a_blocks = detect_blocks(a, m);
    // Detect block structure of B (column blocks)
    let b_blocks = detect_blocks(b, n);

    // Solve A*X + sgn*X*B = scale*C
    // Process columns left to right (L blocks), rows bottom to top (K blocks)
    for l_idx in 0..b_blocks.len() {
        let (l1, l2) = b_blocks[l_idx];

        for k_idx in (0..a_blocks.len()).rev() {
            let (k1, k2) = a_blocks[k_idx];

            // Compute the RHS update: subtract contributions from already-solved blocks
            // R(K,L) = sum_{i>k} A(K,I)*X(I,L) + sgn * sum_{j<l} X(K,J)*B(J,L)

            // Contribution from A blocks below K (already solved, since we go bottom-up)
            if k2 + 1 < m {
                // C(k1..=k2, l1..=l2) -= A(k1..=k2, k2+1..m) * C(k2+1..m, l1..=l2)
                let a_sub = a.submatrix(k1, k2 + 1, k2 - k1 + 1, m - k2 - 1);
                // Copy x_sub to avoid borrow conflict
                let x_rows = m - k2 - 1;
                let x_cols = l2 - l1 + 1;
                let mut x_tmp = Mat::zeros(x_rows, x_cols);
                for jj in 0..x_cols {
                    for ii in 0..x_rows {
                        x_tmp[(ii, jj)] = c[(k2 + 1 + ii, l1 + jj)];
                    }
                }
                let mut c_sub = c.rb_mut().submatrix_mut(k1, l1, k2 - k1 + 1, l2 - l1 + 1);
                faer::linalg::matmul::matmul(
                    c_sub.rb_mut(),
                    Accum::Add,
                    a_sub,
                    x_tmp.as_ref(),
                    -1.0f64,
                    Par::Seq,
                );
            }

            // Contribution from B blocks to the left of L (already solved, since we go left-to-right)
            if l1 > 0 {
                // C(k1..=k2, l1..=l2) -= sgn * C(k1..=k2, 0..l1) * B(0..l1, l1..=l2)
                let b_sub = b.submatrix(0, l1, l1, l2 - l1 + 1);
                // Copy x_sub to avoid borrow conflict
                let x_rows = k2 - k1 + 1;
                let x_cols = l1;
                let mut x_tmp = Mat::zeros(x_rows, x_cols);
                for jj in 0..x_cols {
                    for ii in 0..x_rows {
                        x_tmp[(ii, jj)] = c[(k1 + ii, jj)];
                    }
                }
                let mut c_sub = c.rb_mut().submatrix_mut(k1, l1, k2 - k1 + 1, l2 - l1 + 1);
                faer::linalg::matmul::matmul(
                    c_sub.rb_mut(),
                    Accum::Add,
                    x_tmp.as_ref(),
                    b_sub,
                    -sgn,
                    Par::Seq,
                );
            }

            // Now solve the small block system:
            // A(k1..=k2, k1..=k2) * X + sgn * X * B(l1..=l2, l1..=l2) = C(k1..=k2, l1..=l2)
            let k_size = k2 - k1 + 1;
            let l_size = l2 - l1 + 1;

            let (scaloc, ns) = solve_small_block(
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
                // Rescale entire C by scaloc
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

/// Detects 1x1 and 2x2 diagonal blocks in a quasi-triangular matrix.
///
/// Returns a list of (start, end) index pairs (inclusive) for each block.
pub(crate) fn detect_blocks(mat: MatRef<'_, f64>, n: usize) -> Vec<(usize, usize)> {
    let mut blocks = Vec::new();
    let mut i = 0;
    while i < n {
        if i + 1 < n && mat[(i + 1, i)] != 0.0 {
            // 2x2 block
            blocks.push((i, i + 1));
            i += 2;
        } else {
            // 1x1 block
            blocks.push((i, i));
            i += 1;
        }
    }
    blocks
}

/// Solves the small (1x1 or 2x2) block Sylvester equation:
///   A_kk * X + sgn * X * B_ll = C_kl
///
/// Returns (scale_factor, near_singular).
fn solve_small_block(
    a_kk: MatRef<'_, f64>,
    b_ll: MatRef<'_, f64>,
    mut c_kl: MatMut<'_, f64>,
    sgn: f64,
    smin: f64,
) -> (f64, bool) {
    let k_size = a_kk.nrows();
    let l_size = b_ll.nrows();

    match (k_size, l_size) {
        (1, 1) => solve_1x1(a_kk[(0, 0)], b_ll[(0, 0)], &mut c_kl, sgn, smin),
        (1, 2) => solve_1x2(a_kk[(0, 0)], b_ll, &mut c_kl, sgn, smin),
        (2, 1) => solve_2x1(a_kk, b_ll[(0, 0)], &mut c_kl, sgn, smin),
        (2, 2) => solve_2x2(a_kk, b_ll, &mut c_kl, sgn, smin),
        _ => unreachable!("Block sizes must be 1 or 2"),
    }
}

/// Solve 1x1 block: a11 * x + sgn * x * b11 = c
///   => x = c / (a11 + sgn * b11)
fn solve_1x1(
    a11: f64,
    b11: f64,
    c: &mut MatMut<'_, f64>,
    sgn: f64,
    smin: f64,
) -> (f64, bool) {
    let mut near_singular = false;
    let mut suml = a11 + sgn * b11;

    if suml.abs() <= smin {
        suml = smin;
        near_singular = true;
    }

    let mut scaloc = 1.0;
    let rhs = c[(0, 0)];

    // Overflow check
    if smin != 0.0 && rhs.abs() > 1.0 && suml.abs() < 1.0 {
        if rhs.abs() * smin > suml.abs() {
            scaloc = 1.0 / rhs.abs();
        }
    }

    c[(0, 0)] = (rhs * scaloc) / suml;
    (scaloc, near_singular)
}

/// Solve 1x2 block: a11 * [x1, x2] + sgn * [x1, x2] * B = [c1, c2]
///
/// This becomes a 2x2 linear system when vectorized.
fn solve_1x2(
    a11: f64,
    b: MatRef<'_, f64>,
    c: &mut MatMut<'_, f64>,
    sgn: f64,
    smin: f64,
) -> (f64, bool) {
    // System: (a11 + sgn*b11)*x1 + sgn*b21*x2 = c1
    //         sgn*b12*x1 + (a11 + sgn*b22)*x2 = c2
    let d11 = a11 + sgn * b[(0, 0)];
    let d12 = sgn * b[(1, 0)];
    let d21 = sgn * b[(0, 1)];
    let d22 = a11 + sgn * b[(1, 1)];

    solve_2x2_linear_system(
        d11, d12, d21, d22, c, smin,
    )
}

/// Solve 2x1 block: A * [x1; x2] + sgn * [x1; x2] * b11 = [c1; c2]
///
/// This becomes a 2x2 linear system when vectorized.
fn solve_2x1(
    a: MatRef<'_, f64>,
    b11: f64,
    c: &mut MatMut<'_, f64>,
    sgn: f64,
    smin: f64,
) -> (f64, bool) {
    // System: (a11 + sgn*b11)*x1 + a12*x2 = c1
    //         a21*x1 + (a22 + sgn*b11)*x2 = c2
    let d11 = a[(0, 0)] + sgn * b11;
    let d12 = a[(0, 1)];
    let d21 = a[(1, 0)];
    let d22 = a[(1, 1)] + sgn * b11;

    solve_2x2_linear_system(
        d11, d12, d21, d22, c, smin,
    )
}

/// Solve 2x2 block Sylvester equation using vectorization (lasy2-style).
///
/// A * X + sgn * X * B = C where A is 2x2, B is 2x2, C is 2x2.
///
/// Vectorize: (I ⊗ A + sgn * B^T ⊗ I) * vec(X) = vec(C)
/// This gives a 4x4 linear system.
///
/// # References
///
/// - LAPACK DLASY2 (BSD-3-Clause) for the 2x2 Sylvester solver pattern
fn solve_2x2(
    a: MatRef<'_, f64>,
    b: MatRef<'_, f64>,
    c: &mut MatMut<'_, f64>,
    sgn: f64,
    smin: f64,
) -> (f64, bool) {
    let mut near_singular = false;

    // Build 4x4 system: (I ⊗ A + sgn * B^T ⊗ I) * vec(X) = vec(C)
    // vec(X) = [x11, x21, x12, x22]^T (column-major vectorization)
    let a11 = a[(0, 0)];
    let a12 = a[(0, 1)];
    let a21 = a[(1, 0)];
    let a22 = a[(1, 1)];
    let b11 = b[(0, 0)];
    let b12 = b[(0, 1)];
    let b21 = b[(1, 0)];
    let b22 = b[(1, 1)];

    // 4x4 coefficient matrix T:
    // T = I ⊗ A + sgn * B^T ⊗ I
    // = [[a11+sgn*b11, a12,          sgn*b21,     0          ],
    //    [a21,          a22+sgn*b11,  0,           sgn*b21    ],
    //    [sgn*b12,      0,           a11+sgn*b22, a12         ],
    //    [0,            sgn*b12,     a21,          a22+sgn*b22]]
    let mut t = [[0.0f64; 4]; 4];
    t[0][0] = a11 + sgn * b11;
    t[0][1] = a12;
    t[0][2] = sgn * b21;
    t[0][3] = 0.0;

    t[1][0] = a21;
    t[1][1] = a22 + sgn * b11;
    t[1][2] = 0.0;
    t[1][3] = sgn * b21;

    t[2][0] = sgn * b12;
    t[2][1] = 0.0;
    t[2][2] = a11 + sgn * b22;
    t[2][3] = a12;

    t[3][0] = 0.0;
    t[3][1] = sgn * b12;
    t[3][2] = a21;
    t[3][3] = a22 + sgn * b22;

    // RHS vector (column-major vectorization of C)
    let mut rhs = [c[(0, 0)], c[(1, 0)], c[(0, 1)], c[(1, 1)]];

    // Solve using Gaussian elimination with partial pivoting on the 4x4 system
    let mut ipiv = [0usize, 1, 2, 3];

    for col in 0..4 {
        // Find pivot
        let mut max_val = 0.0f64;
        let mut max_row = col;
        for row in col..4 {
            let val = t[ipiv[row]][col].abs();
            if val > max_val {
                max_val = val;
                max_row = row;
            }
        }

        ipiv.swap(col, max_row);

        let pivot = t[ipiv[col]][col];
        if pivot.abs() <= smin {
            // Near-singular: perturb
            near_singular = true;
            t[ipiv[col]][col] = smin;
            let pivot = smin;
            // Continue with perturbed pivot
            for row in (col + 1)..4 {
                let factor = t[ipiv[row]][col] / pivot;
                for k in (col + 1)..4 {
                    t[ipiv[row]][k] -= factor * t[ipiv[col]][k];
                }
                rhs[ipiv[row]] -= factor * rhs[ipiv[col]];
            }
        } else {
            for row in (col + 1)..4 {
                let factor = t[ipiv[row]][col] / pivot;
                for k in (col + 1)..4 {
                    t[ipiv[row]][k] -= factor * t[ipiv[col]][k];
                }
                rhs[ipiv[row]] -= factor * rhs[ipiv[col]];
            }
        }
    }

    // Back-substitution
    let mut x = [0.0f64; 4];
    for i in (0..4).rev() {
        let mut sum = rhs[ipiv[i]];
        for j in (i + 1)..4 {
            sum -= t[ipiv[i]][j] * x[j];
        }
        let diag = t[ipiv[i]][i];
        if diag.abs() <= smin {
            near_singular = true;
            x[i] = sum / smin;
        } else {
            x[i] = sum / diag;
        }
    }

    // Extract solution back to C (column-major)
    let scaloc = 1.0; // Simplified - full overflow prevention deferred to Phase 6
    c[(0, 0)] = x[0];
    c[(1, 0)] = x[1];
    c[(0, 1)] = x[2];
    c[(1, 1)] = x[3];

    (scaloc, near_singular)
}

/// Solve a 2x2 linear system with partial pivoting:
///   d11*x1 + d12*x2 = c1
///   d21*x1 + d22*x2 = c2
///
/// C is overwritten with the solution.
fn solve_2x2_linear_system(
    d11: f64,
    d12: f64,
    d21: f64,
    d22: f64,
    c: &mut MatMut<'_, f64>,
    smin: f64,
) -> (f64, bool) {
    let mut near_singular = false;

    // Extract RHS values based on whether c is 1×2 (from solve_1x2) or 2×1 (from solve_2x1)
    let (r1_val, r2_val) = if c.nrows() == 1 {
        (c[(0, 0)], c[(0, 1)])
    } else {
        (c[(0, 0)], c[(1, 0)])
    };

    let (a11, a12, a21, a22, r1, r2) = if d11.abs() >= d21.abs() {
        (d11, d12, d21, d22, r1_val, r2_val)
    } else {
        // Swap rows for partial pivoting
        (d21, d22, d11, d12, r2_val, r1_val)
    };

    let pivot = a11;
    let (x1, x2);

    if pivot.abs() <= smin {
        near_singular = true;
        // Perturb and solve anyway
        if a12.abs() > smin {
            x2 = r1 / a12;
            x1 = (r2 - a22 * x2) / if a21.abs() > smin { a21 } else { smin };
        } else {
            x1 = r1 / smin;
            x2 = (r2 - a21 * x1) / if a22.abs() > smin { a22 } else { smin };
        }
    } else {
        let factor = a21 / pivot;
        let new_a22 = a22 - factor * a12;
        let new_r2 = r2 - factor * r1;

        let eff_a22 = if new_a22.abs() <= smin {
            near_singular = true;
            smin
        } else {
            new_a22
        };

        x2 = new_r2 / eff_a22;
        x1 = (r1 - a12 * x2) / pivot;
    }

    // Write solution back. x1 and x2 correspond to the pivoted system.
    // If we swapped (d11.abs() < d21.abs()), x1/x2 correspond to swapped positions.
    let (sol1, sol2) = if d11.abs() >= d21.abs() {
        (x1, x2)
    } else {
        (x2, x1)
    };

    if c.nrows() == 1 {
        c[(0, 0)] = sol1;
        c[(0, 1)] = sol2;
    } else {
        c[(0, 0)] = sol1;
        c[(1, 0)] = sol2;
    }

    (1.0, near_singular)
}

/// Compute the maximum absolute value in a matrix.
pub(crate) fn matrix_max_abs(m: MatRef<'_, f64>) -> f64 {
    let mut max = 0.0f64;
    for j in 0..m.ncols() {
        for i in 0..m.nrows() {
            let v = m[(i, j)].abs();
            if v > max {
                max = v;
            }
        }
    }
    max
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_triangular_1x1() {
        // A = [3], B = [5], C = [16]
        // Solve: 3*x + x*5 = 16  =>  8x = 16  =>  x = 2
        let a = mat![[3.0f64]];
        let b = mat![[5.0f64]];
        let mut c = mat![[16.0f64]];
        let (scale, _ns) = solve_triangular_sylvester(
            a.as_ref(),
            b.as_ref(),
            c.as_mut(),
            1.0,
        );
        assert!((c[(0, 0)] / scale - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_triangular_2x2_diagonal() {
        // Diagonal A and B => each element of X is independent
        // A = diag(1, -1), B = diag(2, 3)
        // AX + XB = C means:
        //   (1+2)*x11 = c11,  (1+3)*x12 = c12
        //   (-1+2)*x21 = c21, (-1+3)*x22 = c22
        let a = mat![[1.0, 0.0], [0.0, -1.0f64]];
        let b = mat![[2.0, 0.0], [0.0, 3.0f64]];
        let mut c = mat![[3.0, 4.0], [1.0, 2.0f64]];
        let (scale, _ns) = solve_triangular_sylvester(
            a.as_ref(),
            b.as_ref(),
            c.as_mut(),
            1.0,
        );
        let s = scale;
        assert!((c[(0, 0)] / s - 1.0).abs() < 1e-12); // 3/3 = 1
        assert!((c[(0, 1)] / s - 1.0).abs() < 1e-12); // 4/4 = 1
        assert!((c[(1, 0)] / s - 1.0).abs() < 1e-12); // 1/1 = 1
        assert!((c[(1, 1)] / s - 1.0).abs() < 1e-12); // 2/2 = 1
    }

    #[test]
    fn test_triangular_upper_triangular() {
        // Upper triangular A and B
        // A = [[1, 2], [0, 3]], B = [[4, 1], [0, 5]]
        // Solve AX + XB = C where C = [[10, 11], [12, 13]]
        let a = mat![[1.0, 2.0], [0.0, 3.0f64]];
        let b = mat![[4.0, 1.0], [0.0, 5.0f64]];
        let c_orig = mat![[10.0, 11.0], [12.0, 13.0f64]];
        let mut c = c_orig.clone();
        let (scale, _ns) = solve_triangular_sylvester(
            a.as_ref(),
            b.as_ref(),
            c.as_mut(),
            1.0,
        );

        // Verify: A*X + X*B should equal scale*C_orig
        let x = &c;
        let ax = &a * x;
        let xb = x * &b;
        let residual = &ax + &xb - scale * &c_orig;

        let mut max_err = 0.0f64;
        for j in 0..2 {
            for i in 0..2 {
                max_err = max_err.max(residual[(i, j)].abs());
            }
        }
        assert!(max_err < 1e-12, "Residual too large: {}", max_err);
    }

    #[test]
    fn test_triangular_3x3() {
        // 3x3 upper triangular system
        let a = mat![
            [2.0, 1.0, 0.5],
            [0.0, 3.0, 0.3],
            [0.0, 0.0, 4.0f64]
        ];
        let b = mat![
            [1.0, 0.2, 0.1],
            [0.0, 2.0, 0.4],
            [0.0, 0.0, 3.0f64]
        ];
        let c_orig = mat![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0f64]
        ];
        let mut c = c_orig.clone();
        let (scale, _ns) = solve_triangular_sylvester(
            a.as_ref(),
            b.as_ref(),
            c.as_mut(),
            1.0,
        );

        // Verify residual
        let x = &c;
        let residual = &a * x + x * &b - scale * &c_orig;
        let mut max_err = 0.0f64;
        for j in 0..3 {
            for i in 0..3 {
                max_err = max_err.max(residual[(i, j)].abs());
            }
        }
        assert!(max_err < 1e-10, "Residual too large: {}", max_err);
    }

    #[test]
    fn test_triangular_with_2x2_block() {
        // Quasi-triangular with a 2x2 block in A
        // A = [[1, 2, 0], [-0.5, 1, 0], [0, 0, 3]]  (2x2 block at top-left)
        // B = [[4, 0], [0, 5]] (diagonal)
        let a = mat![
            [1.0, 2.0, 0.0],
            [-0.5, 1.0, 0.0],
            [0.0, 0.0, 3.0f64]
        ];
        let b = mat![[4.0, 0.0], [0.0, 5.0f64]];
        let c_orig = mat![
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0f64]
        ];
        let mut c = c_orig.clone();
        let (scale, _ns) = solve_triangular_sylvester(
            a.as_ref(),
            b.as_ref(),
            c.as_mut(),
            1.0,
        );

        // Verify residual
        let x = &c;
        let residual = &a * x + x * &b - scale * &c_orig;
        let mut max_err = 0.0f64;
        for j in 0..2 {
            for i in 0..3 {
                max_err = max_err.max(residual[(i, j)].abs());
            }
        }
        assert!(max_err < 1e-10, "Residual too large: {}", max_err);
    }

    #[test]
    fn test_triangular_negative_sign() {
        // Solve A*X - X*B = C (sgn = -1)
        let a = mat![[2.0, 0.0], [0.0, 4.0f64]];
        let b = mat![[1.0, 0.0], [0.0, 3.0f64]];
        let mut c = mat![[5.0, 6.0], [7.0, 8.0f64]];
        let c_orig = c.clone();
        let (scale, _ns) = solve_triangular_sylvester(
            a.as_ref(),
            b.as_ref(),
            c.as_mut(),
            -1.0,
        );

        // Verify: A*X - X*B = scale*C_orig
        let x = &c;
        let residual = &a * x - x * &b - scale * &c_orig;
        let mut max_err = 0.0f64;
        for j in 0..2 {
            for i in 0..2 {
                max_err = max_err.max(residual[(i, j)].abs());
            }
        }
        assert!(max_err < 1e-12, "Residual too large: {}", max_err);
    }

    #[test]
    fn test_triangular_rectangular() {
        // Non-square C: A is 3x3, B is 2x2, C is 3x2
        let a = mat![
            [1.0, 0.5, 0.0],
            [0.0, 2.0, 0.3],
            [0.0, 0.0, 3.0f64]
        ];
        let b = mat![[4.0, 0.5], [0.0, 5.0f64]];
        let c_orig = mat![
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0f64]
        ];
        let mut c = c_orig.clone();
        let (scale, _ns) = solve_triangular_sylvester(
            a.as_ref(),
            b.as_ref(),
            c.as_mut(),
            1.0,
        );

        let x = &c;
        let residual = &a * x + x * &b - scale * &c_orig;
        let mut max_err = 0.0f64;
        for j in 0..2 {
            for i in 0..3 {
                max_err = max_err.max(residual[(i, j)].abs());
            }
        }
        assert!(max_err < 1e-10, "Residual too large: {}", max_err);
    }
}
