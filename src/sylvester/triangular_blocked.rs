//! Blocked triangular Sylvester equation solver (dtrsyl3-style).
//!
//! Solves `A*X + sgn*X*B = scale*C` where A and B are quasi-triangular
//! using a blocked algorithm with Level-3 BLAS (matmul) operations.
//!
//! For large matrices, this is significantly faster than the unblocked
//! version because it makes better use of cache by processing larger
//! blocks of the solution matrix at once.
//!
//! # Algorithm
//!
//! The matrix is partitioned into panels of size `block_size`. Within
//! each block, the unblocked solver is used. Between blocks, updates
//! are performed using matrix-matrix multiplication (Level-3 BLAS).
//!
//! # References
//!
//! - Jonsson & Kågström (2002), "Recursive blocked algorithms for solving
//!   triangular systems — Part II: two-sided and generalized Sylvester
//!   and Lyapunov matrix equations", ACM TOMS 28(4):416-435
//! - LAPACK dtrsyl3 for the blocked Sylvester solver pattern

use faer::prelude::*;
use faer::Accum;

use super::triangular::solve_triangular_sylvester;

/// Default block size for the blocked algorithm.
/// Chosen to fit in L1/L2 cache while being large enough
/// to amortize the overhead of Level-3 BLAS calls.
const DEFAULT_BLOCK_SIZE: usize = 64;

/// Minimum matrix size to use the blocked algorithm.
/// Below this threshold, the unblocked algorithm is faster
/// due to lower overhead.
pub const BLOCKED_THRESHOLD: usize = 64;

/// Solves the triangular Sylvester equation using a blocked algorithm.
///
/// Solves `A*X + sgn*X*B = scale*C` where A (m × m) and B (n × n) are
/// upper quasi-triangular (real Schur form). C (m × n) is overwritten
/// with the solution X. Returns `(scale, near_singular)`.
///
/// This function partitions A and B into blocks and uses Level-3 BLAS
/// (matrix-matrix multiplication) for inter-block updates, which is
/// significantly faster for large matrices.
///
/// # References
///
/// - Jonsson & Kågström (2002), ACM TOMS 28(4):416-435
pub fn solve_triangular_sylvester_blocked(
    a: MatRef<'_, f64>,
    b: MatRef<'_, f64>,
    c: MatMut<'_, f64>,
    sgn: f64,
) -> (f64, bool) {
    let m = a.nrows();
    let n = b.nrows();

    // Fall back to unblocked for small matrices
    if m <= BLOCKED_THRESHOLD && n <= BLOCKED_THRESHOLD {
        return solve_triangular_sylvester(a, b, c, sgn);
    }

    solve_blocked_impl(a, b, c, sgn, DEFAULT_BLOCK_SIZE)
}

/// Internal implementation of the blocked solver.
fn solve_blocked_impl(
    a: MatRef<'_, f64>,
    b: MatRef<'_, f64>,
    mut c: MatMut<'_, f64>,
    sgn: f64,
    block_size: usize,
) -> (f64, bool) {
    let m = a.nrows();
    let n = b.nrows();
    let mut scale = 1.0f64;
    let mut near_singular = false;

    if m == 0 || n == 0 {
        return (scale, near_singular);
    }

    // Compute block boundaries for A (row blocks) and B (column blocks).
    // We need to align block boundaries to not split 2x2 quasi-triangular blocks.
    let a_panels = compute_panel_boundaries(a, m, block_size);
    let b_panels = compute_panel_boundaries(b, n, block_size);

    // Process columns of B left-to-right (L panels), rows of A bottom-to-top (K panels)
    for &(l_start, l_end) in &b_panels {
        let l_size = l_end - l_start;

        for k_panel in (0..a_panels.len()).rev() {
            let (k_start, k_end) = a_panels[k_panel];
            let k_size = k_end - k_start;

            // Update C(K, L) by subtracting contributions from already-solved panels

            // Contribution from A panels below K (already solved going bottom-up)
            if k_end < m {
                let a_sub = a.submatrix(k_start, k_end, k_size, m - k_end);
                let rows_below = m - k_end;

                // Copy X(below K, L) to temporary to avoid borrow conflicts
                let mut x_tmp = Mat::zeros(rows_below, l_size);
                for jj in 0..l_size {
                    for ii in 0..rows_below {
                        x_tmp[(ii, jj)] = c[(k_end + ii, l_start + jj)];
                    }
                }

                // C(K, L) -= A(K, below_K) * X(below_K, L)
                let mut c_sub = c.rb_mut().submatrix_mut(k_start, l_start, k_size, l_size);
                faer::linalg::matmul::matmul(
                    c_sub.rb_mut(),
                    Accum::Add,
                    a_sub,
                    x_tmp.as_ref(),
                    -1.0f64,
                    Par::Seq,
                );
            }

            // Contribution from B panels left of L (already solved going left-to-right)
            if l_start > 0 {
                let b_sub = b.submatrix(0, l_start, l_start, l_size);

                // Copy X(K, left_of_L) to temporary
                let mut x_tmp = Mat::zeros(k_size, l_start);
                for jj in 0..l_start {
                    for ii in 0..k_size {
                        x_tmp[(ii, jj)] = c[(k_start + ii, jj)];
                    }
                }

                // C(K, L) -= sgn * X(K, left_L) * B(left_L, L)
                let mut c_sub = c.rb_mut().submatrix_mut(k_start, l_start, k_size, l_size);
                faer::linalg::matmul::matmul(
                    c_sub.rb_mut(),
                    Accum::Add,
                    x_tmp.as_ref(),
                    b_sub,
                    -sgn,
                    Par::Seq,
                );
            }

            // Solve the diagonal block using unblocked algorithm:
            // A(K,K) * X(K,L) + sgn * X(K,L) * B(L,L) = C(K,L)
            let a_kk = a.submatrix(k_start, k_start, k_size, k_size);
            let b_ll = b.submatrix(l_start, l_start, l_size, l_size);
            let c_kl = c.rb_mut().submatrix_mut(k_start, l_start, k_size, l_size);

            let (scaloc, ns) = solve_triangular_sylvester(a_kk, b_ll, c_kl, sgn);

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

/// Compute panel boundaries that respect 2x2 quasi-triangular block structure.
///
/// Returns a list of (start, end) index pairs (exclusive end) for each panel.
/// The boundaries are adjusted so that no 2x2 block is split across panels.
fn compute_panel_boundaries(
    mat: MatRef<'_, f64>,
    n: usize,
    block_size: usize,
) -> Vec<(usize, usize)> {
    let mut panels = Vec::new();
    let mut i = 0;

    while i < n {
        let mut end = (i + block_size).min(n);

        // If we're splitting a 2x2 block, adjust the boundary
        if end < n && end > 0 && mat[(end, end - 1)] != 0.0 {
            // The (end, end-1) sub-diagonal is nonzero, meaning indices
            // (end-1, end) form a 2x2 block. Include it in this panel.
            end += 1;
        }

        panels.push((i, end));
        i = end;
    }

    panels
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_blocked_small_fallback() {
        // Small matrix should use unblocked path and give correct results
        let a = mat![[2.0, 1.0], [0.0, 3.0f64]];
        let b = mat![[4.0, 0.5], [0.0, 5.0f64]];
        let c_orig = mat![[10.0, 11.0], [12.0, 13.0f64]];
        let mut c = c_orig.clone();

        let (scale, _ns) =
            solve_triangular_sylvester_blocked(a.as_ref(), b.as_ref(), c.as_mut(), 1.0);

        // Verify: A*X + X*B = scale*C_orig
        let x = &c;
        let residual = &a * x + x * &b - scale * &c_orig;
        let mut max_err = 0.0f64;
        for j in 0..2 {
            for i in 0..2 {
                max_err = max_err.max(residual[(i, j)].abs());
            }
        }
        assert!(max_err < 1e-12, "Residual too large: {}", max_err);
    }

    #[test]
    fn test_blocked_matches_unblocked() {
        // Generate a non-trivial upper triangular system and verify
        // blocked and unblocked give the same result
        let n = 8;
        let mut a = Mat::zeros(n, n);
        let mut b = Mat::zeros(n, n);
        for i in 0..n {
            a[(i, i)] = (i + 1) as f64;
            b[(i, i)] = (i + n + 1) as f64;
            if i + 1 < n {
                a[(i, i + 1)] = 0.5;
                b[(i, i + 1)] = 0.3;
            }
        }
        let c_orig = Mat::from_fn(n, n, |i, j| ((i * n + j + 1) as f64) * 0.1);

        // Solve with unblocked
        let mut c_unblocked = c_orig.clone();
        let (scale_u, _) =
            solve_triangular_sylvester(a.as_ref(), b.as_ref(), c_unblocked.as_mut(), 1.0);

        // Solve with blocked (force small block size to exercise the blocked path)
        let mut c_blocked = c_orig.clone();
        let (scale_b, _) = solve_blocked_impl(a.as_ref(), b.as_ref(), c_blocked.as_mut(), 1.0, 3);

        // Results should match
        let mut max_diff = 0.0f64;
        for j in 0..n {
            for i in 0..n {
                let xu = c_unblocked[(i, j)] / scale_u;
                let xb = c_blocked[(i, j)] / scale_b;
                max_diff = max_diff.max((xu - xb).abs());
            }
        }
        assert!(
            max_diff < 1e-10,
            "Blocked and unblocked differ by {}",
            max_diff
        );
    }

    #[test]
    fn test_blocked_with_2x2_blocks() {
        // Quasi-triangular with 2x2 blocks
        let a = mat![
            [1.0, 2.0, 0.0, 0.0, 0.0, 0.0],
            [-0.5, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 3.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 4.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, -0.3, 4.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 5.0f64]
        ];
        let b = mat![
            [6.0, 0.5, 0.0, 0.0, 0.0, 0.0],
            [0.0, 7.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 8.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, -0.4, 8.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 9.0, 0.3],
            [0.0, 0.0, 0.0, 0.0, 0.0, 10.0f64]
        ];
        let c_orig = Mat::from_fn(6, 6, |i, j| ((i + 1) * (j + 1)) as f64);

        let mut c_blocked = c_orig.clone();
        let (scale, _ns) = solve_blocked_impl(a.as_ref(), b.as_ref(), c_blocked.as_mut(), 1.0, 3);

        // Verify residual
        let x = &c_blocked;
        let residual = &a * x + x * &b - scale * &c_orig;
        let mut max_err = 0.0f64;
        for j in 0..6 {
            for i in 0..6 {
                max_err = max_err.max(residual[(i, j)].abs());
            }
        }
        assert!(max_err < 1e-10, "Residual too large: {}", max_err);
    }

    #[test]
    fn test_blocked_negative_sign() {
        // Use well-separated eigenvalues to avoid near-singular issues
        // A eigenvalues: 10, 20, 30, 40, 50, 60 (multiples of 10)
        // B eigenvalues: 1, 2, 3, 4, 5, 6
        // A(i,i) - B(j,j) ranges from 4 to 59, all well-separated
        let n = 6;
        let mut a = Mat::zeros(n, n);
        let mut b = Mat::zeros(n, n);
        for i in 0..n {
            a[(i, i)] = ((i + 1) * 10) as f64;
            b[(i, i)] = (i + 1) as f64;
            if i + 1 < n {
                a[(i, i + 1)] = 0.3;
                b[(i, i + 1)] = 0.2;
            }
        }
        let c_orig = Mat::from_fn(n, n, |i, j| (i + j + 1) as f64);

        // Test blocked
        let mut c_blocked = c_orig.clone();
        let (scale_b, _ns) =
            solve_blocked_impl(a.as_ref(), b.as_ref(), c_blocked.as_mut(), -1.0, 3);

        // Verify: A*X - X*B = scale*C_orig
        let x = &c_blocked;
        let residual = &a * x - x * &b - scale_b * &c_orig;
        let mut max_err = 0.0f64;
        for j in 0..n {
            for i in 0..n {
                max_err = max_err.max(residual[(i, j)].abs());
            }
        }
        assert!(max_err < 1e-10, "Residual too large: {}", max_err);

        // Also verify that blocked matches unblocked
        let mut c_unblocked = c_orig.clone();
        let (scale_u, _) =
            solve_triangular_sylvester(a.as_ref(), b.as_ref(), c_unblocked.as_mut(), -1.0);
        let mut max_diff = 0.0f64;
        for j in 0..n {
            for i in 0..n {
                let xu = c_unblocked[(i, j)] / scale_u;
                let xb = c_blocked[(i, j)] / scale_b;
                max_diff = max_diff.max((xu - xb).abs());
            }
        }
        assert!(
            max_diff < 1e-10,
            "Blocked and unblocked differ by {}",
            max_diff
        );
    }

    #[test]
    fn test_panel_boundaries_simple() {
        // 8x8 upper triangular, block_size=3 → panels: [0,3), [3,6), [6,8)
        let m = Mat::zeros(8, 8);
        let panels = compute_panel_boundaries(m.as_ref(), 8, 3);
        assert_eq!(panels, vec![(0, 3), (3, 6), (6, 8)]);
    }

    #[test]
    fn test_panel_boundaries_with_2x2_block() {
        // 6x6 with 2x2 block at positions (2,3): sub-diagonal (3,2) is nonzero
        // block_size=3 would split at index 3, but (3,2) is nonzero
        // so boundary should shift to 4 to keep the 2x2 block intact
        let mut m = Mat::zeros(6, 6);
        m[(3, 2)] = 1.0; // 2x2 block at (2,3)
        let panels = compute_panel_boundaries(m.as_ref(), 6, 3);
        // First panel [0,4) includes the 2x2 block at (2,3)
        // Second panel [4,6)
        assert_eq!(panels, vec![(0, 4), (4, 6)]);
    }
}
