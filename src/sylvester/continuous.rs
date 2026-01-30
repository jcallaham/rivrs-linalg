//! Continuous-time Sylvester equation solver.
//!
//! Solves `AX + XB = C` for X using the Bartels-Stewart algorithm
//! based on real Schur decomposition.
//!
//! # Algorithm (Bartels-Stewart)
//!
//! 1. Compute real Schur decompositions: `A = U₁ T U₁ᵀ`, `B = U₂ S U₂ᵀ`
//! 2. Transform RHS: `F = U₁ᵀ C U₂`
//! 3. Solve triangular system: `TY + YS = F` (back-substitution)
//! 4. Back-transform solution: `X = U₁ Y U₂ᵀ`
//!
//! # References
//!
//! - Bartels & Stewart (1972), "Solution of the Matrix Equation AX + XB = C",
//!   CACM 15(9):820-826
//! - Golub & Van Loan (2013), "Matrix Computations" (4th Ed), Algorithm 7.6.2

use faer::prelude::*;
use faer::Accum;

use crate::error::SylvesterError;
use super::types::SylvesterSolution;
use super::triangular::solve_triangular_sylvester;
use super::utils::compute_residual;
use super::types::EquationType;
use super::validation::{validate_dimensions, validate_finite};

/// Solves the continuous-time Sylvester equation `AX + XB = C`.
///
/// Uses the Bartels-Stewart algorithm with real Schur decomposition.
/// The algorithm is numerically stable and runs in O(n³ + m³ + n²m) time.
///
/// # Arguments
///
/// - `a`: Square matrix A (n × n)
/// - `b`: Square matrix B (m × m)
/// - `c`: Right-hand side matrix C (n × m)
///
/// # Returns
///
/// A `SylvesterSolution` containing the solution matrix X, a scale factor,
/// residual norm, and near-singular indicator.
///
/// # Errors
///
/// - `DimensionMismatch` if matrix dimensions are incompatible
/// - `NotSquare` if A or B is not square
/// - `InvalidInput` if matrices contain NaN or Inf
/// - `ConvergenceFailure` if Schur decomposition fails to converge
///
/// # Example
///
/// ```rust
/// use csrrs::sylvester::solve_continuous;
/// use faer::prelude::*;
///
/// let a = mat![[1.0, 0.5], [0.0, -1.0f64]];
/// let b = mat![[2.0, 0.0], [0.0, 3.0f64]];
/// let c = mat![[1.0, 2.0], [3.0, 4.0f64]];
///
/// let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
/// let x = &result.solution * (1.0 / result.scale);
/// assert!(result.residual_norm < 1e-10);
/// ```
///
/// # References
///
/// - Bartels & Stewart (1972), CACM 15(9):820-826
/// - Golub & Van Loan (2013), Algorithm 7.6.2
pub fn solve_continuous(
    a: MatRef<'_, f64>,
    b: MatRef<'_, f64>,
    c: MatRef<'_, f64>,
) -> Result<SylvesterSolution, SylvesterError> {
    // Input validation
    validate_dimensions(a, b, c)?;
    validate_finite(a, 'A')?;
    validate_finite(b, 'B')?;
    validate_finite(c, 'C')?;

    let n = a.nrows();
    let m = b.nrows();

    // Handle degenerate cases
    if n == 0 || m == 0 {
        return Ok(SylvesterSolution {
            solution: Mat::zeros(n, m),
            scale: 1.0,
            residual_norm: 0.0,
            near_singular: false,
        });
    }

    // Step 1: Compute real Schur decompositions
    // A = U1 * T * U1^T
    // B = U2 * S * U2^T
    let (schur_t, u1) = compute_real_schur(a)?;
    let (schur_s, u2) = compute_real_schur(b)?;

    // Step 2: Transform RHS: F = U1^T * C * U2
    let mut f = Mat::zeros(n, m);
    {
        // tmp = U1^T * C
        let mut tmp = Mat::zeros(n, m);
        faer::linalg::matmul::matmul(
            tmp.as_mut(),
            Accum::Replace,
            u1.as_ref().transpose(),
            c,
            1.0f64,
            Par::Seq,
        );
        // F = tmp * U2
        faer::linalg::matmul::matmul(
            f.as_mut(),
            Accum::Replace,
            tmp.as_ref(),
            u2.as_ref(),
            1.0f64,
            Par::Seq,
        );
    }

    // Step 3: Solve triangular system: T*Y + Y*S = F
    // sgn = +1 for AX + XB = C
    let (scale, near_singular) = solve_triangular_sylvester(
        schur_t.as_ref(),
        schur_s.as_ref(),
        f.as_mut(),
        1.0,
    );
    // f now contains Y

    // Step 4: Back-transform: X = U1 * Y * U2^T
    let mut solution = Mat::zeros(n, m);
    {
        // tmp = U1 * Y
        let mut tmp = Mat::zeros(n, m);
        faer::linalg::matmul::matmul(
            tmp.as_mut(),
            Accum::Replace,
            u1.as_ref(),
            f.as_ref(),
            1.0f64,
            Par::Seq,
        );
        // X = tmp * U2^T
        faer::linalg::matmul::matmul(
            solution.as_mut(),
            Accum::Replace,
            tmp.as_ref(),
            u2.as_ref().transpose(),
            1.0f64,
            Par::Seq,
        );
    }

    // Step 5: Compute residual norm
    let x_true = &solution * (1.0 / scale);
    let residual_norm = compute_residual(a, b, c, x_true.as_ref(), EquationType::Continuous);

    Ok(SylvesterSolution {
        solution,
        scale,
        residual_norm,
        near_singular,
    })
}

/// Solves the continuous-time Sylvester equation when A and B are already
/// in Schur form.
///
/// This is an advanced API for users who have pre-computed Schur decompositions
/// (e.g., when solving multiple equations with the same A or B).
///
/// # Arguments
///
/// - `schur_a`: Matrix A in real Schur form (quasi-triangular)
/// - `schur_b`: Matrix B in real Schur form (quasi-triangular)
/// - `u`: Orthogonal Schur vectors for A (A = U * schur_a * U^T)
/// - `v`: Orthogonal Schur vectors for B (B = V * schur_b * V^T)
/// - `c`: Right-hand side matrix C
///
/// # References
///
/// - Golub & Van Loan (2013), Algorithm 7.6.2
pub fn solve_continuous_schur(
    schur_a: MatRef<'_, f64>,
    schur_b: MatRef<'_, f64>,
    u: MatRef<'_, f64>,
    v: MatRef<'_, f64>,
    c: MatRef<'_, f64>,
) -> Result<SylvesterSolution, SylvesterError> {
    let n = schur_a.nrows();
    let m = schur_b.nrows();

    // Validate Schur form
    super::validation::validate_quasi_triangular(schur_a, 'T')?;
    super::validation::validate_quasi_triangular(schur_b, 'S')?;

    // Transform RHS: F = U^T * C * V
    let mut f = Mat::zeros(n, m);
    {
        let mut tmp = Mat::zeros(n, m);
        faer::linalg::matmul::matmul(
            tmp.as_mut(), Accum::Replace, u.transpose(), c, 1.0f64, Par::Seq,
        );
        faer::linalg::matmul::matmul(
            f.as_mut(), Accum::Replace, tmp.as_ref(), v, 1.0f64, Par::Seq,
        );
    }

    // Solve triangular system
    let (scale, near_singular) = solve_triangular_sylvester(
        schur_a, schur_b, f.as_mut(), 1.0,
    );

    // Back-transform: X = U * Y * V^T
    let mut solution = Mat::zeros(n, m);
    {
        let mut tmp = Mat::zeros(n, m);
        faer::linalg::matmul::matmul(
            tmp.as_mut(), Accum::Replace, u, f.as_ref(), 1.0f64, Par::Seq,
        );
        faer::linalg::matmul::matmul(
            solution.as_mut(), Accum::Replace, tmp.as_ref(), v.transpose(), 1.0f64, Par::Seq,
        );
    }

    // Reconstruct original A and B for residual computation
    let mut a_orig = Mat::zeros(n, n);
    {
        let mut tmp = Mat::zeros(n, n);
        faer::linalg::matmul::matmul(
            tmp.as_mut(), Accum::Replace, u, schur_a, 1.0f64, Par::Seq,
        );
        faer::linalg::matmul::matmul(
            a_orig.as_mut(), Accum::Replace, tmp.as_ref(), u.transpose(), 1.0f64, Par::Seq,
        );
    }
    let mut b_orig = Mat::zeros(m, m);
    {
        let mut tmp = Mat::zeros(m, m);
        faer::linalg::matmul::matmul(
            tmp.as_mut(), Accum::Replace, v, schur_b, 1.0f64, Par::Seq,
        );
        faer::linalg::matmul::matmul(
            b_orig.as_mut(), Accum::Replace, tmp.as_ref(), v.transpose(), 1.0f64, Par::Seq,
        );
    }

    let x_true = &solution * (1.0 / scale);
    let residual_norm = compute_residual(
        a_orig.as_ref(), b_orig.as_ref(), c, x_true.as_ref(), EquationType::Continuous,
    );

    Ok(SylvesterSolution {
        solution,
        scale,
        residual_norm,
        near_singular,
    })
}

/// Computes the real Schur decomposition of a matrix.
///
/// Returns (schur_form, schur_vectors) where:
/// - `schur_form` is the quasi-triangular Schur form T
/// - `schur_vectors` is the orthogonal matrix U such that A = U * T * U^T
///
/// Uses nalgebra's Schur decomposition internally, converting between
/// faer and nalgebra matrix types.
fn compute_real_schur(
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

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;

    #[test]
    fn test_solve_continuous_2x2_diagonal() {
        // A = diag(1, -1), B = diag(2, 3), C = I
        // Solution: x11 = 1/(1+2) = 1/3, x12 = 0, x21 = 0, x22 = 1/(-1+3) = 1/2
        let a = mat![[1.0, 0.0], [0.0, -1.0f64]];
        let b = mat![[2.0, 0.0], [0.0, 3.0f64]];
        let c = mat![[1.0, 0.0], [0.0, 1.0f64]];

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        let x = &result.solution * (1.0 / result.scale);

        assert!((x[(0, 0)] - 1.0 / 3.0).abs() < 1e-10);
        assert!(x[(0, 1)].abs() < 1e-10);
        assert!(x[(1, 0)].abs() < 1e-10);
        assert!((x[(1, 1)] - 0.5).abs() < 1e-10);
        assert!(result.residual_norm < 1e-10);
    }

    #[test]
    fn test_solve_continuous_2x2_general() {
        let a = mat![[1.0, 2.0], [3.0, 4.0f64]];
        let b = mat![[5.0, 6.0], [7.0, 8.0f64]];
        let c = mat![[1.0, 2.0], [3.0, 4.0f64]];

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(result.residual_norm < 1e-10, "Residual: {}", result.residual_norm);
    }

    #[test]
    fn test_solve_continuous_3x3() {
        let a = mat![
            [1.0, 2.0, 0.0],
            [0.0, 3.0, 1.0],
            [0.0, 0.0, 5.0f64]
        ];
        let b = mat![
            [2.0, 1.0, 0.5],
            [0.0, 4.0, 0.0],
            [0.0, 0.0, 6.0f64]
        ];
        let c = mat![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 9.0f64]
        ];

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(result.residual_norm < 1e-10, "Residual: {}", result.residual_norm);
    }

    #[test]
    fn test_solve_continuous_rectangular_c() {
        // A is 3x3, B is 2x2, C is 3x2
        let a = mat![
            [1.0, 0.5, 0.0],
            [0.0, 2.0, 0.3],
            [0.0, 0.0, 3.0f64]
        ];
        let b = mat![[4.0, 0.5], [0.0, 5.0f64]];
        let c = mat![
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0f64]
        ];

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(result.residual_norm < 1e-10, "Residual: {}", result.residual_norm);
    }

    #[test]
    fn test_solve_continuous_identity_rhs() {
        let a = mat![[1.0, 0.0], [0.0, 2.0f64]];
        let b = mat![[3.0, 0.0], [0.0, 4.0f64]];
        let c = Mat::<f64>::identity(2, 2);

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(result.residual_norm < 1e-10);
    }

    #[test]
    fn test_solve_continuous_empty() {
        let a = Mat::<f64>::zeros(0, 0);
        let b = Mat::zeros(0, 0);
        let c = Mat::zeros(0, 0);

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert_eq!(result.solution.nrows(), 0);
        assert_eq!(result.solution.ncols(), 0);
    }

    #[test]
    fn test_solve_continuous_dimension_mismatch() {
        let a = Mat::<f64>::zeros(2, 2);
        let b = Mat::zeros(3, 3);
        let c = Mat::zeros(2, 2); // Should be 2x3
        assert!(solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).is_err());
    }

    #[test]
    fn test_solve_continuous_nan_input() {
        let a = mat![[f64::NAN, 0.0], [0.0, 1.0]];
        let b = mat![[1.0, 0.0], [0.0, 1.0f64]];
        let c = mat![[1.0, 0.0], [0.0, 1.0f64]];
        assert!(solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).is_err());
    }

    #[test]
    fn test_solve_continuous_4x4() {
        // Larger test with complex eigenvalues
        let a = mat![
            [1.0, 2.0, 0.0, 0.0],
            [-1.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 3.0, 1.0],
            [0.0, 0.0, 0.0, 4.0f64]
        ];
        let b = mat![
            [5.0, 1.0, 0.0],
            [0.0, 6.0, 0.5],
            [0.0, 0.0, 7.0f64]
        ];
        let c = Mat::from_fn(4, 3, |i, j| (i + j + 1) as f64);

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(result.residual_norm < 1e-10, "Residual: {}", result.residual_norm);
    }

    #[test]
    fn test_solve_continuous_1x1() {
        // 3*x + x*5 = 16 => x = 2
        let a = mat![[3.0f64]];
        let b = mat![[5.0f64]];
        let c = mat![[16.0f64]];

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        let x = &result.solution * (1.0 / result.scale);
        assert!((x[(0, 0)] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_solve_continuous_dense_5x5() {
        // Dense 5x5 test
        let a = mat![
            [2.0, 1.0, 0.5, 0.0, 0.0],
            [-1.0, 2.0, 0.0, 0.5, 0.0],
            [0.0, 0.0, 3.0, 1.0, 0.0],
            [0.0, 0.0, -1.0, 3.0, 0.5],
            [0.0, 0.0, 0.0, 0.0, 5.0f64]
        ];
        let b = mat![
            [6.0, 1.0, 0.0, 0.0],
            [0.0, 7.0, 0.5, 0.0],
            [0.0, 0.0, 8.0, 1.0],
            [0.0, 0.0, 0.0, 9.0f64]
        ];
        let c = Mat::from_fn(5, 4, |i, j| ((i * 4 + j + 1) as f64) * 0.1);

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(result.residual_norm < 1e-10, "Residual: {}", result.residual_norm);
    }

    #[test]
    fn test_compute_real_schur() {
        // Test that Schur decomposition is correct: A = U * T * U^T
        let a = mat![
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [7.0, 8.0, 10.0f64]
        ];

        let (t, u) = compute_real_schur(a.as_ref()).unwrap();

        // Reconstruct: A_recon = U * T * U^T
        let mut tmp = Mat::zeros(3, 3);
        faer::linalg::matmul::matmul(
            tmp.as_mut(), Accum::Replace, u.as_ref(), t.as_ref(), 1.0f64, Par::Seq,
        );
        let mut a_recon = Mat::zeros(3, 3);
        faer::linalg::matmul::matmul(
            a_recon.as_mut(), Accum::Replace, tmp.as_ref(), u.as_ref().transpose(), 1.0f64, Par::Seq,
        );

        let mut max_err = 0.0f64;
        for j in 0..3 {
            for i in 0..3 {
                max_err = max_err.max((a_recon[(i, j)] - a[(i, j)]).abs());
            }
        }
        assert!(max_err < 1e-12, "Schur reconstruction error: {}", max_err);
    }
}
