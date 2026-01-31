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

use super::condition::{estimate_separation, SEPARATION_THRESHOLD};
use super::triangular::solve_triangular_sylvester;
use super::triangular_blocked::{solve_triangular_sylvester_blocked, BLOCKED_THRESHOLD};
use super::types::EquationType;
use super::types::SylvesterSolution;
use super::utils::{
    back_transform, compute_real_schur, compute_residual, reconstruct_from_schur, transform_rhs,
};
use super::validation::{validate_dimensions, validate_finite};
use crate::error::SylvesterError;

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
    // A = U1 * T * U1^T, B = U2 * S * U2^T
    let (schur_t, u1) = compute_real_schur(a)?;
    let (schur_s, u2) = compute_real_schur(b)?;

    // Step 1b: Check eigenvalue separation
    let sep_est = estimate_separation(schur_t.as_ref(), schur_s.as_ref(), EquationType::Continuous);
    if sep_est.is_singular {
        return Err(SylvesterError::CommonEigenvalues {
            separation: sep_est.separation,
            threshold: SEPARATION_THRESHOLD,
        });
    }

    // Step 2: Transform RHS: F = U1^T * C * U2
    let mut f = transform_rhs(u1.as_ref(), c, u2.as_ref());

    // Step 3: Solve triangular system: T*Y + Y*S = F
    // Use blocked solver for large matrices for better cache performance
    let (scale, near_singular) = if n > BLOCKED_THRESHOLD || m > BLOCKED_THRESHOLD {
        solve_triangular_sylvester_blocked(schur_t.as_ref(), schur_s.as_ref(), f.as_mut(), 1.0)
    } else {
        solve_triangular_sylvester(schur_t.as_ref(), schur_s.as_ref(), f.as_mut(), 1.0)
    };

    // Step 4: Back-transform: X = U1 * Y * U2^T
    let solution = back_transform(u1.as_ref(), f.as_ref(), u2.as_ref());

    // Step 5: Compute residual norm
    let x_true = &solution * (1.0 / scale);
    let residual_norm = compute_residual(a, b, c, x_true.as_ref(), EquationType::Continuous);

    Ok(SylvesterSolution {
        solution,
        scale,
        residual_norm,
        near_singular: near_singular || sep_est.is_near_singular,
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
    // Validate Schur form
    super::validation::validate_quasi_triangular(schur_a, 'T')?;
    super::validation::validate_quasi_triangular(schur_b, 'S')?;

    // Transform RHS: F = U^T * C * V
    let mut f = transform_rhs(u, c, v);

    // Solve triangular system
    let (scale, near_singular) = solve_triangular_sylvester(schur_a, schur_b, f.as_mut(), 1.0);

    // Back-transform: X = U * Y * V^T
    let solution = back_transform(u, f.as_ref(), v);

    // Reconstruct original A and B for residual computation
    let a_orig = reconstruct_from_schur(u, schur_a);
    let b_orig = reconstruct_from_schur(v, schur_b);

    let x_true = &solution * (1.0 / scale);
    let residual_norm = compute_residual(
        a_orig.as_ref(),
        b_orig.as_ref(),
        c,
        x_true.as_ref(),
        EquationType::Continuous,
    );

    Ok(SylvesterSolution {
        solution,
        scale,
        residual_norm,
        near_singular,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::mat;
    use faer::Mat;

    /// Regression test: continuous solver produces inaccurate solutions for
    /// random dense matrices at n=10 and above. The relative residual
    /// ||AX + XB - C||_F / ||C||_F is O(1) instead of O(machine epsilon).
    ///
    /// Root cause: `solve_2x2_linear_system` in `triangular.rs` incorrectly
    /// swaps the two solution values when partial pivoting reorders rows.
    /// Row pivoting only reorders equations, not unknowns, so the solution
    /// vector should not be permuted. This is triggered whenever the Schur
    /// forms of A or B contain 2x2 blocks (complex eigenvalue pairs) and the
    /// off-diagonal element is larger than the shifted diagonal, which is
    /// common for random matrices.
    ///
    /// The discrete solver is unaffected because it uses a different code path
    /// (Kronecker product vectorization with ipiv-based pivoting) that handles
    /// the permutation correctly.
    #[test]
    fn test_solve_continuous_random_10x10_accuracy() {
        use rand::prelude::*;
        use rand_distr::StandardNormal;

        let mut rng = StdRng::seed_from_u64(42);
        let n = 10;

        // A = randn(n,n) + I, B = randn(n,n) + 5I, C = randn(n,n)
        let mut a = Mat::from_fn(n, n, |_, _| rng.sample::<f64, _>(StandardNormal));
        let mut b = Mat::from_fn(n, n, |_, _| rng.sample::<f64, _>(StandardNormal));
        let c = Mat::from_fn(n, n, |_, _| rng.sample::<f64, _>(StandardNormal));
        for i in 0..n {
            a[(i, i)] += 1.0;
            b[(i, i)] += 5.0;
        }

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        let x = &result.solution * (1.0 / result.scale);

        // Compute relative residual: ||AX + XB - C|| / ||C||
        let residual_norm = super::super::utils::compute_residual(
            a.as_ref(),
            b.as_ref(),
            c.as_ref(),
            x.as_ref(),
            super::super::types::EquationType::Continuous,
        );
        let c_norm = super::super::utils::frobenius_norm(c.as_ref());
        let relative_residual = residual_norm / c_norm;

        assert!(
            relative_residual < 1e-10,
            "Continuous solver 10x10 relative residual too large: {:.2e} (absolute: {:.2e})",
            relative_residual,
            residual_norm,
        );
    }

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
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
    }

    #[test]
    fn test_solve_continuous_3x3() {
        let a = mat![[1.0, 2.0, 0.0], [0.0, 3.0, 1.0], [0.0, 0.0, 5.0f64]];
        let b = mat![[2.0, 1.0, 0.5], [0.0, 4.0, 0.0], [0.0, 0.0, 6.0f64]];
        let c = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0f64]];

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
    }

    #[test]
    fn test_solve_continuous_rectangular_c() {
        // A is 3x3, B is 2x2, C is 3x2
        let a = mat![[1.0, 0.5, 0.0], [0.0, 2.0, 0.3], [0.0, 0.0, 3.0f64]];
        let b = mat![[4.0, 0.5], [0.0, 5.0f64]];
        let c = mat![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0f64]];

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
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
        let b = mat![[5.0, 1.0, 0.0], [0.0, 6.0, 0.5], [0.0, 0.0, 7.0f64]];
        let c = Mat::from_fn(4, 3, |i, j| (i + j + 1) as f64);

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
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
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
    }

    #[test]
    fn test_solve_continuous_common_eigenvalues() {
        // A has eigenvalue 2, B has eigenvalue -2
        // A*X + X*B = C is singular when λ(A) + λ(B) = 0
        let a = mat![[2.0f64]];
        let b = mat![[-2.0f64]];
        let c = mat![[1.0f64]];

        let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref());
        assert!(result.is_err());
        match result.unwrap_err() {
            crate::error::SylvesterError::CommonEigenvalues { .. } => {}
            other => panic!("Expected CommonEigenvalues, got {:?}", other),
        }
    }

    #[test]
    fn test_compute_real_schur() {
        // Test that Schur decomposition is correct: A = U * T * U^T
        let a = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0f64]];

        let (t, u) = compute_real_schur(a.as_ref()).unwrap();

        // Reconstruct: A_recon = U * T * U^T
        let a_recon = reconstruct_from_schur(u.as_ref(), t.as_ref());

        let mut max_err = 0.0f64;
        for j in 0..3 {
            for i in 0..3 {
                max_err = max_err.max((a_recon[(i, j)] - a[(i, j)]).abs());
            }
        }
        assert!(max_err < 1e-12, "Schur reconstruction error: {}", max_err);
    }
}
