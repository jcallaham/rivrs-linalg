//! Discrete-time Sylvester equation solver.
//!
//! Solves `AXB + X = C` for X using a Bartels-Stewart-like algorithm
//! based on real Schur decomposition.
//!
//! # Algorithm
//!
//! 1. Compute real Schur decompositions: `A = U₁ T U₁ᵀ`, `B = U₂ S U₂ᵀ`
//! 2. Transform RHS: `F = U₁ᵀ C U₂`
//! 3. Substitute: `T (U₁ᵀ X U₂) S + (U₁ᵀ X U₂) = F`
//! 4. Let `Y = U₁ᵀ X U₂`, solve: `Y + T Y S = F`
//! 5. Back-transform solution: `X = U₁ Y U₂ᵀ`
//!
//! # References
//!
//! - Golub, Nash & Van Loan (1979), "A Hessenberg-Schur method for AX + XB = C",
//!   IEEE Trans. Auto. Contr. AC-24:909-913
//! - Sima (1996), "Algorithms for Linear-Quadratic Optimization"

use faer::prelude::*;

use super::condition::{estimate_separation, SEPARATION_THRESHOLD};
use super::triangular_discrete::solve_triangular_sylvester_discrete;
use super::types::EquationType;
use super::types::SylvesterSolution;
use super::utils::{
    back_transform, compute_real_schur, compute_residual, reconstruct_from_schur, transform_rhs,
};
use super::validation::{validate_dimensions, validate_finite};
use crate::error::SylvesterError;

/// Solves the discrete-time Sylvester equation `AXB + X = C`.
///
/// Uses a Bartels-Stewart-like algorithm with real Schur decomposition.
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
/// use csrrs::sylvester::solve_discrete;
/// use faer::prelude::*;
///
/// let a = mat![[0.5, 0.1], [0.0, 0.8f64]];
/// let b = mat![[0.6, 0.0], [0.0, 0.9f64]];
/// let c = mat![[1.0, 2.0], [3.0, 4.0f64]];
///
/// let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
/// let x = &result.solution * (1.0 / result.scale);
/// assert!(result.residual_norm < 1e-10);
/// ```
///
/// # References
///
/// - Sima (1996), "Algorithms for Linear-Quadratic Optimization"
pub fn solve_discrete(
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
    let sep_est = estimate_separation(schur_t.as_ref(), schur_s.as_ref(), EquationType::Discrete);
    if sep_est.is_singular {
        return Err(SylvesterError::CommonEigenvalues {
            separation: sep_est.separation,
            threshold: SEPARATION_THRESHOLD,
        });
    }

    // Step 2: Transform RHS: F = U1^T * C * U2
    let mut f = transform_rhs(u1.as_ref(), c, u2.as_ref());

    // Step 3: Solve triangular system: Y + T*Y*S = F
    let (scale, near_singular) =
        solve_triangular_sylvester_discrete(schur_t.as_ref(), schur_s.as_ref(), f.as_mut(), 1.0);

    // Step 4: Back-transform: X = U1 * Y * U2^T
    let solution = back_transform(u1.as_ref(), f.as_ref(), u2.as_ref());

    // Step 5: Compute residual norm
    let x_true = &solution * (1.0 / scale);
    let residual_norm = compute_residual(a, b, c, x_true.as_ref(), EquationType::Discrete);

    Ok(SylvesterSolution {
        solution,
        scale,
        residual_norm,
        near_singular: near_singular || sep_est.is_near_singular,
    })
}

/// Solves the discrete-time Sylvester equation when A and B are already
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
/// - Sima (1996), "Algorithms for Linear-Quadratic Optimization"
pub fn solve_discrete_schur(
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
    let (scale, near_singular) =
        solve_triangular_sylvester_discrete(schur_a, schur_b, f.as_mut(), 1.0);

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
        EquationType::Discrete,
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

    #[test]
    fn test_solve_discrete_2x2_diagonal() {
        // A = diag(0.5, 0.8), B = diag(0.6, 0.9), C = I
        // AXB + X = C => a_ii * x_ij * b_jj + x_ij = c_ij
        // x_ij = c_ij / (1 + a_ii * b_jj)
        let a = mat![[0.5, 0.0], [0.0, 0.8f64]];
        let b = mat![[0.6, 0.0], [0.0, 0.9f64]];
        let c = mat![[1.0, 0.0], [0.0, 1.0f64]];

        let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        let x = &result.solution * (1.0 / result.scale);

        // x11 = 1/(1+0.5*0.6) = 1/1.3
        // x22 = 1/(1+0.8*0.9) = 1/1.72
        assert!((x[(0, 0)] - 1.0 / 1.3).abs() < 1e-10);
        assert!(x[(0, 1)].abs() < 1e-10);
        assert!(x[(1, 0)].abs() < 1e-10);
        assert!((x[(1, 1)] - 1.0 / 1.72).abs() < 1e-10);
        assert!(result.residual_norm < 1e-10);
    }

    #[test]
    fn test_solve_discrete_2x2_general() {
        let a = mat![[0.5, 0.1], [0.2, 0.8f64]];
        let b = mat![[0.6, 0.3], [0.1, 0.9f64]];
        let c = mat![[1.0, 2.0], [3.0, 4.0f64]];

        let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
    }

    #[test]
    fn test_solve_discrete_3x3() {
        let a = mat![[0.5, 0.1, 0.0], [0.0, 0.8, 0.2], [0.0, 0.0, 0.3f64]];
        let b = mat![[0.6, 0.1, 0.0], [0.0, 0.9, 0.05], [0.0, 0.0, 0.4f64]];
        let c = mat![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0f64]];

        let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
    }

    #[test]
    fn test_solve_discrete_rectangular_c() {
        // A is 3x3, B is 2x2, C is 3x2
        let a = mat![[0.5, 0.1, 0.0], [0.0, 0.8, 0.2], [0.0, 0.0, 0.3f64]];
        let b = mat![[0.6, 0.1], [0.0, 0.9f64]];
        let c = mat![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0f64]];

        let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
    }

    #[test]
    fn test_solve_discrete_empty() {
        let a = Mat::<f64>::zeros(0, 0);
        let b = Mat::zeros(0, 0);
        let c = Mat::zeros(0, 0);

        let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert_eq!(result.solution.nrows(), 0);
        assert_eq!(result.solution.ncols(), 0);
    }

    #[test]
    fn test_solve_discrete_dimension_mismatch() {
        let a = Mat::<f64>::zeros(2, 2);
        let b = Mat::zeros(3, 3);
        let c = Mat::zeros(2, 2); // Should be 2x3
        assert!(solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).is_err());
    }

    #[test]
    fn test_solve_discrete_nan_input() {
        let a = mat![[f64::NAN, 0.0], [0.0, 1.0]];
        let b = mat![[1.0, 0.0], [0.0, 1.0f64]];
        let c = mat![[1.0, 0.0], [0.0, 1.0f64]];
        assert!(solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).is_err());
    }

    #[test]
    fn test_solve_discrete_1x1() {
        // 0.5*x*0.6 + x = 7 => x*(1+0.3) = 7 => x = 7/1.3
        let a = mat![[0.5f64]];
        let b = mat![[0.6f64]];
        let c = mat![[7.0f64]];

        let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        let x = &result.solution * (1.0 / result.scale);
        assert!((x[(0, 0)] - 7.0 / 1.3).abs() < 1e-10);
    }

    #[test]
    fn test_solve_discrete_4x4_complex_eigenvalues() {
        // 4x4 with complex eigenvalues
        let a = mat![
            [0.5, 0.3, 0.0, 0.0],
            [-0.2, 0.5, 0.0, 0.0],
            [0.0, 0.0, 0.7, 0.1],
            [0.0, 0.0, 0.0, 0.4f64]
        ];
        let b = mat![[0.6, 0.1, 0.0], [0.0, 0.8, 0.2], [0.0, 0.0, 0.5f64]];
        let c = Mat::from_fn(4, 3, |i, j| (i + j + 1) as f64);

        let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
    }

    #[test]
    fn test_solve_discrete_common_eigenvalues() {
        // A has eigenvalue 2, B has eigenvalue -0.5
        // AXB + X = C is singular when 1 + λ(A)*λ(B) = 0
        // 1 + 2*(-0.5) = 0
        let a = mat![[2.0f64]];
        let b = mat![[-0.5f64]];
        let c = mat![[1.0f64]];

        let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref());
        assert!(result.is_err());
        match result.unwrap_err() {
            crate::error::SylvesterError::CommonEigenvalues { .. } => {}
            other => panic!("Expected CommonEigenvalues, got {:?}", other),
        }
    }

    #[test]
    fn test_solve_discrete_dense_5x5() {
        let a = mat![
            [0.3, 0.1, 0.05, 0.0, 0.0],
            [-0.1, 0.3, 0.0, 0.05, 0.0],
            [0.0, 0.0, 0.5, 0.1, 0.0],
            [0.0, 0.0, -0.1, 0.5, 0.05],
            [0.0, 0.0, 0.0, 0.0, 0.7f64]
        ];
        let b = mat![
            [0.4, 0.1, 0.0, 0.0],
            [0.0, 0.6, 0.05, 0.0],
            [0.0, 0.0, 0.8, 0.1],
            [0.0, 0.0, 0.0, 0.2f64]
        ];
        let c = Mat::from_fn(5, 4, |i, j| ((i * 4 + j + 1) as f64) * 0.1);

        let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
        assert!(
            result.residual_norm < 1e-10,
            "Residual: {}",
            result.residual_norm
        );
    }
}
