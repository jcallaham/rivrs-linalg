//! Integration tests for discrete-time Sylvester equation solver.
//!
//! Tests `solve_discrete` against analytically known solutions and
//! verifies residual norms for various problem sizes and structures.

use csrrs::sylvester::{compute_residual, solve_discrete, EquationType};
use faer::prelude::*;

/// Helper to verify a discrete solution: ||AXB + X - C|| < tol
fn verify_discrete(a: &Mat<f64>, b: &Mat<f64>, c: &Mat<f64>, tol: f64) {
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x = &result.solution * (1.0 / result.scale);
    let residual = compute_residual(
        a.as_ref(),
        b.as_ref(),
        c.as_ref(),
        x.as_ref(),
        EquationType::Discrete,
    );
    assert!(
        residual < tol,
        "Residual {:.2e} exceeds tolerance {:.2e}",
        residual,
        tol,
    );
}

#[test]
fn test_1x1_analytical() {
    // 0.5*x*0.6 + x = 7 => x*(1 + 0.3) = 7 => x = 7/1.3
    let a = mat![[0.5f64]];
    let b = mat![[0.6f64]];
    let c = mat![[7.0f64]];
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x = &result.solution * (1.0 / result.scale);
    assert!((x[(0, 0)] - 7.0 / 1.3).abs() < 1e-12);
}

#[test]
fn test_2x2_diagonal_analytical() {
    // A = diag(0.5, 0.8), B = diag(0.6, 0.9)
    // x_ij = c_ij / (1 + a_ii * b_jj)
    let a = mat![[0.5, 0.0], [0.0, 0.8f64]];
    let b = mat![[0.6, 0.0], [0.0, 0.9f64]];
    let c = mat![[1.3, 1.45], [1.48, 1.72f64]];
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x = &result.solution * (1.0 / result.scale);
    // x11 = 1.3/(1+0.3) = 1.0
    // x12 = 1.45/(1+0.45) = 1.0
    // x21 = 1.48/(1+0.48) = 1.0
    // x22 = 1.72/(1+0.72) = 1.0
    assert!((x[(0, 0)] - 1.0).abs() < 1e-12);
    assert!((x[(0, 1)] - 1.0).abs() < 1e-12);
    assert!((x[(1, 0)] - 1.0).abs() < 1e-12);
    assert!((x[(1, 1)] - 1.0).abs() < 1e-12);
}

#[test]
fn test_3x3_residual() {
    let a = mat![[0.5, 0.1, 0.0], [0.0, 0.8, 0.2], [0.0, 0.0, 0.3f64]];
    let b = mat![[0.6, 0.1, 0.0], [0.0, 0.9, 0.05], [0.0, 0.0, 0.4f64]];
    let c = Mat::from_fn(3, 3, |i, j| (i * 3 + j + 1) as f64);
    verify_discrete(&a, &b, &c, 1e-10);
}

#[test]
fn test_4x4_complex_eigenvalues() {
    let a = mat![
        [0.5, 0.3, 0.0, 0.0],
        [-0.2, 0.5, 0.0, 0.0],
        [0.0, 0.0, 0.7, 0.1],
        [0.0, 0.0, 0.0, 0.4f64]
    ];
    let b = mat![[0.6, 0.1, 0.0], [0.0, 0.8, 0.2], [0.0, 0.0, 0.5f64]];
    let c = Mat::from_fn(4, 3, |i, j| (i + j + 1) as f64);
    verify_discrete(&a, &b, &c, 1e-10);
}

#[test]
fn test_rectangular_3x2() {
    let a = mat![[0.5, 0.1, 0.0], [0.0, 0.8, 0.2], [0.0, 0.0, 0.3f64]];
    let b = mat![[0.6, 0.1], [0.0, 0.9f64]];
    let c = Mat::from_fn(3, 2, |i, j| (i + j + 1) as f64);
    verify_discrete(&a, &b, &c, 1e-10);
}

#[test]
fn test_5x5_dense() {
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
    verify_discrete(&a, &b, &c, 1e-10);
}

#[test]
fn test_zero_rhs() {
    let a = mat![[0.5, 0.0], [0.0, 0.8f64]];
    let b = mat![[0.6, 0.0], [0.0, 0.9f64]];
    let c = Mat::<f64>::zeros(2, 2);
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x = &result.solution * (1.0 / result.scale);
    for j in 0..2 {
        for i in 0..2 {
            assert!(
                x[(i, j)].abs() < 1e-14,
                "X[{},{}] = {} should be zero",
                i,
                j,
                x[(i, j)]
            );
        }
    }
}

#[test]
fn test_small_eigenvalues() {
    // Matrices with small eigenvalues (stable discrete systems)
    let a = mat![[0.1, 0.05], [0.0, 0.2f64]];
    let b = mat![[0.15, 0.0], [0.0, 0.25f64]];
    let c = mat![[1.0, 2.0], [3.0, 4.0f64]];
    verify_discrete(&a, &b, &c, 1e-10);
}
