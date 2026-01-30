//! Integration tests for continuous-time Sylvester equation solver.
//!
//! Tests `solve_continuous` against analytically known solutions and
//! verifies residual norms for various problem sizes and structures.

use csrrs::sylvester::{compute_residual, solve_continuous, EquationType};
use faer::prelude::*;

/// Helper to verify a continuous solution: ||AX + XB - C|| < tol
fn verify_continuous(a: &Mat<f64>, b: &Mat<f64>, c: &Mat<f64>, tol: f64) {
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x = &result.solution * (1.0 / result.scale);
    let residual = compute_residual(
        a.as_ref(),
        b.as_ref(),
        c.as_ref(),
        x.as_ref(),
        EquationType::Continuous,
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
    // 3x + 5x = 16 => x = 2
    let a = mat![[3.0f64]];
    let b = mat![[5.0f64]];
    let c = mat![[16.0f64]];
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x = &result.solution * (1.0 / result.scale);
    assert!((x[(0, 0)] - 2.0).abs() < 1e-12);
}

#[test]
fn test_2x2_diagonal_analytical() {
    // A = diag(1, -1), B = diag(2, 3)
    // x_ij = c_ij / (a_ii + b_jj)
    let a = mat![[1.0, 0.0], [0.0, -1.0f64]];
    let b = mat![[2.0, 0.0], [0.0, 3.0f64]];
    let c = mat![[6.0, 8.0], [2.0, 4.0f64]];
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x = &result.solution * (1.0 / result.scale);
    assert!((x[(0, 0)] - 2.0).abs() < 1e-12); // 6/(1+2)
    assert!((x[(0, 1)] - 2.0).abs() < 1e-12); // 8/(1+3)
    assert!((x[(1, 0)] - 2.0).abs() < 1e-12); // 2/(-1+2)
    assert!((x[(1, 1)] - 2.0).abs() < 1e-12); // 4/(-1+3)
}

#[test]
fn test_3x3_residual() {
    let a = mat![[2.0, 1.0, 0.5], [-1.0, 3.0, 0.0], [0.0, 0.0, 4.0f64]];
    let b = mat![[5.0, 0.5, 0.0], [0.0, 6.0, 1.0], [0.0, 0.0, 7.0f64]];
    let c = Mat::from_fn(3, 3, |i, j| (i * 3 + j + 1) as f64);
    verify_continuous(&a, &b, &c, 1e-10);
}

#[test]
fn test_4x4_complex_eigenvalues() {
    // Matrix with complex eigenvalue pairs
    let a = mat![
        [1.0, 2.0, 0.0, 0.0],
        [-1.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 3.0, 1.0],
        [0.0, 0.0, -0.5, 3.0f64]
    ];
    let b = mat![
        [5.0, 1.0, 0.0, 0.0],
        [0.0, 6.0, 0.5, 0.0],
        [0.0, 0.0, 7.0, 2.0],
        [0.0, 0.0, -1.0, 7.0f64]
    ];
    let c = Mat::from_fn(4, 4, |i, j| ((i + 1) * (j + 1)) as f64);
    verify_continuous(&a, &b, &c, 1e-10);
}

#[test]
fn test_rectangular_3x2() {
    let a = mat![[1.0, 0.5, 0.0], [0.0, 2.0, 0.3], [0.0, 0.0, 3.0f64]];
    let b = mat![[4.0, 0.5], [0.0, 5.0f64]];
    let c = Mat::from_fn(3, 2, |i, j| (i + j + 1) as f64);
    verify_continuous(&a, &b, &c, 1e-10);
}

#[test]
fn test_5x5_dense() {
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
    verify_continuous(&a, &b, &c, 1e-10);
}

#[test]
fn test_identity_rhs() {
    let a = mat![[1.0, 0.0], [0.0, 2.0f64]];
    let b = mat![[3.0, 0.0], [0.0, 4.0f64]];
    let c = Mat::<f64>::identity(2, 2);
    verify_continuous(&a, &b, &c, 1e-10);
}

#[test]
fn test_zero_rhs() {
    let a = mat![[1.0, 0.0], [0.0, 2.0f64]];
    let b = mat![[3.0, 0.0], [0.0, 4.0f64]];
    let c = Mat::<f64>::zeros(2, 2);
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
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
fn test_negative_eigenvalues() {
    // Both A and B have negative eigenvalues (stable systems)
    let a = mat![[-1.0, 0.5], [0.0, -2.0f64]];
    let b = mat![[-3.0, 0.0], [0.0, -4.0f64]];
    let c = mat![[1.0, 2.0], [3.0, 4.0f64]];
    verify_continuous(&a, &b, &c, 1e-10);
}

// === SLICOT Benchmark Tests ===

#[test]
fn test_slicot_sb04md_benchmark() {
    // Test case from SLICOT SB04MD documentation (Release 4.5, NICONET 2002-2005).
    // Solves AX + XB = C (continuous-time, Hessenberg-Schur method).
    // N=3, M=2 (rectangular C).
    //
    // References:
    // - Golub, Nash & Van Loan (1979), IEEE Trans. Auto. Contr., AC-24:909-913
    // - Sima (1996), "Algorithms for Linear-Quadratic Optimization"
    let a = mat![
        [2.0, 1.0, 3.0],
        [0.0, 2.0, 1.0],
        [6.0, 1.0, 2.0f64]
    ];
    let b = mat![
        [2.0, 1.0],
        [1.0, 6.0f64]
    ];
    let c = mat![
        [2.0, 1.0],
        [1.0, 4.0],
        [0.0, 5.0f64]
    ];

    // Known solution from SLICOT SB04MD documentation
    let x_expected = mat![
        [-2.7685, 0.5498],
        [-1.0531, 0.6865],
        [ 4.5257, -0.4389f64]
    ];

    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x = &result.solution * (1.0 / result.scale);

    // Verify solution matches SLICOT reference (documented to 4 decimal places)
    let mut max_err = 0.0f64;
    for j in 0..2 {
        for i in 0..3 {
            max_err = max_err.max((x[(i, j)] - x_expected[(i, j)]).abs());
        }
    }
    assert!(
        max_err < 1e-4,
        "Solution differs from SLICOT SB04MD reference by {:.2e}",
        max_err
    );

    // Verify residual to high precision
    let residual = compute_residual(
        a.as_ref(),
        b.as_ref(),
        c.as_ref(),
        x.as_ref(),
        EquationType::Continuous,
    );
    assert!(
        residual < 1e-11,
        "Residual {:.2e} exceeds tolerance 1e-11",
        residual
    );
}
