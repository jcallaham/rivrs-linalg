//! Integration tests for error handling in Sylvester equation solvers.
//!
//! Verifies that appropriate errors are returned for invalid inputs,
//! singular cases, and other exceptional conditions.

use csrrs::error::SylvesterError;
use csrrs::sylvester::{solve_continuous, solve_discrete};
use faer::prelude::*;

// === Dimension Mismatch Tests ===

#[test]
fn test_continuous_c_wrong_dims() {
    let a = Mat::<f64>::zeros(2, 2);
    let b = Mat::zeros(3, 3);
    let c = Mat::zeros(2, 2); // Should be 2x3
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::DimensionMismatch { .. })));
}

#[test]
fn test_discrete_c_wrong_dims() {
    let a = Mat::<f64>::zeros(3, 3);
    let b = Mat::zeros(2, 2);
    let c = Mat::zeros(3, 3); // Should be 3x2
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::DimensionMismatch { .. })));
}

#[test]
fn test_continuous_a_not_square() {
    let a = Mat::<f64>::zeros(2, 3);
    let b = Mat::zeros(3, 3);
    let c = Mat::zeros(2, 3);
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::NotSquare { matrix: 'A', .. })));
}

#[test]
fn test_discrete_b_not_square() {
    let a = Mat::<f64>::zeros(2, 2);
    let b = Mat::zeros(3, 2);
    let c = Mat::zeros(2, 3);
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::NotSquare { matrix: 'B', .. })));
}

// === Invalid Input Tests ===

#[test]
fn test_continuous_nan_in_a() {
    let a = mat![[f64::NAN, 0.0], [0.0, 1.0]];
    let b = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let c = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::InvalidInput { .. })));
}

#[test]
fn test_continuous_inf_in_b() {
    let a = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let b = mat![[f64::INFINITY, 0.0], [0.0, 1.0]];
    let c = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::InvalidInput { .. })));
}

#[test]
fn test_continuous_nan_in_c() {
    let a = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let b = mat![[2.0, 0.0], [0.0, 3.0f64]];
    let c = mat![[1.0, f64::NAN], [0.0, 1.0]];
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::InvalidInput { .. })));
}

#[test]
fn test_discrete_nan_in_a() {
    let a = mat![[f64::NAN, 0.0], [0.0, 1.0]];
    let b = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let c = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::InvalidInput { .. })));
}

#[test]
fn test_discrete_neg_inf_in_b() {
    let a = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let b = mat![[f64::NEG_INFINITY, 0.0], [0.0, 1.0]];
    let c = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::InvalidInput { .. })));
}

// === Common Eigenvalues Tests ===

#[test]
fn test_continuous_common_eigenvalues_1x1() {
    // λ(A) + λ(B) = 3 + (-3) = 0
    let a = mat![[3.0f64]];
    let b = mat![[-3.0f64]];
    let c = mat![[1.0f64]];
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::CommonEigenvalues { .. })));
}

#[test]
fn test_continuous_common_eigenvalues_2x2() {
    // A has eigenvalue 1, B has eigenvalue -1
    // Make 2x2 matrices with these eigenvalues
    let a = mat![[1.0, 0.0], [0.0, 5.0f64]]; // eigenvalues: 1, 5
    let b = mat![[-1.0, 0.0], [0.0, 3.0f64]]; // eigenvalues: -1, 3
    let c = mat![[1.0, 0.0], [0.0, 1.0f64]];
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::CommonEigenvalues { .. })));
}

#[test]
fn test_discrete_common_eigenvalues_1x1() {
    // 1 + λ(A)*λ(B) = 1 + 2*(-0.5) = 0
    let a = mat![[2.0f64]];
    let b = mat![[-0.5f64]];
    let c = mat![[1.0f64]];
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref());
    assert!(matches!(result, Err(SylvesterError::CommonEigenvalues { .. })));
}

// === Edge Cases ===

#[test]
fn test_continuous_empty() {
    let a = Mat::<f64>::zeros(0, 0);
    let b = Mat::zeros(0, 0);
    let c = Mat::zeros(0, 0);
    let result = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    assert_eq!(result.solution.nrows(), 0);
    assert_eq!(result.solution.ncols(), 0);
}

#[test]
fn test_discrete_empty() {
    let a = Mat::<f64>::zeros(0, 0);
    let b = Mat::zeros(0, 0);
    let c = Mat::zeros(0, 0);
    let result = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    assert_eq!(result.solution.nrows(), 0);
    assert_eq!(result.solution.ncols(), 0);
}

// === Error Message Quality ===

#[test]
fn test_common_eigenvalues_error_message() {
    let a = mat![[2.0f64]];
    let b = mat![[-2.0f64]];
    let c = mat![[1.0f64]];
    let err = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap_err();
    let msg = format!("{}", err);
    assert!(msg.contains("common eigenvalues"), "Error message should mention common eigenvalues: {}", msg);
    assert!(msg.contains("ill-conditioned"), "Error message should mention ill-conditioning: {}", msg);
    assert!(msg.contains("separation"), "Error message should mention separation: {}", msg);
}

#[test]
fn test_dimension_error_message() {
    let a = Mat::<f64>::zeros(2, 2);
    let b = Mat::zeros(3, 3);
    let c = Mat::zeros(2, 2); // Wrong
    let err = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap_err();
    let msg = format!("{}", err);
    assert!(msg.contains("Dimension mismatch"), "Error message: {}", msg);
}
