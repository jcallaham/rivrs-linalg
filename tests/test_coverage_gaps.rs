//! Tests addressing coverage gaps identified in solver analysis.
//!
//! Covers:
//! - Large matrices (>64) to exercise the blocked solver path
//! - solve_continuous_schur / solve_discrete_schur APIs
//! - Explicit quasi-triangular verification through the full pipeline
//! - Moderate-sized matrices (20x20) to stress iteration logic

use csrrs::sylvester::{
    compute_residual, solve_continuous, solve_continuous_schur, solve_discrete,
    solve_discrete_schur, EquationType,
};
use faer::prelude::*;

/// Build a matrix with known eigenvalue structure that guarantees
/// 2x2 quasi-triangular blocks in its real Schur form.
/// Uses 2x2 rotation blocks on the diagonal to produce complex eigenvalue pairs.
fn build_quasi_triangular_test_matrix(n: usize) -> Mat<f64> {
    let mut m = Mat::zeros(n, n);
    let mut i = 0;
    let mut block_idx = 0;
    while i < n {
        if i + 1 < n && block_idx % 3 == 1 {
            // 2x2 rotation block: eigenvalues = center ± freq*i
            let center = (block_idx + 2) as f64;
            let freq = 0.5;
            m[(i, i)] = center;
            m[(i, i + 1)] = freq;
            m[(i + 1, i)] = -freq;
            m[(i + 1, i + 1)] = center;
            i += 2;
        } else {
            // 1x1 block with real eigenvalue
            m[(i, i)] = (block_idx + 1) as f64;
            i += 1;
        }
        block_idx += 1;
    }
    // Add some upper triangular fill to make it non-trivial
    for j in 0..n {
        for ii in 0..j.min(3) {
            m[(ii, j)] += 0.1 * ((ii + j) as f64);
        }
    }
    m
}

/// Compute Schur decomposition using nalgebra (same method as the solver internally).
fn compute_schur(a: MatRef<'_, f64>) -> (Mat<f64>, Mat<f64>) {
    let n = a.nrows();
    let na_matrix = nalgebra::DMatrix::from_fn(n, n, |i, j| a[(i, j)]);
    let schur = nalgebra::Schur::new(na_matrix);
    let (na_u, na_t) = schur.unpack();
    let mut t = Mat::zeros(n, n);
    let mut u = Mat::zeros(n, n);
    for j in 0..n {
        for i in 0..n {
            t[(i, j)] = na_t[(i, j)];
            u[(i, j)] = na_u[(i, j)];
        }
    }
    (t, u)
}

// === Large Matrix Tests (Blocked Solver Path) ===

#[test]
fn test_continuous_large_blocked_path() {
    // Matrices of size 70 > BLOCKED_THRESHOLD (64) to exercise the blocked solver.
    let n = 70;
    let a = build_quasi_triangular_test_matrix(n);
    let b = build_quasi_triangular_test_matrix(n);
    // Shift B eigenvalues to ensure no common eigenvalues with A
    let b = &b + Mat::<f64>::identity(n, n) * 100.0;
    let c = Mat::from_fn(n, n, |i, j| ((i * n + j) % 17) as f64 - 8.0);

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
        residual < 1e-8,
        "Large blocked continuous residual {:.2e} too high",
        residual
    );
}

#[test]
fn test_continuous_large_rectangular_blocked() {
    // A is 80x80, B is 70x70, C is 80x70. Both exceed the blocked threshold.
    let n = 80;
    let m = 70;
    let a = build_quasi_triangular_test_matrix(n);
    let b = build_quasi_triangular_test_matrix(m);
    let b = &b + Mat::<f64>::identity(m, m) * 100.0;
    let c = Mat::from_fn(n, m, |i, j| ((i + j) % 13) as f64 - 6.0);

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
        residual < 1e-8,
        "Large rectangular blocked continuous residual {:.2e} too high",
        residual
    );
}

#[test]
fn test_discrete_large_matrix() {
    // Large discrete test (note: discrete solver does NOT have a blocked path,
    // but this still tests iteration logic at scale).
    let n = 30;
    // Use matrices with eigenvalues < 1 in magnitude for discrete stability.
    let mut a = Mat::zeros(n, n);
    for i in 0..n {
        a[(i, i)] = 0.3 + 0.01 * (i as f64);
        if i + 1 < n {
            a[(i, i + 1)] = 0.05;
        }
    }
    let mut b = Mat::zeros(n, n);
    for i in 0..n {
        b[(i, i)] = 0.2 + 0.015 * (i as f64);
        if i + 1 < n {
            b[(i, i + 1)] = 0.03;
        }
    }
    let c = Mat::from_fn(n, n, |i, j| ((i * n + j) % 11) as f64 - 5.0);

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
        residual < 1e-8,
        "Large discrete residual {:.2e} too high",
        residual
    );
}

// === solve_*_schur API Tests ===

#[test]
fn test_solve_continuous_schur_api() {
    // Pre-compute Schur decompositions and pass them to the Schur API.
    let a = mat![[1.0, 2.0, 0.0], [-1.0, 1.0, 0.0], [0.0, 0.0, 3.0f64]];
    let b = mat![[4.0, 0.5], [0.0, 5.0f64]];
    let c = mat![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0f64]];

    let (schur_a, u) = compute_schur(a.as_ref());
    let (schur_b, v) = compute_schur(b.as_ref());

    let result = solve_continuous_schur(
        schur_a.as_ref(),
        schur_b.as_ref(),
        u.as_ref(),
        v.as_ref(),
        c.as_ref(),
    )
    .unwrap();

    // Compare with the standard API
    let result_std = solve_continuous(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x_schur = &result.solution * (1.0 / result.scale);
    let x_std = &result_std.solution * (1.0 / result_std.scale);

    let mut max_diff = 0.0f64;
    for j in 0..2 {
        for i in 0..3 {
            max_diff = max_diff.max((x_schur[(i, j)] - x_std[(i, j)]).abs());
        }
    }
    assert!(
        max_diff < 1e-10,
        "Schur API and standard API differ by {:.2e}",
        max_diff
    );
    assert!(result.residual_norm < 1e-10);
}

#[test]
fn test_solve_discrete_schur_api() {
    let a = mat![[0.5, 0.1], [0.0, 0.8f64]];
    let b = mat![[0.6, 0.1], [0.0, 0.9f64]];
    let c = mat![[1.0, 2.0], [3.0, 4.0f64]];

    let (schur_a, u) = compute_schur(a.as_ref());
    let (schur_b, v) = compute_schur(b.as_ref());

    let result = solve_discrete_schur(
        schur_a.as_ref(),
        schur_b.as_ref(),
        u.as_ref(),
        v.as_ref(),
        c.as_ref(),
    )
    .unwrap();

    let result_std = solve_discrete(a.as_ref(), b.as_ref(), c.as_ref()).unwrap();
    let x_schur = &result.solution * (1.0 / result.scale);
    let x_std = &result_std.solution * (1.0 / result_std.scale);

    let mut max_diff = 0.0f64;
    for j in 0..2 {
        for i in 0..2 {
            max_diff = max_diff.max((x_schur[(i, j)] - x_std[(i, j)]).abs());
        }
    }
    assert!(
        max_diff < 1e-10,
        "Discrete Schur API and standard API differ by {:.2e}",
        max_diff
    );
    assert!(result.residual_norm < 1e-10);
}

// === Explicit Quasi-Triangular Verification ===

#[test]
fn test_continuous_complex_eigenvalues_verified() {
    // Build a matrix that MUST produce 2x2 Schur blocks (complex eigenvalues)
    // and verify the solver handles them correctly through the full pipeline.
    //
    // A has eigenvalues 1 ± 2i (from 2x2 rotation block) and 5 (real)
    let a = mat![[1.0, 2.0, 0.3], [-2.0, 1.0, 0.1], [0.0, 0.0, 5.0f64]];
    // B has eigenvalues 3 ± i and 7
    let b = mat![[3.0, 1.0, 0.2], [-1.0, 3.0, 0.1], [0.0, 0.0, 7.0f64]];
    let c = Mat::from_fn(3, 3, |i, j| ((i + 1) * (j + 2)) as f64);

    // Verify Schur form actually has 2x2 blocks
    let (schur_a, _) = compute_schur(a.as_ref());
    let has_2x2_block = (0..schur_a.nrows() - 1).any(|i| schur_a[(i + 1, i)].abs() > 1e-10);
    assert!(
        has_2x2_block,
        "Test matrix A should produce 2x2 Schur blocks"
    );

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
        residual < 1e-10,
        "Residual with verified 2x2 Schur blocks: {:.2e}",
        residual
    );
}

#[test]
fn test_discrete_complex_eigenvalues_verified() {
    // A has complex eigenvalues (small magnitude for discrete stability)
    let a = mat![[0.3, 0.4, 0.0], [-0.4, 0.3, 0.0], [0.0, 0.0, 0.5f64]];
    // B has complex eigenvalues
    let b = mat![[0.2, 0.3], [-0.3, 0.2f64]];
    let c = mat![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0f64]];

    // Verify Schur form has 2x2 blocks
    let (schur_a, _) = compute_schur(a.as_ref());
    let has_2x2_block = (0..schur_a.nrows() - 1).any(|i| schur_a[(i + 1, i)].abs() > 1e-10);
    assert!(
        has_2x2_block,
        "Test matrix A should produce 2x2 Schur blocks"
    );

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
        residual < 1e-10,
        "Discrete with verified 2x2 Schur blocks: {:.2e}",
        residual
    );
}

// === Moderate-Size Stress Tests ===

#[test]
fn test_continuous_20x20() {
    let n = 20;
    let a = build_quasi_triangular_test_matrix(n);
    let b = build_quasi_triangular_test_matrix(n);
    let b = &b + Mat::<f64>::identity(n, n) * 50.0;
    let c = Mat::from_fn(n, n, |i, j| (i as f64 - j as f64).sin());

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
        residual < 1e-9,
        "20x20 continuous residual {:.2e} too high",
        residual
    );
}

#[test]
fn test_continuous_50x50() {
    let n = 50;
    let a = build_quasi_triangular_test_matrix(n);
    let b = build_quasi_triangular_test_matrix(n);
    let b = &b + Mat::<f64>::identity(n, n) * 80.0;
    let c = Mat::from_fn(n, n, |i, j| ((i * n + j) % 7) as f64 - 3.0);

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
        residual < 1e-8,
        "50x50 continuous residual {:.2e} too high",
        residual
    );
}
