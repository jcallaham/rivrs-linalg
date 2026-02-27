//! Adversarial matrix perturbation helpers for torture testing.
//!
//! These helpers transform well-formed dense symmetric matrices into adversarial
//! variants that exercise the full pivot decision tree of the APTP kernel.
//!
//! Inspired by the perturbation helpers in SPRAL's `ldlt_app.cxx` (BSD-3),
//! implemented from scratch in Rust (clean room).
//!
//! # Reference
//!
//! - SPRAL `tests/ssids/kernels/ldlt_app.cxx`, lines 78–134 (BSD-3-Clause)
//! - Duff, Hogg & Lopez (2020), Section 4: threshold parameter and pivot failures

use faer::Mat;
use rand::Rng;

/// Configuration for probabilistic torture test generation.
///
/// Controls the perturbation mix applied to randomly generated matrices.
/// Default values match SPRAL's torture test configuration.
#[derive(Debug, Clone)]
pub struct TortureTestConfig {
    /// Number of random instances per configuration.
    pub num_instances: usize,
    /// Probability of applying `cause_delays` (default 0.70).
    pub delay_probability: f64,
    /// Probability of applying `make_singular` (default 0.20).
    pub singular_probability: f64,
    /// Probability of applying `make_dblk_singular` (default 0.10).
    pub dblk_singular_probability: f64,
    /// (m, n) matrix size pairs to test.
    pub matrix_sizes: Vec<(usize, usize)>,
    /// Maximum backward error for non-singular instances.
    pub backward_error_threshold: f64,
    /// Fixed seed for reproducibility.
    pub seed: u64,
}

impl Default for TortureTestConfig {
    fn default() -> Self {
        Self {
            num_instances: 500,
            delay_probability: 0.70,
            singular_probability: 0.20,
            dblk_singular_probability: 0.10,
            matrix_sizes: vec![(32, 32), (64, 64), (128, 128), (128, 48)],
            backward_error_threshold: 5e-11,
            seed: 12345,
        }
    }
}

/// Apply delay-inducing perturbations to a dense symmetric matrix.
///
/// Multiplies n/8 random rows (and their symmetric column counterparts) by 1000,
/// and n/8 random individual entries (symmetrized) by 1000. When `n > block_size`,
/// ensures the first oversized row is past the first block boundary.
///
/// This forces pivot threshold failures in the APTP kernel, triggering delayed columns.
///
/// # Arguments
///
/// * `matrix` - Mutable dense symmetric matrix (modified in place)
/// * `block_size` - The inner block size of the APTP kernel (typically 32)
/// * `rng` - Random number generator
///
/// # Postcondition
///
/// Matrix remains symmetric. At least n/8 rows have been scaled.
pub fn cause_delays(matrix: &mut Mat<f64>, block_size: usize, rng: &mut impl Rng) {
    let n = matrix.nrows();
    if n < 8 {
        return;
    }

    let num_rows = n / 8;
    let num_entries = n / 8;

    // Scale n/8 random rows by 1000 (and symmetric counterpart columns).
    // When n > block_size, ensure first selected row index > block_size.
    let mut scaled_rows = Vec::with_capacity(num_rows);
    for i in 0..num_rows {
        let row = if i == 0 && n > block_size {
            // First oversized row must be past first block boundary
            block_size + rng.gen_range(0..n - block_size)
        } else {
            rng.gen_range(0..n)
        };
        scaled_rows.push(row);

        // Scale row and symmetric column
        for j in 0..n {
            matrix[(row, j)] *= 1000.0;
            matrix[(j, row)] *= 1000.0;
        }
        // Diagonal was scaled twice — correct
        matrix[(row, row)] /= 1000.0;
    }

    // Scale n/8 random individual entries by 1000 (symmetrized)
    for _ in 0..num_entries {
        let i = rng.gen_range(0..n);
        let j = rng.gen_range(0..n);
        matrix[(i, j)] *= 1000.0;
        matrix[(j, i)] *= 1000.0;
        if i == j {
            // Diagonal was scaled twice — correct
            matrix[(i, j)] /= 1000.0;
        }
    }
}

/// Make a matrix rank-deficient by making col2 a scaled copy of col1.
///
/// Reads column `col1`, scales it, writes to column `col2`, then symmetrizes
/// by copying the new column to the corresponding row.
///
/// # Arguments
///
/// * `matrix` - Mutable dense symmetric matrix (modified in place)
/// * `col1` - Source column index
/// * `col2` - Target column index (will become linearly dependent on col1)
///
/// # Postcondition
///
/// Matrix is rank-deficient (rank decreased by 1) and remains symmetric.
///
/// # Panics
///
/// Panics if `col1` or `col2` are out of bounds, or if `col1 == col2`.
pub fn make_singular(matrix: &mut Mat<f64>, col1: usize, col2: usize) {
    assert_ne!(col1, col2, "col1 and col2 must be different");
    let n = matrix.nrows();
    assert!(col1 < n && col2 < n, "column indices out of bounds");

    // Make col2 = scale * col1 in the resulting symmetric matrix.
    //
    // After modification, the matrix must satisfy:
    //   m_new[(i, col2)] = scale * m_new[(i, col1)]  for all i
    //   m_new[(col2, i)] = scale * m_new[(col1, i)]  for all i  (by symmetry)
    //
    // The tricky part: setting row col2 also modifies col1 at row col2.
    // So m_new[(col2, col1)] = scale * m_orig[(col1, col1)].
    // Then m_new[(col2, col2)] must be scale * m_new[(col2, col1)]
    //   = scale * scale * m_orig[(col1, col1)] = scale² * m_orig[(col1, col1)].
    let scale = 2.0;
    let diag_col1 = matrix[(col1, col1)];

    // Read col1 before any modifications
    let col1_vals: Vec<f64> = (0..n).map(|i| matrix[(i, col1)]).collect();

    // Write scaled col1 to col2 column and row col2 (symmetrize)
    for i in 0..n {
        if i != col2 {
            matrix[(i, col2)] = scale * col1_vals[i];
            matrix[(col2, i)] = scale * col1_vals[i];
        }
    }

    // Set cross entry: m[(col2, col1)] = m[(col1, col2)] = scale * diag_col1
    matrix[(col2, col1)] = scale * diag_col1;
    matrix[(col1, col2)] = scale * diag_col1;

    // Set diagonal: m[(col2, col2)] = scale² * diag_col1
    matrix[(col2, col2)] = scale * scale * diag_col1;
}

/// Make a specific diagonal block singular.
///
/// Selects the first and last columns of the specified diagonal block and
/// delegates to [`make_singular`] to create a rank deficiency within the block.
///
/// # Arguments
///
/// * `matrix` - Mutable dense symmetric matrix (modified in place)
/// * `block_row` - Starting row/column index of the block
/// * `block_size` - Size of the diagonal block
///
/// # Postcondition
///
/// The targeted diagonal block is singular. Matrix remains symmetric.
///
/// # Panics
///
/// Panics if the block extends beyond matrix dimensions or if `block_size < 2`.
pub fn make_dblk_singular(matrix: &mut Mat<f64>, block_row: usize, block_size: usize) {
    assert!(block_size >= 2, "block_size must be >= 2");
    let n = matrix.nrows();
    assert!(
        block_row + block_size <= n,
        "block extends beyond matrix dimensions"
    );

    let col1 = block_row;
    let col2 = block_row + block_size - 1;
    make_singular(matrix, col1, col2);
}

/// Generate a random dense symmetric indefinite matrix.
///
/// Creates an n×n symmetric matrix with entries uniformly distributed in [-1, 1],
/// then adds diagonal dominance with mixed signs (half positive, half negative).
///
/// This is the base matrix generator used by torture tests before perturbation.
pub fn generate_dense_symmetric_indefinite(n: usize, rng: &mut impl Rng) -> Mat<f64> {
    let mut a = Mat::zeros(n, n);

    // Fill lower triangle with random values, then symmetrize
    for j in 0..n {
        for i in j..n {
            let v: f64 = rng.gen_range(-1.0..1.0);
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }

    // Add diagonal dominance with mixed signs
    for i in 0..n {
        let row_sum: f64 = (0..n).filter(|&j| j != i).map(|j| a[(i, j)].abs()).sum();
        let margin = 1.0 + rng.r#gen::<f64>();
        if i < n / 2 {
            a[(i, i)] = row_sum + margin;
        } else {
            a[(i, i)] = -(row_sum + margin);
        }
    }

    a
}

/// Generate a random dense symmetric positive definite matrix.
///
/// Creates an n×n symmetric matrix with diagonal dominance ensuring all
/// eigenvalues are positive.
pub fn generate_dense_symmetric_pd(n: usize, rng: &mut impl Rng) -> Mat<f64> {
    let mut a = Mat::zeros(n, n);

    // Fill lower triangle with random values, then symmetrize
    for j in 0..n {
        for i in j..n {
            let v: f64 = rng.gen_range(-1.0..1.0);
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }

    // Add positive diagonal dominance
    for i in 0..n {
        let row_sum: f64 = (0..n).filter(|&j| j != i).map(|j| a[(i, j)].abs()).sum();
        let margin = 1.0 + rng.r#gen::<f64>();
        a[(i, i)] = row_sum + margin;
    }

    a
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn seeded_rng() -> StdRng {
        StdRng::seed_from_u64(42)
    }

    fn is_symmetric(m: &Mat<f64>) -> bool {
        let n = m.nrows();
        for i in 0..n {
            for j in 0..n {
                if (m[(i, j)] - m[(j, i)]).abs() > 1e-15 {
                    return false;
                }
            }
        }
        true
    }

    #[test]
    fn cause_delays_preserves_symmetry() {
        let mut rng = seeded_rng();
        let mut m = generate_dense_symmetric_indefinite(64, &mut rng);
        assert!(is_symmetric(&m));
        cause_delays(&mut m, 32, &mut rng);
        assert!(is_symmetric(&m), "cause_delays broke symmetry");
    }

    #[test]
    fn cause_delays_scales_rows() {
        let mut rng = seeded_rng();
        let n = 64;
        let original = generate_dense_symmetric_indefinite(n, &mut rng);
        let mut perturbed = original.clone();
        cause_delays(&mut perturbed, 32, &mut rng);

        // At least some rows should have been modified
        let mut modified_rows = 0;
        for i in 0..n {
            let orig_norm: f64 = (0..n).map(|j| original[(i, j)].powi(2)).sum::<f64>().sqrt();
            let pert_norm: f64 = (0..n)
                .map(|j| perturbed[(i, j)].powi(2))
                .sum::<f64>()
                .sqrt();
            if (pert_norm - orig_norm).abs() > 1e-10 {
                modified_rows += 1;
            }
        }
        assert!(
            modified_rows >= n / 8,
            "expected at least {} modified rows, got {}",
            n / 8,
            modified_rows
        );
    }

    #[test]
    fn cause_delays_noop_for_small_matrix() {
        let mut rng = seeded_rng();
        let mut m = generate_dense_symmetric_indefinite(4, &mut rng);
        let original = m.clone();
        cause_delays(&mut m, 32, &mut rng);
        // n < 8, so no perturbation should be applied
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(m[(i, j)], original[(i, j)]);
            }
        }
    }

    #[test]
    fn make_singular_creates_rank_deficiency() {
        let mut rng = seeded_rng();
        let n = 16;
        let mut m = generate_dense_symmetric_indefinite(n, &mut rng);
        make_singular(&mut m, 0, 1);

        assert!(is_symmetric(&m), "make_singular broke symmetry");

        // col2 should be exactly 2× col1 in the modified matrix
        for i in 0..n {
            assert!(
                (m[(i, 1)] - 2.0 * m[(i, 0)]).abs() < 1e-12,
                "col1 and col2 not linearly dependent at row {}: m[{},1]={}, 2*m[{},0]={}",
                i,
                i,
                m[(i, 1)],
                i,
                2.0 * m[(i, 0)]
            );
        }
    }

    #[test]
    fn make_singular_preserves_symmetry() {
        let mut rng = seeded_rng();
        let mut m = generate_dense_symmetric_indefinite(32, &mut rng);
        make_singular(&mut m, 5, 10);
        assert!(is_symmetric(&m), "make_singular broke symmetry");
    }

    #[test]
    fn make_dblk_singular_targets_block() {
        let mut rng = seeded_rng();
        let n = 64;
        let mut m = generate_dense_symmetric_indefinite(n, &mut rng);
        make_dblk_singular(&mut m, 16, 8);

        assert!(is_symmetric(&m), "make_dblk_singular broke symmetry");

        // First and last columns of block should be linearly dependent
        let col1 = 16;
        let col2 = 23; // 16 + 8 - 1
        for i in 0..n {
            assert!(
                (m[(i, col2)] - 2.0 * m[(i, col1)]).abs() < 1e-12,
                "block columns not linearly dependent at row {}",
                i
            );
        }
    }

    #[test]
    #[should_panic(expected = "col1 and col2 must be different")]
    fn make_singular_same_col_panics() {
        let mut rng = seeded_rng();
        let mut m = generate_dense_symmetric_indefinite(8, &mut rng);
        make_singular(&mut m, 3, 3);
    }

    #[test]
    #[should_panic(expected = "block_size must be >= 2")]
    fn make_dblk_singular_size_1_panics() {
        let mut rng = seeded_rng();
        let mut m = generate_dense_symmetric_indefinite(8, &mut rng);
        make_dblk_singular(&mut m, 0, 1);
    }

    #[test]
    fn generate_dense_symmetric_indefinite_is_symmetric() {
        let mut rng = seeded_rng();
        let m = generate_dense_symmetric_indefinite(50, &mut rng);
        assert!(is_symmetric(&m));
    }

    #[test]
    fn generate_dense_symmetric_indefinite_has_mixed_signs() {
        let mut rng = seeded_rng();
        let m = generate_dense_symmetric_indefinite(50, &mut rng);
        let mut has_pos = false;
        let mut has_neg = false;
        for i in 0..50 {
            if m[(i, i)] > 0.0 {
                has_pos = true;
            }
            if m[(i, i)] < 0.0 {
                has_neg = true;
            }
        }
        assert!(has_pos && has_neg, "should have mixed diagonal signs");
    }

    #[test]
    fn generate_dense_symmetric_pd_is_symmetric() {
        let mut rng = seeded_rng();
        let m = generate_dense_symmetric_pd(50, &mut rng);
        assert!(is_symmetric(&m));
    }

    #[test]
    fn generate_dense_symmetric_pd_positive_diag() {
        let mut rng = seeded_rng();
        let m = generate_dense_symmetric_pd(50, &mut rng);
        for i in 0..50 {
            assert!(m[(i, i)] > 0.0, "diagonal entry {} should be positive", i);
        }
    }
}
