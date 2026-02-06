//! Random and structured sparse symmetric matrix generators.

use faer::sparse::SparseColMat;
use faer::sparse::Triplet;
use rand::Rng;
use rand_distr::{Distribution, Uniform};

use crate::error::SparseError;

/// Configuration for random sparse symmetric matrix generation.
pub struct RandomMatrixConfig {
    pub size: usize,
    pub target_nnz: usize,
    pub positive_definite: bool,
}

/// Add diagonal entries with diagonal dominance or mixed signs.
fn add_diagonal(
    triplets: &mut Vec<Triplet<usize, usize, f64>>,
    row_abs_sum: &[f64],
    positive_definite: bool,
    rng: &mut impl Rng,
) {
    let n = row_abs_sum.len();
    let half = n / 2;
    for (i, &abs_sum) in row_abs_sum.iter().enumerate() {
        let margin = 1.0 + rng.r#gen::<f64>();
        if positive_definite || i < half {
            triplets.push(Triplet::new(i, i, abs_sum + margin));
        } else {
            triplets.push(Triplet::new(i, i, -(abs_sum + margin)));
        }
    }
}

/// Generate a random sparse symmetric matrix.
///
/// If `positive_definite` is true, uses diagonal dominance to guarantee PD.
/// If false, generates an indefinite matrix with mixed diagonal signs.
/// The `target_nnz` is approximate; actual nnz may differ.
pub fn generate_random_symmetric(
    config: &RandomMatrixConfig,
    rng: &mut impl Rng,
) -> Result<SparseColMat<usize, f64>, SparseError> {
    let n = config.size;

    if n == 0 {
        return Err(SparseError::InvalidInput {
            reason: "matrix size must be > 0".to_string(),
        });
    }

    if n == 1 && !config.positive_definite {
        return Err(SparseError::InvalidInput {
            reason: "cannot generate 1x1 indefinite matrix".to_string(),
        });
    }

    // Off-diagonal entries come in pairs (i,j) and (j,i); plus n diagonal entries
    let target_offdiag_pairs = if config.target_nnz > n {
        (config.target_nnz - n) / 2
    } else {
        0
    };
    let max_offdiag_pairs = n * (n - 1) / 2;
    let actual_pairs = target_offdiag_pairs.min(max_offdiag_pairs);

    let val_dist = Uniform::new(-1.0, 1.0);
    let idx_dist = Uniform::new(0, n);

    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
    let mut placed = std::collections::HashSet::new();

    for _ in 0..actual_pairs * 3 {
        if placed.len() >= actual_pairs {
            break;
        }
        let i = idx_dist.sample(rng);
        let j = idx_dist.sample(rng);
        if i == j {
            continue;
        }
        let (lo, hi) = if i > j { (j, i) } else { (i, j) };
        if placed.contains(&(lo, hi)) {
            continue;
        }
        placed.insert((lo, hi));
        let v = val_dist.sample(rng);
        triplets.push(Triplet::new(lo, hi, v));
        triplets.push(Triplet::new(hi, lo, v));
    }

    // Compute row sums of absolute off-diagonal values
    let mut row_abs_sum = vec![0.0f64; n];
    for t in &triplets {
        if t.row != t.col {
            row_abs_sum[t.row] += t.val.abs();
        }
    }

    add_diagonal(&mut triplets, &row_abs_sum, config.positive_definite, rng);

    SparseColMat::try_new_from_triplets(n, n, &triplets).map_err(|e| SparseError::InvalidInput {
        reason: format!("failed to create sparse matrix from triplets: {:?}", e),
    })
}

/// Generate a sparse symmetric arrow matrix of given size.
///
/// Arrow pattern: dense first row/column, diagonal remainder.
/// If `positive_definite` is true, diagonal dominance is applied.
pub fn generate_arrow(
    size: usize,
    positive_definite: bool,
    rng: &mut impl Rng,
) -> SparseColMat<usize, f64> {
    let n = size;
    let val_dist = Uniform::new(0.1f64, 1.0);
    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
    let mut row_abs_sum = vec![0.0f64; n];

    for j in 1..n {
        let v: f64 = val_dist.sample(rng);
        triplets.push(Triplet::new(0, j, v));
        triplets.push(Triplet::new(j, 0, v));
        row_abs_sum[0] += v.abs();
        row_abs_sum[j] += v.abs();
    }

    add_diagonal(&mut triplets, &row_abs_sum, positive_definite, rng);

    SparseColMat::try_new_from_triplets(n, n, &triplets)
        .expect("arrow matrix construction should not fail")
}

/// Generate a sparse symmetric tridiagonal matrix of given size.
///
/// If `positive_definite` is true, diagonal dominance is applied.
pub fn generate_tridiagonal(
    size: usize,
    positive_definite: bool,
    rng: &mut impl Rng,
) -> SparseColMat<usize, f64> {
    let n = size;
    let val_dist = Uniform::new(0.1f64, 1.0);
    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
    let mut row_abs_sum = vec![0.0f64; n];

    for i in 0..n - 1 {
        let v: f64 = val_dist.sample(rng);
        triplets.push(Triplet::new(i, i + 1, v));
        triplets.push(Triplet::new(i + 1, i, v));
        row_abs_sum[i] += v.abs();
        row_abs_sum[i + 1] += v.abs();
    }

    add_diagonal(&mut triplets, &row_abs_sum, positive_definite, rng);

    SparseColMat::try_new_from_triplets(n, n, &triplets)
        .expect("tridiagonal matrix construction should not fail")
}

/// Generate a sparse symmetric banded matrix of given size and bandwidth.
///
/// If `positive_definite` is true, diagonal dominance is applied.
pub fn generate_banded(
    size: usize,
    bandwidth: usize,
    positive_definite: bool,
    rng: &mut impl Rng,
) -> SparseColMat<usize, f64> {
    let n = size;
    let val_dist = Uniform::new(0.1f64, 1.0);
    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
    let mut row_abs_sum = vec![0.0f64; n];

    for i in 0..n {
        let j_start = i.saturating_sub(bandwidth);
        let j_end = (i + bandwidth + 1).min(n);
        for j in j_start..j_end {
            if j > i {
                let v: f64 = val_dist.sample(rng);
                triplets.push(Triplet::new(i, j, v));
                triplets.push(Triplet::new(j, i, v));
                row_abs_sum[i] += v.abs();
                row_abs_sum[j] += v.abs();
            }
        }
    }

    add_diagonal(&mut triplets, &row_abs_sum, positive_definite, rng);

    SparseColMat::try_new_from_triplets(n, n, &triplets)
        .expect("banded matrix construction should not fail")
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn seeded_rng() -> StdRng {
        StdRng::seed_from_u64(42)
    }

    #[test]
    fn random_pd_matrix_properties() {
        let mut rng = seeded_rng();
        let config = RandomMatrixConfig {
            size: 100,
            target_nnz: 500,
            positive_definite: true,
        };
        let m = generate_random_symmetric(&config, &mut rng).expect("generation failed");
        assert_eq!(m.nrows(), 100);
        assert_eq!(m.ncols(), 100);

        let dense = m.to_dense();
        for i in 0..100 {
            for j in 0..100 {
                assert_eq!(
                    dense[(i, j)],
                    dense[(j, i)],
                    "not symmetric at ({}, {})",
                    i,
                    j
                );
            }
        }

        for i in 0..100 {
            assert!(
                dense[(i, i)] > 0.0,
                "diagonal entry {} should be positive for PD",
                i
            );
        }
    }

    #[test]
    fn random_indefinite_matrix_has_mixed_signs() {
        let mut rng = seeded_rng();
        let config = RandomMatrixConfig {
            size: 50,
            target_nnz: 200,
            positive_definite: false,
        };
        let m = generate_random_symmetric(&config, &mut rng).expect("generation failed");
        let dense = m.to_dense();

        let mut has_positive = false;
        let mut has_negative = false;
        for i in 0..50 {
            if dense[(i, i)] > 0.0 {
                has_positive = true;
            }
            if dense[(i, i)] < 0.0 {
                has_negative = true;
            }
        }
        assert!(has_positive, "should have positive diagonal entries");
        assert!(has_negative, "should have negative diagonal entries");
    }

    #[test]
    fn generate_arrow_pattern() {
        let mut rng = seeded_rng();
        let m = generate_arrow(20, true, &mut rng);
        assert_eq!(m.nrows(), 20);
        let dense = m.to_dense();

        for j in 1..20 {
            assert!(
                dense[(0, j)] != 0.0,
                "first row entry (0, {}) should be nonzero",
                j
            );
            assert!(
                dense[(j, 0)] != 0.0,
                "first col entry ({}, 0) should be nonzero",
                j
            );
        }

        for i in 1..20 {
            for j in 1..20 {
                if i != j {
                    assert_eq!(
                        dense[(i, j)],
                        0.0,
                        "off-diagonal ({}, {}) should be zero in arrow tail",
                        i,
                        j
                    );
                }
            }
        }
    }

    #[test]
    fn generate_tridiagonal_pattern() {
        let mut rng = seeded_rng();
        let m = generate_tridiagonal(30, true, &mut rng);
        assert_eq!(m.nrows(), 30);
        let dense = m.to_dense();

        for i in 0..30 {
            let mut nnz_count = 0;
            for j in 0..30 {
                if dense[(i, j)] != 0.0 {
                    nnz_count += 1;
                    assert!(
                        (i as isize - j as isize).unsigned_abs() <= 1,
                        "nonzero at ({}, {}) outside tridiagonal band",
                        i,
                        j
                    );
                }
            }
            assert!(
                nnz_count <= 3,
                "row {} has {} nonzeros, expected <= 3",
                i,
                nnz_count
            );
        }
    }

    #[test]
    fn generate_banded_pattern() {
        let mut rng = seeded_rng();
        let m = generate_banded(40, 3, true, &mut rng);
        assert_eq!(m.nrows(), 40);
        let dense = m.to_dense();

        for i in 0..40 {
            for j in 0..40 {
                if (i as isize - j as isize).unsigned_abs() > 3 {
                    assert_eq!(
                        dense[(i, j)],
                        0.0,
                        "nonzero at ({}, {}) outside bandwidth 3",
                        i,
                        j
                    );
                }
            }
        }
    }

    #[test]
    fn generation_performance_under_1s() {
        let mut rng = seeded_rng();
        let config = RandomMatrixConfig {
            size: 1000,
            target_nnz: 10000,
            positive_definite: true,
        };
        let start = std::time::Instant::now();
        let _m = generate_random_symmetric(&config, &mut rng).expect("generation failed");
        let elapsed = start.elapsed();
        assert!(
            elapsed.as_secs_f64() < 1.0,
            "generation took {:.3}s, expected < 1s",
            elapsed.as_secs_f64()
        );
    }

    #[test]
    fn infeasible_config_returns_error() {
        let mut rng = seeded_rng();
        let config = RandomMatrixConfig {
            size: 1,
            target_nnz: 1,
            positive_definite: false,
        };
        let result = generate_random_symmetric(&config, &mut rng);
        assert!(result.is_err(), "1x1 indefinite should return error");
    }

    #[test]
    fn excessive_nnz_clamped() {
        let mut rng = seeded_rng();
        let config = RandomMatrixConfig {
            size: 5,
            target_nnz: 1000,
            positive_definite: true,
        };
        let m = generate_random_symmetric(&config, &mut rng).expect("generation failed");
        let dense = m.to_dense();
        let mut actual_nnz = 0;
        for i in 0..5 {
            for j in 0..5 {
                if dense[(i, j)] != 0.0 {
                    actual_nnz += 1;
                }
            }
        }
        assert!(
            actual_nnz <= 25,
            "actual nnz {} exceeds max possible 25",
            actual_nnz
        );
    }
}
