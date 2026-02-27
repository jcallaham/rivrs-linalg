//! Proptest strategies for generating random symmetric matrices.
//!
//! Provides custom [`proptest::strategy::Strategy`] implementations for generating
//! dense and sparse symmetric matrices with automatic shrinking support.
//!
//! Wraps the existing random matrix generators in [`crate::testing::generators`]
//! and [`crate::testing::perturbations`] for use with proptest's infrastructure.

use faer::Mat;
use faer::sparse::{SparseColMat, Triplet};
use proptest::prelude::*;
use rand::RngExt;

use super::perturbations::{generate_dense_symmetric_indefinite, generate_dense_symmetric_pd};

/// Create a proptest strategy that generates random dense symmetric
/// positive definite matrices with sizes in the given range.
///
/// Shrinks toward smaller matrix sizes.
pub fn arb_symmetric_pd(
    size_range: std::ops::RangeInclusive<usize>,
) -> impl Strategy<Value = Mat<f64>> {
    size_range.prop_flat_map(|n| {
        any::<u64>().prop_map(move |seed| {
            use rand::SeedableRng;
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            generate_dense_symmetric_pd(n.max(1), &mut rng)
        })
    })
}

/// Create a proptest strategy that generates random dense symmetric
/// indefinite matrices with sizes in the given range.
///
/// Shrinks toward smaller matrix sizes.
pub fn arb_symmetric_indefinite(
    size_range: std::ops::RangeInclusive<usize>,
) -> impl Strategy<Value = Mat<f64>> {
    size_range.prop_flat_map(|n| {
        any::<u64>().prop_map(move |seed| {
            use rand::SeedableRng;
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            // Need n >= 2 for indefinite (mixed signs need at least 2 diagonal entries)
            generate_dense_symmetric_indefinite(n.max(2), &mut rng)
        })
    })
}

/// Create a proptest strategy that generates random sparse symmetric matrices
/// suitable for `SparseLDLT`, with sizes and densities in the given ranges.
///
/// Generated matrices store full symmetric CSC (both upper and lower triangles),
/// matching the convention used throughout rivrs-sparse.
///
/// Shrinks toward smaller sizes and lower density.
pub fn arb_sparse_symmetric(
    size_range: std::ops::RangeInclusive<usize>,
    density_range: std::ops::RangeInclusive<f64>,
) -> impl Strategy<Value = SparseColMat<usize, f64>> {
    let density_start = *density_range.start();
    let density_end = *density_range.end();
    (size_range, density_start..=density_end).prop_flat_map(move |(n, density)| {
        any::<u64>().prop_map(move |seed| {
            use rand::SeedableRng;
            let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
            let n = n.max(2);
            generate_sparse_symmetric_from_density(n, density, &mut rng)
        })
    })
}

/// Internal helper: generate a sparse symmetric matrix with given size and density.
///
/// Density 0.0 = diagonal only, density 1.0 = fully dense.
/// Uses diagonal dominance with mixed signs (indefinite).
fn generate_sparse_symmetric_from_density(
    n: usize,
    density: f64,
    rng: &mut impl rand::Rng,
) -> SparseColMat<usize, f64> {
    let max_offdiag_pairs = n * (n - 1) / 2;
    let target_pairs = (max_offdiag_pairs as f64 * density.clamp(0.0, 1.0)) as usize;

    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
    let mut placed = std::collections::HashSet::new();
    let mut row_abs_sum = vec![0.0f64; n];

    // Place off-diagonal entries
    for _ in 0..target_pairs.saturating_mul(3) {
        if placed.len() >= target_pairs {
            break;
        }
        let i = rng.random_range(0..n);
        let j = rng.random_range(0..n);
        if i == j {
            continue;
        }
        let (lo, hi) = if i > j { (j, i) } else { (i, j) };
        if placed.contains(&(lo, hi)) {
            continue;
        }
        placed.insert((lo, hi));
        let v: f64 = rng.random_range(-1.0..1.0);
        triplets.push(Triplet::new(lo, hi, v));
        triplets.push(Triplet::new(hi, lo, v));
        row_abs_sum[lo] += v.abs();
        row_abs_sum[hi] += v.abs();
    }

    // Add diagonal entries (diagonal dominance, mixed signs for indefinite)
    let half = n / 2;
    for (i, &abs_sum) in row_abs_sum.iter().enumerate() {
        let margin = 1.0 + rng.random::<f64>();
        let diag = if i < half {
            abs_sum + margin
        } else {
            -(abs_sum + margin)
        };
        triplets.push(Triplet::new(i, i, diag));
    }

    SparseColMat::try_new_from_triplets(n, n, &triplets)
        .expect("sparse matrix generation should not fail")
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::strategy::ValueTree;
    use proptest::test_runner::TestRunner;

    #[test]
    fn arb_symmetric_pd_generates_symmetric_matrices() {
        let mut runner = TestRunner::default();
        let strategy = arb_symmetric_pd(5..=20);
        for _ in 0..10 {
            let value = strategy.new_tree(&mut runner).unwrap().current();
            let n = value.nrows();
            assert_eq!(n, value.ncols());
            for i in 0..n {
                for j in 0..n {
                    assert!(
                        (value[(i, j)] - value[(j, i)]).abs() < 1e-15,
                        "not symmetric at ({}, {})",
                        i,
                        j
                    );
                }
                assert!(
                    value[(i, i)] > 0.0,
                    "diagonal entry {} should be positive",
                    i
                );
            }
        }
    }

    #[test]
    fn arb_symmetric_indefinite_generates_mixed_signs() {
        let mut runner = TestRunner::default();
        let strategy = arb_symmetric_indefinite(10..=50);
        for _ in 0..10 {
            let value = strategy.new_tree(&mut runner).unwrap().current();
            let n = value.nrows();
            let mut has_pos = false;
            let mut has_neg = false;
            for i in 0..n {
                if value[(i, i)] > 0.0 {
                    has_pos = true;
                }
                if value[(i, i)] < 0.0 {
                    has_neg = true;
                }
            }
            assert!(has_pos && has_neg, "should have mixed diagonal signs");
        }
    }

    #[test]
    fn arb_sparse_symmetric_generates_valid_csc() {
        let mut runner = TestRunner::default();
        let strategy = arb_sparse_symmetric(5..=30, 0.1..=0.5);
        for _ in 0..10 {
            let value = strategy.new_tree(&mut runner).unwrap().current();
            let n = value.nrows();
            assert_eq!(n, value.ncols());
            // Verify symmetry in dense form
            let dense = value.to_dense();
            for i in 0..n {
                for j in 0..n {
                    assert!(
                        (dense[(i, j)] - dense[(j, i)]).abs() < 1e-15,
                        "sparse matrix not symmetric at ({}, {})",
                        i,
                        j
                    );
                }
            }
        }
    }
}
