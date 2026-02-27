//! Shared test utilities for integration tests.

#![allow(dead_code)]

use faer::sparse::{SparseColMat, Triplet};

/// Build a SparseColMat from dense lower-triangle triplets (mirrors to upper).
///
/// Input: `(i, j, val)` triplets where `i >= j`.
/// Produces a full symmetric CSC matrix (both triangles stored).
pub fn sparse_from_lower_triplets(
    n: usize,
    entries: &[(usize, usize, f64)],
) -> SparseColMat<usize, f64> {
    let mut triplets = Vec::new();
    for &(i, j, v) in entries {
        triplets.push(Triplet::new(i, j, v));
        if i != j {
            triplets.push(Triplet::new(j, i, v));
        }
    }
    SparseColMat::try_new_from_triplets(n, n, &triplets).unwrap()
}

/// Compute b = A*x for a full symmetric CSC matrix (both triangles stored).
///
/// Uses direct CSC iteration — O(nnz) work, no dense conversion.
pub fn sparse_matvec(a: &SparseColMat<usize, f64>, x: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let symbolic = a.symbolic();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let values = a.val();

    let mut result = vec![0.0f64; n];
    for j in 0..n {
        for idx in col_ptrs[j]..col_ptrs[j + 1] {
            let i = row_indices[idx];
            let v = values[idx];
            result[i] += v * x[j];
        }
    }
    result
}
