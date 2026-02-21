//! Per-supernode triangular solve for multifrontal LDL^T factorization.
//!
//! Provides the core solve function [`aptp_solve`] and per-supernode helper
//! functions that traverse the assembly tree performing gather/solve/scatter
//! operations using faer's dense triangular solve and matmul kernels.
//!
//! The solve operates in **permuted coordinates** — the caller (SparseLDLT)
//! handles permutation (P) and optional scaling (S) around this core.
//!
//! # Algorithm
//!
//! Given P^T A P = L D L^T from multifrontal factorization:
//! 1. Forward solve (postorder): L y = b
//! 2. Diagonal solve: D z = y
//! 3. Backward solve (reverse postorder): L^T x = z
//!
//! # References
//!
//! - Duff, Hogg & Lopez (2020), Section 3: APTP solve algorithm
//! - Liu (1992), Sections 4-5: multifrontal solve traversal
//! - faer `SupernodalIntranodeLbltRef::solve_in_place_no_numeric_permute_with_conj`

use faer::dyn_stack::StackReq;
use faer::linalg::matmul::matmul_with_conj;
use faer::linalg::temp_mat_scratch;
use faer::linalg::triangular_solve::{
    solve_unit_lower_triangular_in_place_with_conj, solve_unit_upper_triangular_in_place_with_conj,
};
use faer::mat::{MatMut, MatRef};
use faer::{Accum, Conj, Par};

use super::numeric::{AptpNumeric, FrontFactors, INTRA_NODE_THRESHOLD};
use super::symbolic::AptpSymbolic;
use crate::error::SparseError;

/// Per-supernode triangular solve through the multifrontal factor structure.
///
/// Solves L D L^T x = b in-place, where b is provided in `rhs` (permuted
/// coordinate system — the caller handles P and S). Traverses the assembly
/// tree performing per-supernode gather/solve/scatter operations.
///
/// # Arguments
///
/// - `symbolic`: Symbolic analysis (supernode structure, traversal order)
/// - `numeric`: Numeric factors (per-supernode L11, D11, L21)
/// - `rhs`: Right-hand side, modified in-place to contain solution
/// - `stack`: Pre-allocated workspace (from `aptp_solve_scratch()`)
///
/// # Errors
///
/// Returns `SparseError::DimensionMismatch` if `rhs.len() != symbolic.nrows()`.
///
/// # References
///
/// - Duff, Hogg & Lopez (2020), Section 3: APTP solve algorithm
/// - Liu (1992), Sections 4-5: multifrontal solve traversal
pub fn aptp_solve(
    symbolic: &AptpSymbolic,
    numeric: &AptpNumeric,
    rhs: &mut [f64],
    _stack: &mut faer::dyn_stack::MemStack,
    par: Par,
) -> Result<(), SparseError> {
    let n = symbolic.nrows();
    if rhs.len() != n {
        return Err(SparseError::DimensionMismatch {
            expected: (n, 1),
            got: (rhs.len(), 1),
            context: "RHS length must match matrix dimension".to_string(),
        });
    }

    if n == 0 {
        return Ok(());
    }

    let factors = numeric.front_factors();
    let n_supernodes = factors.len();

    // Two workspace buffers for per-supernode operations. Each supernode needs
    // at most max_front_size entries. Since supernodes are processed sequentially,
    // we reuse these buffers across all supernodes. This avoids per-supernode
    // heap allocations in the solve hot path.
    let max_size = numeric.stats().max_front_size;
    let mut work = vec![0.0f64; max_size];
    let mut work2 = vec![0.0f64; max_size];

    // 1. Forward solve (postorder = index order)
    for ff in factors {
        forward_solve_supernode(ff, rhs, &mut work, &mut work2, par);
    }

    // 2. Diagonal solve (any order)
    for ff in factors {
        diagonal_solve_supernode(ff, rhs, &mut work);
    }

    // 3. Backward solve (reverse postorder)
    for s in (0..n_supernodes).rev() {
        backward_solve_supernode(&factors[s], rhs, &mut work, &mut work2, par);
    }

    Ok(())
}

/// Workspace requirement for `aptp_solve`.
///
/// Returns `StackReq` for the temporary per-supernode buffer.
/// Size is determined by `max_front_size` from factorization stats.
pub fn aptp_solve_scratch(numeric: &AptpNumeric, rhs_ncols: usize) -> StackReq {
    let max_size = numeric.stats().max_front_size;
    temp_mat_scratch::<f64>(max_size, rhs_ncols)
}

/// Forward solve for a single supernode.
///
/// 1. Gather rhs[col_indices[i]] → work buffer
/// 2. Solve L11 * y = work (unit lower triangular, in-place)
/// 3. Write work → rhs[col_indices[i]]
/// 4. Scatter: rhs[row_indices[j]] -= (L21 * work)[j]
pub(crate) fn forward_solve_supernode(
    ff: &FrontFactors,
    rhs: &mut [f64],
    work: &mut [f64],
    work2: &mut [f64],
    par: Par,
) {
    let ne = ff.num_eliminated();
    if ne == 0 {
        return;
    }

    let col_indices = ff.col_indices();
    let row_indices = ff.row_indices();
    let r = row_indices.len();

    // Apply intra-node threshold: small fronts use Par::Seq
    let effective_par = if ne + r < INTRA_NODE_THRESHOLD {
        Par::Seq
    } else {
        par
    };

    // Gather: work[i] = rhs[col_indices[i]]
    for i in 0..ne {
        work[i] = rhs[col_indices[i]];
    }

    // Solve L11 * y = work (unit lower triangular, in-place)
    if ne > 1 {
        let mut local = MatMut::from_column_major_slice_mut(&mut work[..ne], ne, 1);
        solve_unit_lower_triangular_in_place_with_conj(
            ff.l11().as_ref(),
            Conj::No,
            local.as_mut(),
            effective_par,
        );
    }

    // Write back: rhs[col_indices[i]] = work[i]
    for i in 0..ne {
        rhs[col_indices[i]] = work[i];
    }

    // Scatter: rhs[row_indices[j]] -= (L21 * work[..ne])[j]
    if r > 0 {
        let local = MatRef::from_column_major_slice(&work[..ne], ne, 1);
        let tmp = &mut work2[..r];
        tmp.fill(0.0);
        let mut tmp_mat = MatMut::from_column_major_slice_mut(tmp, r, 1);
        matmul_with_conj(
            tmp_mat.as_mut(),
            Accum::Replace,
            ff.l21().as_ref(),
            Conj::No,
            local,
            Conj::No,
            1.0,
            effective_par,
        );
        for j in 0..r {
            rhs[row_indices[j]] -= tmp[j];
        }
    }
}

/// Diagonal solve for a single supernode.
///
/// 1. Gather rhs[col_indices[i]] → work buffer
/// 2. D11.solve_in_place(work) (handles 1x1, 2x2, and zero pivots)
/// 3. Write work → rhs[col_indices[i]]
pub(crate) fn diagonal_solve_supernode(ff: &FrontFactors, rhs: &mut [f64], work: &mut [f64]) {
    let ne = ff.num_eliminated();
    if ne == 0 {
        return;
    }

    let col_indices = ff.col_indices();

    // Gather
    for i in 0..ne {
        work[i] = rhs[col_indices[i]];
    }

    // D-solve
    ff.d11().solve_in_place(&mut work[..ne]);

    // Write back
    for i in 0..ne {
        rhs[col_indices[i]] = work[i];
    }
}

/// Backward solve for a single supernode.
///
/// 1. Gather rhs[col_indices[i]] → work buffer
/// 2. Gather rhs[row_indices[j]] → work2 buffer (if r > 0)
/// 3. work -= L21^T * work2 (if r > 0)
/// 4. Solve L11^T * z = work (unit upper triangular via transpose)
/// 5. Write work → rhs[col_indices[i]]
pub(crate) fn backward_solve_supernode(
    ff: &FrontFactors,
    rhs: &mut [f64],
    work: &mut [f64],
    work2: &mut [f64],
    par: Par,
) {
    let ne = ff.num_eliminated();
    if ne == 0 {
        return;
    }

    let col_indices = ff.col_indices();
    let row_indices = ff.row_indices();
    let r = row_indices.len();

    // Apply intra-node threshold: small fronts use Par::Seq
    let effective_par = if ne + r < INTRA_NODE_THRESHOLD {
        Par::Seq
    } else {
        par
    };

    // Gather local: work[i] = rhs[col_indices[i]]
    for i in 0..ne {
        work[i] = rhs[col_indices[i]];
    }

    // Gather-update from L21^T
    if r > 0 {
        // Gather row values into work2
        for j in 0..r {
            work2[j] = rhs[row_indices[j]];
        }
        let tmp = MatRef::from_column_major_slice(&work2[..r], r, 1);
        let mut local = MatMut::from_column_major_slice_mut(&mut work[..ne], ne, 1);

        // work -= L21^T * work2
        matmul_with_conj(
            local.as_mut(),
            Accum::Add,
            ff.l21().as_ref().transpose(),
            Conj::No,
            tmp,
            Conj::No,
            -1.0,
            effective_par,
        );
    }

    // Solve L11^T * z = work (unit upper triangular via transpose)
    if ne > 1 {
        let mut local = MatMut::from_column_major_slice_mut(&mut work[..ne], ne, 1);
        solve_unit_upper_triangular_in_place_with_conj(
            ff.l11().as_ref().transpose(),
            Conj::No,
            local.as_mut(),
            effective_par,
        );
    }

    // Write back
    for i in 0..ne {
        rhs[col_indices[i]] = work[i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::dyn_stack::{MemBuffer, MemStack};
    use faer::sparse::linalg::cholesky::SymmetricOrdering;
    use faer::sparse::{SparseColMat, Triplet};

    use crate::aptp::factor::AptpOptions;
    use crate::aptp::numeric::AptpNumeric;
    use crate::aptp::symbolic::AptpSymbolic;

    /// Build a full symmetric CSC matrix from lower-triangular entries.
    fn sparse_from_lower(n: usize, lower: &[(usize, usize, f64)]) -> SparseColMat<usize, f64> {
        let mut triplets = Vec::new();
        for &(i, j, v) in lower {
            triplets.push(Triplet::new(i, j, v));
            if i != j {
                triplets.push(Triplet::new(j, i, v));
            }
        }
        SparseColMat::try_new_from_triplets(n, n, &triplets).expect("triplets")
    }

    /// Factor a matrix with identity ordering, returning (symbolic, numeric).
    fn factor_identity(matrix: &SparseColMat<usize, f64>) -> (AptpSymbolic, AptpNumeric) {
        let symbolic =
            AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Identity).expect("analyze");
        let numeric =
            AptpNumeric::factor(&symbolic, matrix, &AptpOptions::default(), None).expect("factor");
        (symbolic, numeric)
    }

    /// Sparse matvec: y = A * x (full symmetric CSC).
    fn matvec(a: &SparseColMat<usize, f64>, x: &[f64]) -> Vec<f64> {
        let n = a.nrows();
        let col_ptr = a.symbolic().col_ptr();
        let row_idx = a.symbolic().row_idx();
        let values = a.val();
        let mut y = vec![0.0; n];
        for j in 0..n {
            for idx in col_ptr[j]..col_ptr[j + 1] {
                let i = row_idx[idx];
                y[i] += values[idx] * x[j];
            }
        }
        y
    }

    #[test]
    fn forward_solve_gather_scatter() {
        // 3x3 PD: [[4,1,0.5],[1,3,0],[0.5,0,2]]
        let matrix = sparse_from_lower(
            3,
            &[
                (0, 0, 4.0),
                (1, 0, 1.0),
                (1, 1, 3.0),
                (2, 0, 0.5),
                (2, 2, 2.0),
            ],
        );
        let (symbolic, numeric) = factor_identity(&matrix);
        let factors = numeric.front_factors();

        // Known RHS
        let x_exact = [1.0, 2.0, 3.0];
        let b = matvec(&matrix, &x_exact);

        let (perm_fwd, _) = symbolic.perm_vecs();
        let n = 3;
        let mut rhs_perm: Vec<f64> = (0..n).map(|i| b[perm_fwd[i]]).collect();

        // Forward solve only
        let max_size = numeric.stats().max_front_size;
        let mut work = vec![0.0f64; max_size];
        let mut work2 = vec![0.0f64; max_size];
        for ff in factors {
            forward_solve_supernode(ff, &mut rhs_perm, &mut work, &mut work2, Par::Seq);
        }

        // After forward solve, rhs_perm should be L^{-1} b_perm.
        // We can't easily check exact values without reconstructing L,
        // but completing the full solve should recover x.
        for ff in factors {
            diagonal_solve_supernode(ff, &mut rhs_perm, &mut work);
        }
        let n_sn = factors.len();
        for s in (0..n_sn).rev() {
            backward_solve_supernode(&factors[s], &mut rhs_perm, &mut work, &mut work2, Par::Seq);
        }

        let mut x_computed = vec![0.0f64; n];
        for new in 0..n {
            x_computed[perm_fwd[new]] = rhs_perm[new];
        }

        let err: f64 = x_computed
            .iter()
            .zip(x_exact.iter())
            .map(|(c, e)| (c - e).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(err < 1e-12, "full pipeline error: {:.2e}", err);
    }

    #[test]
    fn diagonal_solve_known_values() {
        // 4x4 diagonal matrix: diag = [2, 5, 3, 7]
        let matrix = sparse_from_lower(4, &[(0, 0, 2.0), (1, 1, 5.0), (2, 2, 3.0), (3, 3, 7.0)]);
        let (_, numeric) = factor_identity(&matrix);
        let factors = numeric.front_factors();

        // rhs = [4, 10, 9, 21] → x = [2, 2, 3, 3]
        let mut rhs = vec![4.0, 10.0, 9.0, 21.0];
        let max_size = numeric.stats().max_front_size;
        let mut work = vec![0.0f64; max_size];

        for ff in factors {
            diagonal_solve_supernode(ff, &mut rhs, &mut work);
        }

        let expected = [2.0, 2.0, 3.0, 3.0];
        for i in 0..4 {
            assert!(
                (rhs[i] - expected[i]).abs() < 1e-14,
                "rhs[{}] = {}, expected {}",
                i,
                rhs[i],
                expected[i]
            );
        }
    }

    #[test]
    fn backward_solve_reverses_forward() {
        // 4x4 tridiagonal: diag=[4,4,4,4], off-diag=1
        let matrix = sparse_from_lower(
            4,
            &[
                (0, 0, 4.0),
                (1, 0, 1.0),
                (1, 1, 4.0),
                (2, 1, 1.0),
                (2, 2, 4.0),
                (3, 2, 1.0),
                (3, 3, 4.0),
            ],
        );
        let (symbolic, numeric) = factor_identity(&matrix);
        let factors = numeric.front_factors();

        let x_exact = [1.0, -1.0, 2.0, -2.0];
        let b = matvec(&matrix, &x_exact);

        let (perm_fwd, _) = symbolic.perm_vecs();
        let n = 4;
        let mut rhs: Vec<f64> = (0..n).map(|i| b[perm_fwd[i]]).collect();

        let max_size = numeric.stats().max_front_size;
        let mut work = vec![0.0f64; max_size];
        let mut work2 = vec![0.0f64; max_size];

        // Full pipeline
        for ff in factors {
            forward_solve_supernode(ff, &mut rhs, &mut work, &mut work2, Par::Seq);
        }
        for ff in factors {
            diagonal_solve_supernode(ff, &mut rhs, &mut work);
        }
        let n_sn = factors.len();
        for s in (0..n_sn).rev() {
            backward_solve_supernode(&factors[s], &mut rhs, &mut work, &mut work2, Par::Seq);
        }

        let mut x_computed = vec![0.0f64; n];
        for new in 0..n {
            x_computed[perm_fwd[new]] = rhs[new];
        }

        let err: f64 = x_computed
            .iter()
            .zip(x_exact.iter())
            .map(|(c, e)| (c - e).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(err < 1e-12, "tridiag pipeline error: {:.2e}", err);
    }

    #[test]
    fn per_supernode_solve_indefinite_2x2() {
        // 4x4 indefinite: requires 2x2 pivots
        // [[ 1, 2, 0, 0],
        //  [ 2,-3, 1, 0],
        //  [ 0, 1, 5, 1],
        //  [ 0, 0, 1,-2]]
        let matrix = sparse_from_lower(
            4,
            &[
                (0, 0, 1.0),
                (1, 0, 2.0),
                (1, 1, -3.0),
                (2, 1, 1.0),
                (2, 2, 5.0),
                (3, 2, 1.0),
                (3, 3, -2.0),
            ],
        );
        let (symbolic, numeric) = factor_identity(&matrix);
        let factors = numeric.front_factors();

        let x_exact = [1.0, -1.0, 2.0, 0.5];
        let b = matvec(&matrix, &x_exact);

        let (perm_fwd, _) = symbolic.perm_vecs();
        let n = 4;
        let mut rhs: Vec<f64> = (0..n).map(|i| b[perm_fwd[i]]).collect();

        let max_size = numeric.stats().max_front_size;
        let mut work = vec![0.0f64; max_size];
        let mut work2 = vec![0.0f64; max_size];

        for ff in factors {
            forward_solve_supernode(ff, &mut rhs, &mut work, &mut work2, Par::Seq);
        }
        for ff in factors {
            diagonal_solve_supernode(ff, &mut rhs, &mut work);
        }
        let n_sn = factors.len();
        for s in (0..n_sn).rev() {
            backward_solve_supernode(&factors[s], &mut rhs, &mut work, &mut work2, Par::Seq);
        }

        let mut x_computed = vec![0.0f64; n];
        for new in 0..n {
            x_computed[perm_fwd[new]] = rhs[new];
        }

        let err: f64 = x_computed
            .iter()
            .zip(x_exact.iter())
            .map(|(c, e)| (c - e).powi(2))
            .sum::<f64>()
            .sqrt();
        assert!(err < 1e-12, "indefinite 2x2 pipeline error: {:.2e}", err);
    }

    #[test]
    fn front_factors_index_consistency() {
        // 8x8 arrow matrix: tests multi-supernode factorization
        let mut entries = Vec::new();
        for i in 0..8 {
            entries.push((i, i, if i < 7 { 10.0 } else { 20.0 }));
            if i < 7 {
                entries.push((7, i, 1.0));
            }
        }
        let matrix = sparse_from_lower(8, &entries);
        let (symbolic, numeric) = factor_identity(&matrix);
        let factors = numeric.front_factors();

        // Verify structural properties
        let n = symbolic.nrows();
        let mut total_eliminated = 0;
        for ff in factors {
            let ne = ff.num_eliminated();
            total_eliminated += ne;

            // L11 dimensions must match num_eliminated
            assert_eq!(ff.l11().nrows(), ne);
            assert_eq!(ff.l11().ncols(), ne);

            // L21 rows = row_indices length
            assert_eq!(ff.l21().nrows(), ff.row_indices().len());
            assert_eq!(ff.l21().ncols(), ne);

            // All col_indices must be valid
            for &ci in ff.col_indices() {
                assert!(ci < n, "col_index {} >= n={}", ci, n);
            }
            for &ri in ff.row_indices() {
                assert!(ri < n, "row_index {} >= n={}", ri, n);
            }
        }
        assert_eq!(total_eliminated, n, "total eliminated must equal n");
    }

    #[test]
    fn aptp_solve_dimension_mismatch() {
        let matrix = sparse_from_lower(3, &[(0, 0, 1.0), (1, 1, 2.0), (2, 2, 3.0)]);
        let (symbolic, numeric) = factor_identity(&matrix);

        let mut rhs = vec![1.0, 2.0]; // wrong size
        let scratch = aptp_solve_scratch(&numeric, 1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);

        let result = aptp_solve(&symbolic, &numeric, &mut rhs, stack, Par::Seq);
        assert!(
            matches!(result, Err(SparseError::DimensionMismatch { .. })),
            "expected DimensionMismatch, got {:?}",
            result
        );
    }

    #[test]
    fn aptp_solve_empty() {
        let triplets: Vec<Triplet<usize, usize, f64>> = vec![];
        let matrix = SparseColMat::try_new_from_triplets(0, 0, &triplets).expect("empty matrix");
        let symbolic =
            AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Identity).expect("analyze");
        let numeric =
            AptpNumeric::factor(&symbolic, &matrix, &AptpOptions::default(), None).expect("factor");

        let mut rhs: Vec<f64> = vec![];
        let scratch = aptp_solve_scratch(&numeric, 1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);

        let result = aptp_solve(&symbolic, &numeric, &mut rhs, stack, Par::Seq);
        assert!(result.is_ok(), "0x0 solve should succeed");
    }
}
