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
use faer::{Accum, Conj, Mat, Par};

use super::numeric::{AptpNumeric, FrontFactors};
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

    // Single allocation for per-supernode workspace. Each supernode needs at most
    // max_front_size entries. Since supernodes are processed sequentially, we reuse
    // this buffer across all supernodes.
    let max_size = numeric.stats().max_front_size;
    let mut work_buf = vec![0.0f64; max_size];

    // 1. Forward solve (postorder = index order)
    for ff in factors {
        forward_solve_supernode(ff, rhs, &mut work_buf);
    }

    // 2. Diagonal solve (any order)
    for ff in factors {
        diagonal_solve_supernode(ff, rhs, &mut work_buf);
    }

    // 3. Backward solve (reverse postorder)
    for s in (0..n_supernodes).rev() {
        backward_solve_supernode(&factors[s], rhs, &mut work_buf);
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
/// 1. Gather rhs[col_indices[i]] → local buffer
/// 2. Solve L11 * y = local (unit lower triangular, in-place)
/// 3. Write local → rhs[col_indices[i]]
/// 4. Scatter: rhs[row_indices[j]] -= (L21 * local)[j]
pub(crate) fn forward_solve_supernode(ff: &FrontFactors, rhs: &mut [f64], work: &mut [f64]) {
    let ne = ff.num_eliminated();
    if ne == 0 {
        return;
    }

    let col_indices = ff.col_indices();
    let row_indices = ff.row_indices();
    let r = row_indices.len();

    // Gather: local[i] = rhs[col_indices[i]]
    for i in 0..ne {
        work[i] = rhs[col_indices[i]];
    }

    // Solve L11 * y = local (unit lower triangular, in-place)
    if ne > 1 {
        let mut local_mat = Mat::<f64>::from_fn(ne, 1, |i, _| work[i]);
        solve_unit_lower_triangular_in_place_with_conj(
            ff.l11().as_ref(),
            Conj::No,
            local_mat.as_mut(),
            Par::Seq,
        );
        for i in 0..ne {
            work[i] = local_mat[(i, 0)];
        }
    }

    // Write back: rhs[col_indices[i]] = local[i]
    for i in 0..ne {
        rhs[col_indices[i]] = work[i];
    }

    // Scatter: rhs[row_indices[j]] -= (L21 * local)[j]
    if r > 0 {
        let local_mat = Mat::<f64>::from_fn(ne, 1, |i, _| work[i]);
        let mut tmp = Mat::<f64>::zeros(r, 1);
        matmul_with_conj(
            tmp.as_mut(),
            Accum::Replace,
            ff.l21().as_ref(),
            Conj::No,
            local_mat.as_ref(),
            Conj::No,
            1.0,
            Par::Seq,
        );
        for j in 0..r {
            rhs[row_indices[j]] -= tmp[(j, 0)];
        }
    }
}

/// Diagonal solve for a single supernode.
///
/// 1. Gather rhs[col_indices[i]] → local buffer
/// 2. D11.solve_in_place(local) (handles 1x1, 2x2, and zero pivots)
/// 3. Write local → rhs[col_indices[i]]
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
/// 1. Gather rhs[col_indices[i]] → local buffer
/// 2. Gather rhs[row_indices[j]] → tmp buffer (if r > 0)
/// 3. local -= L21^T * tmp (if r > 0)
/// 4. Solve L11^T * z = local (unit upper triangular via transpose)
/// 5. Write local → rhs[col_indices[i]]
pub(crate) fn backward_solve_supernode(ff: &FrontFactors, rhs: &mut [f64], work: &mut [f64]) {
    let ne = ff.num_eliminated();
    if ne == 0 {
        return;
    }

    let col_indices = ff.col_indices();
    let row_indices = ff.row_indices();
    let r = row_indices.len();

    // Gather local: local[i] = rhs[col_indices[i]]
    for i in 0..ne {
        work[i] = rhs[col_indices[i]];
    }

    // Gather-update from L21^T
    if r > 0 {
        let tmp = Mat::<f64>::from_fn(r, 1, |j, _| rhs[row_indices[j]]);
        let mut local_mat = Mat::<f64>::from_fn(ne, 1, |i, _| work[i]);

        // local -= L21^T * tmp
        matmul_with_conj(
            local_mat.as_mut(),
            Accum::Add,
            ff.l21().as_ref().transpose(),
            Conj::No,
            tmp.as_ref(),
            Conj::No,
            -1.0,
            Par::Seq,
        );

        for i in 0..ne {
            work[i] = local_mat[(i, 0)];
        }
    }

    // Solve L11^T * z = local (unit upper triangular via transpose)
    if ne > 1 {
        let mut local_mat = Mat::<f64>::from_fn(ne, 1, |i, _| work[i]);
        solve_unit_upper_triangular_in_place_with_conj(
            ff.l11().as_ref().transpose(),
            Conj::No,
            local_mat.as_mut(),
            Par::Seq,
        );
        for i in 0..ne {
            work[i] = local_mat[(i, 0)];
        }
    }

    // Write back
    for i in 0..ne {
        rhs[col_indices[i]] = work[i];
    }
}
