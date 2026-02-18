//! Diagnostic test: isolate whether backward error is in factorization or solve.
//!
//! For CI matrices that fail backward error thresholds, this test:
//! 1. Runs the full solver and records backward error + forward error
//! 2. Performs one step of iterative refinement (solve A*dx = r, x += dx)
//! 3. Compares improvement — if refinement helps dramatically, factorization
//!    is accurate and the solve is correct; the issue is the factorization
//!    quality itself (accumulated rounding). If refinement barely helps,
//!    something more fundamental is wrong.
//! 4. Per-supernode reconstruction error for the largest supernodes.

mod common;

use faer::Col;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::SparseColMat;

use rivrs_sparse::aptp::{
    AnalyzeOptions, AptpNumeric, AptpOptions, AptpSymbolic, FactorOptions,
    OrderingStrategy, SparseLDLT,
};
use rivrs_sparse::validate::sparse_backward_error;

/// Sparse matrix-vector product: y = A*x using CSC iteration.
fn csc_matvec(a: &SparseColMat<usize, f64>, x: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let symbolic = a.symbolic();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let values = a.val();

    let mut result = vec![0.0f64; n];
    for j in 0..n {
        for idx in col_ptrs[j]..col_ptrs[j + 1] {
            let i = row_indices[idx];
            result[i] += values[idx] * x[j];
        }
    }
    result
}

/// Compute ||v||_2.
fn norm2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// Run the full diagnostic for a single matrix.
fn diagnose_matrix(
    name: &str,
    a: &SparseColMat<usize, f64>,
    ordering: OrderingStrategy,
) {
    let n = a.nrows();

    // Build a known exact solution and RHS
    let x_exact: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
    let b_vec = csc_matvec(a, &x_exact);
    let b = Col::from_fn(n, |i| b_vec[i]);

    eprintln!("\n{}", "=".repeat(72));
    eprintln!("DIAGNOSING: {} (n={}, ordering={})", name, n, ordering_name(&ordering));
    eprintln!("{}", "=".repeat(72));

    // 1. Analyze + factor
    let analyze_opts = AnalyzeOptions {
        ordering: ordering.clone(),
    };
    let factor_opts = FactorOptions::default();

    let mut solver = SparseLDLT::analyze_with_matrix(a, &analyze_opts)
        .unwrap_or_else(|e| panic!("analyze '{}': {}", name, e));
    solver.factor(a, &factor_opts)
        .unwrap_or_else(|e| panic!("factor '{}': {}", name, e));

    let stats = solver.stats().expect("stats after factor");
    eprintln!("  Factorization stats:");
    eprintln!("    1x1 pivots:    {}", stats.total_1x1_pivots);
    eprintln!("    2x2 pivots:    {}", stats.total_2x2_pivots);
    eprintln!("    delayed:       {}", stats.total_delayed);
    eprintln!("    zero pivots:   {}", stats.zero_pivots);
    eprintln!("    max front:     {}", stats.max_front_size);

    // 2. Solve
    let scratch_req = solver.solve_scratch(1);
    let mut mem = MemBuffer::new(scratch_req);
    let stack = MemStack::new(&mut mem);
    let x_computed = solver.solve(&b, stack)
        .unwrap_or_else(|e| panic!("solve '{}': {}", name, e));

    // 3. Backward error
    let be = sparse_backward_error(a, &x_computed, &b);
    eprintln!("\n  Backward error:       {:.4e}", be);

    // 4. Forward error
    let forward_err: f64 = (0..n)
        .map(|i| (x_computed[i] - x_exact[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    let x_exact_norm = norm2(&x_exact);
    let rel_forward_err = forward_err / x_exact_norm;
    eprintln!("  Relative forward err: {:.4e}", rel_forward_err);

    // 5. Compute residual: r = b - A*x_computed
    let ax = csc_matvec(a, &(0..n).map(|i| x_computed[i]).collect::<Vec<_>>());
    let residual: Vec<f64> = (0..n).map(|i| b_vec[i] - ax[i]).collect();
    let residual_norm = norm2(&residual);
    let b_norm = norm2(&b_vec);
    eprintln!("  ||r|| / ||b||:        {:.4e}", residual_norm / b_norm);

    // 6. Iterative refinement: solve A * dx = r
    let r_col = Col::from_fn(n, |i| residual[i]);
    let scratch_req2 = solver.solve_scratch(1);
    let mut mem2 = MemBuffer::new(scratch_req2);
    let stack2 = MemStack::new(&mut mem2);
    let dx = solver.solve(&r_col, stack2)
        .unwrap_or_else(|e| panic!("refinement solve '{}': {}", name, e));

    // 7. Refined solution: x_refined = x_computed + dx
    let x_refined = Col::from_fn(n, |i| x_computed[i] + dx[i]);
    let be_refined = sparse_backward_error(a, &x_refined, &b);
    let forward_refined: f64 = (0..n)
        .map(|i| (x_refined[i] - x_exact[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    let rel_forward_refined = forward_refined / x_exact_norm;

    eprintln!("\n  After 1 refinement step:");
    eprintln!("    Backward error:     {:.4e} (was {:.4e}, ratio={:.1}x)",
        be_refined, be, be / be_refined.max(1e-20));
    eprintln!("    Forward error:      {:.4e} (was {:.4e}, ratio={:.1}x)",
        rel_forward_refined, rel_forward_err, rel_forward_err / rel_forward_refined.max(1e-20));

    // 8. Second refinement step
    let ax2 = csc_matvec(a, &(0..n).map(|i| x_refined[i]).collect::<Vec<_>>());
    let residual2: Vec<f64> = (0..n).map(|i| b_vec[i] - ax2[i]).collect();
    let r2_col = Col::from_fn(n, |i| residual2[i]);
    let scratch_req3 = solver.solve_scratch(1);
    let mut mem3 = MemBuffer::new(scratch_req3);
    let stack3 = MemStack::new(&mut mem3);
    let dx2 = solver.solve(&r2_col, stack3)
        .unwrap_or_else(|e| panic!("refinement step 2 '{}': {}", name, e));
    let x_refined2 = Col::from_fn(n, |i| x_refined[i] + dx2[i]);
    let be_refined2 = sparse_backward_error(a, &x_refined2, &b);
    let forward_refined2: f64 = (0..n)
        .map(|i| (x_refined2[i] - x_exact[i]).powi(2))
        .sum::<f64>()
        .sqrt();
    let rel_forward_refined2 = forward_refined2 / x_exact_norm;
    eprintln!("\n  After 2 refinement steps:");
    eprintln!("    Backward error:     {:.4e}", be_refined2);
    eprintln!("    Forward error:      {:.4e}", rel_forward_refined2);

    // 9. Interpretation
    eprintln!("\n  DIAGNOSIS:");
    if be < 1e-10 {
        eprintln!("    Matrix passes backward error threshold. No bug.");
    } else if be_refined < be * 0.01 {
        // One refinement step improved by 100x+ → factorization quality is the issue
        // (the solve is propagating correctly, the factors are just not precise enough)
        eprintln!("    Refinement improved BE by {:.0}x.", be / be_refined.max(1e-20));
        eprintln!("    This indicates the FACTORIZATION has accumulated rounding error");
        eprintln!("    (large fronts, many BLAS-2 updates) but the factors are");
        eprintln!("    still useful for iterative refinement. The SOLVE is correct.");
        eprintln!("    FIX: Restore inner block processing (ib=32 BLAS-3 sub-blocks)");
        eprintln!("    and/or add iterative refinement at the SparseLDLT level.");
    } else if be_refined < be * 0.1 {
        eprintln!("    Refinement improved BE by {:.1}x (modest).", be / be_refined.max(1e-20));
        eprintln!("    The factorization is somewhat imprecise but the factors contain");
        eprintln!("    useful information. Iterative refinement can help.");
    } else {
        eprintln!("    Refinement barely helped (ratio={:.1}x).", be / be_refined.max(1e-20));
        eprintln!("    This suggests the factorization is severely inaccurate,");
        eprintln!("    OR there is a bug in the SOLVE phase.");
        eprintln!("    FURTHER INVESTIGATION NEEDED.");
    }
    eprintln!();
}

fn ordering_name(o: &OrderingStrategy) -> &'static str {
    match o {
        OrderingStrategy::Amd => "AMD",
        OrderingStrategy::Metis => "METIS",
        OrderingStrategy::MatchOrderMetis => "MatchOrderMetis",
        OrderingStrategy::UserSupplied(_) => "UserSupplied",
    }
}

/// Per-supernode reconstruction diagnostic. For each supernode with
/// front_size > threshold, reconstruct PAP^T locally and compare with L*D*L^T.
fn diagnose_per_supernode_reconstruction(
    name: &str,
    a: &SparseColMat<usize, f64>,
    ordering: OrderingStrategy,
    front_threshold: usize,
) {
    let n = a.nrows();
    eprintln!("\n  Per-supernode reconstruction (front > {}):", front_threshold);

    // We need the internal factors. Use the three-phase API to get them.
    let analyze_opts = AnalyzeOptions {
        ordering: ordering.clone(),
    };
    let factor_opts = FactorOptions::default();

    let mut solver = SparseLDLT::analyze_with_matrix(a, &analyze_opts)
        .unwrap_or_else(|e| panic!("analyze '{}': {}", name, e));
    solver.factor(a, &factor_opts)
        .unwrap_or_else(|e| panic!("factor '{}': {}", name, e));

    // Access internals: we need to build the AptpSymbolic + AptpNumeric directly
    // to get per-supernode data. The SparseLDLT wraps them, so let's factor
    // directly via AptpNumeric.
    // We can't easily extract internals from SparseLDLT. Instead, factor directly:
    // Build the ordering permutation first.
    use rivrs_sparse::aptp::ordering::{match_order_metis, metis_ordering};
    use faer::sparse::linalg::cholesky::SymmetricOrdering;

    let (aptp_symbolic, scaling) = match &ordering {
        OrderingStrategy::MatchOrderMetis => {
            let result = match_order_metis(a).unwrap();
            let sym = AptpSymbolic::analyze(
                a.symbolic(),
                SymmetricOrdering::Custom(result.ordering.as_ref()),
            ).unwrap();
            let (perm_fwd, _) = sym.perm_vecs();
            let elim_scaling: Vec<f64> = (0..n).map(|i| result.scaling[perm_fwd[i]]).collect();
            (sym, Some(elim_scaling))
        }
        OrderingStrategy::Metis => {
            let perm = metis_ordering(a.symbolic()).unwrap();
            let sym = AptpSymbolic::analyze(
                a.symbolic(),
                SymmetricOrdering::Custom(perm.as_ref()),
            ).unwrap();
            (sym, None)
        }
        _ => {
            let sym = AptpSymbolic::analyze(
                a.symbolic(),
                SymmetricOrdering::Amd,
            ).unwrap();
            (sym, None)
        }
    };

    let aptp_options = AptpOptions {
        threshold: factor_opts.threshold,
        fallback: factor_opts.fallback,
        outer_block_size: factor_opts.outer_block_size,
        inner_block_size: factor_opts.inner_block_size,
        ..AptpOptions::default()
    };

    let numeric = AptpNumeric::factor(
        &aptp_symbolic,
        a,
        &aptp_options,
        scaling.as_deref(),
    ).unwrap();

    let factors = numeric.front_factors();
    let (perm_fwd, _perm_inv) = aptp_symbolic.perm_vecs();

    // For each large supernode, reconstruct PAP^T locally
    let mut large_count = 0;
    let mut max_recon_err = 0.0f64;
    let mut worst_sn = 0;
    let mut worst_front = 0;

    for (s, ff) in factors.iter().enumerate() {
        let ne = ff.num_eliminated();
        let r = ff.row_indices().len();
        let front_size = ne + r;

        if front_size < front_threshold || ne == 0 {
            continue;
        }
        large_count += 1;

        // Collect all global indices for this front:
        // - col_indices[0..ne]: eliminated columns
        // - row_indices[0..r]: L21 rows
        let col_indices = ff.col_indices();
        let row_indices = ff.row_indices();

        let m = ne + r;
        let mut global_indices: Vec<usize> = Vec::with_capacity(m);
        global_indices.extend_from_slice(col_indices);
        global_indices.extend_from_slice(row_indices);

        // Build PAP^T submatrix from original matrix.
        // global_indices are in permuted order. For PAP^T[i,j]:
        // PAP^T[i,j] = A[perm_fwd[gi], perm_fwd[gj]]
        // We need to extract these entries from the sparse matrix.
        let a_dense = extract_submatrix(a, &global_indices, &perm_fwd, scaling.as_deref());

        // Build L*D*L^T for the eliminated portion.
        // L = [L11; L21] is m x ne, D = ne x ne
        let l11 = ff.l11();
        let l21 = ff.l21();
        let d11 = ff.d11();

        // Build full L (m x ne)
        let mut l_full = faer::Mat::<f64>::zeros(m, ne);
        for i in 0..ne {
            for j in 0..ne {
                l_full[(i, j)] = l11[(i, j)];
            }
        }
        for i in 0..r {
            for j in 0..ne {
                l_full[(ne + i, j)] = l21[(i, j)];
            }
        }

        // Build D (ne x ne)
        let mut d_mat = faer::Mat::<f64>::zeros(ne, ne);
        let mut col = 0;
        while col < ne {
            match d11.get_pivot_type(col) {
                rivrs_sparse::aptp::PivotType::OneByOne => {
                    d_mat[(col, col)] = d11.get_1x1(col);
                    col += 1;
                }
                rivrs_sparse::aptp::PivotType::TwoByTwo { partner } => {
                    if partner > col {
                        let blk = d11.get_2x2(col);
                        d_mat[(col, col)] = blk.a;
                        d_mat[(col, col + 1)] = blk.b;
                        d_mat[(col + 1, col)] = blk.b;
                        d_mat[(col + 1, col + 1)] = blk.c;
                        col += 2;
                    } else {
                        col += 1;
                    }
                }
                _ => { col += 1; }
            }
        }

        // Compute L*D
        let mut ld = faer::Mat::<f64>::zeros(m, ne);
        faer::linalg::matmul::matmul(
            ld.as_mut(),
            faer::Accum::Replace,
            l_full.as_ref(),
            d_mat.as_ref(),
            1.0,
            faer::Par::Seq,
        );

        // Compute L*D*L^T (m x m)
        let mut ldlt = faer::Mat::<f64>::zeros(m, m);
        faer::linalg::matmul::matmul(
            ldlt.as_mut(),
            faer::Accum::Replace,
            ld.as_ref(),
            l_full.as_ref().transpose(),
            1.0,
            faer::Par::Seq,
        );

        // Reconstruction error: ||PAP^T[0..m, 0..m] - L*D*L^T||_F / ||PAP^T[0..m, 0..m]||_F
        // BUT: L*D*L^T only covers the eliminated portion. The trailing (m-ne) x (m-ne) block
        // is the Schur complement that gets passed to the parent. So we only compare
        // the upper-left ne x ne block and the L21 area.
        // Actually, LDL^T covers the full m x m front: L is m x ne, D is ne x ne,
        // so LDL^T is m x m. This accounts for the entire frontal contribution from
        // this supernode's elimination. The difference PAP^T - LDL^T should be the
        // contribution block (Schur complement) that goes to the parent.
        // For the eliminated rows/cols (0..ne x 0..ne), the match should be exact.
        // For the mixed part (ne..m x 0..ne) also exact.
        // The Schur complement is F22 - L21 * D11 * L21^T, which is what goes up.

        // So check the first ne columns of PAP^T vs LDL^T:
        let mut diff_sq = 0.0f64;
        let mut ref_sq = 0.0f64;
        for i in 0..m {
            for j in 0..ne {
                let pap_val = a_dense[(i, j)];
                let ldlt_val = ldlt[(i, j)];
                diff_sq += (pap_val - ldlt_val).powi(2);
                ref_sq += pap_val.powi(2);
            }
        }
        // Also check the symmetric part (j > ne not applicable — LDL^T has these)
        // Actually, also check columns ne..m for the LDL^T contribution:
        // LDL^T[i, j] for i,j >= ne comes from L21 * D * L21^T, which is part of the Schur complement
        // The frontal matrix is PAP^T = LDL^T + contribution, so PAP^T[i,j] for i,j >= ne
        // = LDL^T[i,j] + Schur[i-ne, j-ne]. We can't check that without the Schur complement.
        // Just check the eliminated portion:

        let recon_err = if ref_sq > 0.0 {
            diff_sq.sqrt() / ref_sq.sqrt()
        } else {
            0.0
        };

        if large_count <= 10 || recon_err > 1e-10 {
            eprintln!(
                "    SN {:>4}: front={:>5}, ne={:>5}, r={:>5}, recon_err={:.4e}",
                s, front_size, ne, r, recon_err
            );
        }

        if recon_err > max_recon_err {
            max_recon_err = recon_err;
            worst_sn = s;
            worst_front = front_size;
        }
    }

    if large_count == 0 {
        eprintln!("    No supernodes with front > {}", front_threshold);
    } else {
        eprintln!("    {} large supernodes checked", large_count);
        eprintln!(
            "    Worst reconstruction: SN {} (front={}): {:.4e}",
            worst_sn, worst_front, max_recon_err
        );
        if max_recon_err < 1e-10 {
            eprintln!("    => Per-supernode factorization is ACCURATE.");
            eprintln!("    => If backward error is bad, the bug is in ASSEMBLY or SOLVE.");
        } else {
            eprintln!("    => Per-supernode factorization has ERRORS.");
            eprintln!("    => The dense APTP kernel is producing inaccurate factors.");
        }
    }
}

/// Extract a dense submatrix from a sparse matrix.
/// global_indices are in permuted space. PAP^T[i,j] = A[perm_fwd[gi], perm_fwd[gj]].
/// If scaling is present, also apply S: entry *= scaling[gi] * scaling[gj].
fn extract_submatrix(
    a: &SparseColMat<usize, f64>,
    global_indices: &[usize],
    perm_fwd: &[usize],
    scaling: Option<&[f64]>,
) -> faer::Mat<f64> {
    let m = global_indices.len();
    let n = a.nrows();
    let symbolic = a.symbolic();
    let col_ptrs = symbolic.col_ptr();
    let row_indices = symbolic.row_idx();
    let values = a.val();

    // Build reverse map: original index -> local position (if present)
    let mut orig_to_local = vec![usize::MAX; n];
    for (local, &gi) in global_indices.iter().enumerate() {
        let orig = perm_fwd[gi];
        orig_to_local[orig] = local;
    }

    let mut result = faer::Mat::<f64>::zeros(m, m);

    // For each original column that corresponds to a global index in our set:
    for (local_col, &gi_col) in global_indices.iter().enumerate() {
        let orig_col = perm_fwd[gi_col];
        let start = col_ptrs[orig_col];
        let end = col_ptrs[orig_col + 1];
        for idx in start..end {
            let orig_row = row_indices[idx];
            let local_row = orig_to_local[orig_row];
            if local_row == usize::MAX {
                continue;
            }
            let mut val = values[idx];
            if let Some(s) = scaling {
                // scaling is in elimination order (indexed by permuted col index)
                // gi_col and the permuted index of orig_row
                // Actually, we need the permuted index of orig_row to get the scaling.
                // perm_inv[orig_row] would give the permuted index.
                // But we don't have perm_inv here. We have perm_fwd.
                // For orig_row that maps to local_row, global_indices[local_row] is
                // the permuted index. So scaling[global_indices[local_row]] and
                // scaling[gi_col].
                val *= s[global_indices[local_row]] * s[gi_col];
            }
            // Place symmetrically
            result[(local_row, local_col)] += val;
            // Note: full symmetric CSC stores both triangles, so we get
            // both (i,j) and (j,i) from the CSC. We just accumulate.
            // But we need to be careful not to double-count diagonal.
            // Actually, since we're iterating CSC columns and each entry
            // (i,j) in column j is placed at (local_row, local_col),
            // if the matrix stores both triangles, the entry (j,i)
            // will be found when we iterate column i. So the above
            // is correct: just add each entry once.
        }
    }

    result
}

// ---------------------------------------------------------------------------
// Actual test cases
// ---------------------------------------------------------------------------

#[test]
#[ignore = "requires CI SuiteSparse matrices — run with --ignored --nocapture"]
fn diagnose_bratu3d() {
    use rivrs_sparse::io::registry;

    let test = registry::load_test_matrix("GHS_indef/bratu3d")
        .expect("registry error")
        .expect("bratu3d not found — is test-data/suitesparse-ci/hard-indefinite/bratu3d.mtx present?");
    let a = &test.matrix;

    // Test with MatchOrderMetis (recommended for hard-indefinite)
    diagnose_matrix("bratu3d", a, OrderingStrategy::MatchOrderMetis);
    diagnose_per_supernode_reconstruction("bratu3d", a, OrderingStrategy::MatchOrderMetis, 100);

    // Also test with plain METIS (known to be worse on this matrix)
    diagnose_matrix("bratu3d", a, OrderingStrategy::Metis);
    diagnose_per_supernode_reconstruction("bratu3d", a, OrderingStrategy::Metis, 100);
}

#[test]
#[ignore = "requires CI SuiteSparse matrices — run with --ignored --nocapture"]
fn diagnose_sparsine() {
    use rivrs_sparse::io::registry;

    let test = registry::load_test_matrix("GHS_indef/sparsine")
        .expect("registry error")
        .expect("sparsine not found — is test-data/suitesparse-ci/easy-indefinite/sparsine.mtx present?");
    let a = &test.matrix;

    // Test with plain METIS (recommended for easy-indefinite)
    diagnose_matrix("sparsine", a, OrderingStrategy::Metis);
    diagnose_per_supernode_reconstruction("sparsine", a, OrderingStrategy::Metis, 100);
}

#[test]
#[ignore = "requires CI SuiteSparse matrices — run with --ignored --nocapture"]
fn diagnose_cvxqp3() {
    use rivrs_sparse::io::registry;

    let test = registry::load_test_matrix("GHS_indef/cvxqp3")
        .expect("registry error")
        .expect("cvxqp3 not found");
    let a = &test.matrix;

    diagnose_matrix("cvxqp3", a, OrderingStrategy::MatchOrderMetis);
    diagnose_per_supernode_reconstruction("cvxqp3", a, OrderingStrategy::MatchOrderMetis, 100);
}

/// Quick test on a small hand-constructed matrix that should pass.
/// This validates the diagnostic machinery works correctly.
#[test]
fn diagnose_small_pd_sanity_check() {
    // 5x5 PD tridiagonal
    let entries: &[(usize, usize, f64)] = &[
        (0, 0, 4.0), (1, 0, 1.0),
        (1, 1, 4.0), (2, 1, 1.0),
        (2, 2, 4.0), (3, 2, 1.0),
        (3, 3, 4.0), (4, 3, 1.0),
        (4, 4, 4.0),
    ];
    let a = common::sparse_from_lower_triplets(5, entries);
    diagnose_matrix("small-pd-5", &a, OrderingStrategy::Metis);
}

/// Diagnose ALL CI matrices to get a complete picture.
#[test]
#[ignore = "requires CI SuiteSparse matrices — run with --ignored --nocapture --test-threads=1"]
fn diagnose_all_ci_matrices() {
    use rivrs_sparse::io::registry;

    let all_meta = registry::load_registry().expect("load registry");
    let ci_matrices: Vec<_> = all_meta.iter().filter(|m| m.ci_subset).collect();

    for meta in &ci_matrices {
        let test = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(tm)) => tm,
            Ok(None) => {
                eprintln!("SKIP: {} (not found on disk)", meta.name);
                continue;
            }
            Err(e) => {
                eprintln!("ERROR loading {}: {}", meta.name, e);
                continue;
            }
        };

        let ordering = if meta.category == "hard-indefinite" {
            OrderingStrategy::MatchOrderMetis
        } else {
            OrderingStrategy::Metis
        };

        diagnose_matrix(&meta.name, &test.matrix, ordering.clone());
        diagnose_per_supernode_reconstruction(
            &meta.name,
            &test.matrix,
            ordering,
            100,
        );
    }
}
