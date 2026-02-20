//! MC64 pipeline isolation: matching vs scaling vs ordering.
//!
//! Systematically isolates which part of the MC64→condensation→METIS pipeline
//! affects backward error for specific SuiteSparse matrices.
//!
//! Experiments:
//! 1. Baseline: Plain METIS vs MatchOrderMetis backward error
//! 2. MC64 quality diagnostics: matching stats, scaling distribution, cycle structure
//! 3. Scaling isolation: MatchOrderMetis ordering WITHOUT scaling
//! 4. Ordering isolation: plain METIS ordering WITH MC64 scaling
//!
//! Usage:
//!   cargo run --example mc64_isolation --release [matrix_names...]
//!
//! If no matrix names are given, runs on a default set of interesting matrices.
//! Matrix names use the SuiteSparse format: "GHS_indef/dawson5", "Newman/astro-ph", etc.

use faer::Col;
use faer::dyn_stack::{MemBuffer, MemStack};
use faer::sparse::linalg::cholesky::SymmetricOrdering;

use rivrs_sparse::aptp::factor::AptpOptions;
use rivrs_sparse::aptp::matching::count_cycles;
use rivrs_sparse::aptp::numeric::AptpNumeric;
use rivrs_sparse::aptp::solve::{aptp_solve, aptp_solve_scratch};
use rivrs_sparse::aptp::{
    AnalyzeOptions, AptpSymbolic, FactorOptions, Mc64Job, OrderingStrategy, SparseLDLT,
    match_order_metis, mc64_matching, metis_ordering,
};
use rivrs_sparse::io::registry;
use rivrs_sparse::validate::sparse_backward_error;

/// Sparse matrix-vector product using CSC storage.
fn sparse_matvec(a: &faer::sparse::SparseColMat<usize, f64>, x: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let sym = a.symbolic();
    let cp = sym.col_ptr();
    let ri = sym.row_idx();
    let vals = a.val();
    let mut result = vec![0.0f64; n];
    for j in 0..n {
        for idx in cp[j]..cp[j + 1] {
            result[ri[idx]] += vals[idx] * x[j];
        }
    }
    result
}

/// Print a separator line.
fn sep() {
    println!("{}", "-".repeat(100));
}

/// Experiment 1: Baseline backward error comparison.
fn experiment_baseline(_name: &str, a: &faer::sparse::SparseColMat<usize, f64>, b: &Col<f64>) {
    println!("\n=== EXPERIMENT 1: Baseline Backward Error ===");
    println!(
        "{:<15} {:<18} {:<12} {:<8} {:<8} {:<8} {:<8}",
        "Ordering", "BE", "Status", "1x1", "2x2", "delays", "max_frt"
    );
    sep();

    for (ordering_name, ordering) in [
        ("Metis", OrderingStrategy::Metis),
        ("MatchOrderMetis", OrderingStrategy::MatchOrderMetis),
    ] {
        let opts = AnalyzeOptions { ordering };
        let mut solver = match SparseLDLT::analyze_with_matrix(a, &opts) {
            Ok(s) => s,
            Err(e) => {
                println!("{:<15} analyze error: {}", ordering_name, e);
                continue;
            }
        };

        let factor_opts = FactorOptions::default();
        if let Err(e) = solver.factor(a, &factor_opts) {
            println!("{:<15} factor error: {}", ordering_name, e);
            continue;
        }

        let scratch = solver.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);
        match solver.solve(b, stack) {
            Ok(x) => {
                let be = sparse_backward_error(a, &x, b);
                let stats = solver.stats().unwrap();
                let status = if be < 5e-11 { "PASS" } else { "FAIL" };
                println!(
                    "{:<15} {:<18.2e} {:<12} {:<8} {:<8} {:<8} {:<8}",
                    ordering_name,
                    be,
                    status,
                    stats.total_1x1_pivots,
                    stats.total_2x2_pivots,
                    stats.total_delayed,
                    stats.max_front_size,
                );
            }
            Err(e) => {
                println!("{:<15} solve error: {}", ordering_name, e);
            }
        }
    }
}

/// Experiment 2: MC64 quality diagnostics.
fn experiment_mc64_diagnostics(_name: &str, a: &faer::sparse::SparseColMat<usize, f64>) {
    println!("\n=== EXPERIMENT 2: MC64 Quality Diagnostics ===");

    let mc64_result = match mc64_matching(a, Mc64Job::MaximumProduct) {
        Ok(r) => r,
        Err(e) => {
            println!("  MC64 failed: {}", e);
            return;
        }
    };

    let n = a.nrows();
    let (fwd, _) = mc64_result.matching.as_ref().arrays();

    // Matching cardinality
    println!(
        "  Matched: {} / {} ({:.1}%)",
        mc64_result.matched,
        n,
        100.0 * mc64_result.matched as f64 / n as f64
    );
    let unmatched = mc64_result.is_matched.iter().filter(|&&m| !m).count();
    println!("  Unmatched rows: {}", unmatched);

    // Cycle structure
    let (singletons, two_cycles, longer_cycles) = count_cycles(fwd);
    println!(
        "  Cycle decomposition: {} singletons, {} two-cycles, {} longer-cycles",
        singletons, two_cycles, longer_cycles
    );

    // Matching weight: sum of |a_{i,match[i]}|
    let sym = a.symbolic();
    let cp = sym.col_ptr();
    let ri = sym.row_idx();
    let vals = a.val();

    let mut matching_weight = 0.0f64;
    let mut matched_diag = 0usize;
    let mut matched_offdiag = 0usize;
    for i in 0..n {
        let j = fwd[i];
        if i == j {
            matched_diag += 1;
        } else {
            matched_offdiag += 1;
        }
        // Find |a_{i,j}| in the CSC
        let mut found = false;
        for idx in cp[j]..cp[j + 1] {
            if ri[idx] == i {
                matching_weight += vals[idx].abs();
                found = true;
                break;
            }
        }
        // Try symmetric (j,i) if not found directly
        if !found && i != j {
            for idx in cp[i]..cp[i + 1] {
                if ri[idx] == j {
                    matching_weight += vals[idx].abs();
                    break;
                }
            }
        }
    }
    println!(
        "  Matching weight (sum |a_{{i,sigma(i)}}|): {:.6e}",
        matching_weight
    );
    println!(
        "  Matched diagonal: {}, off-diagonal: {}",
        matched_diag, matched_offdiag
    );

    // Scaling statistics
    let scaling = &mc64_result.scaling;
    let min_s = scaling.iter().copied().fold(f64::INFINITY, f64::min);
    let max_s = scaling.iter().copied().fold(0.0f64, f64::max);
    let mean_s: f64 = scaling.iter().sum::<f64>() / n as f64;
    let log_scaling: Vec<f64> = scaling.iter().map(|&s| s.ln()).collect();
    let mean_log = log_scaling.iter().sum::<f64>() / n as f64;
    let std_log = (log_scaling
        .iter()
        .map(|&l| (l - mean_log).powi(2))
        .sum::<f64>()
        / n as f64)
        .sqrt();

    println!("  Scaling factors:");
    println!(
        "    min={:.6e}  max={:.6e}  ratio={:.6e}",
        min_s,
        max_s,
        max_s / min_s
    );
    println!(
        "    mean={:.6e}  geometric_mean={:.6e}",
        mean_s,
        mean_log.exp()
    );
    println!("    log-scale: mean={:.4}  std={:.4}", mean_log, std_log);

    // Extreme scaling factors
    let extreme_large: Vec<(usize, f64)> = scaling
        .iter()
        .enumerate()
        .filter(|&(_, s)| *s > 1e10)
        .map(|(i, s)| (i, *s))
        .collect();
    let extreme_small: Vec<(usize, f64)> = scaling
        .iter()
        .enumerate()
        .filter(|&(_, s)| *s < 1e-10)
        .map(|(i, s)| (i, *s))
        .collect();
    println!("    Extreme large (>1e10): {} entries", extreme_large.len());
    if !extreme_large.is_empty() {
        for &(i, s) in extreme_large.iter().take(5) {
            println!("      scaling[{}] = {:.6e}", i, s);
        }
        if extreme_large.len() > 5 {
            println!("      ... and {} more", extreme_large.len() - 5);
        }
    }
    println!(
        "    Extreme small (<1e-10): {} entries",
        extreme_small.len()
    );
    if !extreme_small.is_empty() {
        for &(i, s) in extreme_small.iter().take(5) {
            println!("      scaling[{}] = {:.6e}", i, s);
        }
        if extreme_small.len() > 5 {
            println!("      ... and {} more", extreme_small.len() - 5);
        }
    }

    // Scaling bound check: max |s_i * a_ij * s_j| over all entries
    let mut global_max_scaled = 0.0f64;
    let mut violations = 0usize;
    for j in 0..n {
        for idx in cp[j]..cp[j + 1] {
            let i = ri[idx];
            let scaled = (scaling[i] * vals[idx] * scaling[j]).abs();
            if scaled > global_max_scaled {
                global_max_scaled = scaled;
            }
            if scaled > 1.0 + 1e-8 {
                violations += 1;
            }
        }
    }
    println!(
        "  Scaling bound: max |s_i * a_ij * s_j| = {:.6e}",
        global_max_scaled
    );
    println!("  Scaling violations (>1+1e-8): {}", violations);

    // Scaled diagonal dominance check
    let mut diag_dominant = 0usize;
    let mut diag_weak = 0usize;
    let mut diag_missing = 0usize;
    for j in 0..n {
        let mut diag_val = 0.0f64;
        let mut offdiag_max = 0.0f64;
        let mut has_diag = false;
        for idx in cp[j]..cp[j + 1] {
            let i = ri[idx];
            let scaled = (scaling[i] * vals[idx] * scaling[j]).abs();
            if i == j {
                diag_val = scaled;
                has_diag = true;
            } else {
                if scaled > offdiag_max {
                    offdiag_max = scaled;
                }
            }
        }
        if !has_diag {
            diag_missing += 1;
        } else if diag_val >= offdiag_max {
            diag_dominant += 1;
        } else {
            diag_weak += 1;
        }
    }
    println!(
        "  Scaled diagonal: {} dominant, {} weak, {} missing",
        diag_dominant, diag_weak, diag_missing
    );

    // Condensation diagnostics from match_order_metis
    let mor_result = match match_order_metis(a) {
        Ok(r) => r,
        Err(e) => {
            println!("  match_order_metis failed: {}", e);
            return;
        }
    };
    println!(
        "  Condensation: dim={} (from {}), singletons={}, two_cycles={}",
        mor_result.condensed_dim, n, mor_result.singletons, mor_result.two_cycles
    );
    println!(
        "  Condensation ratio: {:.1}%",
        100.0 * mor_result.condensed_dim as f64 / n as f64
    );
}

/// Experiment 3: Scaling isolation.
/// Tests: (a) MatchOrderMetis ordering with scaling (normal),
///        (b) MatchOrderMetis ordering WITHOUT scaling,
///        (c) Plain METIS ordering WITH MC64 scaling.
fn experiment_scaling_isolation(
    _name: &str,
    a: &faer::sparse::SparseColMat<usize, f64>,
    b: &Col<f64>,
) {
    println!("\n=== EXPERIMENT 3: Scaling Isolation ===");
    println!(
        "{:<40} {:<18} {:<12} {:<8} {:<8}",
        "Configuration", "BE", "Status", "delays", "max_frt"
    );
    sep();

    let n = a.nrows();

    // 3a: MatchOrderMetis with scaling (baseline)
    {
        let opts = AnalyzeOptions {
            ordering: OrderingStrategy::MatchOrderMetis,
        };
        let mut solver = match SparseLDLT::analyze_with_matrix(a, &opts) {
            Ok(s) => s,
            Err(e) => {
                println!("{:<40} analyze error: {}", "MatchOrderMetis (ord+scale)", e);
                return;
            }
        };
        let factor_opts = FactorOptions::default();
        solver.factor(a, &factor_opts).unwrap();
        let scratch = solver.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);
        let x = solver.solve(b, stack).unwrap();
        let be = sparse_backward_error(a, &x, b);
        let stats = solver.stats().unwrap();
        let status = if be < 5e-11 { "PASS" } else { "FAIL" };
        println!(
            "{:<40} {:<18.2e} {:<12} {:<8} {:<8}",
            "MatchOrderMetis (ord+scale)", be, status, stats.total_delayed, stats.max_front_size
        );
    }

    // 3b: MatchOrderMetis ordering WITHOUT scaling
    // Use match_order_metis to get ordering, then pass as UserSupplied (no scaling)
    {
        let mor_result = match match_order_metis(a) {
            Ok(r) => r,
            Err(e) => {
                println!(
                    "{:<40} match_order_metis error: {}",
                    "MatchOrder ord, NO scale", e
                );
                return;
            }
        };
        let perm = mor_result.ordering;

        let opts = AnalyzeOptions {
            ordering: OrderingStrategy::UserSupplied(perm),
        };
        let mut solver = SparseLDLT::analyze_with_matrix(a, &opts).unwrap();
        solver.factor(a, &FactorOptions::default()).unwrap();
        let scratch = solver.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);
        let x = solver.solve(b, stack).unwrap();
        let be = sparse_backward_error(a, &x, b);
        let stats = solver.stats().unwrap();
        let status = if be < 5e-11 { "PASS" } else { "FAIL" };
        println!(
            "{:<40} {:<18.2e} {:<12} {:<8} {:<8}",
            "MatchOrder ord, NO scale", be, status, stats.total_delayed, stats.max_front_size
        );
    }

    // 3c: Plain METIS ordering WITH MC64 scaling
    // Use metis_ordering for the ordering, but manually inject MC64 scaling
    // This requires building SparseLDLT manually at the lower level.
    {
        let mc64_result = match mc64_matching(a, Mc64Job::MaximumProduct) {
            Ok(r) => r,
            Err(e) => {
                println!("{:<40} mc64 error: {}", "Metis ord + MC64 scale", e);
                return;
            }
        };
        let scaling = mc64_result.scaling;

        let perm = match metis_ordering(a.symbolic()) {
            Ok(p) => p,
            Err(e) => {
                println!("{:<40} metis error: {}", "Metis ord + MC64 scale", e);
                return;
            }
        };

        let symbolic =
            match AptpSymbolic::analyze(a.symbolic(), SymmetricOrdering::Custom(perm.as_ref())) {
                Ok(s) => s,
                Err(e) => {
                    println!("{:<40} symbolic error: {}", "Metis ord + MC64 scale", e);
                    return;
                }
            };

        // Transform scaling to elimination order
        let (perm_fwd, _) = symbolic.perm_vecs();
        let elim_scaling: Vec<f64> = (0..n).map(|i| scaling[perm_fwd[i]]).collect();

        let aptp_options = AptpOptions::default();
        let numeric = match AptpNumeric::factor(&symbolic, a, &aptp_options, Some(&elim_scaling)) {
            Ok(num) => num,
            Err(e) => {
                println!("{:<40} factor error: {}", "Metis ord + MC64 scale", e);
                return;
            }
        };

        // Solve manually with scaling
        let mut rhs_perm = vec![0.0f64; n];
        for new in 0..n {
            rhs_perm[new] = b[perm_fwd[new]];
        }
        // Apply scaling
        for i in 0..n {
            rhs_perm[i] *= elim_scaling[i];
        }

        let scratch_req = aptp_solve_scratch(&numeric, 1);
        let mut mem = MemBuffer::new(scratch_req);
        let stack = MemStack::new(&mut mem);
        aptp_solve(&symbolic, &numeric, &mut rhs_perm, stack).unwrap();

        // Unscale
        for i in 0..n {
            rhs_perm[i] *= elim_scaling[i];
        }

        // Unpermute
        let mut x = Col::zeros(n);
        for new in 0..n {
            x[perm_fwd[new]] = rhs_perm[new];
        }

        let be = sparse_backward_error(a, &x, b);
        let stats = numeric.stats();
        let status = if be < 5e-11 { "PASS" } else { "FAIL" };
        println!(
            "{:<40} {:<18.2e} {:<12} {:<8} {:<8}",
            "Metis ord + MC64 scale", be, status, stats.total_delayed, stats.max_front_size
        );
    }

    // 3d: Plain Metis, no scaling (reference)
    {
        let opts = AnalyzeOptions {
            ordering: OrderingStrategy::Metis,
        };
        let mut solver = SparseLDLT::analyze_with_matrix(a, &opts).unwrap();
        solver.factor(a, &FactorOptions::default()).unwrap();
        let scratch = solver.solve_scratch(1);
        let mut mem = MemBuffer::new(scratch);
        let stack = MemStack::new(&mut mem);
        let x = solver.solve(b, stack).unwrap();
        let be = sparse_backward_error(a, &x, b);
        let stats = solver.stats().unwrap();
        let status = if be < 5e-11 { "PASS" } else { "FAIL" };
        println!(
            "{:<40} {:<18.2e} {:<12} {:<8} {:<8}",
            "Metis ord, NO scale (reference)",
            be,
            status,
            stats.total_delayed,
            stats.max_front_size
        );
    }
}

fn main() {
    let all = registry::load_registry().expect("registry");

    let cli_args: Vec<String> = std::env::args().skip(1).collect();
    let default_targets = vec![
        "GHS_indef/dawson5".to_string(),
        "Newman/astro-ph".to_string(),
        "GHS_indef/copter2".to_string(),
        "GHS_indef/helm3d01".to_string(),
        "GHS_indef/sparsine".to_string(),
        "TSOPF/TSOPF_FS_b162_c1".to_string(),
        "TSOPF/TSOPF_FS_b39_c7".to_string(),
        "Cote/vibrobox".to_string(),
        "Boeing/crystk02".to_string(),
    ];
    let targets = if cli_args.is_empty() {
        &default_targets
    } else {
        &cli_args
    };

    for target in targets {
        let meta = all.iter().find(|m| m.name == *target);
        let meta = match meta {
            Some(m) => m,
            None => {
                println!("\n{}: not found in registry", target);
                continue;
            }
        };

        let test = match registry::load_test_matrix_from_entry(meta) {
            Ok(Some(t)) => t,
            Ok(None) => {
                println!("\n{}: not found on disk (gitignored?)", target);
                continue;
            }
            Err(e) => {
                println!("\n{}: load error: {}", target, e);
                continue;
            }
        };

        let a = &test.matrix;
        let n = a.nrows();
        let nnz = a.compute_nnz();

        println!("\n{}", "=".repeat(100));
        println!("MATRIX: {} (n={}, nnz={})", target, n, nnz);
        println!("{}", "=".repeat(100));

        // Build RHS
        let x_exact: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
        let b_vec = sparse_matvec(a, &x_exact);
        let b = Col::from_fn(n, |i| b_vec[i]);

        experiment_baseline(target, a, &b);
        experiment_mc64_diagnostics(target, a);
        experiment_scaling_isolation(target, a, &b);
    }

    println!("\n{}", "=".repeat(100));
    println!("INVESTIGATION COMPLETE");
    println!("{}", "=".repeat(100));
}
