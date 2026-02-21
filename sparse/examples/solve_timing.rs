//! End-to-end solve timing and accuracy for CI SuiteSparse matrices.
//!
//! Uses MatchOrderMetis ordering (default). Reports factor/solve timing,
//! backward error, and factorization statistics.
//!
//! Usage: cargo run --example solve_timing --release

use faer::Col;
use faer::Par;
use rivrs_sparse::aptp::SparseLDLT;
use rivrs_sparse::io::registry;
use rivrs_sparse::validate::sparse_backward_error;
use std::time::Instant;

fn main() {
    let all = registry::load_registry().expect("registry");
    let ci: Vec<_> = all.iter().filter(|m| m.ci_subset).collect();

    println!(
        "{:<30} {:>8} {:>10} {:>8} {:>8} {:>12} {:>7} {:>7} {:>7} {:>6}",
        "Matrix", "n", "nnz", "fac_ms", "slv_ms", "backward_err", "1x1", "2x2", "delay", "status"
    );
    println!("{}", "-".repeat(120));

    let mut pass = 0;
    let mut fail = 0;

    for meta in &ci {
        let test = match registry::load_test_matrix(&meta.name) {
            Ok(Some(t)) => t,
            _ => {
                println!("{:<30} SKIP (not found)", meta.name);
                continue;
            }
        };
        let a = &test.matrix;
        let n = a.nrows();

        // Create RHS: b = A * x_true
        let x_true: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
        let sym = a.symbolic();
        let cp = sym.col_ptr();
        let ri = sym.row_idx();
        let vals = a.val();
        let mut b_vec = vec![0.0f64; n];
        for j in 0..n {
            for idx in cp[j]..cp[j + 1] {
                b_vec[ri[idx]] += vals[idx] * x_true[j];
            }
        }
        let b = Col::from_fn(n, |i| b_vec[i]);

        // Analyze + factor (MatchOrderMetis default)
        let t0 = Instant::now();
        let opts = rivrs_sparse::aptp::AnalyzeOptions::default();
        let mut solver = match SparseLDLT::analyze_with_matrix(a, &opts) {
            Ok(s) => s,
            Err(e) => {
                println!(
                    "{:<30} {:>8} {:>10} ANALYZE FAILED: {}",
                    meta.name, n, meta.nnz, e
                );
                fail += 1;
                continue;
            }
        };
        let factor_opts = rivrs_sparse::aptp::FactorOptions::default();
        match solver.factor(a, &factor_opts) {
            Ok(()) => {}
            Err(e) => {
                println!(
                    "{:<30} {:>8} {:>10} FACTOR FAILED: {}",
                    meta.name, n, meta.nnz, e
                );
                fail += 1;
                continue;
            }
        };
        let factor_ms = t0.elapsed().as_millis();

        // Solve
        let t1 = Instant::now();
        let scratch = solver.solve_scratch(1);
        let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
        let stack = faer::dyn_stack::MemStack::new(&mut mem);
        let x = match solver.solve(&b, stack, Par::Seq) {
            Ok(x) => x,
            Err(e) => {
                println!(
                    "{:<30} {:>8} {:>10} {:>8} SOLVE FAILED: {}",
                    meta.name, n, meta.nnz, factor_ms, e
                );
                fail += 1;
                continue;
            }
        };
        let solve_ms = t1.elapsed().as_millis();

        let be = sparse_backward_error(a, &x, &b);
        let stats = solver.stats().unwrap();
        let status = if be < 5e-11 {
            "PASS"
        } else if be < 1e-6 {
            "WARN"
        } else {
            "FAIL"
        };
        if be < 5e-11 {
            pass += 1;
        } else {
            fail += 1;
        }

        println!(
            "{:<30} {:>8} {:>10} {:>8} {:>8} {:>12.2e} {:>7} {:>7} {:>7} {:>6}",
            meta.name,
            n,
            meta.nnz,
            factor_ms,
            solve_ms,
            be,
            stats.total_1x1_pivots,
            stats.total_2x2_pivots,
            stats.total_delayed,
            status
        );
    }

    println!("\n{}/{} passed (threshold: 5e-11)", pass, pass + fail);
}
