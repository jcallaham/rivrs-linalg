//! Isolate whether MC64 scaling or ordering causes backward error degradation.
//!
//! For linverse: compare MatchOrderMetis (ordering+scaling) vs the same ordering
//! WITHOUT scaling, to determine if the ordering or the scaling is the problem.
//!
//! Usage: cargo run --example scaling_diagnosis --release

use faer::Col;
use rivrs_sparse::aptp::{
    match_order_metis, AnalyzeOptions, FactorOptions, OrderingStrategy, SparseLDLT,
};
use rivrs_sparse::io::registry;
use rivrs_sparse::validate::sparse_backward_error;

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

fn main() {
    let all = registry::load_registry().expect("registry");

    let targets = ["GHS_indef/linverse"];

    for target in &targets {
        let meta = all.iter().find(|m| m.name == *target).unwrap();
        let test = registry::load_test_matrix_from_entry(meta)
            .unwrap()
            .unwrap();
        let a = &test.matrix;
        let n = a.nrows();

        let x_exact: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
        let b_vec = sparse_matvec(a, &x_exact);
        let b = Col::from_fn(n, |i| b_vec[i]);

        println!("\n{} (n={})", target, n);

        // 1. MatchOrderMetis (ordering + scaling)
        {
            let opts = AnalyzeOptions::default(); // MatchOrderMetis
            let mut solver = SparseLDLT::analyze_with_matrix(a, &opts).unwrap();
            solver.factor(a, &FactorOptions::default()).unwrap();
            let scratch = solver.solve_scratch(1);
            let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
            let stack = faer::dyn_stack::MemStack::new(&mut mem);
            let x = solver.solve(&b, stack).unwrap();
            let be = sparse_backward_error(a, &x, &b);
            let stats = solver.stats().unwrap();
            println!(
                "  MatchOrderMetis (ordering+scaling): BE={:.2e}  delays={} max_front={}",
                be, stats.total_delayed, stats.max_front_size
            );
        }

        // 2. Same ordering, NO scaling (UserSupplied with match_order_metis perm)
        {
            let result = match_order_metis(a).unwrap();
            let perm = result.ordering;

            // Print scaling stats
            let min_s = result
                .scaling
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min);
            let max_s = result
                .scaling
                .iter()
                .copied()
                .fold(0.0f64, f64::max);
            let mean_s: f64 = result.scaling.iter().sum::<f64>() / n as f64;
            println!(
                "  MC64 scaling: min={:.4e} max={:.4e} mean={:.4e} ratio={:.4e}",
                min_s,
                max_s,
                mean_s,
                max_s / min_s
            );
            println!(
                "  MC64 stats: matched={} condensed_dim={} two_cycles={}",
                result.matched, result.condensed_dim, result.two_cycles
            );

            let opts = AnalyzeOptions {
                ordering: OrderingStrategy::UserSupplied(perm),
            };
            let mut solver = SparseLDLT::analyze_with_matrix(a, &opts).unwrap();
            solver.factor(a, &FactorOptions::default()).unwrap();
            let scratch = solver.solve_scratch(1);
            let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
            let stack = faer::dyn_stack::MemStack::new(&mut mem);
            let x = solver.solve(&b, stack).unwrap();
            let be = sparse_backward_error(a, &x, &b);
            let stats = solver.stats().unwrap();
            println!(
                "  Same ordering, NO scaling:          BE={:.2e}  delays={} max_front={}",
                be, stats.total_delayed, stats.max_front_size
            );
        }

        // 3. Plain Metis (different ordering, no scaling)
        {
            let opts = AnalyzeOptions {
                ordering: OrderingStrategy::Metis,
            };
            let mut solver = SparseLDLT::analyze_with_matrix(a, &opts).unwrap();
            solver.factor(a, &FactorOptions::default()).unwrap();
            let scratch = solver.solve_scratch(1);
            let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
            let stack = faer::dyn_stack::MemStack::new(&mut mem);
            let x = solver.solve(&b, stack).unwrap();
            let be = sparse_backward_error(a, &x, &b);
            let stats = solver.stats().unwrap();
            println!(
                "  Plain Metis (no scaling):            BE={:.2e}  delays={} max_front={}",
                be, stats.total_delayed, stats.max_front_size
            );
        }
    }
}
