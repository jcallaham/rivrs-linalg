//! Random matrix accuracy benchmark for the APTP solver.
//!
//! Generates random symmetric indefinite matrices at various sizes and densities,
//! solves Ax = b with a known x, and checks backward error. This isolates accuracy
//! issues from sparsity-pattern or ordering effects.
//!
//! Also includes bratu3d from the CI subset as an intermediate real-world check.
//!
//! Usage: cargo run --example accuracy_benchmark --release

use faer::Col;
use faer::Par;
use faer::sparse::{SparseColMat, Triplet};
use rivrs_sparse::aptp::SparseLDLT;
use rivrs_sparse::io::registry;
use rivrs_sparse::validate::sparse_backward_error;
use std::collections::HashSet;

/// Simple LCG RNG for reproducibility without needing the rand crate in examples.
struct Lcg(u64);

impl Lcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next_f64(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 33) as f64) / (u32::MAX as f64) * 2.0 - 1.0 // [-1, 1]
    }
    fn next_usize(&mut self, bound: usize) -> usize {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 33) as usize) % bound
    }
}

/// Generate a random symmetric indefinite sparse matrix.
/// Half the diagonal entries are positive, half negative (diagonally dominant).
fn random_symmetric_indefinite(n: usize, density: f64, seed: u64) -> SparseColMat<usize, f64> {
    let mut rng = Lcg::new(seed);

    let target_offdiag = ((n as f64 * (n as f64 - 1.0) / 2.0) * density) as usize;
    let target_offdiag = target_offdiag.max(n); // at least n off-diagonal pairs

    let mut triplets: Vec<Triplet<usize, usize, f64>> = Vec::new();
    let mut placed = HashSet::new();
    let mut row_abs_sum = vec![0.0f64; n];

    // Place off-diagonal entries
    for _ in 0..target_offdiag * 4 {
        if placed.len() >= target_offdiag {
            break;
        }
        let i = rng.next_usize(n);
        let j = rng.next_usize(n);
        if i == j {
            continue;
        }
        let (lo, hi) = if i < j { (i, j) } else { (j, i) };
        if placed.contains(&(lo, hi)) {
            continue;
        }
        placed.insert((lo, hi));
        let v = rng.next_f64();
        triplets.push(Triplet::new(lo, hi, v));
        triplets.push(Triplet::new(hi, lo, v));
        row_abs_sum[lo] += v.abs();
        row_abs_sum[hi] += v.abs();
    }

    // Diagonal: first half positive, second half negative (diagonally dominant)
    let half = n / 2;
    for (i, &abs_sum) in row_abs_sum.iter().enumerate() {
        let margin = 1.0 + rng.next_f64().abs();
        let d = if i < half {
            abs_sum + margin
        } else {
            -(abs_sum + margin)
        };
        triplets.push(Triplet::new(i, i, d));
    }

    SparseColMat::try_new_from_triplets(n, n, &triplets).expect("triplet creation failed")
}

/// Compute b = A * x_true using CSC access.
fn matvec(a: &SparseColMat<usize, f64>, x: &[f64]) -> Vec<f64> {
    let n = a.nrows();
    let sym = a.symbolic();
    let col_ptrs = sym.col_ptr();
    let row_indices = sym.row_idx();
    let values = a.val();
    let mut b = vec![0.0f64; n];
    for j in 0..n {
        for idx in col_ptrs[j]..col_ptrs[j + 1] {
            let i = row_indices[idx];
            b[i] += values[idx] * x[j];
        }
    }
    b
}

fn run_solve(
    label: &str,
    a: &SparseColMat<usize, f64>,
    threshold: f64,
    ordering: rivrs_sparse::aptp::OrderingStrategy,
) {
    let n = a.nrows();
    let nnz = a.compute_nnz();

    // Known solution
    let x_true: Vec<f64> = (0..n).map(|i| ((i % 7) as f64 - 3.0) / 3.0).collect();
    let b_vec = matvec(a, &x_true);
    let b = Col::from_fn(n, |i| b_vec[i]);

    // Analyze + factor
    let opts = rivrs_sparse::aptp::AnalyzeOptions {
        ordering,
        ..Default::default()
    };
    let mut solver = match SparseLDLT::analyze_with_matrix(a, &opts) {
        Ok(s) => s,
        Err(e) => {
            println!("{:<40} {:>6} {:>8} ANALYZE FAILED: {}", label, n, nnz, e);
            return;
        }
    };

    let factor_opts = rivrs_sparse::aptp::FactorOptions {
        threshold,
        ..Default::default()
    };
    let t0 = std::time::Instant::now();
    if let Err(e) = solver.factor(a, &factor_opts) {
        println!("{:<40} {:>6} {:>8} FACTOR FAILED: {}", label, n, nnz, e);
        return;
    }
    let factor_ms = t0.elapsed().as_millis();

    // Solve
    let scratch = solver.solve_scratch(1);
    let mut mem = faer::dyn_stack::MemBuffer::new(scratch);
    let stack = faer::dyn_stack::MemStack::new(&mut mem);
    let x = match solver.solve(&b, stack, Par::Seq) {
        Ok(x) => x,
        Err(e) => {
            println!(
                "{:<40} {:>6} {:>8} {:>8} SOLVE FAILED: {}",
                label, n, nnz, factor_ms, e
            );
            return;
        }
    };

    let be = sparse_backward_error(a, &x, &b);
    let inertia = solver
        .inertia()
        .map(|i| format!("({},{},{})", i.positive, i.negative, i.zero));
    let status = if be < 1e-10 {
        "PASS"
    } else if be < 1e-6 {
        "WARN"
    } else {
        "FAIL"
    };

    println!(
        "{:<40} {:>6} {:>8} {:>8} {:>12.2e}  {}  {}",
        label,
        n,
        nnz,
        factor_ms,
        be,
        status,
        inertia.unwrap_or_default()
    );
}

fn main() {
    println!(
        "{:<40} {:>6} {:>8} {:>8} {:>12}  {:6}  {:7}",
        "Test", "n", "nnz", "fac_ms", "bwd_err", "status", "inertia"
    );
    println!("{}", "-".repeat(105));

    use rivrs_sparse::aptp::OrderingStrategy;

    // ===== Part 1: Random symmetric indefinite matrices (MatchOrderMetis) =====
    println!(
        "\n--- Random symmetric indefinite (MatchOrderMetis, density=0.05, threshold=0.01) ---"
    );
    for &n in &[50, 100, 200, 500, 1000] {
        for seed in 0..3u64 {
            let label = format!("random n={} seed={}", n, seed);
            let a = random_symmetric_indefinite(n, 0.05, seed * 1000 + n as u64);
            run_solve(&label, &a, 0.01, OrderingStrategy::MatchOrderMetis);
        }
    }

    // ===== Part 2: Varying density (MatchOrderMetis) =====
    println!("\n--- Random n=200, varying density (MatchOrderMetis, threshold=0.01) ---");
    for &density in &[0.01, 0.05, 0.10, 0.20, 0.50] {
        for seed in 0..3u64 {
            let label = format!("random n=200 dens={:.2} seed={}", density, seed);
            let a =
                random_symmetric_indefinite(200, density, seed * 100 + (density * 100.0) as u64);
            run_solve(&label, &a, 0.01, OrderingStrategy::MatchOrderMetis);
        }
    }

    // ===== Part 3: bratu3d comparison =====
    println!("\n--- SuiteSparse: bratu3d — METIS vs MatchOrderMetis ---");
    if let Ok(all) = registry::load_registry() {
        let bratu = all.iter().find(|m| m.name.contains("bratu3d"));
        if let Some(meta) = bratu {
            if let Ok(Some(test)) = registry::load_test_matrix(&meta.name) {
                run_solve("bratu3d METIS", &test.matrix, 0.01, OrderingStrategy::Metis);
                run_solve(
                    "bratu3d MatchOrderMetis",
                    &test.matrix,
                    0.01,
                    OrderingStrategy::MatchOrderMetis,
                );
            } else {
                println!("bratu3d: not found in test-data");
            }
        }
    }
}
