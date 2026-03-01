# rivrs-sparse

Sparse linear algebra functionality for Rust.

This library builds on [`faer`](https://crates.io/crates/faer) for foundational data types, dense solvers, and some sparse linear algebra routines.
Algorithms are implemented from academic literature and permissively-licensed reference code.

This library is written as numerical implementations for the [`rivrs`](https://crates.io/crates/rivrs) symbolic-numeric framework but functions as a standalone crate for scientific computing applications.
It is also re-exported by the [`rivrs-linalg`](https://crates.io/crates/rivrs-linalg) crate, containing more general numerical linear algebra functionality.

**AI Disclaimer**: Large amounts of the codebase were written using Claude Code.
Every effort was made to adhere to reference literature and codebases, to use structured development processes to maintain quality, and to extensively test the generated code using comparable test suites and ground truth data (e.g. [SuiteSparse](https://sparse.tamu.edu/) matrices).
However, as with any numerical code that hasn't existed for very long, it is recommended that you add your own correctness checks (e.g. backwards residual calculations) when using these solvers.
**Please report any issues**

## Features

### Sparse Symmetric Indefinite Direct Solver

Based on [SPRAL](https://github.com/ralna/spral). Implemented using [faer](https://github.com/sarah-ek/faer-rs) for dense linear algebra.

Use this for KKT matrices, saddle-point problems, and optimization formulations where the diagonal has negative or zero entries.
Compared to faer's built-in sparse LDL^T, it handles difficult indefinite matrices via APTP pivoting, MC64 matching/scaling, and METIS ordering.
Compared to SPRAL and MUMPS, it is native Rust, 100% safe code (except for calls to METIS), and competitive on sequential and parallel performance.

- Multifrontal LDL^T factorization with A Posteriori Threshold Pivoting (APTP)
- Two-level pivoting: TPP for small fronts, complete pivoting with BLAS-3 blocking for large fronts
- MC64 weighted bipartite matching & scaling + METIS nested dissection ordering
- Parallel factorization & solve via rayon (tree-level) and faer `Par` (intra-node BLAS)
- Three-stage API: analyze → factor → solve (symbolic reuse across re-factorizations)

#### Performance

Benchmarked on 65 SuiteSparse matrices (factorization time, backward error < 1e-12 on all, typically ~1e-16):

| Configuration | vs SPRAL | vs MUMPS |
|---------------|----------|----------|
| Sequential (1 thread) | median 5% speedup (36% faster, 38% comparable, 24% slower) | median 0.46x (~2x faster) |
| Parallel (8 threads) | median 10% speedup (53% faster, 23% comparable, 23% slower) | — |

See [`comparisons/README.md`](comparisons/README.md) for full results per matrix.

## API

```rust
use faer::sparse::{SparseColMat, Triplet};
use faer::Col;
use rivrs_sparse::symmetric::{SparseLDLT, SolverOptions};

// Build a symmetric matrix (full storage, both triangles)
let triplets = vec![
    Triplet::new(0, 0, 4.0),
    Triplet::new(1, 0, 1.0), Triplet::new(0, 1, 1.0),
    Triplet::new(1, 1, 3.0),
];
let a = SparseColMat::try_new_from_triplets(2, 2, &triplets).unwrap();
let b = Col::from_fn(2, |i| [5.0, 4.0][i]);

// One-shot solve
let x = SparseLDLT::solve_full(&a, &b, &SolverOptions::default()).unwrap();
```

For the three-stage API with parallel factorization:

```rust
use faer::{Par, dyn_stack::{MemBuffer, MemStack}};
use rivrs_sparse::symmetric::{AnalyzeOptions, FactorOptions, SparseLDLT};

let mut solver = SparseLDLT::analyze_with_matrix(&a, &AnalyzeOptions::default()).unwrap();
solver.factor(&a, &FactorOptions { par: Par::rayon(4), ..Default::default() }).unwrap();

let mut mem = MemBuffer::new(solver.solve_scratch(1));
let x = solver.solve(&b, MemStack::new(&mut mem), Par::rayon(4)).unwrap();
```

## Building & Testing

```bash
cargo build                    # Standard build (opt-level=3 in dev profile)
cargo build --release          # Optimized build (needed for large matrices)

# Unit tests (~360 tests) + CI subset (10 SuiteSparse matrices)
cargo test

# Full SuiteSparse collection (65 matrices, requires test-data/suitesparse/)
cargo test -- --ignored --test-threads=1

# Lint
cargo clippy --all-targets
cargo clippy --all-targets --features diagnostic
```

### Test Data

| Tier | Path | Size | Git status |
|------|------|------|------------|
| Hand-constructed | `test-data/hand-constructed/` | ~144KB | Tracked |
| CI subset | `test-data/suitesparse-ci/` | ~19MB | Tracked |
| Full suite | `test-data/suitesparse/` | ~500MB | Gitignored |

The full test suite is enumerated in `test-data/metadata.json` and can be downloaded from SuiteSparse (recommended) or using the download script (throttled downloads - very slow).

## Examples

See [`examples/README.md`](examples/README.md) for full documentation.

```bash
# Self-contained hello world (no external data needed)
cargo run --example basic_usage

# Multiple right-hand sides with workspace reuse
cargo run --example multiple_rhs

# Refactorization: same sparsity, different values
cargo run --example refactorization

# End-to-end solve timing on CI SuiteSparse matrices
cargo run --example solve_timing --release

# Single-matrix profiling with Chrome Trace export
cargo run --example profile_matrix --features diagnostic --release -- d_pretok --trace /tmp/trace.json

# Parallel scaling across thread counts
cargo run --example parallel_scaling --release -- --ci-only --threads 1,2,4,8

# Performance baseline collection for regression tracking
cargo run --example baseline_collection --features diagnostic --release -- --ci-only
```

## External Solver Comparisons

The `comparisons/` directory contains tools for benchmarking against reference
solvers (SPRAL and MUMPS). See
[`comparisons/README.md`](comparisons/README.md) for build instructions.

```bash
# Build SPRAL driver (prerequisite)
comparisons/drivers/build_spral.sh

# Side-by-side rivrs vs SPRAL on CI subset
cargo run --bin spral-comparison --release -- --ci-only --rivrs

# Control SPRAL thread count
cargo run --bin spral-comparison --release -- --ci-only --rivrs --threads 4
```

## License

Apache-2.0. See [LICENSE](../LICENSE) and [NOTICE](../NOTICE) for attribution.
