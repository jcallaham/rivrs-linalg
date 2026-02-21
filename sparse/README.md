# rivrs-sparse

Sparse Symmetric Indefinite Direct Solver (SSIDS) for Rust, inspired by [SPRAL](https://github.com/ralna/spral). Clean-room implementation using [faer](https://github.com/sarah-ek/faer-rs) for dense linear algebra.

**Status**: Feature-complete through Phase 8.2 (parallel factorization & solve). 65/65 SuiteSparse test matrices passing. Not yet published to crates.io.

## Features

- Multifrontal LDL^T factorization with A Posteriori Threshold Pivoting (APTP)
- Two-level pivoting: TPP for small fronts, complete pivoting with BLAS-3 blocking for large fronts
- MC64 weighted bipartite matching & scaling + METIS nested dissection ordering
- Parallel factorization & solve via rayon (tree-level) and faer `Par` (intra-node BLAS)
- Three-phase API: analyze → factor → solve (symbolic reuse across re-factorizations)

## API

```rust
use faer::sparse::{SparseColMat, Triplet};
use faer::Col;
use rivrs_sparse::aptp::{SparseLDLT, SolverOptions};

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

For the three-phase API with parallel factorization:

```rust
use faer::{Par, dyn_stack::{MemBuffer, MemStack}};
use rivrs_sparse::aptp::{AnalyzeOptions, FactorOptions, SparseLDLT};

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

## Benchmarking

### SPRAL comparison

Requires building the SPRAL driver first (`tools/build_spral.sh`).

```bash
# Side-by-side rivrs vs SPRAL on CI subset
cargo run --example spral_benchmark --release -- --ci-only --rivrs

# Full collection, 4 threads
cargo run --example spral_benchmark --release -- --rivrs --threads 4

# Filter by category
cargo run --example spral_benchmark --release -- --rivrs --category hard-indefinite
```

### Parallel scaling

```bash
# CI subset, default thread counts
cargo run --example parallel_scaling --features diagnostic --release -- --ci-only

# Custom thread counts
cargo run --example parallel_scaling --features diagnostic --release -- --threads 1,2,4,8

# Full collection
cargo run --example parallel_scaling --features diagnostic --release
```

### Baseline collection

Structured JSON with per-phase timing, per-supernode stats, and backward error.

```bash
cargo run --example baseline_collection --features diagnostic --release -- --ci-only
cargo run --example baseline_collection --features diagnostic --release -- \
  --compare target/benchmarks/baselines/prev.json
```

## Documentation

- [CLAUDE.md](CLAUDE.md) — Development guidance, algorithm architecture, code layout
- [docs/ssids-plan.md](docs/ssids-plan.md) — Development plan and phase roadmap
- [docs/ssids-log.md](docs/ssids-log.md) — Development changelog

## License

Apache-2.0. See [LICENSE](../LICENSE) and [NOTICE](../NOTICE) for attribution.
