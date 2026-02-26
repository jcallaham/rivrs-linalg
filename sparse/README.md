# rivrs-sparse

Sparse Symmetric Indefinite Direct Solver (SSIDS) for Rust, inspired by [SPRAL](https://github.com/ralna/spral). Clean-room implementation using [faer](https://github.com/sarah-ek/faer-rs) for dense linear algebra.

**Status**: Feature-complete (parallel multifrontal LDL^T with APTP pivoting). 65/65 SuiteSparse test matrices passing. Not yet published to crates.io.

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

For the three-phase API with parallel factorization:

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

## Examples

See [`examples/README.md`](examples/README.md) for full documentation.

```bash
# Self-contained hello world (no external data needed)
cargo run --example basic_usage

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
solvers (SPRAL, with MUMPS and HSL MA27 planned). See
[`comparisons/README.md`](comparisons/README.md) for build instructions.

```bash
# Build SPRAL driver (prerequisite)
comparisons/drivers/build_spral.sh

# Side-by-side rivrs vs SPRAL on CI subset
cargo run --bin spral-comparison --release -- --ci-only --rivrs

# Control SPRAL thread count
cargo run --bin spral-comparison --release -- --ci-only --rivrs --threads 4
```

## Documentation

- [CLAUDE.md](CLAUDE.md) — Development guidance, algorithm architecture, code layout
- [docs/ssids-plan.md](docs/ssids-plan.md) — Development plan and phase roadmap
- [docs/ssids-log.md](docs/ssids-log.md) — Development changelog

## License

Apache-2.0. See [LICENSE](../LICENSE) and [NOTICE](../NOTICE) for attribution.
