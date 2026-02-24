# Quickstart: Small Leaf Subtree Fast Path

**Feature**: 025-small-leaf-fastpath

## Build

```bash
cd /workspace/rivrs-linalg/sparse
cargo build                          # Standard build
cargo build --features diagnostic    # With per-supernode timing
cargo build --release                # Optimized (needed for benchmarking)
```

## Test

```bash
# Unit tests (includes new small-leaf classification and fast-path tests)
cargo test

# With diagnostic feature (verifies cfg-gated code compiles)
cargo test --features diagnostic

# Full SuiteSparse regression check (requires test-data/suitesparse/)
cargo test -- --ignored --test-threads=1
```

## Benchmark

```bash
# Profile simplicial matrices (target workload)
cargo run --example profile_matrix --features diagnostic --release -- dixmaanl
cargo run --example profile_matrix --features diagnostic --release -- bloweybq
cargo run --example profile_matrix --features diagnostic --release -- mario001

# Full baseline comparison
cargo run --example baseline_collection --features diagnostic --release -- --ci-only

# Compare against previous baseline
cargo run --example baseline_collection --features diagnostic --release -- \
  --compare target/benchmarks/baselines/prev.json
```

## Configuration

The small-leaf threshold is configurable via `FactorOptions`:

```rust
use rivrs_sparse::aptp::solver::{SparseLDLT, FactorOptions};

let mut factor_opts = FactorOptions::default();
factor_opts.small_leaf_threshold = 256;  // default: 256, 0 = disabled

solver.factor(&matrix, &factor_opts)?;
```

## Verify Correctness

After changes, verify no regression:

```bash
# All 65 SuiteSparse matrices must pass (backward error < 5e-11)
cargo test -- --ignored --test-threads=1

# Check simplicial matrices specifically
cargo test suitesparse -- --ignored --test-threads=1 | grep -E 'dixmaanl|bloweybq|mario001|linverse|spmsrtls'
```

## Key Files

| File | Role |
|------|------|
| `src/aptp/numeric.rs` | Subtree classification, fast-path factorization, level-set integration |
| `src/aptp/solver.rs` | FactorOptions/SolverOptions threshold field |
| `tests/` | Classification unit tests, fast-path correctness tests, regression tests |
