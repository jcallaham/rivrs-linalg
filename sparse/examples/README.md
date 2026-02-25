# Examples

rivrs-sparse implements a multifrontal LDL^T solver for sparse symmetric
indefinite systems using A Posteriori Threshold Pivoting (APTP). The solver
exposes a three-phase API:

1. **Analyze** — symbolic factorization (reusable for matrices with the same sparsity pattern)
2. **Factor** — numeric LDL^T factorization with APTP pivoting
3. **Solve** — forward/diagonal/backward substitution

## Quick Start

```sh
cargo run --example basic_usage
```

No external data files or feature flags required.

## Examples

### basic_usage

Self-contained example that constructs a small symmetric indefinite matrix from
triplets, solves `Ax = b` for a known solution, and validates the result with
backward error analysis. Demonstrates the core `SparseLDLT` API without any
external dependencies.

```sh
cargo run --example basic_usage
```

### solve_timing

Solves all CI-subset SuiteSparse matrices and reports factorization/solve timing,
backward error, and pivot statistics. Good for a quick end-to-end validation.

Requires the CI test matrices in `test-data/suitesparse-ci/` (checked into git).

```sh
cargo run --example solve_timing --release
```

### profile_matrix

Deep profiling tool for a single matrix. Reports per-phase timing, factorization
statistics, and per-supernode breakdowns. Optionally exports
[Chrome Trace](https://ui.perfetto.dev) JSON for visual analysis.

Requires `--features diagnostic` for per-supernode timing instrumentation.

```sh
# Profile by registry name
cargo run --example profile_matrix --features diagnostic --release -- d_pretok

# Profile by file path
cargo run --example profile_matrix --features diagnostic --release -- --path matrix.mtx

# Export Chrome Trace
cargo run --example profile_matrix --features diagnostic --release -- d_pretok --trace /tmp/trace.json

# List available matrices
cargo run --example profile_matrix --features diagnostic --release -- --list
```

### baseline_collection

Collects structured performance baselines (JSON) across the SuiteSparse test
suite for regression tracking. Records per-phase timing, per-supernode
statistics, peak RSS, and backward error. Supports comparison against a
previously saved baseline to detect regressions.

Requires `--features diagnostic` for per-supernode timing fields.

```sh
# Collect baselines for CI subset
cargo run --example baseline_collection --features diagnostic --release -- --ci-only

# Collect for all available matrices
cargo run --example baseline_collection --features diagnostic --release

# Compare against a previous baseline
cargo run --example baseline_collection --features diagnostic --release -- --compare target/benchmarks/baselines/prev.json
```

JSON output is written to `target/benchmarks/baselines/`.

### parallel_scaling

Measures factorization and solve performance across multiple thread counts.
Produces a structured report showing speedup and efficiency per matrix.

```sh
# Default thread counts (1, 2, 4) on CI subset
cargo run --example parallel_scaling --release -- --ci-only

# Custom thread counts
cargo run --example parallel_scaling --release -- --ci-only --threads 1,2,4,8

# All available matrices
cargo run --example parallel_scaling --release
```

JSON output is written to `target/benchmarks/parallel/`.

## Feature Flags

| Flag | Purpose | Used by |
|------|---------|---------|
| `diagnostic` | Per-supernode timing instrumentation (zero overhead when disabled) | `profile_matrix`, `baseline_collection` |

The `test-util` feature (random matrix generators, benchmarking utilities) is
automatically enabled for examples via the dev-dependency in `Cargo.toml`.

## Test Data

Some examples require SuiteSparse test matrices:

| Tier | Path | Size | Git status |
|------|------|------|------------|
| CI subset | `test-data/suitesparse-ci/` | ~19 MB | Tracked |
| Full suite | `test-data/suitesparse/` | ~500 MB | Gitignored |

`basic_usage` requires no external data. All other examples use the CI subset
by default.
