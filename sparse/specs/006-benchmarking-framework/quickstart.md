# Quickstart: Benchmarking Framework

**Feature**: 006-benchmarking-framework

## Running Benchmarks

### All benchmarks
```bash
cargo bench
```

### Component benchmarks only (e.g., factorize phase)
```bash
cargo bench -- "ssids/factorize"
```

### Specific matrix
```bash
cargo bench -- "bcsstk14"
```

### With a saved baseline
```bash
# Save current results as a named baseline
cargo bench -- --save-baseline main

# Compare a feature branch against the baseline
cargo bench -- --baseline main
```

## Benchmark Organization

Benchmarks are organized as Criterion groups:

```
ssids/analyze/{matrix_name}      — Symbolic analysis per matrix
ssids/factorize/{matrix_name}    — Numeric factorization per matrix
ssids/solve/{matrix_name}        — Triangular solve per matrix
ssids/roundtrip/{matrix_name}    — Full pipeline per matrix
```

## Matrix Selection

Benchmarks use the same `TestCaseFilter` as the test infrastructure:

- **Hand-constructed matrices** (15 matrices, ~144KB): Always available, fast benchmarks for development iteration
- **CI subset** (10 SuiteSparse matrices, ~73MB): Representative real-world matrices for CI
- **Full SuiteSparse** (67 matrices): For comprehensive performance evaluation (requires extracted archive)

Matrix selection is configured in the benchmark binary. Missing matrices are silently skipped.

## Output Locations

| Output | Location |
|--------|----------|
| Criterion HTML reports | `target/criterion/{group}/{matrix}/report/index.html` |
| Raw timing data | `target/criterion/{group}/{matrix}/new/raw.csv` |
| Statistical estimates | `target/criterion/{group}/{matrix}/new/estimates.json` |
| Saved baselines | `target/benchmarks/baselines/{name}.json` |
| Exported CSV | `target/benchmarks/results.csv` |
| Peak RSS | Printed to stderr at end of benchmark run |

## Adding a New Solver to Benchmarks

1. Implement the `Benchmarkable` trait for your solver
2. Return `None` from any phase method not yet implemented (the harness will skip it)
3. The solver appears in benchmark output alongside any other registered solvers

## Regression Detection Workflow

```bash
# 1. Establish baseline on main branch
git checkout main
cargo bench -- --save-baseline main

# 2. Switch to feature branch and compare
git checkout feature-branch
cargo bench -- --baseline main

# 3. Check Criterion output for "Performance has regressed" warnings
```

Criterion flags regressions exceeding the noise threshold (default 1%) with statistical significance testing. The custom regression detector provides a configurable threshold (default 5%) and structured reporting.
