# Data Model: Benchmarking Framework

**Feature**: 006-benchmarking-framework
**Date**: 2026-02-07

## Entities

### BenchmarkConfig

Describes what to benchmark and how to measure it.

**Fields**:
- `filter: TestCaseFilter` — Selects which matrices to include (reuses existing type)
- `phases: Vec<BenchmarkPhase>` — Which solver phases to benchmark
- `sample_size: Option<usize>` — Criterion sample count override (default: Criterion's default of 100)
- `measurement_time: Option<Duration>` — Criterion measurement time override (default: 5s)
- `warm_up_time: Option<Duration>` — Criterion warm-up time override (default: 3s)
- `timeout_per_matrix: Option<Duration>` — Maximum time for any single matrix benchmark before skipping

**Relationships**: References `TestCaseFilter` (from `testing::cases`). Used by benchmark harness functions.

### BenchmarkPhase

Enumeration of measurable solver operations.

**Variants**:
- `Analyze` — Symbolic analysis
- `Factor` — Numeric factorization
- `Solve` — Triangular solve (forward/backward substitution)
- `Roundtrip` — Full pipeline (analyze → factor → solve)

**Notes**: Mirrors `TestKind` from `testing::harness` but is a separate type to allow independent evolution.

### BenchmarkResult

A single measurement outcome from one (matrix, phase) pair.

**Fields**:
- `matrix_name: String` — Which test matrix was benchmarked
- `phase: BenchmarkPhase` — Which solver operation
- `mean_ns: f64` — Mean execution time in nanoseconds
- `std_dev_ns: f64` — Standard deviation in nanoseconds
- `median_ns: f64` — Median execution time in nanoseconds
- `iterations: u64` — Total iterations Criterion ran
- `throughput_nnz_per_sec: Option<f64>` — Operations/sec normalized by matrix nnz (if throughput set)
- `matrix_size: usize` — Matrix dimension (n)
- `matrix_nnz: usize` — Number of nonzeros

**Lifecycle**: Created after a Criterion benchmark group completes, populated from Criterion's output files.

### BenchmarkSuiteResult

Aggregated results from a complete benchmark run.

**Fields**:
- `results: Vec<BenchmarkResult>` — All individual measurements
- `peak_rss_kb: Option<u64>` — Process-wide peak RSS in kilobytes (from VmHWM)
- `skipped: Vec<SkippedBenchmark>` — Matrices/phases that were skipped and why
- `timestamp: String` — ISO 8601 timestamp of the run
- `config: BenchmarkConfig` — The configuration used

**Relationships**: Contains `BenchmarkResult` entries. Serializable to JSON for baseline storage and export.

### SkippedBenchmark

Records a benchmark that was not executed.

**Fields**:
- `matrix_name: String` — Which matrix
- `phase: BenchmarkPhase` — Which phase
- `reason: String` — Why skipped (e.g., "matrix file not found", "phase not implemented")

### Baseline

A saved suite result used for regression detection.

**Fields**:
- `name: String` — Baseline identifier (e.g., "main", "v0.1.0")
- `suite_result: BenchmarkSuiteResult` — The saved results
- `created: String` — ISO 8601 timestamp

**Storage**: JSON file at a configurable path (default: `target/benchmarks/baselines/{name}.json`).

### RegressionReport

Comparison of current results against a baseline.

**Fields**:
- `baseline_name: String` — Which baseline was compared against
- `regressions: Vec<Regression>` — Detected regressions
- `improvements: Vec<Improvement>` — Detected improvements
- `unchanged: Vec<String>` — Benchmark IDs with no significant change
- `threshold_pct: f64` — The threshold used (default: 5.0)

### Regression / Improvement

A single detected performance change.

**Fields**:
- `matrix_name: String`
- `phase: BenchmarkPhase`
- `baseline_mean_ns: f64`
- `current_mean_ns: f64`
- `change_pct: f64` — Percentage change (positive = regression, negative = improvement)

## Entity Relationships

```
BenchmarkConfig ──uses──▶ TestCaseFilter (existing)
       │
       ▼
BenchmarkSuiteResult ──contains──▶ BenchmarkResult (many)
       │                           SkippedBenchmark (many)
       │
       ▼
    Baseline ──saved as──▶ JSON file
       │
       ▼
RegressionReport ──compares──▶ BenchmarkSuiteResult (current vs baseline)
       │
       ├──contains──▶ Regression (many)
       └──contains──▶ Improvement (many)
```

## Reused Types (from existing codebase)

- `TestCaseFilter` — Matrix selection (from `testing::cases`)
- `SolverTestCase` — Matrix + metadata + reference (from `testing::cases`)
- `TestMatrixProperties` — Matrix metadata (from `testing::cases`)
- `SparseColMat<usize, f64>` — Matrix data (from `faer`)

## Serialization

All benchmark data entities (`BenchmarkResult`, `BenchmarkSuiteResult`, `Baseline`, `RegressionReport`) are serializable to JSON via serde. This supports:
- Baseline save/load (FR-008)
- Machine-readable export (FR-007)
- Markdown table generation from deserialized data (FR-012)
