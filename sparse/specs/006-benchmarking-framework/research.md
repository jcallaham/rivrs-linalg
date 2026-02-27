# Research: Benchmarking Framework

**Feature**: 006-benchmarking-framework
**Date**: 2026-02-07

## R-001: Criterion.rs Integration Patterns

**Decision**: Use Criterion 0.5 (current dep) with `BenchmarkGroup` per phase, parameterized by matrix name via `BenchmarkId::from_parameter`.

**Rationale**: The `BenchmarkGroup` + `BenchmarkId` pattern maps directly to our (phase, matrix) parameter space. Grouping by phase (one group for "analyze", one for "factorize", etc.) lets Criterion generate within-group comparison plots across matrices. Criterion remains actively maintained (0.8.2 as of Feb 2026, repo moved to `criterion-rs/criterion.rs`), but updating from 0.5 is not required for this feature — the API is stable across versions.

**Alternatives considered**:
- **Divan**: Simpler attribute-based API (`#[divan::bench]`), built-in allocation tracking. Rejected: lacks baseline management, HTML reports, and `BenchmarkId` parameterization.
- **CodSpeed**: CI-as-a-service wrapping Criterion. Rejected: adds external dependency for something achievable locally.
- **Group-per-matrix** (alternative Criterion layout): Rejected: group-per-phase is more natural for comparing solver performance across problem sizes.

## R-002: Multi-Matrix Benchmark Structure

**Decision**: One `BenchmarkGroup` per solver phase, parameterized by matrix name. Benchmarks conditionally registered (skip missing matrices with `continue`).

**Rationale**: Criterion has no built-in `skip()` mechanism; the idiomatic approach is to not register the benchmark. This aligns with FR-009 (graceful handling of missing matrices). The benchmark IDs will look like `ssids/factorize/bcsstk14`, enabling CLI filtering via `cargo bench -- "ssids/factorize"`.

**Per-group configuration**:
- `sample_size`, `measurement_time`, `warm_up_time` configurable per group
- `SamplingMode::Flat` available for slow benchmarks (large matrices)
- `Throughput::Elements(nnz)` to report operations/sec normalized by matrix size

## R-003: Baseline and Regression Detection

**Decision**: Use Criterion's built-in baseline management (`--save-baseline`, `--baseline` flags). Custom regression reporting layered on top by parsing Criterion output files.

**Rationale**: Criterion automatically compares each run against `base/` (previous run). Named baselines (`--save-baseline main`) persist in `target/criterion/.../main/`. The `change/estimates.json` file contains percentage change with confidence intervals. Criterion's two-stage filter (hypothesis test at `significance_level` + `noise_threshold`) handles statistical rigor.

**CI workflow**:
1. On `main`: `cargo bench -- --save-baseline main`
2. On feature branch: `cargo bench -- --baseline main`
3. Parse `change/estimates.json` or terminal output for regression flags

**Alternatives considered**:
- **cargo-criterion with `--message-format=json`**: Provides structured JSON output per benchmark. Worth considering for custom reporting but adds a separate binary dependency.
- **Custom baseline storage**: Rejected: Criterion's built-in mechanism is sufficient; no need to reimplement.

## R-004: Machine-Readable Output

**Decision**: Parse Criterion's `raw.csv` files (stable format) for result export. Supplement with a custom JSON sidecar file for metadata (peak RSS, matrix properties, run configuration).

**Rationale**: Criterion produces per-benchmark `raw.csv` files with columns: group, function, value, throughput_num, throughput_type, sample_measured_value, unit, iteration_count. This is documented as a stable format. The `estimates.json` files contain richer statistical data but are not guaranteed stable across versions.

**Alternatives considered**:
- **cargo-criterion `--message-format=json`**: Structured JSON to stdout. Good for CI pipelines but requires installing a separate binary. Can be added later.
- **Custom `Measurement` trait**: Would allow embedding metadata in Criterion's output. Rejected: Criterion only supports one measurement type per group, and wall-clock time is the primary metric.

## R-005: Peak RSS Measurement

**Decision**: Read `/proc/self/status` and parse `VmHWM` line. Zero dependencies, no `unsafe`. Record once per benchmark suite run as a whole-run metric.

**Rationale**: VmHWM (High Water Mark) is the kernel's peak RSS counter. It requires no crates, no `unsafe`, and is available on all Linux systems. The value is in kB.

**Key finding**: VmHWM *can* be reset on Linux 4.0+ by writing `"5"` to `/proc/self/clear_refs`. This means per-phase measurement is technically possible, but the spec decision (whole-run only, defer per-phase to 1.4) is the right call for simplicity.

**Alternatives considered**:
- **`libc::getrusage`**: Also gives peak RSS (`ru_maxrss` in KB on Linux). Requires `libc` crate and `unsafe`. `ru_maxrss` is truly monotonic (no reset). No advantage over `/proc/self/status`.
- **`peak_alloc` crate**: Measures allocator-level peak, not OS RSS. Misses mmap, stack, allocator fragmentation. Wrong metric for solver benchmarking.
- **`sysinfo` / `procfs` crates**: Heavy dependencies for a single metric. Overkill.

## R-006: Benchmark Trait vs Existing SolverTest Trait

**Decision**: Create a separate `Benchmarkable` trait for benchmarking, distinct from `SolverTest`. Reuse `SolverTestCase` and `TestCaseFilter` for input selection.

**Rationale**: `SolverTest` is validation-focused — its methods return `TestResult` with correctness metrics (reconstruction error, inertia checks). Benchmarking needs a different interface: methods that execute a solver phase and return nothing (Criterion handles timing externally via closure wrapping). The traits serve different purposes:
- `SolverTest`: "Is the output correct?"
- `Benchmarkable`: "How fast does it run?"

However, the input layer is shared: both use `SolverTestCase` for matrix data and `TestCaseFilter` for selection. The benchmark framework should depend on the `testing` module's case/filter infrastructure.

**Interface sketch** (notional — describes *what*, not *how*):
- `bench_analyze(case)`: Run symbolic analysis phase
- `bench_factor(case)`: Run numeric factorization phase
- `bench_solve(case)`: Run triangular solve phase
- `bench_roundtrip(case)`: Run full pipeline

Each returns an opaque result (to prevent Criterion from optimizing away the computation).

## R-007: Unimplemented Phase Handling

**Decision**: The `Benchmarkable` trait methods return `Option` or `Result` to signal whether a phase is implemented. The benchmark harness skips `None`/`Err` results with a diagnostic message.

**Rationale**: During iterative development (Phases 2-8), only some solver phases will be implemented at any time. The framework must not fail when a phase is unavailable. This aligns with FR-009 and the edge case spec.

## R-008: Markdown Report Generation

**Decision**: Post-process Criterion's `raw.csv` + `estimates.json` into a Markdown summary table. This is a standalone utility function, not integrated into Criterion's measurement loop.

**Rationale**: Criterion produces HTML reports natively but not Markdown. A simple post-processing step can read the estimates and format them into a table like:

```
| Matrix       | Phase     | Mean (ms) | Std Dev | Throughput (nnz/s) |
|--------------|-----------|-----------|---------|---------------------|
| bcsstk14     | factorize | 12.3      | 0.4     | 1.2M               |
```

This satisfies FR-012 without modifying Criterion's internals.
