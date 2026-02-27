# Feature Specification: Benchmarking Framework

**Feature Branch**: `006-benchmarking-framework`
**Created**: 2026-02-07
**Status**: Draft
**Input**: User description: "implement phase 1.2 in dev/ssids-plan.md. Treat any code snippets as notional APIs that help to explain *what* not *how*"

## User Scenarios & Testing *(mandatory)*

### User Story 1 - Benchmark Individual Solver Phases (Priority: P1)

A developer working on the SSIDS solver wants to measure the wall-clock time and memory usage of individual solver phases (symbolic analysis, numeric factorization, triangular solve) in isolation. They select one or more test matrices, run the benchmark, and receive reproducible timing results via Criterion's statistical framework.

**Why this priority**: Phase-level benchmarks are the fundamental building block. Every subsequent benchmark capability (end-to-end, regression detection) depends on the ability to time individual operations against individual matrices. This also directly supports the iterative development workflow of Phases 2-8 where each component is built and optimized incrementally.

**Independent Test**: Can be fully tested by defining a trivial no-op solver stub, running a single-matrix benchmark, and confirming that Criterion produces valid statistical output (mean, std-dev, iterations).

**Acceptance Scenarios**:

1. **Given** a registered test matrix and a solver implementing the benchmark interface, **When** a component benchmark is executed for the "factor" phase, **Then** Criterion reports timing statistics (mean, standard deviation, throughput) for that phase on that matrix.
2. **Given** the benchmark suite, **When** a developer runs `cargo bench`, **Then** all configured component benchmarks execute without manual intervention and produce results in Criterion's standard output formats (terminal summary and HTML report).
3. **Given** a benchmark run, **When** results are collected, **Then** process-wide peak RSS is recorded alongside timing data as a whole-run metric.

---

### User Story 2 - End-to-End Solve Benchmarks (Priority: P2)

A developer wants to measure the full solve pipeline (analyze → factor → solve) as a single timed operation across multiple test matrices. This captures the real-world cost a user would experience, including overhead from phase transitions and data passing.

**Why this priority**: End-to-end benchmarks are the primary metric users and publications care about. However, they require the component benchmarks from P1 to exist first, and they only become meaningful once at least one solver phase is implemented.

**Independent Test**: Can be tested by running the full-roundtrip benchmark on a hand-constructed matrix with a mock solver and verifying aggregated timing is reported.

**Acceptance Scenarios**:

1. **Given** a set of test matrices spanning multiple sizes and difficulty levels, **When** end-to-end benchmarks execute, **Then** results report total wall-clock time for the full pipeline for each matrix. Per-phase breakdown is available separately via component benchmarks (US1).
2. **Given** a benchmark configuration specifying matrix subsets (e.g., "hand-constructed only" or "CI subset"), **When** run, **Then** only the specified matrices are benchmarked.

---

### User Story 3 - Performance Regression Detection (Priority: P3)

A developer wants to know if a code change has introduced a performance regression exceeding a configurable threshold. After running benchmarks, the system compares results against a saved baseline and flags regressions.

**Why this priority**: Regression detection closes the feedback loop — developers need to know promptly when an optimization attempt backfires or a refactor degrades throughput. It becomes useful as soon as a baseline exists from P1/P2 benchmarks.

**Independent Test**: Can be tested by saving a baseline result, artificially degrading a benchmark (e.g., adding a sleep), re-running, and confirming the regression is flagged.

**Acceptance Scenarios**:

1. **Given** a saved baseline and new benchmark results, **When** compared, **Then** any phase whose mean time increased beyond a configurable threshold (default 5%) is flagged as a regression.
2. **Given** a CI pipeline, **When** benchmarks run on a pull request, **Then** regressions are reported in a human-readable summary suitable for inclusion in PR comments or CI output.

---

### Edge Cases

- What happens when a test matrix is missing from disk? The benchmark should skip the matrix with a warning rather than aborting the entire suite.
- What happens when a solver phase is not yet implemented (returns an error or is a no-op)? The framework should gracefully skip that benchmark and report it as "not available."
- What happens when a benchmark runs for an unexpectedly long time (e.g., large matrix on slow hardware)? Criterion's configurable measurement time and sample count provide natural bounds; the framework should also support per-matrix timeout configuration.
- What happens when system load causes noisy timing? Criterion's statistical framework handles this through sufficient sampling and outlier detection; the framework should surface Criterion's confidence intervals.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: The system MUST provide a trait-based interface that any solver can implement to participate in benchmarking, covering the three solver phases: symbolic analysis, numeric factorization, and triangular solve.
- **FR-002**: The system MUST support benchmarking individual solver phases in isolation (component benchmarks) on any registered test matrix.
- **FR-003**: The system MUST support benchmarking the full solve pipeline (analyze → factor → solve) as a single operation (end-to-end benchmarks).
- **FR-004**: The system MUST integrate with Criterion.rs for statistical measurement, producing Criterion's standard outputs: terminal summary, HTML reports, and machine-readable JSON.
- **FR-005**: The system MUST support configurable matrix selection using the existing `TestCaseFilter` mechanism so benchmarks can target specific matrix subsets (hand-constructed, CI subset, by size, by difficulty).
- **FR-006**: The system MUST record process-wide peak RSS as a single metric per benchmark run. Fine-grained per-phase memory measurement is deferred to Phase 1.4 (profiling tools).
- **FR-007**: The system MUST support exporting benchmark results to at least one machine-readable format (CSV or JSON) for historical tracking.
- **FR-008**: The system MUST detect performance regressions by comparing new results against a saved baseline, flagging regressions exceeding a configurable threshold (default 5%).
- **FR-009**: The system MUST gracefully handle missing matrices (skip with warning) and unimplemented solver phases (report as unavailable) without aborting the suite.
- **FR-010**: The system MUST support per-matrix timeout configuration to prevent runaway benchmarks on unexpectedly large or slow inputs.
- **FR-011**: The system MUST be usable through the standard `cargo bench` command for component and end-to-end benchmarks.
- **FR-012**: The system MUST produce Markdown-formatted summary tables suitable for documentation and publication.

### Key Entities

- **Benchmark Configuration**: Specifies which matrices to benchmark, which solver phases to measure, and measurement parameters (sample size, warm-up time, timeout).
- **Benchmark Result**: A single measurement outcome containing the matrix name, operation phase, timing statistics (mean, std-dev, median, confidence interval), and arbitrary metadata. Whole-run peak RSS is recorded at the suite level, not per-result.
- **Benchmark Suite**: A collection of benchmark configurations and their aggregated results, supporting filtering and export.
- **Baseline**: A saved set of benchmark results used as a reference point for regression detection.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: Any solver component can be benchmarked on any registered test matrix with a single `cargo bench` invocation, producing statistically valid timing results within 5 minutes for the CI matrix subset.
- **SC-002**: Benchmark results are reproducible: re-running the same benchmark on the same hardware produces timing means within 10% of each other (excluding system noise outliers identified by Criterion).
- **SC-003**: Performance regressions exceeding 5% are detectable by comparing against a stored baseline, with clear identification of which component and matrix regressed.
- **SC-004**: Benchmark results for the full CI matrix subset can be exported to a machine-readable format and used to generate a Markdown comparison table without manual intervention.
- **SC-005**: The benchmarking framework adds zero overhead to production builds (gated behind feature flags or dev-dependency boundaries).
- **SC-006**: The benchmark trait interface is extensible enough that a future external solver adapter (e.g., SPRAL) can participate in benchmarks without modifying the framework itself.

## Assumptions

- Criterion.rs will remain the benchmarking framework; its statistical methodology (bootstrap confidence intervals, outlier detection) is sufficient for the project's needs.
- The existing `TestCaseFilter` and matrix registry infrastructure from Phase 0.5 is stable and will be reused for matrix selection in benchmarks.
- Memory measurement will use OS-level peak RSS as a whole-run metric. Per-phase memory instrumentation is deferred to Phase 1.4 (profiling tools with allocator tracking).
- The Criterion HTML report and built-in regression detection satisfy most reporting needs; custom reporting is layered on top rather than replacing Criterion.
- Per-matrix timeouts will be implemented at the benchmark harness level, not within Criterion itself, since Criterion does not natively support per-input timeouts.

## Scope Boundaries

### In Scope

- Criterion.rs integration for component and end-to-end benchmarks
- Matrix selection using existing `TestCaseFilter`
- Process-wide peak RSS measurement per benchmark run
- Machine-readable result export (CSV or JSON)
- Baseline save/load and regression detection
- Markdown summary table generation
- Graceful handling of missing matrices and unimplemented phases

### Out of Scope

- Parallel/multi-threaded benchmarking configurations (deferred to when the solver supports parallelism)
- Comparison benchmarks against external solvers (SPRAL, MUMPS, PARDISO); deferred until the SSIDS implementation is substantially complete
- GUI-based benchmark dashboards or interactive visualization
- Continuous benchmarking infrastructure (e.g., Bencher.dev integration); CI integration is limited to running benchmarks and comparing against baselines
- Custom allocator-based memory profiling (covered by Phase 1.4 profiling tools)

## Clarifications

### Session 2026-02-07

- Q: How should per-phase memory measurement work given that peak RSS is process-wide? → A: Whole-run peak RSS only; per-phase memory tracking deferred to Phase 1.4 (profiling tools).

## Dependencies

- Phase 0.5 test infrastructure (complete): `TestCaseFilter`, `SolverTest` trait, matrix registry
- Criterion.rs (already a dev-dependency)
- Existing `benches/matrix_loading.rs` will be superseded or refactored into the new framework
