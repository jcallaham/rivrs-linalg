# Tasks: Benchmarking Framework

**Input**: Design documents from `/specs/006-benchmarking-framework/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/

**Tests**: Tests are included per Constitution Principle III (TDD). Data structures tested via unit tests; Criterion integration verified with mock solver.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Create module structure, data types, and trait definitions that all user stories depend on

- [X] T001 Create `src/benchmarking/mod.rs` with module declarations and public re-exports for: `config`, `traits`, `results`, `rss`, `baseline`, `report`
- [X] T002 Register `benchmarking` module in `src/lib.rs` behind `#[cfg(feature = "test-util")]` gate (same as `testing` module)
- [X] T003 Register new benchmark binary in `Cargo.toml`: add `[[bench]] name = "solver_benchmarks"` with `harness = false`

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core data types, trait definition, and RSS measurement that MUST be complete before any user story benchmarks can run

**CRITICAL**: No user story work can begin until this phase is complete

- [X] T004 [P] Implement `BenchmarkPhase` enum (`Analyze`, `Factor`, `Solve`, `Roundtrip`) with `Display` impl and serde `Serialize`/`Deserialize`, and `BenchmarkConfig` struct with fields (`filter: TestCaseFilter`, `phases: Vec<BenchmarkPhase>`, `sample_size: Option<usize>`, `measurement_time: Option<Duration>`, `warm_up_time: Option<Duration>`, `timeout_per_matrix: Option<Duration>`) and builder methods, in `src/benchmarking/config.rs`
- [X] T005 [P] Implement `Benchmarkable` trait in `src/benchmarking/traits.rs` with methods: `bench_analyze`, `bench_factor`, `bench_solve`, `bench_roundtrip` (with default impl chaining the first three). Each returns `Option<Box<dyn Any>>` or `Option<Vec<f64>>` per contract
- [X] T006 [P] Implement `BenchmarkResult`, `BenchmarkSuiteResult`, and `SkippedBenchmark` structs with serde `Serialize`/`Deserialize` in `src/benchmarking/results.rs`. Include `Display` impls for human-readable output
- [X] T007 [P] Implement `read_peak_rss_kb()` function in `src/benchmarking/rss.rs` that parses `VmHWM` from `/proc/self/status`. Return `Option<u64>` (None on non-Linux). Include unit test with `#[cfg(target_os = "linux")]` guard
- [X] T008 Implement `MockBenchmarkable` struct in `src/benchmarking/traits.rs` (or a test submodule) that implements `Benchmarkable` with configurable no-op or sleep-based phases for testing the harness. Uses existing `SolverTestCase` matrix data
- [X] T009 Add unit tests for `BenchmarkConfig` builder, `BenchmarkPhase` display/serde round-trip, and `BenchmarkResult` serde round-trip in their respective modules

**Checkpoint**: Foundation ready — data types compile, trait defined, RSS measurable, mock solver available. User story implementation can now begin.

---

## Phase 3: User Story 1 — Benchmark Individual Solver Phases (Priority: P1) MVP

**Goal**: A developer can run `cargo bench` and get Criterion statistical output for individual solver phases (analyze, factor, solve) across selected test matrices.

**Independent Test**: Run `cargo bench -- "ssids/"` with the `MockBenchmarkable` and verify Criterion produces timing statistics and HTML reports in `target/criterion/`.

### Tests for User Story 1

- [X] T010 [US1] Write test in `src/benchmarking/config.rs` (or `tests/`) that loads hand-constructed test cases via `TestCaseFilter::hand_constructed()` and verifies they can be iterated for benchmarking (non-empty, matrices loadable)
- [X] T011 [US1] Write test that constructs a `MockBenchmarkable`, calls each phase method with a hand-constructed matrix, and verifies `Some(...)` is returned for implemented phases and `None` for unimplemented ones

### Implementation for User Story 1

- [X] T012 [US1] Implement `run_component_benchmarks` function in `benches/solver_benchmarks.rs` that: loads test cases via `TestCaseFilter`, creates one `BenchmarkGroup` per phase (named `ssids/{phase}`), parameterizes by matrix name via `BenchmarkId::from_parameter`, skips missing matrices with `eprintln!` warning, skips unimplemented phases (when `Benchmarkable` method returns `None`), enforces `timeout_per_matrix` from `BenchmarkConfig` by recording a `SkippedBenchmark` when a matrix exceeds the deadline, and sets `Throughput::Elements(nnz)` per benchmark
- [X] T013 [US1] Wire up `criterion_group!` and `criterion_main!` in `benches/solver_benchmarks.rs` with the `MockBenchmarkable` as default solver, using `BenchmarkConfig` with `TestCaseFilter::hand_constructed()` as default matrix set
- [X] T014 [US1] Add peak RSS measurement: call `read_peak_rss_kb()` before and after the benchmark suite in `benches/solver_benchmarks.rs`, print the result to stderr
- [X] T015 [US1] Verify `cargo bench -- "ssids/"` executes successfully, producing Criterion terminal output and HTML reports in `target/criterion/ssids/`

**Checkpoint**: `cargo bench` runs component benchmarks on hand-constructed matrices with mock solver. Criterion HTML reports generated. Peak RSS printed. US1 is independently functional.

---

## Phase 4: User Story 2 — End-to-End Solve Benchmarks (Priority: P2)

**Goal**: A developer can benchmark the full solve pipeline (analyze → factor → solve) as a single timed operation, with configurable matrix subset selection.

**Independent Test**: Run `cargo bench -- "ssids/roundtrip"` and verify roundtrip benchmarks appear with per-matrix timing.

### Tests for User Story 2

- [X] T016 [US2] Write test that calls `MockBenchmarkable::bench_roundtrip` with a hand-constructed matrix and verifies the full pipeline produces a result

### Implementation for User Story 2

- [X] T017 [US2] Implement `run_e2e_benchmarks` function in `benches/solver_benchmarks.rs` that creates a `BenchmarkGroup` named `ssids/roundtrip`, benchmarks the full pipeline per matrix via `bench_roundtrip`, skips missing matrices, enforces timeout, and reports throughput
- [X] T018 [US2] Add matrix subset selection support: accept a `BenchmarkConfig` that controls which `TestCaseFilter` preset is used (hand-constructed, CI subset, or all). Demonstrate with at least two filter presets wired in the benchmark binary
- [X] T019 [US2] Register `run_e2e_benchmarks` in the `criterion_group!` alongside component benchmarks so both run under `cargo bench`

**Checkpoint**: `cargo bench` runs both component and roundtrip benchmarks. Matrix filtering works. US1 and US2 both function independently.

---

## Phase 5: User Story 3 — Performance Regression Detection (Priority: P3)

**Goal**: A developer can save benchmark results as a baseline, re-run benchmarks, and detect regressions exceeding a configurable threshold (default 5%).

**Independent Test**: Save a baseline with the mock solver, introduce an artificial delay, re-run, and verify the regression report correctly flags the degraded phase.

### Tests for User Story 3

- [X] T020 [P] [US3] Write unit tests for `save_baseline` and `load_baseline` in `src/benchmarking/baseline.rs`: create a `BenchmarkSuiteResult` with synthetic data, save to a temp file, load back, and verify round-trip equality
- [X] T021 [P] [US3] Write unit tests for `detect_regressions` in `src/benchmarking/baseline.rs`: construct two `BenchmarkSuiteResult` instances with known timing differences, run regression detection at 5% threshold, and verify correct classification (regression, improvement, unchanged)

### Implementation for User Story 3

- [X] T022 [US3] Implement `Baseline` struct, `save_baseline`, and `load_baseline` functions in `src/benchmarking/baseline.rs`. Use serde JSON serialization to `target/benchmarks/baselines/{name}.json`. Create parent directories if needed
- [X] T023 [US3] Implement `RegressionReport`, `Regression`, and `Improvement` structs in `src/benchmarking/baseline.rs` with serde derives and `Display` impl for human-readable summary
- [X] T024 [US3] Implement `detect_regressions` function in `src/benchmarking/baseline.rs` that matches (matrix_name, phase) pairs between current and baseline results, computes percentage change, and classifies each as regression/improvement/unchanged based on configurable threshold
- [X] T025 [US3] Implement `collect_results` function in `src/benchmarking/results.rs` (or `baseline.rs`) that parses Criterion's `raw.csv` and `estimates.json` from `target/criterion/` directory tree into a `BenchmarkSuiteResult`

**Checkpoint**: Baselines can be saved/loaded. Regressions are detected and reported. US3 is independently testable with synthetic data.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Export, reporting, documentation, and integration improvements that span all user stories

- [X] T026 [P] Implement `export_csv` function in `src/benchmarking/report.rs` that writes `BenchmarkSuiteResult` to a CSV file with columns: matrix_name, phase, mean_ns, std_dev_ns, median_ns, iterations, matrix_size, matrix_nnz, throughput_nnz_per_sec
- [X] T027 [P] Implement `export_json` function in `src/benchmarking/report.rs` that serializes `BenchmarkSuiteResult` to a JSON file via serde
- [X] T028 [P] Implement `generate_markdown_table` function in `src/benchmarking/report.rs` that formats `BenchmarkSuiteResult` into an aligned Markdown table (matrix, phase, mean, std_dev, throughput columns). Add a variant for `RegressionReport` that includes change percentage
- [X] T029 Write unit tests for `export_csv`, `export_json`, and `generate_markdown_table` in `src/benchmarking/report.rs` using synthetic `BenchmarkSuiteResult` data: verify CSV column count, JSON round-trip, and Markdown table structure
- [X] T030 Add rustdoc comments with `# Examples` sections to all public items in `src/benchmarking/` modules: trait, config builder, results, baseline, report, RSS
- [X] T031 Run `cargo test` and `cargo bench` end-to-end to verify all modules integrate: tests pass, benchmarks produce output, no warnings from `cargo clippy`
- [X] T032 Update `docs/ssids-log.md` with Phase 1.2 completion entry documenting what was built, key design decisions (separate `Benchmarkable` trait, whole-run RSS, Criterion group-per-phase), and any deviations from the original plan

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion — BLOCKS all user stories
- **US1 (Phase 3)**: Depends on Foundational phase completion
- **US2 (Phase 4)**: Depends on Foundational phase completion. Shares benchmark binary with US1 but is independently testable via `cargo bench -- "ssids/roundtrip"`
- **US3 (Phase 5)**: Depends on Foundational phase completion. Independent of US1/US2 at the data-structure level, but most useful after benchmarks exist to generate baselines
- **Polish (Phase 6)**: Depends on all user stories being complete

### User Story Dependencies

- **User Story 1 (P1)**: Can start after Foundational (Phase 2) — No dependencies on other stories
- **User Story 2 (P2)**: Can start after Foundational (Phase 2) — Shares `benches/solver_benchmarks.rs` with US1 but adds its own function. Best done after US1 to extend the same benchmark binary
- **User Story 3 (P3)**: Can start after Foundational (Phase 2) — Baseline and regression types are independent. `collect_results` depends on Criterion output existing (from US1/US2 runs)

### Within Each User Story

- Tests written first, verified to fail (or compile with mock), then implementation
- Data types before functions that use them
- Core logic before integration/wiring
- Story complete before moving to next priority

### Parallel Opportunities

- T004, T005, T006, T007 can all run in parallel (different files, no dependencies)
- T010, T011 can run in parallel with each other
- T020, T021 can run in parallel with each other
- T026, T027, T028 can all run in parallel (different functions in same file, but independent)
- US3 implementation (T022-T025) is parallelizable with US1/US2 at the type level

---

## Parallel Example: Foundational Phase

```
# Launch all data type tasks together (different files):
T004: BenchmarkPhase + BenchmarkConfig in src/benchmarking/config.rs
T005: Benchmarkable trait in src/benchmarking/traits.rs
T006: Result structs in src/benchmarking/results.rs
T007: Peak RSS in src/benchmarking/rss.rs
```

## Parallel Example: User Story 3

```
# Launch both test tasks together:
T020: Unit tests for save/load baseline
T021: Unit tests for detect_regressions
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T003)
2. Complete Phase 2: Foundational (T004-T009)
3. Complete Phase 3: User Story 1 (T010-T015)
4. **STOP and VALIDATE**: Run `cargo bench -- "ssids/"` and verify Criterion output
5. This delivers: component benchmarks for all solver phases, peak RSS, mock solver demonstration

### Incremental Delivery

1. Setup + Foundational → Foundation ready
2. Add User Story 1 → `cargo bench` works with component benchmarks (MVP!)
3. Add User Story 2 → End-to-end roundtrip benchmarks added
4. Add User Story 3 → Baseline management and regression detection
5. Polish → Export, Markdown reports, documentation
6. Each story adds value without breaking previous stories

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Each user story should be independently completable and testable
- Commit after each task or logical group
- Stop at any checkpoint to validate story independently
- The `MockBenchmarkable` from T008 is the primary test vehicle until real solver phases exist (Phase 2+)
- Existing `benches/matrix_loading.rs` is retained unchanged; `solver_benchmarks.rs` is a new, separate benchmark binary
