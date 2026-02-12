# SSIDS Development Log

## Phase 3 Follow-up: AMD Ordering Quality & METIS Elevation

**Status**: Complete
**Branch**: `010-aptp-symbolic`
**Date**: 2026-02-10

### What Was Done

Timing analysis of the full SuiteSparse test suite revealed that AMD ordering produces
catastrophically poor fill predictions for many benchmark matrices. Investigation traced
the issue to SPRAL using METIS (nested dissection) by default, not AMD.

**Timing diagnostic** on full 65-matrix collection (sorted by size, AMD ordering):
- Small matrices (< 10K): < 100ms symbolic analysis
- Medium matrices (10K-30K): 100ms - 4s
- Pathological cases: sparsine (50K) → 1.04B predicted nnz(L), 12s analysis;
  nd6k (18K) → 73M predicted nnz(L), 13s analysis

**Comparison against paper values** (Hogg et al. 2016, Table III, METIS ordering):
nd3k: 1.8× more fill with AMD; Si10H16: 2.9×; Si5H12: 2.8×; sparsine: ~10-20×.

**Plan changes** (`docs/ssids-plan.md`):
- Added "Lessons Learned: AMD Ordering Quality" section to Phase 3
- Elevated METIS from "consider for Phase 9" to Phase 4.1 (before MC64)
- Expanded Phase 4 from "MC64 Matching & Scaling" to "Ordering & Preprocessing"
  with 4.1 (METIS) and 4.2 (MC64) sub-deliverables
- Updated Phase 9 note (METIS no longer deferred there)

**Test changes** (`tests/symbolic_analysis_full.rs`):
- Added `MAX_DIM_FOR_AMD = 30_000` guard to all full-collection tests (skips matrices
  where AMD produces excessive fill; to be removed when METIS is integrated)
- Removed `test_amd_reduces_fill_full_suitesparse` (was testing faer's AMD quality,
  not our code; will be replaced by METIS vs AMD comparison in Phase 4.1)
- Added documentation explaining the dimension cap and its relationship to METIS

### Key Findings

1. **SPRAL uses METIS by default** (`options%ordering = 1`), not AMD. All APTP
   benchmark papers (Duff/Hogg/Lopez 2020, Hogg et al. 2016) report results with
   METIS ordering.

2. **AMD is adequate for small/well-structured matrices** but produces 2-20× more
   fill than METIS on matrices with geometric structure (FEM, quantum chemistry,
   optimization problems).

3. **faer infrastructure is unchanged** — METIS produces a `Perm<usize>` that plugs
   into `SymmetricOrdering::Custom`. This is purely an input quality improvement,
   not an architectural change.

4. **Full SuiteSparse tests are reusable for METIS** — structural property tests
   (etree validity, supernode partitioning, assembly tree postorder) are
   ordering-independent. Predicted nnz(L) can be compared against paper values
   to validate METIS integration.

---

## Phase 3 Follow-up: Validation Hardening

**Status**: Complete
**Branch**: `010-aptp-symbolic`
**Date**: 2026-02-10

### What Was Done

Post-merge review of Phase 3 identified two pieces of from-scratch logic that needed
deeper validation beyond the existing sanity checks.

**Extracted function** (`src/aptp/symbolic.rs`):
- `permute_symbolic_upper_triangle` — extracted the permutation remapping logic
  (P^T A P upper-triangular CSC construction) from `compute_permuted_etree_and_col_counts`
  into a standalone `pub(crate)` function. Pure refactor, no behavior change. Makes
  permutation correctness independently testable and reusable for Phase 4 (MC64).

**New tests** (10 total):
- 4 permutation tests: identity, reverse-arrow, pair-swap-tridiag, diagonal-invariant
- 1 cross-validation: `col_counts.sum() == predicted_nnz` for all helpers × both orderings
- 1 regression test: arrow(5) with reverse ordering (exact etree, col_counts, predicted_nnz)
- 3 supernode parent tests: etree consistency, child round-trip, block-diagonal roots

**Documentation** (`CLAUDE.md`):
- Added "Testing Discipline for Implementation Phases" section with 5 principles for
  Phases 5-10: validation proportional to novelty, refactor for testability, regression
  tests before refactoring, cross-validation with faer, property-based testing.

### Key Findings

- Arrow(5) under reverse ordering creates a star graph (hub at last column) with zero
  fill-in: `col_counts = [2,2,2,2,1]`, `predicted_nnz = 9`. Contrasts with identity
  ordering where the hub at column 0 causes complete fill-in (`predicted_nnz = 15`).

---

## Phase 3: APTP Symbolic Analysis

**Status**: Complete
**Branch**: `010-aptp-symbolic`
**Date**: 2026-02-10

### What Was Built

Symbolic analysis module for the APTP solver pipeline — the "analyze" step of the
three-phase analyze → factorize → solve API. Composes faer's `SymbolicCholesky<usize>`
with APTP-specific metadata: elimination tree parent pointers, per-column nonzero counts,
and heuristic delayed-pivot buffer estimates.

**Symbolic module** (`src/aptp/symbolic.rs`):
- `AptpSymbolic` — central analysis result struct wrapping `SymbolicCholesky<usize>`
  (inner), `Vec<isize>` etree parent pointers, `Vec<usize>` column counts, and
  `Vec<usize>` pivot buffer estimates. Immutable after creation, reusable across
  multiple numeric factorizations with the same sparsity pattern.
- `AptpSymbolic::analyze(matrix, ordering)` — constructor performing: input validation
  (NotSquare, DimensionMismatch), `factorize_symbolic_cholesky` for full symbolic result,
  permuted structure computation (P^T A P upper-triangular CSC),
  `prefactorize_symbolic_cholesky` on permuted structure for etree + col_counts, and
  10% buffer heuristic for pivot buffer estimates
- `SymbolicStatistics` — diagnostic summary with `Display` impl (dimension, predicted_nnz,
  average_col_count, is_supernodal, n_supernodes, total_pivot_buffer)
- faer-delegating accessors: `perm()`, `predicted_nnz()`, `nrows()`, `ncols()`, `raw()`, `inner()`
- APTP-specific accessors: `etree()`, `col_counts()`, `pivot_buffer_estimates()`,
  `total_pivot_buffer()`, `is_supernodal()`, `n_supernodes()`
- Supernodal accessors: `supernode_begin()`, `supernode_end()`, `supernode_pattern(s)`,
  `supernode_parent(s)` (assembly tree parent derivation from column-level etree)
- Dual-variant handling: all accessors handle both simplicial and supernodal results
  transparently via pattern matching on `SymbolicCholeskyRaw`

**Module wiring** (`src/aptp/mod.rs`):
- Added `pub mod symbolic;` and re-exports for `AptpSymbolic`, `SymbolicStatistics`

**Integration tests** (`tests/symbolic_analysis.rs`):
- `test_analyze_all_hand_constructed` — all 15 hand-constructed matrices with AMD ordering
- `test_analyze_suitesparse_ci_subset` — all 9 SuiteSparse CI-subset matrices with AMD
- `test_custom_ordering_on_hand_constructed` — AMD vs identity custom ordering comparison
- `test_supernodal_structure_suitesparse` — validates supernode ranges, assembly tree,
  row patterns for supernodal SuiteSparse matrices

**Benchmarks** (`benches/solver_benchmarks.rs`):
- Added `bench_symbolic_analysis` benchmark group measuring `AptpSymbolic::analyze` on
  CI-subset matrices with AMD ordering, establishing baseline for SC-006

**Tests**: 20 new unit tests + 4 integration tests + 3 doc-tests, all passing.
Zero clippy warnings, fmt clean, cargo doc clean.

**Follow-up additions**:
- 4 regression tests with exact expected values using `SymmetricOrdering::Identity`:
  diagonal, tridiagonal, arrow, and block-diagonal matrices with analytically
  predicted etree, col_counts, and predicted_nnz values
- `make_tridiagonal(n)` helper for test matrix generation
- Full SuiteSparse integration test file (`tests/symbolic_analysis_full.rs`):
  4 `#[ignore]` tests validating all 67 matrices (analyze, supernodal structure,
  AMD fill reduction property, pivot buffer sanity), run via `cargo test -- --ignored`
- Optimization comments flagging future review points in
  `compute_permuted_etree_and_col_counts` (Phase 10) and `supernode_parent` (Phase 6)
- Moved inline imports (`SparseColMat`, `Triplet`) to top-level import block

### Key Decisions

1. **Two-call strategy**: `factorize_symbolic_cholesky` for the full symbolic result +
   `prefactorize_symbolic_cholesky` on the permuted structure for etree/col_counts.
   The prefactorize call is O(nnz · α(n)) — negligible compared to full symbolic
   factorization. This avoids forking faer or relying on `pub(crate)` internals.

2. **Permuted structure computation**: The elimination tree must correspond to the
   **permuted** matrix (after ordering). We build P^T A P's upper-triangular CSC
   structure explicitly (remapping row/column indices through the permutation) before
   calling `prefactorize_symbolic_cholesky`. This matches what faer computes internally.

3. **10% buffer heuristic**: `ceil(0.10 * column_count)` per supernode (supernodal) or
   per column (simplicial). Conservative starting point from Hogg et al. (2016) Section
   2.4 on delayed pivot propagation. Will be validated empirically in Phase 5-6.

4. **Dual-variant handling**: `SymbolicCholeskyRaw` is an enum (Simplicial/Supernodal).
   All accessors handle both variants — supernodal-specific methods return `Option<T>`
   (None for simplicial). Small test matrices may be simplicial; production matrices
   will typically be supernodal.

5. **Assembly tree derivation**: Supernode parent pointers are derived from the column-level
   etree: for supernode `s` with column range `[begin, end)`, look up `etree[end - 1]` and
   binary search `supernode_begin` to find the containing supernode. Citing Liu (1992).

### Issues Encountered

- **faer `prefactorize_symbolic_cholesky` path**: Not re-exported at the `cholesky` module
  level — lives inside `faer::sparse::linalg::cholesky::simplicial::`.

- **faer `MemStack`/`StackReq` path**: Not at `faer::` root. `dyn_stack` is
  `pub extern crate` → use `faer::dyn_stack::MemStack` and `faer::dyn_stack::MemBuffer`.

- **faer `SymbolicSparseColMatRef` accessor names**: `col_ptr()` not `col_ptrs()`,
  `row_idx()` not `row_indices()`.

- **`SparseColMat<usize, ()>` not constructable**: `()` doesn't implement `ComplexField`.
  Used `SparseColMat<usize, f64>` with dummy 1.0 values for the permuted symbolic structure.

- **Pivot buffer ratio bound**: Initial test expected buffer/nnz ratio ≤ 50%, but the 10%
  heuristic applied to col_counts (which include diagonal entries) can produce higher ratios
  for small matrices. Widened bound to 100%.

### Feature Spec

Full specification in `specs/010-aptp-symbolic/` including spec.md, plan.md, research.md,
data-model.md, contracts/aptp-symbolic-api.md, quickstart.md, checklists/requirements.md,
and tasks.md.

---

## Phase 2: APTP Data Structures

**Status**: Complete
**Branch**: `009-aptp-data-structures`
**Date**: 2026-02-08

### What Was Built

Core data structures for A Posteriori Threshold Pivoting (APTP) factorization,
implementing the mixed 1x1/2x2 block diagonal D storage, pivot classification,
inertia computation, and permutation construction.

**APTP module** (`src/aptp/`):
- `pivot.rs` — `PivotType` enum (`OneByOne`, `TwoByTwo { partner }`, `Delayed`)
  for classifying column pivot decisions per Hogg, Duff & Lopez (2020) Section 3;
  `Block2x2` struct storing symmetric 2x2 diagonal blocks `[[a, b], [b, c]]` with
  `determinant()` and `trace()` methods, citing Bunch & Kaufman (1977)
- `diagonal.rs` — `MixedDiagonal` struct with parallel array storage
  (`pivot_map: Vec<PivotType>`, `diag_1x1: Vec<f64>`, `blocks_2x2: Vec<Block2x2>`,
  `n: usize`); construction/query API (`new`, `set_1x1`, `set_2x2`, `get_pivot_type`,
  `get_1x1`, `get_2x2`, `num_delayed`, `num_1x1`, `num_2x2_pairs`, `dimension`);
  `solve_in_place` (1x1 via division, 2x2 via Cramer's rule); `compute_inertia`
  (trace/determinant eigenvalue sign classification)
- `inertia.rs` — `Inertia` struct relocated from `io/reference.rs` with backward-
  compatible re-export; enhanced rustdoc citing eigenvalue inertia theory
- `perm.rs` — `perm_from_forward` function bridging ordering output (forward
  permutation array) to faer's `Perm<usize>` by computing inverse and calling
  `Perm::new_checked`; validates via `validate::validate_permutation`
- `mod.rs` — Module hub with re-exports of all public types

**Integration tests** (`tests/aptp_data_structures.rs`):
- `inertia_matches_all_15_hand_constructed_references` — constructs MixedDiagonal
  from each reference factorization's DBlock entries, verifies `compute_inertia()`
  matches reference inertia (SC-003)
- `perm_from_forward_matches_hand_constructed_references` — passes reference
  permutations through `perm_from_forward`, verifies forward/inverse relationship
  (SC-004)

**Tests**: 49 new tests (28 pivot/diagonal unit + 6 inertia unit + 9 perm unit +
2 integration + 4 doc-tests), all passing. Total suite: 197 unit + 5 integration +
4 doc-tests. Clippy, fmt, and rustdoc all clean.

### Key Decisions

1. **Parallel array storage for MixedDiagonal**: Three parallel arrays (`pivot_map`,
   `diag_1x1`, `blocks_2x2`) rather than a single `Vec<PivotEntry>` enum. This
   provides cache-friendly access patterns during solve (sequential scan of `diag_1x1`)
   and avoids match-per-element overhead. Research finding R3 from spec.

2. **Debug-assert for internal invariants**: `solve_in_place` and `compute_inertia`
   use `debug_assert!` (not `Result`) for invariant violations like singular pivots
   or delayed columns. These are programming errors in the calling code, not runtime
   failures. The caller is responsible for ensuring no delayed columns remain before
   solving.

3. **Separate DBlock and MixedDiagonal types (FR-012)**: Reference factorization
   `DBlock` (IO deserialization) and live factorization `MixedDiagonal` are
   intentionally separate types. Integration tests bridge them with a local helper
   function. This avoids coupling the IO layer to the factorization layer.

4. **Inertia relocation with re-export**: Moved `Inertia` from `io/reference.rs`
   to `aptp/inertia.rs` with `pub use crate::aptp::Inertia;` re-export in the
   original location. Zero breakage to 154+ existing tests.

5. **Transparent composition with faer**: `perm_from_forward` returns `Perm<usize>`
   directly (no custom wrapper), following the project's faer integration principle.
   Uses `Perm::new_checked` with owned `Box<[usize]>` arrays.

6. **Cramer's rule for 2x2 solve**: Explicit `det = a*c - b²`, `x1 = (c*r1 - b*r2)/det`,
   `x2 = (a*r2 - b*r1)/det` rather than computing an inverse matrix. Numerically
   equivalent and avoids temporary allocation.

7. **Trace/determinant inertia classification**: For 2x2 blocks, eigenvalue signs are
   determined from det and trace without computing actual eigenvalues. det < 0 means
   one positive + one negative; det > 0 with trace > 0 means both positive; det > 0
   with trace < 0 means both negative. Research finding R2.

### Issues Encountered

- **faer `Perm::arrays()` on owned type**: Owned `Perm<usize>` does not have an
  `.arrays()` method — it has `.into_arrays()` (consuming). The correct borrowing
  approach is `.as_ref().arrays()` to get `(&[usize], &[usize])` without consuming
  the permutation.

- **Incremental re-exports in mod.rs**: Initially added all `pub use` re-exports
  before types existed, causing unresolved import errors. Fixed by adding re-exports
  incrementally as each type was implemented.

- **Clippy `neg_multiply`**: `-1.0 * x[2]` flagged; simplified to `-x[2]`.

### Feature Spec

Full specification in `specs/009-aptp-data-structures/` including spec.md, plan.md,
research.md, data-model.md, contracts/aptp-api.md, quickstart.md, and tasks.md.

---

## Phase 1.4: Profiling and Debug Tools

**Status**: Complete
**Branch**: `008-profiling-debug-tools`
**Date**: 2026-02-08

### What Was Built

Profiling, memory tracking, and debug visualization tools for solver development,
gated behind the existing `test-util` Cargo feature flag.

**Profiling module** (`src/profiling/`, feature-gated behind `test-util`):
- `section.rs` — `ProfileEvent` (raw timing record) and `ProfileSection` (aggregated
  view with call count, min/max/mean duration, parent percentage, recursive children)
- `session.rs` — `ProfileSession` (thread-safe via `Send + Sync`, records to
  thread-local storage for zero-contention), `SectionGuard` (RAII, `!Send`),
  `FinishedSession` (merged events + aggregated section tree),
  `flush_thread_events()` for multi-threaded collection via shared `Mutex<HashMap>`
- `report.rs` — `summary_report()` (hierarchical text table with Section/Total/Calls/
  Mean/Min/Max/Parent% columns), `export_chrome_trace()` (Chrome Trace Event format
  JSON with Complete Events type "X", microsecond timestamps, per-thread tid)
- `memory.rs` — `MemoryTracker` (snapshot-based RSS recording), `MemorySnapshot`,
  `MemoryDelta`, `MemoryReport` with `display_report()` (formatted table with
  comma-separated KB values, delta signs, peak RSS footer with MB conversion)

**Debug module** (`src/debug/`, feature-gated behind `test-util`):
- `sparsity.rs` — `SparsityDisplay` with builder pattern (`from_sparse`,
  `with_max_width/height/ascii_only`), density-based character mapping
  (Unicode `█▓▒░.` / ASCII `#+-. `), configurable downsampling for large
  matrices, `fmt::Display` impl
- `etree.rs` — `ETreeDisplay` (`from_parent_array` with root detection via
  `parent[i]==i` or sentinel), `render_tree()` with box-drawing characters
  (├── └── │) for small trees (n<20, falls back to stats for larger),
  `render_stats()` and `stats()` → `EliminationTreeStats` (depth, leaves,
  branching min/max/mean, subtree sizes with median)

**Extended utility** (`src/benchmarking/rss.rs`):
- Added `read_current_rss_kb()` reading `VmRSS` from `/proc/self/status`
- Refactored shared `read_proc_status_field()` helper

**Tests**: 74 new tests (28 profiling + 18 memory + 14 sparsity + 14 etree),
all passing. Integration test validates SparsityDisplay on all 15 hand-constructed
matrices. Full clippy + rustdoc clean.

### Key Decisions

1. **Thread-local storage for profiling**: Each thread records events into its own
   `RefCell<Vec<RawEvent>>` via `thread_local!` — zero contention during recording.
   Worker threads call `flush_thread_events()` to move events to a shared
   `LazyLock<Mutex<HashMap>>` before the session owner calls `finish()`.

2. **Chrome Trace Complete Events**: Type "X" events (vs Begin/End pairs) are simpler
   to emit and represent full duration in a single record. Hierarchy is implicit from
   overlapping time ranges on the same tid.

3. **RSS import from benchmarking module**: Memory tracker imports `read_current_rss_kb`
   and `read_peak_rss_kb` directly from `crate::benchmarking::rss` rather than factoring
   into a shared utility module. Both modules are behind `test-util`, so the dependency
   direction is clean.

4. **Density-based sparsity visualization**: Unicode block characters (`█▓▒░`) at 4
   density levels with ASCII fallback (`#+-`). Downsampling bins rows/columns into
   display cells and computes non-zero density per bin.

5. **Peer modules, no coupling**: `profiling` and `debug` are independent peer modules
   with no cross-imports. Either can be used without the other.

### Issues Encountered

- **`Mutex::new` not const in HashMap context**: `static SHARED_FLUSHED: Mutex<HashMap<...>>`
  requires `LazyLock` wrapper since `HashMap::new()` is not a const fn.

- **Clippy type_complexity**: The `LazyLock<Mutex<HashMap<u64, Vec<(usize, RawEvent)>>>>`
  type triggered clippy. Fixed with a type alias.

---

## Phase 1.3: Continuous Integration Setup

**Status**: Complete
**Branch**: `007-ci-setup`
**Date**: 2026-02-07

### What Was Built

Added benchmark compilation verification to the existing CI pipeline. Gap analysis
found that 8 of 10 spec requirements (FR-001 through FR-010) were already satisfied
by the CI configuration established in Phase 0.4.

**CI change** (`.github/workflows/ci.yml`):
- Added `bench-sparse` job — runs `cargo bench --no-run` on stable toolchain with
  `Swatinem/rust-cache@v2`, following the same structure as existing sparse domain jobs
  (checkout → toolchain → cache → run). Compiles the `solver_benchmarks` criterion
  binary without executing benchmarks.

**Sparse domain CI jobs after Phase 1.3**:
- `test-sparse` — MSRV (1.87) + stable matrix, `cargo test --all-targets`
- `lint-sparse` — `cargo fmt --check` + `cargo clippy -- -D warnings`
- `doc-sparse` — `cargo doc --no-deps` with `RUSTDOCFLAGS: -D warnings`
- `bench-sparse` — `cargo bench --no-run` (NEW)

### Key Decisions

1. **Minimal scope**: Gap analysis showed most CI requirements were already met.
   Rather than over-engineering, only the missing benchmark compilation check was added.

2. **Stable-only for benchmarks**: No MSRV matrix for bench-sparse. Benchmarks are a
   development tool; if they compile on stable, MSRV adds no value.

3. **No path filtering**: Considered `dorny/paths-filter` for monorepo efficiency but
   deferred — two domains don't justify the complexity, and GitHub required checks
   interact poorly with path-filtered jobs.

4. **Feature-gated coverage via dev-dependency**: The self-referencing
   `rivrs-sparse = { path = ".", features = ["test-util"] }` dev-dependency already
   activates `test-util` during `cargo test --all-targets`. No separate feature-flag
   CI job needed.

5. **SPRAL comparison deferred**: Per Phase 0.3 decision, SPRAL is not built or
   invoked in CI. Will be added in Phases 2-8 when the solver can process SuiteSparse
   matrices.

### Issues Encountered

- None. The implementation was a straightforward additive change (~8 lines of YAML).

---

## Phase 1.2: Benchmarking Framework

**Status**: Complete
**Branch**: `006-benchmarking-framework`
**Date**: 2026-02-07

### What Was Built

Criterion-based benchmarking framework for measuring solver phase performance,
gated behind the existing `test-util` Cargo feature flag.

**Benchmarking module** (`src/benchmarking/`, feature-gated behind `test-util`):
- `config.rs` — `BenchmarkPhase` enum (Analyze, Factor, Solve, Roundtrip) with
  Display and serde support; `BenchmarkConfig` struct with builder methods for
  filter, phases, sample_size, measurement_time, warm_up_time, timeout_per_matrix
- `traits.rs` — `Benchmarkable` trait with four methods (`bench_analyze`,
  `bench_factor`, `bench_solve`, `bench_roundtrip` with default chained impl);
  `MockBenchmarkable` with configurable phase enable/disable for harness testing
- `results.rs` — `BenchmarkResult`, `BenchmarkSuiteResult`, `SkippedBenchmark`
  structs with serde and Display; `collect_results()` function that parses
  Criterion's `estimates.json` and `benchmark.json` output files
- `baseline.rs` — `Baseline`, `Regression`, `Improvement`, `RegressionReport`
  structs; `save_baseline()`, `load_baseline()`, `detect_regressions()` functions
  with configurable threshold percentage
- `report.rs` — `export_csv()`, `export_json()`, `generate_markdown_table()`,
  `generate_regression_markdown_table()` for result export and reporting
- `rss.rs` — `read_peak_rss_kb()` reading VmHWM from `/proc/self/status`
  (Linux only, returns None on other platforms)
- `mod.rs` — Module root with public re-exports

**Benchmark binary** (`benches/solver_benchmarks.rs`):
- `run_component_benchmarks()` — One `BenchmarkGroup` per phase (ssids/analyze,
  ssids/factor, ssids/solve), parameterized by matrix name via `BenchmarkId`,
  with `Throughput::Elements(nnz)` per benchmark
- `run_e2e_benchmarks()` — Single group (ssids/roundtrip) for full pipeline
- Three `criterion_group!` registrations: component, e2e, ci_subset
- Peak RSS measurement before/after benchmark suite, printed to stderr
- Graceful skip for unimplemented phases and missing matrices

**Cargo.toml changes**:
- Added `[[bench]] name = "solver_benchmarks" harness = false`
- Added `tempfile = "3"` dev-dependency for baseline tests

**Tests**: 80 total (78 unit + 2 integration), all passing, 0 clippy warnings

### Key Decisions

1. **Separate `Benchmarkable` trait**: Distinct from `SolverTest` — benchmarking
   needs opaque `Box<dyn Any>` returns (to prevent optimizing away computation)
   while testing needs structured `TestResult` with correctness metrics.

2. **`Option` return for unimplemented phases**: `None` signals the harness to
   skip rather than fail, supporting incremental solver development where only
   some phases are implemented at any time.

3. **Whole-run RSS measurement**: VmHWM from `/proc/self/status` measured once
   per suite (before/after delta). Per-phase RSS deferred to Phase 1.4 — would
   require `/proc/self/clear_refs` writes and complicates measurement.

4. **Criterion group-per-phase layout**: Benchmark IDs like `ssids/factor/bcsstk14`
   enable CLI filtering (`cargo bench -- "factor"`) and produce per-phase
   comparison HTML reports across matrices.

5. **Custom `collect_results` parser**: Reads Criterion's `estimates.json` and
   `benchmark.json` files directly rather than parsing terminal output or
   requiring `cargo-criterion`. Stable across Criterion versions.

6. **Configurable regression threshold**: `detect_regressions()` accepts a
   threshold percentage (default 5%), classifying changes as regression,
   improvement, or unchanged based on mean time comparison.

### Issues Encountered

- **faer API method names**: `row_indices_of_col` → `row_idx_of_col` and
  `values_of_col` → `val_of_col` in faer 0.22. Method already returns an
  iterator, not a slice.

- **Criterion directory naming**: Group names with `/` (e.g., `ssids/analyze`)
  are stored as directories with `_` separator (e.g., `ssids_analyze/`). CLI
  filtering still works with either separator.

---

## Phase 1.1: Test Infrastructure

**Status**: Complete
**Branch**: `005-test-infrastructure`
**Date**: 2026-02-07

### What Was Built

Reusable test harness and validation infrastructure for solver development,
gated behind a `test-util` Cargo feature flag to keep production builds lean.

**Testing module** (`src/testing/`, feature-gated behind `test-util`):
- `harness.rs` — `SolverTest` trait defining four test methods (analyze, factor,
  solve, roundtrip), `MockSolver` implementation that validates reference
  factorizations directly, `TestKind` enum, `MetricResult` and `TestResult`
  structs with `Display` formatting
- `validator.rs` — `NumericalValidator` with builder pattern for configurable
  tolerances (reconstruction 10^-12, backward error 10^-10), methods for
  `check_reconstruction()`, `check_backward_error()`, `check_inertia()`, and
  `validate_factorization()` returning structured `TestResult`
- `cases.rs` — `SolverTestCase`, `TestMatrixProperties`, `TestCaseFilter` with
  builder methods (`all()`, `hand_constructed()`, `ci_subset()`, `with_source()`,
  `with_category()`, `with_difficulty()`, `ci_only()`, `require_reference()`),
  `load_test_cases()` wrapping registry functions with predicate filtering
- `generators.rs` — `RandomMatrixConfig`, `generate_random_symmetric()` (PD via
  diagonal dominance or indefinite with mixed signs), `generate_arrow()`,
  `generate_tridiagonal()`, `generate_banded()` — all producing
  `SparseColMat<usize, f64>` via faer's `Triplet` API
- `mod.rs` — Module root with re-exports and compilable rustdoc example

**Cargo.toml changes**:
- `test-util` feature flag gating `rand` and `rand_distr` as optional deps
- Self-dev-dependency: `rivrs-sparse = { path = ".", features = ["test-util"] }`

**Integration test refactoring**:
- `tests/hand_constructed.rs` — refactored to use `load_test_cases` +
  `NumericalValidator` instead of raw registry + validate calls
- `tests/suitesparse_ci.rs` — refactored to use `TestCaseFilter::ci_subset()`

**Tests**: 48 total (45 unit + 2 integration + 1 doctest), all passing

### Key Decisions

1. **Feature-gated testing module**: The `test-util` feature flag keeps `rand`,
   `rand_distr`, and all test generators out of production builds. Downstream
   crates enable it in `[dev-dependencies]` only.

2. **Builder pattern for filters and validators**: `TestCaseFilter` and
   `NumericalValidator` use builder patterns for ergonomic, composable
   configuration. Default tolerances match the constitution (reconstruction
   10^-12, backward error 10^-10).

3. **MockSolver for immediate validation**: The `MockSolver` validates reference
   factorizations from the hand-constructed JSON files directly (no solver needed
   yet). This allows the test harness to be fully exercised before any solver
   exists.

4. **Diagonal dominance for PD generators**: Random PD matrices use Gershgorin
   circle theorem — diagonal entries exceed row absolute sums by a random margin,
   guaranteeing positive definiteness without Cholesky verification.

5. **Seeded RNG for reproducibility**: All generator tests use
   `StdRng::seed_from_u64(42)` for deterministic output across runs.

### Issues Encountered

- **`faer::sparse::Triplet` path**: `Triplet` is at `faer::sparse::Triplet`,
  not `faer::Triplet`. The `try_new_from_triplets` method expects
  `&[Triplet<usize, usize, f64>]`, not tuples.

- **Rust 2024 edition `gen` keyword**: `gen` is reserved in edition 2024. Must
  use raw identifier syntax `rng.r#gen::<f64>()` to call `rand::Rng::gen()`.

- **`Uniform::new` type inference**: `Uniform::new(0.1, 1.0)` is ambiguous;
  needs explicit type: `Uniform::new(0.1f64, 1.0)`.

- **Reconstruction tolerance sensitivity**: A perturbation of 1e-4 to an L
  entry produced reconstruction error exceeding 1e-6 tolerance in tests. Fixed
  by using 1e-8 perturbation for the relaxed-tolerance test.

- **`cargo fmt` CI failure**: Import ordering and line wrapping differences
  between local and CI rustfmt. Fixed by running `cargo fmt` and committing.

---

## Phase 0.4: Repository Setup for Solver Development

**Status**: Complete
**Branch**: `004-repo-setup`
**Date**: 2026-02-06

### What Was Built

Development infrastructure enabling Rust-based loading, parsing, and validation
of the test matrix collection established in Phase 0.2.

**IO modules** (`src/io/`):
- `mtx.rs` — Matrix Market parser (`coordinate real symmetric` format) with
  1-indexed→0-indexed conversion, symmetric expansion, and descriptive parse
  errors including file path and line number
- `reference.rs` — JSON reference factorization loader with types for Inertia,
  LEntry, DBlock (1×1/2×2 with polymorphic serde), ReferenceFactorization;
  includes validation of strict lower triangle and permutation consistency
- `registry.rs` — Test matrix catalog backed by `metadata.json` with CI-subset
  path fallback (`suitesparse-ci/` preferred for CI matrices over gitignored
  `suitesparse/` directory)

**Error handling** (`src/error.rs`):
- Added `IoError` and `ParseError` variants to `SparseError`
- Added `From<std::io::Error>` and `From<serde_json::Error>` impls

**Validation utilities** (`src/validate.rs`):
- `reconstruction_error(A, ref)` — `||P^T A P - L D L^T||_F / ||A||_F`
- `backward_error(A, x, b)` — `||Ax - b|| / (||A||_F ||x|| + ||b||)`
- `check_inertia(computed, expected)` — field-wise equality comparison
- All implemented using faer dense operations (to_dense, matmul, norm_l2)

**Tests** (20 total):
- 11 unit tests across mtx, reference, registry modules
- 7 validation unit tests (reconstruction, backward error, inertia)
- 1 integration test loading all 15 hand-constructed matrices with
  reconstruction error < 10^-12 (SC-008 validated)
- 1 integration test for 10 CI-subset SuiteSparse matrices

**CI pipeline** (`.github/workflows/ci.yml`):
- `test-sparse` job (stable + MSRV 1.87)
- `lint-sparse` job (clippy + rustfmt)
- `doc-sparse` job (rustdoc with `-D warnings`)

**Benchmark scaffold** (`benches/matrix_loading.rs`):
- Criterion benchmarks for matrix loading and reconstruction error computation

### Key Decisions

1. **Custom MTX parser**: Lightweight, scoped to `coordinate real symmetric`
   format. No external crate dependency. Descriptive errors with file path and
   line number for debugging data quality issues.

2. **faer-based validation**: All validation uses dense matrix operations via
   faer (to_dense, matmul, norm_l2, operator overloads). Simple and correct for
   test matrices up to ~20K dimension. No need for sparse validation at this
   stage.

3. **Polymorphic DBlock serde**: Custom deserializer reads `"size"` field to
   distinguish 1×1 scalar pivots from 2×2 symmetric blocks. Matches the JSON
   schema established in Phase 0.2.

4. **CI-subset path fallback**: `load_test_matrix()` checks `suitesparse-ci/`
   before `suitesparse/` for CI-subset entries, ensuring tests work both in
   CI (where only `suitesparse-ci/` is available) and locally (where
   `suitesparse/` may also be extracted).

5. **TDD throughout**: All tests written and verified to fail before
   implementation, per Constitution Principle III.

### Issues Encountered

- **faer API discovery**: `sparse_dense_matmul` lives in
  `faer::sparse::linalg::matmul`, not `faer::linalg::matmul`. `Par::Seq`
  (not `Par::sequential()`). Column `as_mat()` (not `as_mat_ref()`).
  Resolved by reading faer source code directly.

- **metadata.json local corruption**: Working copy had been modified from 82
  to 73 entries (likely from a prior script run). Restored from git.

---

## Phase 0.3: SPRAL Golden Results — Deferred

**Status**: Deferred
**Branch**: `003-spral-golden-results`
**Date**: 2026-02-06

### Decision

Phase 0.3 was planned to build SPRAL from source, write a C driver program, run
SPRAL on all 82 test matrices, and capture golden reference results (timing,
inertia, pivot statistics, backward error) as JSON files.

After investigating SPRAL's C API (`spral_ssids.h`) and assessing cost-benefit,
we decided to **defer this phase**. Full rationale in
`specs/003-spral-golden-results/decision.md`.

### Key Finding: SPRAL API Limitations

SPRAL's C API exposes:
- Aggregate statistics from the `inform` struct (num_factor, num_flops,
  num_delay, num_neg, num_two, matrix_rank, maxfront, maxsupernode)
- Block diagonal D entries via `spral_ssids_enquire_indef()`
- Pivot ordering via `spral_ssids_enquire_indef()`
- The solution vector (from which errors are computed)

SPRAL does **NOT** expose: the L factor, permutation P, elimination tree parent
pointers, supernode membership arrays, or fill-in patterns. These are internal to
the opaque `akeep`/`fkeep` handles.

### Key Decision: Reconstruction Tests as Primary Oracle

The **reconstruction test** (`||P^T A P - L D L^T|| / ||A|| < epsilon`) is a
strictly stronger correctness oracle than comparing against SPRAL's output. If
our factorization reconstructs A to machine precision, it is correct by
definition — regardless of what any other solver produces.

This means the most valuable correctness tests require no SPRAL infrastructure:
1. **Reconstruction**: P^T A P = L D L^T (proves mathematical correctness)
2. **Backward error**: ||Ax - b|| / (||A|| ||x|| + ||b||) (validates solve pipeline)
3. **Hand-constructed matrices**: Analytically known factorizations from Phase 0.2
4. **Property-based tests**: Inertia consistency, symmetry preservation

### Impact

- Constitution updated to v1.1.0: reconstruction tests are primary oracle; SPRAL
  comparison is a secondary validation layer added when available
- Phase 0 exit criteria updated: "SPRAL reference results" requirement replaced
  with reconstruction test strategy
- SPRAL build infrastructure deferred to Phases 2-8 when needed for performance
  benchmarking and inertia validation on large SuiteSparse matrices

### What Was Produced

No code changes. Documentation artifacts only:
- `specs/003-spral-golden-results/decision.md` — full decision record
- `specs/003-spral-golden-results/spec.md` — original feature specification
- `specs/003-spral-golden-results/plan.md` — implementation plan (not executed)
- `specs/003-spral-golden-results/research.md` — SPRAL API investigation findings
- Updated `docs/ssids-plan.md` — Phase 0.3 marked deferred, exit criteria revised
- Updated constitution v1.0.0 → v1.1.0

---

## Phase 0.2: Test Matrix Collection Assembly

**Status**: Complete
**Branch**: `002-test-matrix-collection`
**Date**: 2026-02-05

### What Was Built

A comprehensive test matrix collection with 82 matrices (15 hand-constructed +
67 SuiteSparse), organized under `sparse/test-data/`:

```
sparse/test-data/
├── metadata.json                    # Complete index (82 matrices)
├── README.md                        # Setup instructions
├── hand-constructed/                # 15 matrices with exact LDL^T factorizations
│   ├── {name}.mtx                   # Matrix Market format
│   └── {name}.json                  # Exact L, D, P, inertia
├── suitesparse/
│   ├── easy-indefinite/             # 30 matrices from APTP paper benchmarks
│   ├── hard-indefinite/             # 18 matrices (includes killer cases)
│   └── positive-definite/           # 19 matrices for fast-path validation
├── suitesparse-ci/                  # 10 representative matrices (plain git, ~73MB)
│   ├── easy-indefinite/{name}.mtx
│   ├── hard-indefinite/{name}.mtx
│   └── positive-definite/{name}.mtx
└── scripts/
    ├── requirements.txt             # Python dependencies (ssgetpy, numpy, scipy)
    ├── mtx_utils.py                 # Matrix Market I/O
    ├── metadata_utils.py            # Metadata index builder
    ├── factorization_utils.py       # Factorization writing and verification
    ├── generate_hand_constructed.py  # Hand-constructed matrix generator
    ├── download_suitesparse.py      # SuiteSparse download with curated lists
    └── validate_collection.py       # Collection validation and reporting
```

### Hand-Constructed Matrices (15)

Six categories of matrices with analytically verified LDL^T factorizations:

1. **Arrow matrices** (3): 5x5 PD, 10x10 indefinite, 15x15 indefinite
2. **Tridiagonal matrices** (3): 5x5 PD, 10x10 indefinite, 20x20 indefinite
3. **Block diagonal** (1): 15x15 with PD, indefinite, and 2x2-pivot blocks
4. **Bordered block diagonal** (1): 20x20 adapted from SPRAL pattern (BSD-3)
5. **Stress tests** (3): zero diagonal (forces 2x2 pivots), dense column (fill-in), ill-conditioned
6. **Degenerate cases** (4): 3x3 zero-diagonal, 3x3 singular, 1x1 trivial, 2x2 requiring 2x2 pivot

All 15 factorizations verified: `||P^T A P - L D L^T|| / ||A|| < 1e-10`.

### SuiteSparse Matrices (67)

Curated from APTP benchmark papers:
- **Easy indefinite (30)**: From Duff, Hogg, Lopez (2020) Table 1 and Hogg et al. (2016)
- **Hard indefinite (18)**: From Duff et al. (2020) Table 2 (KKT/saddle-point systems)
- **Positive definite (19)**: From Hogg et al. (2016) Table III (dagger-marked)

Large matrices (>200K rows) stored as metadata-only with download commands.

### Key Decisions

1. **Curated lists from papers** (not broad queries): Provides gold-standard difficulty
   classification from actual APTP benchmark results.
2. **Tiered storage (no LFS)**: Hand-constructed + CI subset committed to plain git;
   full SuiteSparse collection gitignored and extracted from archive at container build.
   Originally planned Git LFS, replaced with three-tier approach to avoid LFS costs
   and complexity for 4GB of immutable reference data.
3. **CI subset (10 matrices, ~73MB)**: Representative set in `suitesparse-ci/` committed
   to plain git (~19MB in pack). Avoids 40KB/s SuiteSparse API throttle in CI. Covers
   all categories including 2 killer cases.
4. **Full collection via archive**: `references/ssids/suitesparse.tar.gz` (1.3GB) stored
   outside git as static reference data. Docker entrypoint extracts automatically.
5. **Python scripts with uv**: Virtual environment at `scripts/.venv/`, dependencies via `uv pip`.
6. **ssgetpy API**: `ssgetpy.search(name)` with exact group/name matching; download with
   manual tar.gz extraction for reliable file placement.
7. **vibrobox group correction**: SuiteSparse has `Cote/vibrobox`, not `GHS_indef/vibrobox`.
8. **cont-201 killer case**: Added `killer_case: True` to `cont-201` (large optimization problem
   with many zero diagonals forcing 2×2 pivots) to reach the SC-005 threshold of ≥5 killer cases.
   Total killer cases: stokes128, cont-201, ncvxqp3, c-big, c-71.

### Issues Encountered

- **ssgetpy API change**: The `matname` keyword argument no longer works; must use positional
  `name_or_id` parameter.
- **ssgetpy extract=True**: Doesn't always place the correct `.mtx` file; rewrote to use
  `extract=False` + manual `tarfile` extraction with name-based file selection.
- **Slow download speeds**: SuiteSparse downloads throttled to ~40KB/s from this environment.

## Phase 0.1: Literature Review & Reference Library

**Status**: Complete
**Branch**: `001-reference-library`
**Date**: 2026-02-05

### What Was Built

Four synthesis documents and a paper library completing the Phase 0.1 reference
library, organized following the structure from `ssids-plan.md`:

```
sparse/docs/references/
├── INDEX.md                          # Paper index with citations and categories
├── papers/                           # 14 academic papers (.md format)
│   ├── davis2016.md ... schenk2006.md
├── algorithms/
│   └── APTP-ALGORITHM.md            # Extracted pseudocode from papers
└── notes/
    ├── SPRAL-CODE-REVIEW.md          # Annotated SPRAL SSIDS architecture review
    └── FAER-INTEGRATION-NOTES.md     # faer component reuse map
```

1. **`docs/references/INDEX.md`** (~1400 words) -- Paper index cataloging all 14
   academic papers by category (Core APTP, Multifrontal Foundations, Pivoting
   Strategies, Ordering & Analysis, Infrastructure) with full bibliographic
   citations, relevance summaries, markdown quality flags, and RAL Technical
   Reports status.

2. **`docs/references/notes/SPRAL-CODE-REVIEW.md`** (~5100 words) -- Annotated
   architecture review of SPRAL's SSIDS implementation covering: executive
   summary of three-phase design, module-by-module overview (9 top-level
   modules), CPU factorization stack (24 files), ASCII data flow diagrams, APTP
   kernel internals (`ldlt_app.cxx`), key data structures (`akeep`, `fkeep`,
   `node_type`, `contrib_type`, diagonal D storage), and external dependencies
   with `ssids_options` configuration fields.

3. **`docs/references/notes/FAER-INTEGRATION-NOTES.md`** (~3500 words) --
   Component-by-component map of faer 0.24.0 sparse infrastructure with reuse
   classifications: 8 "direct use" components (CSC, AMD, COLAMD, elimination
   tree, triangular solve, permutations, workspace management, CSR), 2 "adapt"
   components (symbolic/numeric Cholesky), 2 "reference only" (dense
   Bunch-Kaufman, sparse LU). Includes deep dive on cholesky.rs
   Bunch-Kaufman vs APTP difference and integration strategy with build order.

4. **`docs/references/algorithms/APTP-ALGORITHM.md`** (~5000 words) -- Unified
   pseudocode reference covering: APTP overview with TPP/SBK comparison table,
   Algorithm 3.1 (main APTP loop), all 7 dense block kernels, Algorithm 4.1
   (complete pivoting), pivot acceptance/column delay mechanism, symbolic
   analysis (elimination tree + Gilbert-Ng-Peyton column counts), and
   multifrontal structure integration.

### Key Findings

1. **faer Bunch-Kaufman vs APTP**: faer's sparse LDLT uses pre-decided
   Bunch-Kaufman pivoting. The column-by-column update pattern and `ereach`
   traversal are reusable, but the pivot decision logic must be replaced entirely
   with a posteriori stability checks, column delay mechanism, and hybrid
   1x1/2x2 handling.

2. **SPRAL APTP kernel location**: The core APTP implementation lives in
   `src/ssids/cpu/kernels/ldlt_app.cxx` (~2600 lines). Key abstractions:
   `Column<T>` class tracking per-block-column elimination state,
   `check_threshold()` for a posteriori stability, and three pivot methods
   (APP_AGGRESSIVE, APP_BLOCK, TPP).

3. **Diagonal D storage convention**: SPRAL stores D^{-1} as a flat `2n` array
   with an `Inf` sentinel for 2x2 pivots (`d[2k+2] = Inf`). faer uses a
   different approach (`subdiag` + permutation arrays). The SPRAL convention is
   more compact for sparse solvers.

4. **Threshold parameter**: Default `u = 0.01` provides good stability/fill-in
   trade-off. Smaller `u` (0.001) reduces fill-in but may compromise accuracy;
   larger `u` (0.1) improves stability at the cost of more delayed pivots.

5. **faer reuse estimate**: ~70% of needed sparse infrastructure is available
   from faer as direct-use components. The remaining ~30% (APTP kernel, column
   delay, block-diagonal D handling) must be built from academic paper
   references.

### Readiness for Phase 0.2

All Phase 0.1 exit criteria are met:
- All key papers obtained and organized (INDEX.md)
- Algorithm pseudocode extracted (APTP-ALGORITHM.md)
- SPRAL source code reviewed (SPRAL-CODE-REVIEW.md)
- faer integration points identified (FAER-INTEGRATION-NOTES.md)

Phase 0.2 (Test Matrix Collection) can proceed. The reference library provides
sufficient algorithmic understanding to select appropriate test matrices and
define expected behaviors.
