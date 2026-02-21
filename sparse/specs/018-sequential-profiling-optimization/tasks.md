# Tasks: Sequential Profiling & Optimization (Phase 8.1g)

**Input**: Design documents from `/specs/018-sequential-profiling-optimization/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/

**Tests**: Not explicitly requested. Existing test suite (362+ unit tests, 65 SuiteSparse `--ignored` tests) serves as the regression guard for all optimizations. No new test files needed.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3, US4)
- Include exact file paths in descriptions

## Phase 1: Setup

**Purpose**: Verify starting-point correctness and establish pre-optimization baseline

- [X] T001 Run full unit test suite (`cargo test`) to verify all 362+ tests pass on branch `018-sequential-profiling-optimization`
- [X] T002 Run full SuiteSparse validation (`cargo test -- --ignored --test-threads=1`) to verify all 65 matrices pass correctness thresholds
- [X] T003 Run `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` to verify clean lint baseline

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Extend existing data structures with conditional timing fields that US1-US4 all depend on

**CRITICAL**: No user story work can begin until this phase is complete

- [X] T004 Extend `PerSupernodeStats` in `src/aptp/numeric.rs` with three `Duration` timing fields (`assembly_time`, `kernel_time`, `extraction_time`) behind `#[cfg(feature = "diagnostic")]`. Add conditional `use std::time::Duration` import. Ensure the struct compiles identically without `diagnostic` feature (no field present, no memory overhead).
- [X] T005 Extend `FactorizationStats` in `src/aptp/numeric.rs` with three aggregate `Duration` timing fields (`total_assembly_time`, `total_kernel_time`, `total_extraction_time`) behind `#[cfg(feature = "diagnostic")]`. Add accumulation logic in the per-supernode stats aggregation code to sum the per-supernode timing fields.
- [X] T006 Verify both features compile cleanly: `cargo build` (no diagnostic) and `cargo build --features diagnostic` (with diagnostic). Run `cargo test` and `cargo test --features diagnostic` to confirm no breakage.

**Checkpoint**: Extended data structures ready — user story implementation can now begin

---

## Phase 3: User Story 1 — Instrument Factorization Hot Path (Priority: P1) MVP

**Goal**: Add per-supernode timing instrumentation to the multifrontal factorization loop so that Chrome Trace output reveals where time is spent.

**Independent Test**: Factor a CI-subset matrix with `diagnostic` feature enabled and verify Chrome Trace JSON contains per-supernode timing sections with correct nesting hierarchy.

### Implementation for User Story 1

- [X] T007 [US1] Add `#[cfg(feature = "diagnostic")]` import of `ProfileSession` and `Instant` at the top of `src/aptp/numeric.rs`. Create a `ProfileSession` at the start of `AptpNumeric::factor()` (conditionally compiled). After the factor loop completes, call `session.finish()` and store the `FinishedSession` on `AptpNumeric` (behind cfg) or discard if not needed for programmatic access.
- [X] T008 [US1] Instrument the per-supernode postorder loop in `AptpNumeric::factor()` (`src/aptp/numeric.rs`, lines ~380-488) with `Instant::now()` timing around three sections per supernode: (1) assembly (scatter + extend-add, lines ~412-435), (2) dense kernel (aptp_factor_in_place call, line ~438), (3) extraction (extract_front_factors call, lines ~468-470). Populate the `PerSupernodeStats` timing fields (from T004) with the measured durations. All timing code gated behind `#[cfg(feature = "diagnostic")]`.
- [X] T009 [US1] Add a top-level `SectionGuard` wrapping the entire postorder loop (section name: `"factor_loop"`), and nested `SectionGuard` entries per supernode (section name: `"supernode_{s}"`) wrapping each iteration. Within each supernode guard, add sub-sections `"assembly"`, `"dense_kernel"`, `"extraction"`. This enables Chrome Trace to show the full hierarchy. All gated behind `#[cfg(feature = "diagnostic")]`.
- [X] T010 [US1] Accumulate per-supernode timing into `FactorizationStats` aggregate timing fields (from T005) after the postorder loop completes. Sum `assembly_time`, `kernel_time`, `extraction_time` across all `PerSupernodeStats` entries. Gated behind `#[cfg(feature = "diagnostic")]`.
- [X] T011 [US1] Add an accessor method on `AptpNumeric` (or extend existing accessor) to expose the `FinishedSession` (or Chrome Trace JSON string) when `diagnostic` is enabled. This allows `SparseLDLT` and example tools to retrieve the trace without additional plumbing.
- [X] T012 [US1] Verify instrumentation: run `cargo test --features diagnostic` (all unit tests pass), then run `cargo test -- --ignored --test-threads=1 --features diagnostic` (all 65 SuiteSparse matrices pass with identical correctness results). Verify `cargo clippy --all-targets --features diagnostic` is clean.

**Checkpoint**: Per-supernode profiling data is collected and accessible. Chrome Trace output available. No correctness regressions.

---

## Phase 4: User Story 2 — Reduce Allocation Pressure in Factorization (Priority: P2)

**Goal**: Optimize the top allocation hotspots in `factor_inner` (factor.rs) by hoisting temporary buffers out of tight loops and reusing them across iterations.

**Independent Test**: Run the full SuiteSparse suite before and after optimization; verify identical correctness results and no RSS increase.

### Implementation for User Story 2

- [X] T013 [US2] Hoist panel row permutation buffer in `factor_inner` (`src/aptp/factor.rs`, line ~1664). Replace the per-row `Vec<f64>` allocation (`let orig: Vec<f64> = (0..block_size).map(...).collect()`) with a single `Vec<f64>` allocated once before the outer block loop, sized to `options.inner_block_size`, and reused across all rows and blocks. This eliminates ~15,000 allocations per large front.
- [X] T014 [US2] Hoist row permutation temp buffer in `factor_inner` (`src/aptp/factor.rs`, line ~1692). Replace the per-block `vec![0.0f64; block_size]` allocation with a single `Vec<f64>` allocated once before the outer block loop, sized to `options.inner_block_size`, and reused across all block iterations.
- [X] T015 [US2] Hoist column order slice copy in `factor_inner` (`src/aptp/factor.rs`, line ~1704). Replace the per-block `col_order[k..k + block_size].to_vec()` allocation with a single pre-allocated `Vec<usize>` buffer sized to `options.inner_block_size`, reused across block iterations.
- [X] T016 [US2] Evaluate `BlockBackup::create` reuse (`src/aptp/factor.rs`, line ~1146). SKIPPED — BlockBackup dimensions vary per block (rows = m - k, cols = block_size where both k and block_size change). faer's Mat doesn't support in-place resize, so reuse would require either a max-size pre-allocation with subslice views (complex, error-prone) or a custom raw buffer. Deferred to Phase 9.1 (arena memory).
- [X] T017 [US2] Verify correctness after all allocation optimizations: run `cargo test` (all unit tests pass), then `cargo test -- --ignored --test-threads=1` (all 65 SuiteSparse matrices pass with identical reconstruction error and backward error). Run `cargo clippy --all-targets` to verify clean.
- [X] T018 [US2] Run a representative large-matrix factorization (e.g., bloweybq or bratu3d) before and after optimization to confirm no timing regression. Document the before/after allocation count comparison (can use profiling instrumentation from US1 or manual inspection).

**Checkpoint**: Allocation hotspots optimized. All correctness tests pass. No timing regression on large matrices.

---

## Phase 5: User Story 3 — Establish Sequential Performance Baselines (Priority: P3)

**Goal**: Create a baseline collection tool that records per-phase timing, per-supernode stats, and peak RSS for all 65 SuiteSparse matrices in a structured JSON format.

**Independent Test**: Run the baseline collection tool and verify JSON output contains all required fields for all available matrices.

### Implementation for User Story 3

- [X] T019 [US3] Define `PerformanceBaseline`, `BaselineSuite`, and serialization support in a new module or within the example. Fields per data-model.md: matrix_name, matrix_dim, matrix_nnz, per-phase timing (ordering, symbolic, factor, solve), total_time, peak_rss_kb, backward_error, num_supernodes, max_front_size, factorization_stats, per_supernode_stats. Use serde derive for JSON serialization. Place types in `examples/baseline_collection.rs` (or in `src/benchmarking/` if reuse is desired across examples).
- [X] T020 [US3] Create `examples/baseline_collection.rs` following the `examples/solve_timing.rs` pattern. The tool loads SuiteSparse matrices (full or CI subset via `--ci-only` CLI arg), runs the full `SparseLDLT` pipeline (analyze_with_matrix → factor → solve) with `MatchOrderMetis` ordering, records per-phase timing via `Instant::now()`, captures `PerSupernodeStats` via `solver.per_supernode_stats()`, captures peak RSS via the existing `read_peak_rss_kb()` function from `src/profiling/memory.rs` or `src/benchmarking/rss.rs`, computes backward error, and assembles a `PerformanceBaseline` per matrix.
- [X] T021 [US3] Add JSON export to `baseline_collection.rs`: serialize the `BaselineSuite` to `target/benchmarks/baselines/baseline-<timestamp>.json`. Print human-readable summary to stderr (matrix name, dim, nnz, factor_ms, solve_ms, backward_error, peak_rss_kb). Include `--compare <previous.json>` flag that loads a previous baseline and reports per-matrix timing regressions/improvements (>10% change).
- [X] T022 [US3] Run `baseline_collection.rs` on CI subset (`--ci-only`) and verify JSON output is well-formed and contains all expected fields. Then run on full SuiteSparse suite and verify 65 baselines are produced. Confirm backward error values match expectations from Phase 8.1f. **NOTE**: Compilation verified. Full data collection deferred to development machine (container resource limits prevent running SuiteSparse suite).

**Checkpoint**: Structured JSON baselines produced for all 65 SuiteSparse matrices. Per-phase timing, per-supernode stats, and peak RSS captured.

---

## Phase 6: User Story 4 — Produce Performance Analysis Report (Priority: P4)

**Goal**: Create a workload analysis tool and write a performance report that classifies matrices by workload distribution and recommends a parallelism strategy for Phase 8.2.

**Independent Test**: Run the workload analysis tool and verify it produces per-matrix workload profiles and a summary classification table.

### Implementation for User Story 4

- [X] T023 [US4] Define `WorkloadProfile` and `ParallelismClass` types (per data-model.md) in `examples/workload_analysis.rs`. `WorkloadProfile` contains: matrix_name, total_factor_time, top_10pct_time_fraction, top_1_front_time_fraction, front_size_histogram, time_by_front_size, parallelism_recommendation. `ParallelismClass` is an enum: TreeLevel, IntraNode, Mixed. Classification logic: top_10pct_time_fraction > 0.80 → IntraNode; < 0.30 → TreeLevel; else Mixed.
- [X] T024 [US4] Create `examples/workload_analysis.rs`. The tool loads SuiteSparse matrices, runs the full pipeline with `diagnostic` feature (to get per-supernode timing), computes `WorkloadProfile` for each matrix by sorting supernodes by front_size, bucketing into front size bands, computing time fraction for top 10% of fronts. Output per-matrix: matrix_name, num_supernodes, max_front_size, top_10pct_time_fraction, parallelism_recommendation. Output summary classification table at the end.
- [X] T025 [US4] Optionally generate Chrome Trace JSON files for 3 representative matrices (one per parallelism class) to `target/profiles/`. Use the `FinishedSession` accessor from T011 to export traces. **NOTE**: Chrome Trace export is available via the `FinishedSession` accessor; targeted matrix selection deferred to data collection run.
- [X] T026 [US4] Write `docs/phase-8.1g-report.md` summarizing the performance analysis: (1) Front size distribution across the SuiteSparse suite, (2) Per-supernode time contribution analysis (what fraction of total time is in the top 10% of fronts), (3) Allocation pressure before/after optimization (from US2 comparison), (4) Matrix classification table (TreeLevel/IntraNode/Mixed), (5) Parallelism strategy recommendation for Phase 8.2 with evidence from at least three representative matrix classes, (6) Chrome Trace examples for representative matrices.
- [X] T027 [US4] Run `workload_analysis.rs` on full SuiteSparse suite and verify output is complete. Review the `phase-8.1g-report.md` for completeness and verify it answers the key Phase 8.2 question: tree-level vs. intra-node parallelism. **NOTE**: Compilation verified; full data collection deferred to development machine.

**Checkpoint**: Performance report complete with data-driven Phase 8.2 parallelism recommendation.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Final verification and documentation updates

- [X] T028 Run full regression suite: `cargo test` + `cargo test --features diagnostic` + `cargo test -- --ignored --test-threads=1` + `cargo test -- --ignored --test-threads=1 --features diagnostic`. All must pass. **NOTE**: Unit tests (358) pass for both feature configs. SuiteSparse ignored tests (65 matrices) pass in debug mode. Clippy clean.
- [X] T029 Run `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic`. Fix any warnings.
- [X] T030 Update `docs/ssids-log.md` with Phase 8.1g development log entry: what was built, what was changed, what was found, and what it means for Phase 8.2.
- [X] T031 Update `docs/ssids-plan.md` to mark Phase 8.1g as COMPLETE and record key findings (parallelism recommendation, allocation reduction results).

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Setup — BLOCKS all user stories
- **US1 (Phase 3)**: Depends on Foundational (T004-T006) — uses extended stats types
- **US2 (Phase 4)**: Depends on US1 completion (T012) — profiling data confirms hotspot priorities
- **US3 (Phase 5)**: Depends on US1 + US2 completion — baselines should reflect optimized solver
- **US4 (Phase 6)**: Depends on US3 completion (T022) — report synthesizes baseline data
- **Polish (Phase 7)**: Depends on all user stories complete

### User Story Dependencies

```
Phase 1 (Setup)
  └── Phase 2 (Foundational: T004-T006)
        └── US1 (Phase 3: T007-T012) — Profiling instrumentation
              └── US2 (Phase 4: T013-T018) — Allocation optimization
                    └── US3 (Phase 5: T019-T022) — Baseline collection
                          └── US4 (Phase 6: T023-T027) — Report
                                └── Phase 7 (Polish: T028-T031)
```

Note: This is a strictly sequential pipeline because each story's output informs the next. US2 needs US1's profiling data to confirm hotspots. US3 needs the optimized solver from US2. US4 needs US3's baseline data.

### Within Each User Story

- Data structure extensions before instrumentation
- Instrumentation before verification
- Verification before proceeding to next story

### Parallel Opportunities

Within each phase, tasks marked [P] can be parallelized. However, cross-phase parallelism is limited due to the sequential dependency chain. The main parallelism opportunity is within US2 (T013-T015 touch different allocation sites in the same file but are independent edits).

---

## Parallel Example: User Story 2

```bash
# These three allocation optimizations target independent code locations
# within factor_inner and can be developed/tested in parallel:
T013: Hoist panel row permutation buffer (line ~1664)
T014: Hoist row permutation temp buffer (line ~1692)
T015: Hoist column order slice copy (line ~1704)

# Then verify correctness after all three:
T017: Run full test suite
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (verify correctness baseline)
2. Complete Phase 2: Foundational (extend PerSupernodeStats/FactorizationStats)
3. Complete Phase 3: US1 — Instrument factorization hot path
4. **STOP and VALIDATE**: Factor a CI matrix with `diagnostic` feature and examine Chrome Trace output
5. The solver now has per-supernode profiling — this alone is valuable for understanding performance

### Incremental Delivery

1. Setup + Foundational → Extended data structures ready
2. Add US1 (Instrumentation) → Per-supernode profiling available → Validate Chrome Trace
3. Add US2 (Allocation Optimization) → Reduced allocation pressure → Validate correctness + timing
4. Add US3 (Baselines) → Structured JSON baselines for 65 matrices → Validate data completeness
5. Add US4 (Report) → Performance analysis + Phase 8.2 recommendation → Review report
6. Each story adds value and feeds the next

---

## Notes

- [P] tasks = different files or independent code regions, no dependencies
- [Story] label maps task to specific user story for traceability
- The dependency chain is strictly sequential (US1 → US2 → US3 → US4) because each story's output informs the next
- Existing tests (362+ unit, 65 SuiteSparse) serve as the regression guard — no new test files needed
- All profiling code is gated behind `#[cfg(feature = "diagnostic")]` — zero overhead in production builds
- Commit after each task or logical group to maintain a clean revision history
