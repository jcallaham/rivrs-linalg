# Tasks: Profiling and Debug Tools

**Input**: Design documents from `/specs/008-profiling-debug-tools/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, contracts/

**Tests**: Included — the project constitution (Principle III) requires TDD for all components.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3, US4)
- Include exact file paths in descriptions

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Create module scaffolding and shared utilities

- [x] T001 Register `profiling` and `debug` modules in `src/lib.rs` behind `#[cfg(feature = "test-util")]`, following the existing pattern used by `benchmarking` and `testing` modules
- [x] T002 [P] Create `src/profiling/mod.rs` with module declarations for `session`, `section`, `report`, and `memory` submodules, plus public re-exports
- [x] T003 [P] Create `src/debug/mod.rs` with module declarations for `sparsity` and `etree` submodules, plus public re-exports
- [x] T004 Add `read_current_rss_kb() -> Option<u64>` to `src/benchmarking/rss.rs` reading `VmRSS` from `/proc/self/status`, with `#[cfg(target_os = "linux")]` gating and tests matching the existing `read_peak_rss_kb` pattern

**Checkpoint**: `cargo test` and `cargo clippy` pass. New modules compile (empty). RSS utility extended.

---

## Phase 2: User Story 1 — Profile Solver Component Execution Times (Priority: P1) 🎯 MVP

**Goal**: Developers can instrument code with profiling sections, get hierarchical timing breakdowns, and export Chrome Trace files for visual analysis. Profiling compiles to zero-cost no-ops when `test-util` is disabled.

**Independent Test**: Instrument a mock function with nested profiling sections, verify timing data is collected, hierarchy is correct, summary report shows expected structure, and Chrome Trace JSON is valid.

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T005 [P] [US1] Write tests for `ProfileEvent` and `ProfileSection` data types in `src/profiling/section.rs`: construction, field access, aggregation of repeated labels (total time, call count, min/max/mean), and serde roundtrip for `ProfileSection`
- [ ] T006 [P] [US1] Write tests for `ProfileSession` in `src/profiling/session.rs`: new session records start time; `enter_section` returns guard; dropping guard records event; nested sections produce correct depth; `finish()` returns `FinishedSession` with merged sections; repeated section names are aggregated with correct call count, min, max, mean; panic in a profiled section (use `std::panic::catch_unwind`) still records partial timing data for sections completed before the panic
- [ ] T007 [P] [US1] Write tests for `FinishedSession` reporting in `src/profiling/report.rs`: `summary_report()` produces hierarchical text with correct columns (name, total, calls, mean, min, max, parent%); `export_chrome_trace()` produces valid JSON with `traceEvents` array containing Complete Events (ph="X") with correct `name`, `ts`, `dur`, `pid`, `tid` fields; nested events have overlapping time ranges
- [ ] T008 [US1] Write thread-safety test in `src/profiling/session.rs`: spawn multiple threads each calling `enter_section` on a shared `ProfileSession` (via `Arc`), verify `finish()` collects events from all threads without data races or lost events

### Implementation for User Story 1

- [ ] T009 [P] [US1] Implement `ProfileEvent` and `ProfileSection` structs in `src/profiling/section.rs` per data-model.md: `ProfileEvent` with name, start Instant, duration, thread_id, depth; `ProfileSection` with name, total_duration, call_count, min/max/mean_duration, parent_pct, children Vec. Include `From<Vec<ProfileEvent>>` aggregation logic that groups by (name, depth, parent) and computes min/max/mean
- [ ] T010 [US1] Implement `ProfileSession`, `SectionGuard`, and `FinishedSession` in `src/profiling/session.rs` per contracts/profiling-api.md: `ProfileSession::new()` records start Instant; `enter_section()` returns RAII `SectionGuard` that records a `ProfileEvent` into thread-local storage on drop; `finish()` collects all thread-local buffers, reconstructs hierarchy from timestamps and depth, aggregates into `ProfileSection` tree, returns `FinishedSession`. Use `thread_local!` with `RefCell<Vec<ProfileEvent>>` for zero-contention recording. `ProfileSession` must be `Send + Sync`; `SectionGuard` must NOT be `Send`
- [ ] T011 [US1] Implement `summary_report()` in `src/profiling/report.rs`: generate hierarchical text table with columns (Section, Total, Calls, Mean, Min, Max, Parent%) using indentation for nesting depth. Include session total duration header
- [ ] T012 [US1] Implement `export_chrome_trace()` in `src/profiling/report.rs`: serialize `FinishedSession` events to Chrome Trace Event format JSON with `traceEvents` array of Complete Events (ph="X"), timestamps in microseconds relative to session start, pid from `std::process::id()`, tid from thread index, `args` containing call_count. Individual invocations (not aggregated) in trace export per FR-005a
- [ ] T013 [US1] Verify all T005-T008 tests pass. Run `cargo test --features test-util` and `cargo clippy --all-targets`

**Checkpoint**: Profiling module fully functional. Can instrument code, get summary reports, export Chrome Trace files. Thread-safe.

---

## Phase 3: User Story 2 — Track Memory Usage During Solver Operations (Priority: P2)

**Goal**: Developers can record RSS snapshots at labeled points during execution and receive a memory report showing peak RSS, per-section deltas, and a formatted terminal display.

**Independent Test**: Create a `MemoryTracker`, take snapshots around allocation-heavy operations, verify report shows plausible deltas and correct peak RSS identification.

### Tests for User Story 2

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T014 [P] [US2] Write tests for `MemorySnapshot` and `MemoryDelta` in `src/profiling/memory.rs`: construction, field access, `None` handling for RSS values on non-Linux
- [ ] T015 [P] [US2] Write tests for `MemoryTracker` in `src/profiling/memory.rs`: new tracker starts empty; `snapshot()` adds entries; multiple snapshots produce correct ordering; `report()` computes deltas between consecutive snapshots; `peak_rss_kb()` returns maximum across all snapshots; empty tracker produces empty report
- [ ] T016 [US2] Write tests for `MemoryReport::display_report()` in `src/profiling/memory.rs`: output contains header, column labels, one row per snapshot with RSS and delta values, peak RSS footer; handles `None` RSS gracefully with "N/A" display

### Implementation for User Story 2

- [ ] T017 [P] [US2] Implement `MemorySnapshot`, `MemoryDelta`, and `MemoryReport` structs in `src/profiling/memory.rs` per data-model.md: `MemorySnapshot` with label, timestamp, rss_kb Option, peak_rss_kb Option; `MemoryDelta` with from_label, to_label, rss_delta_kb Option<i64>; `MemoryReport` with snapshots Vec, peak_rss_kb Option, deltas Vec
- [ ] T018 [US2] Implement `MemoryTracker` in `src/profiling/memory.rs` per contracts/memory-api.md: `new()` creates empty tracker; `snapshot(label)` calls `read_current_rss_kb()` and `read_peak_rss_kb()` from `crate::benchmarking::rss` and appends a `MemorySnapshot`; `report()` computes deltas between consecutive snapshots and finds peak RSS
- [ ] T019 [US2] Implement `MemoryReport::display_report()` in `src/profiling/memory.rs`: formatted text table with columns (Label, RSS KB, Peak KB, Delta KB), separator lines, peak RSS footer with MB conversion. Handle `None` values with "N/A" display
- [ ] T020 [US2] Verify all T014-T016 tests pass. Run `cargo test --features test-util` and `cargo clippy --all-targets`

**Checkpoint**: Memory tracking module fully functional. Can snapshot RSS at labeled points, compute deltas, display formatted report.

---

## Phase 4: User Story 3 — Visualize Sparse Matrix Structure (Priority: P3)

**Goal**: Developers can generate text-based sparsity pattern visualizations with configurable dimensions, density-based characters, and downsampling for large matrices.

**Independent Test**: Generate visualizations for hand-constructed test matrices and verify patterns match known non-zero structure. Test downsampling on larger matrices.

### Tests for User Story 3

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T021 [P] [US3] Write tests for `SparsityDisplay` construction and rendering in `src/debug/sparsity.rs`: `from_sparse()` on a small known matrix (e.g., 5x5 tridiagonal) produces correct grid; non-zero positions show filled characters; zero positions show dots; header contains correct dimensions and nnz count
- [ ] T022 [P] [US3] Write tests for downsampling in `src/debug/sparsity.rs`: matrix larger than max_width/max_height is downsampled; density values are correct (count non-zeros in each bin / bin area); display dimensions match configured limits; header includes "[WxH view]" annotation
- [ ] T023 [P] [US3] Write tests for ASCII mode and Display impl in `src/debug/sparsity.rs`: `with_ascii_only(true)` uses only ASCII characters (`#`, `+`, `-`, `.`); `fmt::Display` produces same output as `render()`; empty matrix (0x0 or all-zeros) produces meaningful representation

### Implementation for User Story 3

- [ ] T024 [US3] Implement `SparsityDisplay` struct and builder methods in `src/debug/sparsity.rs` per contracts/debug-viz-api.md: `from_sparse()` computes bin dimensions based on matrix size and max_width/max_height defaults (80x40); iterate CSC column pointers to count non-zeros per bin cell; store density as f64 (0.0-1.0) per cell. Builder methods `with_max_width()`, `with_max_height()`, `with_ascii_only()` return `Self`
- [ ] T025 [US3] Implement `render()` and `fmt::Display` in `src/debug/sparsity.rs`: generate header line with matrix name/dimensions/nnz/density/downsampling note; map each cell density to character (Unicode: `.`/`░`/`▒`/`▓`/`█`; ASCII: `.`/`-`/`+`/`#`); join rows with newlines. `Display` delegates to `render()`
- [ ] T026 [US3] Verify all T021-T023 tests pass. Run `cargo test --features test-util` and `cargo clippy --all-targets`

**Checkpoint**: Sparsity visualization fully functional. Can render any sparse matrix as text with configurable dimensions.

---

## Phase 5: User Story 4 — Visualize Elimination Tree Structure (Priority: P4)

**Goal**: Developers can inspect elimination trees as text tree diagrams (small trees) or summary statistics (any size), supporting Phase 3 (Symbolic Analysis) debugging.

**Independent Test**: Construct elimination trees from known parent arrays and verify tree rendering matches expected structure and statistics are correct.

### Tests for User Story 4

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T027 [P] [US4] Write tests for `EliminationTreeStats` in `src/debug/etree.rs`: from a known parent array (e.g., linear chain of 5 nodes), verify depth, leaf count, branching min/max/mean, and subtree sizes are correct. Test a balanced binary tree and a star topology for diverse structures
- [ ] T028 [P] [US4] Write tests for `ETreeDisplay` rendering in `src/debug/etree.rs`: `render_tree()` on a small tree (n<20) produces text with `├──`, `└──`, `│` box-drawing characters showing correct parent-child relationships; `render_tree()` on a large tree (n>=20) falls back to stats view; `render_stats()` produces formatted statistics block

### Implementation for User Story 4

- [ ] T029 [US4] Implement `EliminationTreeStats` and `ETreeDisplay::from_parent_array()` in `src/debug/etree.rs` per contracts/debug-viz-api.md: parse parent array to build children lists; compute depth via DFS/BFS; count leaves; compute branching factor stats; compute subtree sizes. Handle root detection (parent[i]==i or sentinel)
- [ ] T030 [US4] Implement `render_tree()` and `render_stats()` in `src/debug/etree.rs`: `render_tree()` uses recursive traversal with box-drawing characters for small trees (n<20), delegates to `render_stats()` for larger trees; `render_stats()` formats depth, leaves (with percentage), branching factor range and mean, subtree size range and median
- [ ] T031 [US4] Verify all T027-T028 tests pass. Run `cargo test --features test-util` and `cargo clippy --all-targets`

**Checkpoint**: Elimination tree visualization fully functional. Can inspect tree structure and statistics.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Verify cross-cutting requirements, documentation, and integration

- [ ] T032 Verify feature gating (FR-004, FR-016, FR-017, SC-006): build without `test-util` feature (`cargo build`) and confirm profiling/debug modules are excluded; compare binary sizes with and without feature; write a compile-test (or integration test behind a `cfg` check) that verifies the no-op API compiles with identical signatures when `test-util` is disabled (ProfileSession::new, enter_section, finish, summary_report all return empty/ZST values); verify `profiling` and `debug` modules can each be used independently without importing the other
- [ ] T033 [P] Add rustdoc documentation with `# Examples` sections to all public types and functions in `src/profiling/mod.rs`, `src/debug/mod.rs`, and all submodules. Verify `cargo doc --no-deps` passes with `RUSTDOCFLAGS="-D warnings"`
- [ ] T034 [P] Run all hand-constructed test matrices through `SparsityDisplay` to verify SC-004: all 15 matrices produce correct visualizations matching their known sparsity structures
- [ ] T035 Update `docs/ssids-log.md` with Phase 1.4 entry documenting what was built, key decisions, and issues encountered
- [ ] T036 Update `docs/ssids-plan.md` Phase 1.4 section to mark as complete and note any deviations from original plan (similar to how Phase 1.3 was updated)
- [ ] T037 Run full test suite (`cargo test --all-targets`), clippy (`cargo clippy --all-targets -- -D warnings`), and benchmark compilation (`cargo bench --no-run`) to verify nothing is broken

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **US1 Profiling (Phase 2)**: Depends on Phase 1 (T001-T003 for module scaffolding)
- **US2 Memory (Phase 3)**: Depends on Phase 1 (T004 for RSS extension) and module scaffolding (T001-T002)
- **US3 Sparsity (Phase 4)**: Depends on Phase 1 (T001, T003 for debug module scaffolding)
- **US4 ETree (Phase 5)**: Depends on Phase 1 (T001, T003 for debug module scaffolding)
- **Polish (Phase 6)**: Depends on all user stories being complete

### User Story Dependencies

- **US1 (P1)**: Depends only on Phase 1 setup. No dependency on other stories.
- **US2 (P2)**: Depends on Phase 1 setup (T004 RSS extension). No dependency on US1.
- **US3 (P3)**: Depends on Phase 1 setup (T003 debug module). No dependency on US1/US2.
- **US4 (P4)**: Depends on Phase 1 setup (T003 debug module). No dependency on US1/US2/US3.

### Within Each User Story

- Tests MUST be written and FAIL before implementation (TDD per constitution)
- Data types before session/tracker logic
- Core logic before reporting/export
- Verify tests pass as final step

### Parallel Opportunities

- T002 and T003 can run in parallel (different modules)
- T005, T006, T007 can run in parallel (different test files within US1)
- T009 can run in parallel with T010 initial scaffolding (different files)
- T014, T015 can run in parallel (different test aspects within US2)
- T021, T022, T023 can run in parallel (different test aspects within US3)
- T027, T028 can run in parallel (different test aspects within US4)
- US3 and US4 can run in parallel after Phase 1 (both use debug module, different files)

---

## Parallel Example: User Story 1

```
# Phase 1 setup (sequential):
T001: Register modules in lib.rs
T002: Create profiling/mod.rs   ─┐
T003: Create debug/mod.rs       ─┘ (parallel)

# US1 tests (parallel):
T005: ProfileEvent/ProfileSection tests  ─┐
T006: ProfileSession tests               ─┤ (parallel)
T007: Report/Chrome Trace tests          ─┘

# US1 implementation:
T009: ProfileEvent + ProfileSection types  ─┐
T010: ProfileSession + SectionGuard        ─┘ (T009 first, T010 depends on T009)
T011: summary_report()                     ─┐
T012: export_chrome_trace()                ─┘ (parallel, both depend on T010)
T008: Thread-safety test (after T010)
T013: Verify all tests pass
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T004)
2. Complete Phase 2: US1 Profiling (T005-T013)
3. **STOP and VALIDATE**: Test profiling independently — instrument a mock function, verify summary report and Chrome Trace export
4. This alone delivers the most valuable development tool for all subsequent solver phases

### Incremental Delivery

1. Setup → Module scaffolding ready
2. US1 Profiling → Timing breakdowns available (MVP!)
3. US2 Memory → RSS tracking available
4. US3 Sparsity → Matrix visualization available
5. US4 ETree → Tree visualization available
6. Each story adds independent value without breaking previous stories

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Constitution requires TDD: all test tasks must run and fail before corresponding implementation tasks
- Commit after each task or logical group per CLAUDE.md guidelines
- The profiler's thread-local storage design means T010 is the most complex task — allow extra time
- Chrome Trace export (T012) uses existing serde_json dependency — no new deps needed
- RSS extension (T004) is a small addition to existing `benchmarking/rss.rs` — no code move needed
