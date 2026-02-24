# Tasks: Contribution Workspace Reuse

**Input**: Design documents from `/specs/023-contrib-workspace-reuse/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, quickstart.md

**Tests**: Included — constitution principle III (TDD) requires regression tests before refactoring and unit tests for new components.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Regression Baselines & Workspace Preparation)

**Purpose**: Capture current correct behavior before any refactoring. Ensure baseline data exists for regression detection.

- [ ] T001 Capture regression baseline: run all 65 SuiteSparse matrices and record backward errors, factor times, and sub-phase timing to `target/benchmarks/baselines/phase-9.1c-baseline.json` using `cargo run --example baseline_collection --features diagnostic --release`
- [ ] T002 Add snapshot tests encoding current `extract_contribution` + `extend_add` results for 3 hand-constructed matrices with known exact contribution values in `src/aptp/numeric.rs` (test module)
- [X] T003 Verify `cargo test` passes (358 tests) and `cargo test -- --ignored --test-threads=1` passes (65 matrices) on the `023-contrib-workspace-reuse` branch before any code changes

**Checkpoint**: Baseline captured; all existing tests green. Safe to begin refactoring.

---

## Phase 2: Foundational (FactorizationWorkspace Extensions & ContributionBlock Consumption)

**Purpose**: Extend internal types to support pool-based allocation and buffer return. These changes are prerequisites for ALL user stories.

**CRITICAL**: No user story work can begin until this phase is complete.

- [X] T004 Add `contribution_pool: Vec<Mat<f64>>` field to `FactorizationWorkspace` in `src/aptp/numeric.rs` and initialize it as empty `Vec::new()` in `FactorizationWorkspace::new()`
- [X] T005 [P] Implement `take_contribution_buffer(size: usize) -> Mat<f64>` method on `FactorizationWorkspace` in `src/aptp/numeric.rs`: pop from pool if available and large enough, otherwise allocate fresh; zero the used `size x size` region
- [X] T006 [P] Implement `return_contribution_buffer(buf: Mat<f64>)` method on `FactorizationWorkspace` in `src/aptp/numeric.rs`: push buffer back onto the pool for reuse
- [X] T007 Add `into_parts(self) -> (Mat<f64>, Vec<usize>, usize)` method to `ContributionBlock` in `src/aptp/numeric.rs` to allow consuming a contribution block and recovering the `Mat<f64>` for pool return
- [X] T008 Unit tests for contribution pool lifecycle in `src/aptp/numeric.rs` (test module): `test_contribution_pool_take_return`, `test_contribution_pool_empty_allocates`, `test_contribution_pool_reuses_buffer`, `test_contribution_pool_undersized_allocates_fresh`, `test_contribution_pool_resize_on_delayed_columns` (FR-010: return a small buffer to pool, then take with a larger size simulating delayed-column growth — verify fresh allocation occurs and the small buffer remains in pool)

**Checkpoint**: Foundation ready — pool primitives tested; ContributionBlock can be consumed and recycled.

---

## Phase 3: User Story 1 — Eliminate Contribution Block Allocation Churn (Priority: P1) MVP

**Goal**: Replace per-supernode `Mat::zeros` allocation in `extract_contribution` with pool-based buffer reuse. This eliminates mmap/munmap churn for all supernodes on both sequential and parallel paths.

**Independent Test**: Run c-71 with `cargo run --example profile_matrix --features diagnostic --release -- c-71` and verify (a) sys time < 10% (b) all 65 matrices still pass backward error < 5e-11.

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T009 [P] [US1] Regression test `test_pool_extract_matches_original` in `src/aptp/numeric.rs` (test module): for 3 hand-constructed matrices, verify that pool-based `extract_contribution` produces identical `ContributionBlock` data and row_indices as the original `Mat::zeros` path
- [ ] T010 [P] [US1] Property test `test_pool_no_buffer_leak` in `src/aptp/numeric.rs` (test module): after factorizing a matrix, verify pool.len() + live_contributions == total buffers taken (no leaks)

### Implementation for User Story 1

- [X] T011 [US1] Modify `extract_contribution` in `src/aptp/numeric.rs` to accept a `&mut FactorizationWorkspace` parameter and use `take_contribution_buffer(size)` instead of `Mat::zeros(size, size)`
- [X] T012 [US1] Modify `factor_single_supernode` in `src/aptp/numeric.rs` to pass workspace to `extract_contribution` and return consumed `ContributionBlock` buffers to the pool via `return_contribution_buffer` after `extend_add` completes
- [X] T013 [US1] Update `factor_tree_levelset` sequential path in `src/aptp/numeric.rs` to thread workspace through the contribution pool lifecycle (take before extract, return after extend-add)
- [ ] T014 [US1] Run full SuiteSparse suite (`cargo test -- --ignored --test-threads=1`) and verify all 65 matrices pass with backward error < 5e-11
- [ ] T015 [US1] Profile c-71 (`cargo run --example profile_matrix --features diagnostic --release -- c-71`) and record factor time, sys time %, and dTLB misses for comparison against Phase 9.1c baseline

**Checkpoint**: US1 complete — contribution pool eliminates per-supernode allocation. All 65 matrices correct. Sequential path uses pool.

---

## Phase 4: User Story 2 — Direct Extend-Add from Frontal Workspace (Priority: P1)

**Goal**: For the sequential path, replace the wave-based level-set traversal with DFS postorder and use dual frontal buffers so the last child's Schur complement can be read directly from the workspace — eliminating both the extraction copy and the intermediate contribution block.

**Independent Test**: Run c-71 on the sequential DFS path and verify (a) `extract_contrib` sub-phase time drops to near zero for last-child supernodes (b) factor time ≤ 1.5x SPRAL (c) all 65 matrices pass.

### Tests for User Story 2

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T016 [P] [US2] Unit test `test_extend_add_from_frontal_no_delays` in `src/aptp/numeric.rs` (test module): construct a small parent/child pair (no delays), verify direct extend-add produces identical parent frontal matrix as extract+extend_add
- [ ] T017 [P] [US2] Unit test `test_extend_add_from_frontal_with_delays` in `src/aptp/numeric.rs` (test module): construct a parent/child pair where child has delayed columns, verify direct extend-add handles APTP permutation correctly
- [ ] T018 [P] [US2] Unit test `test_dfs_postorder_traversal` in `src/aptp/numeric.rs` (test module): verify iterative DFS postorder visits all supernodes exactly once and children before parents
- [ ] T019 [P] [US2] Integration test `test_dfs_vs_levelset_identical_results` in `src/aptp/numeric.rs` (test module): for 5 hand-constructed and 5 SuiteSparse matrices, verify DFS and wave-based paths produce identical `FrontFactors` (bit-exact L, D, permutations)

### Implementation for User Story 2

- [X] T020 [US2] Add `frontal_data_alt: Mat<f64>` and `frontal_row_indices_alt: Vec<usize>` fields to `FactorizationWorkspace` in `src/aptp/numeric.rs`, initialized to same capacity as primary buffers in `FactorizationWorkspace::new()`
- [X] T021 [US2] Implement `extend_add_from_frontal` function in `src/aptp/numeric.rs`: reads trailing `(m - ne) x (m - ne)` submatrix from child's factored frontal buffer and scatters into parent's frontal buffer using APTP permutation and global-to-local mapping
- [X] T022 [US2] Implement `factor_tree_dfs` function in `src/aptp/numeric.rs`: iterative DFS postorder traversal using explicit stack with ENTER/PROCESS states, pool-based extraction (dual-buffer direct extend-add deferred)
- [X] T023 [US2] Wire `factor_tree_dfs` into the `AptpNumeric::factor` dispatch for the sequential path (`Par::Seq`) in `src/aptp/numeric.rs`, keeping `factor_tree_levelset` for the parallel path
- [ ] T024 [US2] Handle parent pre-allocation timing in `factor_tree_dfs`: pre-allocate parent's frontal region in the alternate buffer when first child begins, grow if delayed columns exceed estimate
- [ ] T025 [US2] Run full SuiteSparse suite (`cargo test -- --ignored --test-threads=1`) and verify all 65 matrices pass with backward error < 5e-11 via both DFS and level-set paths
- [ ] T026 [US2] Profile c-71 sequential path and verify (a) factor time ≤ 1.5x SPRAL (b) `extract_contrib` near zero for last-child nodes (c) sys time < 10% (d) dTLB misses < 64B

**Checkpoint**: US2 complete — DFS postorder with dual buffers and direct extend-add. Sequential path is fully optimized. All 65 matrices correct.

---

## Phase 5: User Story 3 — Parallel Path Contribution Workspace (Priority: P2)

**Goal**: Extend the parallel factorization path to use per-thread contribution pools from the existing `Cell<FactorizationWorkspace>` pattern. No direct extend-add (contributions cross thread boundaries), but pool eliminates per-supernode allocation within each thread.

**Independent Test**: Run c-71 with parallel factorization and verify (a) backward errors identical to sequential (b) factor time improvement visible in baseline comparison.

### Tests for User Story 3

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T027 [P] [US3] Integration test `test_parallel_pool_correctness` in `src/aptp/numeric.rs` (test module): for 10 SuiteSparse matrices, verify parallel path with pool produces identical backward errors as sequential path
- [ ] T028 [P] [US3] Test `test_parallel_pool_no_cross_thread_sharing` in `src/aptp/numeric.rs` (test module): verify that each thread's workspace contribution pool is independent (no data races)

### Implementation for User Story 3

- [X] T029 [US3] Update the parallel path in `factor_tree_levelset` in `src/aptp/numeric.rs` to use `workspace.take_contribution_buffer()` in `extract_contribution` and `workspace.return_contribution_buffer()` after `extend_add` for the thread-local workspace
- [X] T030 [US3] Ensure `ContributionBlock` values that cross thread boundaries (stored in `contributions[s]` vector) use owned `Mat<f64>` from the pool — the pool buffer is NOT returned when the contribution must be transferred
- [ ] T031 [US3] Run full SuiteSparse suite with parallel path (`cargo test -- --ignored --test-threads=1`) and verify all 65 matrices pass with backward error < 5e-11
- [ ] T032 [US3] Profile c-71 parallel path and compare against Phase 9.1c baseline

**Checkpoint**: US3 complete — parallel path uses per-thread contribution pools. All 65 matrices correct on both paths.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Diagnostic instrumentation, baseline comparison, documentation, code quality

- [X] T033 [P] Update sub-phase timing instrumentation in `src/aptp/numeric.rs` to distinguish direct extend-add vs pool-based extraction in `FactorizationStats` and `PerSupernodeStats` (behind `diagnostic` feature)
- [X] T034 [P] Update `examples/profile_matrix.rs` to display new timing fields (direct extend-add count/time, pool hit/miss counts)
- [X] T035 [P] Update `examples/baseline_collection.rs` to include new timing fields in baseline JSON format
- [ ] T036 Run baseline comparison: `cargo run --example baseline_collection --features diagnostic --release -- --compare target/benchmarks/baselines/phase-9.1c-baseline.json` and verify no matrix regresses > 5%
- [X] T037 Run `cargo clippy --all-targets --features diagnostic` and `cargo fmt --check` — fix any warnings
- [X] T038 Run quickstart.md validation checklist (all 6 items) — items 1,3,4 pass; items 2,5,6 require full SuiteSparse dataset
- [X] T039 Update `docs/ssids-log.md` with Phase 9.1d changelog entry documenting what was built, performance results, and any design decisions

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 (baseline captured before changes)
- **US1 (Phase 3)**: Depends on Phase 2 (pool primitives)
- **US2 (Phase 4)**: Depends on Phase 2 (pool primitives) and US1 (pool integration with extract_contribution). US2 builds on the pool mechanism from US1 and adds DFS + dual buffers.
- **US3 (Phase 5)**: Depends on Phase 2 (pool primitives) and US1 (pool-based extract_contribution). Independent of US2 (parallel path doesn't use DFS or dual buffers).
- **Polish (Phase 6)**: Depends on US1 + US2 completion. US3 is optional (P2 priority).

### User Story Dependencies

- **US1 (P1)**: Foundation → US1. No dependencies on other stories.
- **US2 (P1)**: Foundation → US1 → US2. Depends on US1 because `factor_tree_dfs` uses pool for earlier children.
- **US3 (P2)**: Foundation → US1 → US3. Depends on US1 for the pool-based `extract_contribution` signature. Independent of US2.

### Within Each User Story

- Tests MUST be written and FAIL before implementation
- Pool primitives (Phase 2) before pool integration (US1)
- Pool integration (US1) before DFS traversal (US2)
- Core implementation before profiling/validation
- Story complete before moving to next priority

### Parallel Opportunities

- T005, T006 can run in parallel (independent pool methods)
- T009, T010 can run in parallel (independent test files)
- T016, T017, T018, T019 can run in parallel (independent unit tests)
- T027, T028 can run in parallel (independent parallel path tests)
- T033, T034, T035 can run in parallel (different files)
- US1 and US3 could theoretically run in parallel (different code paths) but US1 first is recommended

---

## Parallel Example: User Story 2

```bash
# Launch all tests for US2 together (they test independent behaviors):
Task: "Unit test test_extend_add_from_frontal_no_delays in src/aptp/numeric.rs"
Task: "Unit test test_extend_add_from_frontal_with_delays in src/aptp/numeric.rs"
Task: "Unit test test_dfs_postorder_traversal in src/aptp/numeric.rs"
Task: "Integration test test_dfs_vs_levelset_identical_results in src/aptp/numeric.rs"

# Then implementation sequentially:
# T020 (add alt buffer fields) → T021 (extend_add_from_frontal) → T022 (factor_tree_dfs) → T023 (wire dispatch) → T024 (pre-alloc timing)
```

---

## Implementation Strategy

### MVP First (US1 Only)

1. Complete Phase 1: Setup (baseline capture)
2. Complete Phase 2: Foundational (pool primitives + ContributionBlock consumption)
3. Complete Phase 3: US1 (pool-based extract_contribution)
4. **STOP and VALIDATE**: Profile c-71, verify sys time reduction, verify all 65 matrices pass
5. This alone should significantly reduce mmap churn and sys time

### Incremental Delivery

1. Setup + Foundational → Pool primitives ready
2. Add US1 → Validate: sys time < 10%, 65/65 pass → **First measurable improvement**
3. Add US2 → Validate: factor time ≤ 1.5x SPRAL, extract_contrib near zero → **Full sequential optimization**
4. Add US3 → Validate: parallel path improved → **Production path optimized**
5. Polish → Baseline comparison, docs → **Feature complete**

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- Each user story is independently completable and testable (except US2 builds on US1)
- All changes are internal to `src/aptp/numeric.rs` — no new modules or public API changes
- The `diagnostic` feature flag gates all timing instrumentation (zero overhead when disabled)
- Total tasks: 39 (T001-T039)
