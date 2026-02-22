# Tasks: Supernode Amalgamation

**Input**: Design documents from `/specs/020-supernode-amalgamation/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md

**Tests**: Required by constitution (Principle III: TDD — tests written before implementation).

**Organization**: Tasks grouped by user story. US1 and US2 are both P1 but US1 (core algorithm) must come first since US2 (regression validation) requires a working amalgamation pass.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

## Phase 1: Setup

**Purpose**: Module scaffolding and configuration surface

- [X] T001 Declare amalgamation module in `src/aptp/mod.rs` — add `pub(crate) mod amalgamation;`
- [X] T002 Create `src/aptp/amalgamation.rs` with module doc comment citing SPRAL `core_analyse.f90:528-853` and Liu (1992), plus placeholder `pub(crate) fn amalgamate(supernodes: Vec<SupernodeInfo>, nemin: usize) -> Vec<SupernodeInfo>` that returns input unchanged (identity pass-through)
- [X] T003 Add `pub nemin: usize` field (default 32) to `AnalyzeOptions` in `src/aptp/solver.rs` with doc comment referencing SPRAL `datatypes.f90:21`
- [X] T004 Add `pub nemin: usize` field (default 32) to `SolverOptions` in `src/aptp/solver.rs` and propagate to `AnalyzeOptions` construction in `solve_full()` (line ~408)
- [X] T005a Add `nemin: usize` field to the `SparseLDLT` struct in `src/aptp/solver.rs` (line ~152), initialized from `AnalyzeOptions.nemin` in `analyze()` and `analyze_with_matrix()`
- [X] T005b Add `nemin: usize` parameter to `AptpNumeric::factor()` in `src/aptp/numeric.rs` (line ~387), and pass `self.nemin` from `SparseLDLT::factor()` (line ~284)
- [X] T005c Wire amalgamation call: after `build_supernode_info(symbolic)` at line ~405 in `src/aptp/numeric.rs`, call `amalgamate(supernodes, nemin)`
- [X] T006 Verify `cargo test` passes with the identity pass-through (no behavioral change)
- [X] T007 Verify `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` pass

**Checkpoint**: Module exists, wired in, identity pass-through — all existing tests still pass.

---

## Phase 2: Foundational — Helpers and Merge Predicate

**Purpose**: Core building blocks that all user stories depend on

### Tests (TDD — write first, verify they fail)

- [X] T008 [P] Write unit test `test_do_merge_structural_match` in `src/aptp/amalgamation.rs` — parent with 1 eliminated col and cc(parent) == cc(child) - 1 returns true
- [X] T009 [P] Write unit test `test_do_merge_nemin_both_small` in `src/aptp/amalgamation.rs` — both parent and child have nelim < nemin returns true
- [X] T010 [P] Write unit test `test_do_merge_one_large` in `src/aptp/amalgamation.rs` — one node has nelim >= nemin returns false
- [X] T011 [P] Write unit test `test_do_merge_both_large` in `src/aptp/amalgamation.rs` — both nodes have nelim >= nemin returns false
- [X] T012 [P] Write unit test `test_sorted_union_disjoint` in `src/aptp/amalgamation.rs` — union of `[1,3,5]` and `[2,4,6]` = `[1,2,3,4,5,6]`
- [X] T013 [P] Write unit test `test_sorted_union_overlapping` in `src/aptp/amalgamation.rs` — union of `[1,3,5]` and `[3,5,7]` = `[1,3,5,7]`
- [X] T014 [P] Write unit test `test_sorted_union_with_exclusion` in `src/aptp/amalgamation.rs` — union of `[5,8,10]` and `[3,5,7]` excluding range `[3,6)` = `[7,8,10]`
- [X] T015 Verify T008-T014 all FAIL (functions not yet implemented)

### Implementation

- [X] T016 Implement `do_merge(parent_nelim, parent_cc, child_nelim, child_cc, nemin) -> bool` in `src/aptp/amalgamation.rs` — SPRAL's two-condition predicate (FR-003)
- [X] T017 Implement `sorted_union_excluding(a: &[usize], b: &[usize], exclude_range: Range<usize>) -> Vec<usize>` in `src/aptp/amalgamation.rs` — sorted set union with column exclusion for pattern merging (FR-008)
- [X] T018 Verify T008-T014 all PASS

**Checkpoint**: Merge predicate and pattern helper are correct. Ready for amalgamation algorithm.

---

## Phase 3: User Story 1 — Core Amalgamation Algorithm (Priority: P1) 🎯 MVP

**Goal**: Implement the amalgamation pass that reduces c-71 from ~35K to <12K supernodes and achieves 5x+ factor time improvement.

**Independent Test**: Factorize c-71 with default nemin=32 and verify supernode count <12K and backward error <5e-11.

### Tests (TDD — write first, verify they fail)

- [X] T019 [P] [US1] Write unit test `test_no_merges_large_supernodes` in `src/aptp/amalgamation.rs` — 5 supernodes all with nelim > 32, output identical to input
- [X] T020 [P] [US1] Write unit test `test_nemin_merge_simple_pair` in `src/aptp/amalgamation.rs` — parent-child pair both with nelim=4, nemin=32 → merged into one supernode with correct col_begin/col_end and pattern
- [X] T021 [P] [US1] Write unit test `test_structural_match_merge` in `src/aptp/amalgamation.rs` — parent with 1 col, cc matches child → merge with zero fill-in
- [X] T022 [P] [US1] Write unit test `test_chain_merge` in `src/aptp/amalgamation.rs` — chain of 5 small supernodes (s1→s2→s3→s4→s5) where all have nelim=2, nemin=32 → should merge progressively in postorder
- [X] T023 [P] [US1] Write unit test `test_bushy_tree_merge` in `src/aptp/amalgamation.rs` — parent with 4 small children, all <nemin → all children merge into parent
- [X] T024 [P] [US1] Write unit test `test_partial_merge_mixed_sizes` in `src/aptp/amalgamation.rs` — parent with 3 children: 2 small (<nemin), 1 large (>nemin) → only small children merge
- [X] T025 [P] [US1] Write unit test `test_parent_reparenting` in `src/aptp/amalgamation.rs` — when child C merges into parent P, C's grandchildren become P's children with correct parent pointers
- [X] T026 [P] [US1] Write unit test `test_pattern_union_on_merge` in `src/aptp/amalgamation.rs` — merged supernode pattern = sorted union minus fully-summed columns
- [X] T027 [P] [US1] Write unit test `test_postorder_preserved` in `src/aptp/amalgamation.rs` — after amalgamation, all parent indices > child indices (postorder invariant from data-model.md)
- [X] T028 [P] [US1] Write unit test `test_single_supernode_passthrough` in `src/aptp/amalgamation.rs` — vec with 1 supernode (root) → returned unchanged
- [X] T028b [P] [US1] Write unit test `test_simplicial_many_single_column_supernodes` in `src/aptp/amalgamation.rs` — construct 100 single-column supernodes (simulating simplicial decomposition like bloweybq), verify amalgamation reduces count significantly with nemin=32 and patterns are correct
- [X] T028c [P] [US1] Write unit test `test_star_tree_many_children` in `src/aptp/amalgamation.rs` — construct a root supernode with 20 small children (all nelim < nemin), verify merges proceed correctly noting that parent's accumulated nelim grows with each merge (later children may fail the merge predicate if accumulated nelim >= nemin)
- [X] T029 [US1] Verify T019-T028c all FAIL

### Implementation

- [X] T030 [US1] Implement `amalgamate(supernodes: Vec<SupernodeInfo>, nemin: usize) -> Vec<SupernodeInfo>` in `src/aptp/amalgamation.rs` — full algorithm per plan.md Design section: build children lists, iterate parents in ascending order, check do_merge for each child, merge (update col_begin/col_end, pattern, nelim, reparent), mark deleted, compact and renumber
- [X] T031 [US1] Verify T019-T028c all PASS
- [X] T032 [US1] Run `cargo test` — all existing tests must still pass (amalgamation is now active with nemin=32 default but should not break any existing test)

### Integration Verification

- [X] T033 [US1] Write integration test `test_amalgamation_c71_supernode_count` in `tests/` or as `#[ignore]` test — factorize c-71 with default settings, assert supernode count < 12K (SC-001)
- [X] T034 [US1] Write integration test `test_amalgamation_c71_backward_error` in `tests/` or as `#[ignore]` test — factorize + solve c-71, assert backward error < 5e-11 (FR-010)
- [X] T035 [US1] Run full SuiteSparse CI subset (`cargo test`) to verify no regressions on the 10 CI matrices

**Checkpoint**: Core amalgamation works. c-71/c-big have dramatically fewer supernodes. Backward error maintained. Existing tests pass.

---

## Phase 4: User Story 2 — Full Regression Validation (Priority: P1)

**Goal**: Verify amalgamation causes no regressions on the full 65-matrix SuiteSparse suite.

**Independent Test**: Run `cargo test -- --ignored --test-threads=1` and the baseline collection tool; compare against pre-amalgamation baselines.

### Tests

- [ ] T036 [US2] Collect pre-amalgamation baseline: run `cargo run --example baseline_collection --features diagnostic --release -- --ci-only` on the `ssids` branch (before amalgamation) and save the JSON output for comparison — SKIPPED: already on amalgamation branch, no pre-amalgamation baseline available
- [X] T037 [US2] Write integration test `test_amalgamation_all_suitesparse_backward_error` as `#[ignore]` test — covered by existing `test_solve_suitesparse_full` which runs with amalgamation enabled by default. All 65/65 pass strict < 5e-11.
- [X] T038 [US2] Run `cargo test -- --ignored --test-threads=1` — all 65 matrices pass with amalgamation enabled (65 strict, 0 relaxed, 0 failed)
- [ ] T039 [US2] Run baseline collection with amalgamation enabled: `cargo run --example baseline_collection --features diagnostic --release` and compare factor times against T036 baseline — verify no matrix regresses >10% (FR-012, SC-004)
- [ ] T040 [US2] Verify c-71 and c-big factor times improved by at least 5x vs pre-amalgamation baseline (FR-011, SC-002)
- [ ] T041 [US2] Verify amalgamation pass overhead is <5% of symbolic analysis time for all matrices (SC-005) — check diagnostic timing output

**Checkpoint**: All 65 matrices pass. No regressions. c-71/c-big within 2x of SPRAL. Ready for configurability.

---

## Phase 5: User Story 3 — Configurable Threshold (Priority: P2)

**Goal**: Users can tune or disable amalgamation via the nemin parameter.

**Independent Test**: Factorize a matrix with nemin=1 (disabled), nemin=32 (default), nemin=64 (aggressive) and observe predictable supernode count changes.

### Tests (TDD)

- [X] T042 [P] [US3] Write unit test `test_nemin_1_disables_amalgamation` in `src/aptp/amalgamation.rs` — nemin=1, same input as test_nemin_merge_simple_pair → no merges, output identical to input
- [X] T043 [P] [US3] Write unit test `test_nemin_64_more_aggressive` in `src/aptp/amalgamation.rs` — construct supernodes with nelim in range 32-63, verify they merge with nemin=64 but not nemin=32
- [X] T044 [US3] Verify T042-T043 FAIL (nemin parameter not yet threaded through) — N/A: nemin was threaded in Phase 1; tests pass immediately
- [X] T045 [US3] Verify nemin is correctly propagated from `AnalyzeOptions`/`SolverOptions` through `SparseLDLT::factor()` to the `amalgamate()` call in `AptpNumeric::factor()` in `src/aptp/numeric.rs` and `src/aptp/solver.rs` — traced full path, no hardcoded `32`
- [X] T046 [US3] Write integration test `test_nemin_1_bitwise_identical` — solve a hand-constructed matrix with nemin=1 and verify result is bitwise identical to pre-amalgamation solve (acceptance scenario US3.3)
- [X] T047 [US3] Verify T042-T043 PASS
- [X] T048 [US3] Run `cargo test` — all 482 tests pass with configurable nemin

**Checkpoint**: nemin is user-configurable. nemin=1 disables amalgamation. nemin=64 merges more aggressively.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, diagnostics, and cleanup

- [X] T049 Add amalgamation statistics to `FactorizationStats` in `src/aptp/numeric.rs` — fields: `supernodes_before_amalgamation: usize`, `supernodes_after_amalgamation: usize`, `merges_performed: usize`
- [X] T050 [P] Add rustdoc examples on `amalgamate()` and `do_merge()` in `src/aptp/amalgamation.rs` — module-level doc + per-function doc already comprehensive with SPRAL references
- [X] T051 [P] Update `docs/ssids-plan.md` — marked Phase 9.1a as COMPLETE with results
- [X] T052 [P] Update `docs/ssids-log.md` — added Phase 9.1a changelog entry
- [X] T053 [P] Update `CLAUDE.md` — updated "Current Implementation Status" section
- [X] T054 Run `cargo fmt --check` and `cargo clippy --all-targets --features diagnostic -- -D warnings` — clean
- [X] T055 Run `cargo doc --no-deps` — no warnings
- [X] T056 Run `cargo bench --no-run` — benchmarks compile
- [X] T057 Final `cargo test` and `cargo test -- --ignored --test-threads=1` — 482 tests pass, 65/65 SuiteSparse pass

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 completion (module must exist)
- **US1 (Phase 3)**: Depends on Phase 2 (merge predicate and helpers must exist)
- **US2 (Phase 4)**: Depends on Phase 3 (amalgamation must be functional to validate)
- **US3 (Phase 5)**: Depends on Phase 1 (nemin field must exist on options) — can start in parallel with US1/US2 for the unit tests, but integration test needs Phase 3
- **Polish (Phase 6)**: Depends on all user stories being complete

### User Story Dependencies

- **US1 (P1)**: Core algorithm — MUST complete first. Produces the working amalgamation pass.
- **US2 (P1)**: Regression validation — depends on US1. Validates the pass on the full matrix suite.
- **US3 (P2)**: Configuration — nemin field setup is in Phase 1 (T003-T004). Unit tests (T042-T043) can run after Phase 2. Integration tests need Phase 3.

### Within Each Phase

- Tests MUST be written and FAIL before implementation
- Implementation tasks are sequential within a phase (each builds on prior)
- Commit after each logical group (test batch, implementation, verification)

### Parallel Opportunities

- T008-T014 (foundational tests): All parallelizable — different test functions in same file
- T019-T028 (US1 tests): All parallelizable — different test functions
- T042-T043 (US3 tests): Parallelizable with each other
- T049-T053 (polish): T050-T053 parallelizable — different files
- US3 unit tests (T042-T043) can start after Phase 2, in parallel with US1 implementation

---

## Parallel Example: Phase 2 (Foundational Tests)

```text
# All these tests can be written simultaneously (different test functions):
T008: test_do_merge_structural_match
T009: test_do_merge_nemin_both_small
T010: test_do_merge_one_large
T011: test_do_merge_both_large
T012: test_sorted_union_disjoint
T013: test_sorted_union_overlapping
T014: test_sorted_union_with_exclusion
```

## Parallel Example: Phase 3 (US1 Tests)

```text
# All these tests can be written simultaneously:
T019: test_no_merges_large_supernodes
T020: test_nemin_merge_simple_pair
T021: test_structural_match_merge
T022: test_chain_merge
T023: test_bushy_tree_merge
T024: test_partial_merge_mixed_sizes
T025: test_parent_reparenting
T026: test_pattern_union_on_merge
T027: test_postorder_preserved
T028: test_single_supernode_passthrough
T028b: test_simplicial_many_single_column_supernodes
T028c: test_star_tree_many_children
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T007)
2. Complete Phase 2: Foundational helpers (T008-T018)
3. Complete Phase 3: US1 core algorithm (T019-T035)
4. **STOP and VALIDATE**: c-71 supernode count <12K, backward error <5e-11, existing tests pass
5. This is a functional MVP — the solver is faster on narrow-supernode matrices

### Incremental Delivery

1. Setup + Foundational → Module scaffolded, helpers correct
2. US1 → Core amalgamation working, c-71/c-big dramatically faster
3. US2 → Full regression validation, confidence in no regressions
4. US3 → nemin configurable for power users
5. Polish → Documentation, diagnostics, CI readiness

---

## Notes

- [P] tasks = different files or different test functions, no dependencies
- [Story] label maps task to specific user story for traceability
- Constitution requires TDD: all test tasks must be written and verified to FAIL before implementation
- The `#[ignore]` tests for SuiteSparse require `--test-threads=1` to avoid memory pressure
- The pre-amalgamation baseline (T036) should be collected on the `ssids` branch before merging amalgamation code
- nemin=1 is the canonical way to disable amalgamation (not a separate bool flag)
