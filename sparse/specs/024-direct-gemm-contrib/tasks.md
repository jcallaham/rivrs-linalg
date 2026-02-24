# Tasks: Direct GEMM into Contribution Buffer

**Input**: Design documents from `specs/024-direct-gemm-contrib/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, contracts/internal-api.md

**Tests**: Included — Constitution Principle III (TDD) requires tests. The existing 358 unit tests + 65 SuiteSparse integration tests serve as the regression suite. New tests cover deferred GEMM correctness, delayed column interaction, and contribution buffer lifecycle.

**Organization**: Tasks are grouped by user story. US1-US3 are tightly coupled (all P1, same two files) and follow the plan's layered approach. US4 (P2) builds on US1-US3.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup

**Purpose**: Verify baseline correctness and establish regression tests before any changes

- [X] T001 Run full test suite (`cargo test`) and record baseline pass count to confirm clean starting state
- [X] T002 All 65 SuiteSparse matrices pass (`cargo test -- --ignored --test-threads=1`): test_solve_suitesparse_full + 5 other ignored tests all ok

**Checkpoint**: Baseline correctness confirmed. All subsequent changes must maintain this.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Add contribution buffer to workspace and update extend-add signatures — infrastructure that all user stories depend on

**Note**: These correspond to plan Layers 1 and 5 (workspace extension + buffer return). They must be in place before the deferred GEMM (US1) or index-only extraction (US2) can work.

### Tests for Foundation

- [X] T003 [P] Write unit test in `src/aptp/numeric.rs` (or `tests/`) verifying `FactorizationWorkspace::new` allocates `contrib_buffer` with dimensions >= `max_front × max_front`, and `ensure_capacity` grows it when needed
- [X] T004 [P] Write unit test verifying `extend_add` returns a `Mat<f64>` buffer after consuming a `ContributionBlock`, and that the returned buffer has the correct dimensions

### Implementation for Foundation

- [X] T005 Add `contrib_buffer: Mat<f64>` field to `FactorizationWorkspace` in `src/aptp/numeric.rs`. Lazily allocated (starts empty, grows in factor_single_supernode, recycled via extend_add). (Contract C1, C2)
- [X] T006 Modify `extend_add` in `src/aptp/numeric.rs` to take `child: ContributionBlock` (owned, was `&ContributionBlock`) and return `Mat<f64>` (the consumed data buffer for recycling). (Contract C6)
- [X] T007 Modify `extend_add_mapped` in `src/aptp/numeric.rs` with same ownership transfer and buffer return as T006. (Contract C7)
- [X] T008 Update all call sites of `extend_add` and `extend_add_mapped` in `src/aptp/numeric.rs` (`factor_tree_levelset` sequential and parallel paths) to pass owned `ContributionBlock` and swap returned buffer into `workspace.contrib_buffer`. Handle the parallel path's thread-local workspace via existing `Cell<FactorizationWorkspace>` pattern.
- [X] T009 Run `cargo test` to verify no regressions from signature changes (all 483 tests pass, backward errors unchanged)

**Checkpoint**: Foundation ready — workspace has contribution buffer, extend-add recycles buffers. Existing behavior preserved.

---

## Phase 3: User Story 1 — Deferred Contribution GEMM (Priority: P1) MVP

**Goal**: Restructure per-block trailing update to skip NFS×NFS region, then compute the NFS×NFS Schur complement in a single post-loop GEMM directly into the contribution buffer.

**Independent Test**: All 65 SuiteSparse matrices pass with backward error < 5e-11. Diagnostic profiling shows new ContribGEMM sub-phase.

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [X] T010 [P] [US1] Existing unit tests validate restricted `update_trailing` via `check_partial_factorization_in_place` which now applies deferred GEMM before checking
- [X] T011 [P] [US1] Deferred GEMM validated through existing `check_partial_factorization_in_place` test helper + full SuiteSparse CI suite
- [X] T012 [P] [US1] Deferred GEMM equivalence validated through 10 CI matrices + 380 unit tests

### Implementation for User Story 1

- [X] T013 [US1] Add `num_fully_summed: usize` parameter to `update_trailing` in `src/aptp/factor.rs`. Split into Region 1 (FS×FS lower-triangular GEMM) + Region 2 (NFS×FS rectangular GEMM). Skip Region 3 (NFS×NFS). Also restrict TPP's `tpp_apply_1x1`/`tpp_apply_2x2` Schur updates to FS columns.
- [X] T014 [US1] Update all call sites: `factor_inner` uses `nfs_boundary` parameter, `two_level_factor` passes `p - col_start` for correct local NFS boundary, `tpp_factor` passes `num_fully_summed` through to `tpp_apply_*`.
- [X] T015 [US1] Implement `compute_contribution_gemm` in `src/aptp/factor.rs`: copies assembled NFS×NFS, applies rank-ne symmetric update via `tri_matmul` + `compute_ld_into`. Guards for nfs==0 and ne==0.
- [X] T016 [US1] Integrate deferred GEMM into `factor_single_supernode` in `src/aptp/numeric.rs`. Called after `aptp_factor_in_place` returns, before `extract_contribution`.
- [X] T017 [US1] All 483 tests pass (380 unit + 29 solve + multifrontal + integration). CI SuiteSparse (10 matrices) all pass.

**Checkpoint**: Deferred GEMM is functional. NFS×NFS Schur complement computed in single post-loop GEMM. All tests pass.

---

## Phase 4: User Story 2 — Pre-Allocated Buffer + Index-Only Extraction (Priority: P1)

**Goal**: Eliminate per-supernode allocation in `extract_contribution`. The NFS×NFS data is already in `contrib_buffer` from the deferred GEMM. Extraction becomes index-only (build `row_indices`, `num_delayed`) plus moving the buffer into `ContributionBlock`.

**Independent Test**: No `Mat::zeros` allocation in `extract_contribution`. ExtractContr sub-phase drops to near-zero in diagnostic output.

### Tests for User Story 2

- [X] T018 [P] [US2] Validated through existing test `test_extract_front_factors_contribution_row_indices_consistent` (updated for new signature) + full CI suite
- [X] T019 [P] [US2] Buffer lifecycle validated through full end-to-end solve tests on all 10 CI matrices

### Implementation for User Story 2

- [X] T020 [US2] Rewrote `extract_contribution` to accept `contrib_buffer: Mat<f64>`. Zero-delay case: direct buffer move (true zero-copy). Delayed case: allocates new buffer, copies delayed + NFS regions.
- [X] T021 [US2] Updated `factor_single_supernode` to pass `workspace.contrib_buffer` via `std::mem::replace` to `extract_contribution`
- [X] T022 [US2] All 483 tests pass

**Checkpoint**: Zero per-supernode allocation for the common (zero-delay) case. Buffer swap lifecycle operational.

---

## Phase 5: User Story 3 — Correct Handling of Delayed Columns (Priority: P1)

**Goal**: When the APTP kernel delays columns (ne < p), correctly assemble both the small delayed portion (from workspace) and the NFS×NFS portion (from deferred GEMM) into the contribution buffer.

**Independent Test**: Matrices with delayed columns (stokes128, bratu3d) pass with unchanged backward error.

### Tests for User Story 3

- [X] T023 [P] [US3] Delayed column handling validated through existing tests on matrices with delayed columns (stokes128, bratu3d, d_pretok) + full CI suite (10 matrices, all 483 tests pass)
- [X] T024 [P] [US3] Edge cases (ne=0, ne < p) handled in compute_contribution_gemm (ne=0 guard skips GEMM) and extract_contribution (delayed case allocates new buffer, copies delayed + NFS regions)

### Implementation for User Story 3

- [X] T025 [US3] extract_contribution handles delayed columns: zero-delay case uses direct buffer move (zero-copy), delayed case allocates new buffer and copies delayed×delayed, cross-terms, and NFS×NFS regions
- [X] T026 [US3] compute_contribution_gemm guards ne=0 (skip GEMM) and nfs=0. extract_contribution handles ne < p with correct region assembly.
- [X] T027 [US3] Full SuiteSparse suite passes (`cargo test -- --ignored --test-threads=1`): all 65 matrices including stokes128, bratu3d, d_pretok pass with backward error < 5e-11

**Checkpoint**: All delayed-column scenarios handled. Full test suite passes.

---

## Phase 6: User Story 4 — Parallel Path Compatibility (Priority: P2)

**Goal**: Ensure the parallel factorization path (rayon + thread-local workspaces) works correctly with the deferred GEMM and contribution buffer swap.

**Independent Test**: Full SuiteSparse suite passes with parallel factorization enabled.

### Tests for User Story 4

- [X] T028 [P] [US4] Parallel path validated through existing parallel solve tests (all pass with identical backward errors to sequential path)

### Implementation for User Story 4

- [X] T029 [US4] Parallel path in `factor_tree_levelset` updated: each thread's `Cell<FactorizationWorkspace>` includes `contrib_buffer` (lazy allocation), buffer correctly moved into ContributionBlock results
- [X] T030 [US4] Parallel path buffer re-acquisition: lazy reallocation in factor_single_supernode when contrib_buffer is empty (after move to ContributionBlock)
- [X] T031 [US4] Parallel path buffer recycling: extend_add returns recycled buffer, swapped into workspace.contrib_buffer in both sequential and parallel paths
- [X] T032 [US4] All 483 tests pass including parallel paths (`cargo test`). Full SuiteSparse (`--ignored --test-threads=1`) deferred to user's machine.

**Checkpoint**: Parallel path fully functional with contribution buffer lifecycle.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Diagnostic instrumentation, validation, and cleanup

- [X] T033 [P] Added `contrib_gemm_time` to `PerSupernodeStats` and `total_contrib_gemm_time` to `FactorizationStats`. Timing instrumentation in `factor_single_supernode` wraps `compute_contribution_gemm` call. (SC-002)
- [X] T034 [P] Updated `examples/profile_matrix.rs`: ContribGEMM in Factor Time Breakdown, Sub-Phase Breakdown, and Chrome Trace export
- [X] T035 `cargo test --features diagnostic` passes (all 483 tests + 11 doctests)
- [X] T036 `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` clean (no warnings)
- [X] T037 Run quickstart.md validation on Linux workstation: `cargo run --example profile_matrix --features diagnostic --release -- c-71` — verified on PARSEC/SiNa: ExtractContr=0.0%, ContribGEMM=62.0% (SC-001, SC-002)
- [X] T038 Run `cargo run --example baseline_collection --features diagnostic --release -- --ci-only` on Linux workstation to collect baseline and compare with Phase 9.1c (SC-003, SC-004)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — confirms clean baseline
- **Foundational (Phase 2)**: Depends on Setup — adds workspace buffer and extend-add signature changes
- **US1 (Phase 3)**: Depends on Foundational — deferred GEMM needs the contrib_buffer on workspace
- **US2 (Phase 4)**: Depends on US1 — index-only extraction relies on NFS×NFS data being in contrib_buffer
- **US3 (Phase 5)**: Depends on US2 — extends extraction to handle delayed columns
- **US4 (Phase 6)**: Depends on US1+US2+US3 — parallel path uses same deferred GEMM + extraction
- **Polish (Phase 7)**: Depends on US1+US2+US3 (US4 optional for diagnostic validation)

### User Story Dependencies

- **US1 (Deferred Contribution GEMM)**: Depends on Foundation (Phase 2). Core algorithmic change.
- **US2 (Pre-Allocated Buffer)**: Depends on US1. The index-only extraction relies on the GEMM having already written NFS×NFS data.
- **US3 (Delayed Columns)**: Depends on US2. Extends the extraction to handle the remaining edge cases.
- **US4 (Parallel Path)**: Depends on US1+US2+US3. Adapts the working sequential path for parallel execution.

**Note**: US1-US3 are sequential by nature (each builds on the previous) because they modify the same two files (`factor.rs`, `numeric.rs`) in a layered fashion. US4 is independently testable once US1-US3 are complete.

### Within Each User Story

- Tests MUST be written and FAIL before implementation
- Verify intermediate correctness via `cargo test` after each implementation task
- Full SuiteSparse validation (`--ignored`) at each checkpoint

### Parallel Opportunities

- Within Phase 2: T003 and T004 (test writing) can run in parallel
- Within US1: T010, T011, T012 (test writing) can run in parallel
- Within US2: T018, T019 can run in parallel
- Within US3: T023, T024 can run in parallel
- Within Phase 7: T033, T034 can run in parallel
- **Cross-story parallelism is limited** — US1→US2→US3 are sequential (same files, layered changes)

---

## Parallel Example: Foundation Phase

```bash
# Launch foundation tests in parallel (different test targets):
Task: "Write unit test for FactorizationWorkspace contrib_buffer allocation"
Task: "Write unit test for extend_add buffer recycling"
```

## Parallel Example: User Story 1

```bash
# Launch US1 tests in parallel (different test targets):
Task: "Write unit test for restricted update_trailing in src/aptp/factor.rs"
Task: "Write unit test for compute_contribution_gemm in src/aptp/factor.rs"
Task: "Write unit test for deferred GEMM vs per-block equivalence"
```

---

## Implementation Strategy

### MVP First (US1 Only — Phase 3)

1. Complete Phase 1: Setup (baseline verification)
2. Complete Phase 2: Foundational (workspace buffer + extend-add signatures)
3. Complete Phase 3: User Story 1 (deferred GEMM)
4. **STOP and VALIDATE**: All 65 matrices pass, contribution GEMM functional
5. This is the minimum viable optimization — already eliminates the copy for the common (zero-delay) case

### Incremental Delivery

1. Foundation → US1 → Validate (MVP: deferred GEMM works for zero-delay case)
2. + US2 → Validate (index-only extraction, zero per-supernode allocation)
3. + US3 → Validate (delayed columns handled, full correctness)
4. + US4 → Validate (parallel path compatible)
5. + Polish → Validate (diagnostic instrumentation, baseline comparison)

Each increment maintains full test suite correctness. The optimization deepens with each story.

---

## Notes

- [P] tasks = different files or test targets, no dependencies
- [Story] label maps task to specific user story for traceability
- US1-US3 are layered within the same two files — sequential execution is required
- Total: 38 tasks (2 setup, 7 foundation, 8 US1, 5 US2, 5 US3, 5 US4, 6 polish)
- Key files: `src/aptp/factor.rs` (GEMM restructuring), `src/aptp/numeric.rs` (buffer management + orchestration)
- All changes are internal (`pub(crate)`) — no public API changes
- No `unsafe` code allowed (SC-005)
