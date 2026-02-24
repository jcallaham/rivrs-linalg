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

- [ ] T001 Run full test suite (`cargo test`) and record baseline pass count to confirm clean starting state
- [ ] T002 Run `cargo test -- --ignored --test-threads=1` to confirm all 65 SuiteSparse matrices pass with backward error < 5e-11

**Checkpoint**: Baseline correctness confirmed. All subsequent changes must maintain this.

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Add contribution buffer to workspace and update extend-add signatures — infrastructure that all user stories depend on

**Note**: These correspond to plan Layers 1 and 5 (workspace extension + buffer return). They must be in place before the deferred GEMM (US1) or index-only extraction (US2) can work.

### Tests for Foundation

- [ ] T003 [P] Write unit test in `src/aptp/numeric.rs` (or `tests/`) verifying `FactorizationWorkspace::new` allocates `contrib_buffer` with dimensions >= `max_front × max_front`, and `ensure_capacity` grows it when needed
- [ ] T004 [P] Write unit test verifying `extend_add` returns a `Mat<f64>` buffer after consuming a `ContributionBlock`, and that the returned buffer has the correct dimensions

### Implementation for Foundation

- [ ] T005 Add `contrib_buffer: Mat<f64>` field to `FactorizationWorkspace` in `src/aptp/numeric.rs`. Initialize in `new()` to `Mat::zeros(max_front, max_front)`. Grow in `ensure_capacity()` if needed. Also handle the runtime fallback (FR-009): if `contrib_buffer` is too small for a specific supernode's contribution (e.g., after delayed-column cascades inflate the contribution beyond `max_front`), resize before the deferred GEMM — matching the frontal workspace's existing `ensure_capacity` pattern. (Contract C1, C2)
- [ ] T006 Modify `extend_add` in `src/aptp/numeric.rs` to take `child: ContributionBlock` (owned, was `&ContributionBlock`) and return `Mat<f64>` (the consumed data buffer for recycling). (Contract C6)
- [ ] T007 Modify `extend_add_mapped` in `src/aptp/numeric.rs` with same ownership transfer and buffer return as T006. (Contract C7)
- [ ] T008 Update all call sites of `extend_add` and `extend_add_mapped` in `src/aptp/numeric.rs` (`factor_tree_levelset` sequential and parallel paths) to pass owned `ContributionBlock` and swap returned buffer into `workspace.contrib_buffer`. Handle the parallel path's thread-local workspace via existing `Cell<FactorizationWorkspace>` pattern.
- [ ] T009 Run `cargo test` to verify no regressions from signature changes (all 358 tests pass, backward errors unchanged)

**Checkpoint**: Foundation ready — workspace has contribution buffer, extend-add recycles buffers. Existing behavior preserved.

---

## Phase 3: User Story 1 — Deferred Contribution GEMM (Priority: P1) MVP

**Goal**: Restructure per-block trailing update to skip NFS×NFS region, then compute the NFS×NFS Schur complement in a single post-loop GEMM directly into the contribution buffer.

**Independent Test**: All 65 SuiteSparse matrices pass with backward error < 5e-11. Diagnostic profiling shows new ContribGEMM sub-phase.

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T010 [P] [US1] Write unit test for restricted `update_trailing` in `src/aptp/factor.rs`: construct a small frontal matrix (e.g., 8×8 with p=3), run one block's trailing update, verify NFS×NFS region (`A[p..m, p..m]`) is unchanged while FS×FS and cross-terms are updated correctly
- [ ] T011 [P] [US1] Write unit test for `compute_contribution_gemm` in `src/aptp/factor.rs`: construct a small frontal matrix with known L21, D, and assembled NFS×NFS values, invoke the function, verify `contrib_buffer` contains `assembled - L21 * D * L21^T` (lower triangle) with exact expected values
- [ ] T012 [P] [US1] Write unit test verifying the deferred GEMM produces the same NFS×NFS Schur complement as the current per-block approach on a hand-constructed 12×12 frontal matrix with 4 fully-summed columns and 2 blocks

### Implementation for User Story 1

- [ ] T013 [US1] Add `num_fully_summed: usize` parameter to `update_trailing` in `src/aptp/factor.rs`. Restrict the lower-triangular GEMM to `A[ts..p, ts..p]` (FS×FS, region 1) and add rectangular GEMM on `A[p..m, ts..p]` (NFS×FS cross-term, region 2). Skip `A[p..m, p..m]` (region 3). (Contract C3)
- [ ] T014 [US1] Update all call sites of `update_trailing` in `src/aptp/factor.rs` (`factor_inner`, `two_level_factor`, `tpp_factor_as_primary`) to pass `num_fully_summed`
- [ ] T015 [US1] Implement `compute_contribution_gemm` function in `src/aptp/factor.rs` per Contract C4: copies assembled NFS×NFS from `frontal_data[p..m, p..m]` into `contrib_buffer`, then applies rank-`ne` update in-place (`contrib -= L21_NFS * D * L21_NFS^T`). Handle D's mixed 1×1/2×2 structure via `MixedDiagonal` multiply. Use `faer::linalg::matmul::triangular::matmul` for symmetric output. Guard: skip entirely when `nfs == 0` (no contribution, FR-010) or `ne == 0` (copy only, no GEMM — handled in US3/T026).
- [ ] T016 [US1] Integrate the deferred GEMM into `factor_single_supernode` in `src/aptp/numeric.rs` (Contract C8): after `aptp_factor_in_place` returns, call `compute_contribution_gemm` passing the workspace frontal data and `workspace.contrib_buffer` as the target. The function is self-contained — it copies assembled NFS×NFS from workspace into `contrib_buffer` and applies the rank-`ne` update in-place.
- [ ] T017 [US1] Run `cargo test` to verify all unit tests pass, then `cargo test -- --ignored --test-threads=1` for full SuiteSparse validation

**Checkpoint**: Deferred GEMM is functional. NFS×NFS Schur complement computed in single post-loop GEMM. All tests pass.

---

## Phase 4: User Story 2 — Pre-Allocated Buffer + Index-Only Extraction (Priority: P1)

**Goal**: Eliminate per-supernode allocation in `extract_contribution`. The NFS×NFS data is already in `contrib_buffer` from the deferred GEMM. Extraction becomes index-only (build `row_indices`, `num_delayed`) plus moving the buffer into `ContributionBlock`.

**Independent Test**: No `Mat::zeros` allocation in `extract_contribution`. ExtractContr sub-phase drops to near-zero in diagnostic output.

### Tests for User Story 2

- [ ] T018 [P] [US2] Write unit test for index-only `extract_contribution` in `src/aptp/numeric.rs`: given a workspace with NFS×NFS data already in `contrib_buffer` (zero delayed columns case), verify `extract_contribution` returns a `ContributionBlock` whose `data` is the moved-in buffer, `row_indices` is correct, and `num_delayed` is 0 — with no new allocation
- [ ] T019 [P] [US2] Write unit test verifying the contribution buffer swap lifecycle: workspace starts with buffer → deferred GEMM fills it → `extract_contribution` moves it into `ContributionBlock` → `extend_add` recycles it back to workspace → workspace has buffer for next supernode

### Implementation for User Story 2

- [ ] T020 [US2] Rewrite `extract_contribution` in `src/aptp/numeric.rs` per Contract C5: accept `contrib_buffer: Mat<f64>` as parameter (moved in, already containing NFS×NFS data from deferred GEMM). Build `row_indices` and compute `num_delayed`. For zero-delay case, move buffer directly into `ContributionBlock`. Remove the `Mat::zeros` allocation and column-by-column copy loop.
- [ ] T021 [US2] Update `factor_single_supernode` in `src/aptp/numeric.rs` to pass `workspace.contrib_buffer` (via `std::mem::take` or similar) to the rewritten `extract_contribution`
- [ ] T022 [US2] Run `cargo test` and `cargo test -- --ignored --test-threads=1` to verify all tests pass

**Checkpoint**: Zero per-supernode allocation for the common (zero-delay) case. Buffer swap lifecycle operational.

---

## Phase 5: User Story 3 — Correct Handling of Delayed Columns (Priority: P1)

**Goal**: When the APTP kernel delays columns (ne < p), correctly assemble both the small delayed portion (from workspace) and the NFS×NFS portion (from deferred GEMM) into the contribution buffer.

**Independent Test**: Matrices with delayed columns (stokes128, bratu3d) pass with unchanged backward error.

### Tests for User Story 3

- [ ] T023 [P] [US3] Write unit test for `extract_contribution` with delayed columns in `src/aptp/numeric.rs`: construct a scenario where ne < p (e.g., ne=2, p=4, m=10), verify the contribution buffer contains delayed×delayed from workspace in `[0..2, 0..2]`, NFS×delayed cross-terms in `[2..8, 0..2]`, and NFS×NFS from deferred GEMM in `[2..8, 2..8]`
- [ ] T024 [P] [US3] Write unit test for the all-columns-delayed edge case (ne=0): verify the deferred GEMM is a no-op (rank-0 update), and the entire contribution is copied from the workspace

### Implementation for User Story 3

- [ ] T025 [US3] Extend `extract_contribution` in `src/aptp/numeric.rs` to handle delayed columns: when `ne < p`, copy the delayed-column region (`A[ne..p, ne..p]`) and cross-terms (`A[p..m, ne..p]`) from the workspace into the contribution buffer at the correct positions (`[0..num_delayed, 0..num_delayed]` and `[num_delayed..size, 0..num_delayed]`)
- [ ] T026 [US3] Handle the ne=0 edge case in `compute_contribution_gemm` (skip GEMM when ne=0) and in `extract_contribution` (copy entire trailing submatrix from workspace)
- [ ] T027 [US3] Run full SuiteSparse suite (`cargo test -- --ignored --test-threads=1`) specifically checking matrices known to produce delayed columns (stokes128, bratu3d, d_pretok) — verify backward errors unchanged

**Checkpoint**: All delayed-column scenarios handled. Full test suite passes.

---

## Phase 6: User Story 4 — Parallel Path Compatibility (Priority: P2)

**Goal**: Ensure the parallel factorization path (rayon + thread-local workspaces) works correctly with the deferred GEMM and contribution buffer swap.

**Independent Test**: Full SuiteSparse suite passes with parallel factorization enabled.

### Tests for User Story 4

- [ ] T028 [P] [US4] Write integration test verifying parallel factorization produces identical backward errors to sequential path on a matrix that exercises both tree-level and intra-node parallelism (e.g., c-71 or bcsstk18)

### Implementation for User Story 4

- [ ] T029 [US4] Review and update the parallel path in `factor_tree_levelset` in `src/aptp/numeric.rs`: ensure each thread's `Cell<FactorizationWorkspace>` includes `contrib_buffer`, and that the buffer is correctly moved into `ContributionBlock` results before cross-thread transfer
- [ ] T030 [US4] Handle the parallel path buffer re-acquisition: when a thread's `contrib_buffer` has been moved into a `ContributionBlock` for cross-thread transfer, allocate a new buffer for the next supernode (lazy reallocation, matching R4 design)
- [ ] T031 [US4] Handle the parallel path buffer recycling: when the parent thread's `extend_add` returns a recycled buffer from a child contribution, swap it into the parent's workspace `contrib_buffer` (reusing the allocation)
- [ ] T032 [US4] Run full test suite including parallel paths: `cargo test` and `cargo test -- --ignored --test-threads=1`

**Checkpoint**: Parallel path fully functional with contribution buffer lifecycle.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Diagnostic instrumentation, validation, and cleanup

- [ ] T033 [P] Add `ContribGEMM` sub-phase timing to diagnostic instrumentation in `src/aptp/numeric.rs` (behind `diagnostic` feature flag). Time the `compute_contribution_gemm` call in `factor_single_supernode`. (SC-002)
- [ ] T034 [P] Update `examples/profile_matrix.rs` to display the new `ContribGEMM` sub-phase in diagnostic output
- [ ] T035 Run `cargo test --features diagnostic` to verify diagnostic builds pass
- [ ] T036 Run `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` — fix any warnings
- [ ] T037 Run quickstart.md validation: `cargo run --example profile_matrix --features diagnostic --release -- c-71` and verify ExtractContr near zero and ContribGEMM visible (SC-001, SC-002)
- [ ] T038 Run `cargo run --example baseline_collection --features diagnostic --release -- --ci-only` to collect baseline and compare with Phase 9.1c (SC-003, SC-004)

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
