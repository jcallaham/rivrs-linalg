# Tasks: Workspace Reuse & Per-Supernode Allocation Optimization

**Input**: Design documents from `/specs/021-workspace-reuse/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md

**Tests**: Not explicitly requested in the spec — existing test suite provides regression coverage. No new test tasks generated. All implementation tasks must pass `cargo test` and `cargo test -- --ignored --test-threads=1` after completion.

**Organization**: Tasks grouped by user story (US1 = workspace reuse, US2 = contribution block optimization), with shared setup phase and final benchmark validation.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2)
- Include exact file paths in descriptions

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Collect pre-optimization baseline measurements for comparison after optimization

- [ ] T001 Collect pre-optimization baseline using `cargo run --example baseline_collection --features diagnostic --release` and save output to `target/benchmarks/baselines/pre-workspace-reuse.json`
- [ ] T002 Record current allocation counts per supernode by adding temporary instrumentation or noting current allocation sites in `src/aptp/numeric.rs` and `src/aptp/factor.rs` for comparison

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Define workspace types that both US1 and US2 depend on

**CRITICAL**: No user story work can begin until this phase is complete

- [ ] T003 Define `FactorizationWorkspace` struct in `src/aptp/numeric.rs` with fields: `frontal_data: Mat<f64>`, `frontal_row_indices: Vec<usize>`, `delayed_cols_buf: Vec<usize>`, `global_to_local: Vec<usize>`. Add `pub(crate)` visibility. Include rustdoc with algorithm reference (Duff et al. 2020, Liu 1992 — workspace reuse strategy). Document the workspace invariant: never shared between concurrent supernodes.
- [ ] T004 Implement `FactorizationWorkspace::new(max_front: usize, n: usize)` constructor in `src/aptp/numeric.rs`. Allocate `frontal_data` as `Mat::zeros(max_front, max_front)`, `frontal_row_indices` with capacity `max_front`, `delayed_cols_buf` with capacity `max_front`, `global_to_local` as `vec![NOT_IN_FRONT; n]`. Handle degenerate case `max_front == 0` gracefully (empty buffers).
- [ ] T005 Implement `FactorizationWorkspace::prepare_for_supernode(&mut self, m: usize)` in `src/aptp/numeric.rs`. Zero the lower triangle of the m x m subregion of `frontal_data` using column-wise `fill_zero()` on the submatrix view. Clear `frontal_row_indices` and `delayed_cols_buf` (truncate to length 0, preserving capacity). Assert `m <= self.frontal_data.nrows()`.
- [ ] T006 Verify foundational types compile: run `cargo build` and `cargo clippy --all-targets` with no errors or warnings in `src/aptp/numeric.rs`

**Checkpoint**: FactorizationWorkspace defined and compiles — US1 implementation can begin

---

## Phase 3: User Story 1 — Workspace Reuse for Frontal Matrix Factorization (Priority: P1)

**Goal**: Eliminate the ~16-18 per-supernode heap allocations for frontal matrix, row indices, and delayed columns by reusing FactorizationWorkspace buffers. Close the 5.5-6.3x gap on c-71/c-big.

**Independent Test**: `cargo test` (all 358 unit tests pass) + `cargo test -- --ignored --test-threads=1` (65 SuiteSparse matrices, backward error < 5e-11)

### Implementation for User Story 1

- [ ] T007 [US1] Modify `factor_single_supernode` signature in `src/aptp/numeric.rs` to accept `workspace: &mut FactorizationWorkspace` instead of `global_to_local: &mut Vec<usize>`. Use `workspace.frontal_data.as_mut().submatrix_mut(0, 0, m, m)` as the frontal matrix data instead of `Mat::zeros(m, m)`. Use `workspace.frontal_row_indices` instead of allocating a new Vec for row indices. Use `workspace.delayed_cols_buf` instead of allocating a new Vec for delayed columns. Use `workspace.global_to_local` instead of the separate g2l parameter. Call `workspace.prepare_for_supernode(m)` at the start. Reset `global_to_local` entries at the end (existing cleanup logic).
- [ ] T008 [US1] Update `FrontalMatrix` usage in `factor_single_supernode` in `src/aptp/numeric.rs` so that `FrontalMatrix.data` references the workspace submatrix view rather than owning a fresh `Mat`. This may require changing `FrontalMatrix` to use `MatMut<'_, f64>` instead of `Mat<f64>`, or restructuring the code to pass the workspace's `MatMut` directly to `scatter_original_entries_multi`, `extend_add`, and `aptp_factor_in_place`. Ensure all three callees continue to work correctly with the workspace buffer.
- [ ] T009 [US1] Modify the sequential path in `factor_tree_levelset` in `src/aptp/numeric.rs` to allocate `FactorizationWorkspace::new(max_front_size, n)` once before the level-set loop and pass `&mut workspace` to each `factor_single_supernode` call. Remove the separate `shared_g2l` variable — it is now part of the workspace. Obtain `max_front_size` from `AptpSymbolic` (already available as a method or computable from supernode metadata).
- [ ] T010 [US1] Run `cargo test` to verify all unit tests pass after sequential workspace reuse refactoring in `src/aptp/numeric.rs`. Fix any compilation errors or test failures. Pay special attention to: hand-constructed matrix tests, reconstruction tests, and backward error tests.
- [ ] T011 [US1] Modify the parallel path in `factor_tree_levelset` in `src/aptp/numeric.rs` to use thread-local workspace. Add `thread_local! { static FACTOR_WORKSPACE: Cell<FactorizationWorkspace> = const { Cell::new(FactorizationWorkspace::empty()) }; }`. In the `par_iter` closure: take workspace via `FACTOR_WORKSPACE.take()`, resize if needed (lazy initialization matching the existing g2l pattern at line 948-950), pass to `factor_single_supernode`, then set back via `FACTOR_WORKSPACE.set(workspace)`. Add `FactorizationWorkspace::empty()` const constructor for the Cell default. Remove the separate `G2L_BUF` thread-local — fold its buffer into the workspace.
- [ ] T012 [US1] Run `cargo test` (including parallel factorization tests) to verify all tests pass with thread-local workspace. Verify with `cargo test --features diagnostic` to confirm diagnostic instrumentation still works (FR-010).
- [ ] T013 [US1] Run `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` to verify no new warnings. Fix any clippy issues in modified code in `src/aptp/numeric.rs`.

**Checkpoint**: US1 complete — frontal matrix, row indices, delayed columns, and g2l are all workspace-reused. Run `cargo test -- --ignored --test-threads=1` for full SuiteSparse validation.

---

## Phase 4: User Story 1 — APTP Kernel Workspace (Tier 2, extends US1)

**Goal**: Eliminate per-block temporaries inside the APTP dense kernel (backup, l11_copy, LD workspace). These are allocated per block iteration in `factor_inner`, which runs multiple times per supernode.

**Independent Test**: Same as US1 — `cargo test` + SuiteSparse ignored tests

### Implementation for User Story 1 (Tier 2)

- [ ] T014 [P] [US1] Define `AptpKernelWorkspace` struct in `src/aptp/factor.rs` with fields: `backup_data: Mat<f64>` (max_front x inner_block_size), `l11_temp: Mat<f64>` (inner_block_size x inner_block_size), `ld_workspace: Mat<f64>` (max_front x inner_block_size), `copy_workspace: Mat<f64>` (max_front x inner_block_size). Add `pub(crate)` visibility. Add `AptpKernelWorkspace::new(max_front: usize, inner_block_size: usize)` constructor. Include rustdoc referencing SPRAL's per-thread workspace pattern (NumericSubtree.hxx:75-81).
- [ ] T015 [US1] Modify `factor_inner` in `src/aptp/factor.rs` to accept `kernel_ws: &mut AptpKernelWorkspace`. Replace `BlockBackup::create` allocation (line ~1654) with subview into `kernel_ws.backup_data`. Replace `l11_copy` allocation in `apply_and_check` (line ~1287) with subview into `kernel_ws.l11_temp`. Replace `w` allocation in `update_trailing` (line ~1397) with subview into `kernel_ws.ld_workspace`. Replace `l21_copy` allocation in `update_trailing` (line ~1426) with subview into `kernel_ws.copy_workspace`.
- [ ] T016 [US1] Modify `apply_and_check` in `src/aptp/factor.rs` to accept a workspace `MatMut` parameter for the L11 copy buffer instead of allocating via `.to_owned()`. Use a submatrix view of the provided buffer sized to `block_nelim x block_nelim`.
- [ ] T017 [US1] Modify `update_trailing` in `src/aptp/factor.rs` to accept workspace `MatMut` parameters for the LD product workspace and the L21 copy buffer instead of allocating via `Mat::zeros` and `.to_owned()`. Use submatrix views of the provided buffers sized to `trailing_size x nelim`.
- [ ] T018 [US1] Modify `update_cross_terms` in `src/aptp/factor.rs` to accept workspace `MatMut` parameters for the L block copy and LD product workspace instead of allocating via `.to_owned()` and `compute_ld`. Use submatrix views of the provided buffers.
- [ ] T019 [US1] Modify `two_level_factor` in `src/aptp/factor.rs` to accept and pass through `&mut AptpKernelWorkspace` to `factor_inner`. The per-outer-block `temp` buffer (line ~1934) can also use workspace storage.
- [ ] T020 [US1] Modify `aptp_factor_in_place` dispatch in `src/aptp/factor.rs` to allocate `AptpKernelWorkspace` once at the start of the call and pass to `factor_inner` / `two_level_factor`. Alternatively, accept `Option<&mut AptpKernelWorkspace>` from the caller (numeric.rs could pass workspace from FactorizationWorkspace). The `tpp_factor_as_primary` path does not need kernel workspace (tpp_factor has no per-block temporaries).
- [ ] T021 [US1] Run `cargo test` to verify all unit tests pass after kernel workspace refactoring in `src/aptp/factor.rs`. Pay special attention to: factor_inner block loop tests, two_level_factor tests, BlockBackup restore-on-failure tests, and BLAS-3 pipeline correctness.
- [ ] T022 [US1] Run `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` to verify no new warnings in `src/aptp/factor.rs`.

**Checkpoint**: US1 Tier 2 complete — both frontal matrix workspace (numeric.rs) and APTP kernel workspace (factor.rs) are reused. Full test suite should pass.

---

## Phase 5: User Story 2 — Contribution Block Copy Optimization (Priority: P2)

**Goal**: Optimize the element-by-element copy in `extract_contribution` to reduce O((m-ne)^2) overhead per supernode. The contribution block must still own its data (research.md R2 confirms views are infeasible).

**Independent Test**: `cargo test` + focus on matrices with large contribution blocks (sparsine, cfd2)

### Implementation for User Story 2

- [ ] T023 [US2] Replace the element-by-element copy loop in `extract_contribution` in `src/aptp/numeric.rs` (lines ~1331-1335) with a bulk column-wise copy. Use faer's submatrix view: extract the trailing (m-ne) x (m-ne) lower-triangular submatrix from `frontal.data` and copy via column iteration with `copy_from` or `zipped` operations. Verify that only the lower triangle is copied (matching current behavior). If faer supports `submatrix(ne, ne, size, size).to_owned()` with correct lower-triangle semantics, use that instead.
- [ ] T024 [US2] Run `cargo test` to verify contribution block extraction produces identical results after the copy optimization. Pay special attention to: extend_add correctness (contribution data must be read correctly by parent), hand-constructed matrices with known contribution block values, and SuiteSparse CI subset matrices.
- [ ] T025 [US2] Run `cargo clippy --all-targets` to verify no new warnings in `src/aptp/numeric.rs`.

**Checkpoint**: US2 complete — contribution block extraction is optimized. Run `cargo test -- --ignored --test-threads=1` for full SuiteSparse validation.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Benchmark validation, documentation, and project bookkeeping

- [ ] T026 Collect post-optimization baseline using `cargo run --example baseline_collection --features diagnostic --release` and save to `target/benchmarks/baselines/post-workspace-reuse.json`. Compare against pre-optimization baseline from T001. Report per-matrix factor time ratios vs SPRAL.
- [ ] T027 Verify success criteria SC-001 through SC-006 against post-optimization baseline: (SC-001) c-71 and c-big within 2x SPRAL, (SC-002) simplicial matrices within 1.5x, (SC-003) median at or below 1.0x, (SC-004) all 65 backward errors < 5e-11, (SC-005) no regressions, (SC-006) no memory increase. Document results.
- [ ] T028 If any success criteria are not met after T027, profile remaining hotspots using `cargo run --example workload_analysis --features diagnostic --release` and document findings for Phase 9.1c. Identify which allocations remain as bottlenecks.
- [ ] T029 Update `docs/ssids-plan.md` Phase 9.1b section with actual results: supernode counts, factor times, SPRAL ratios, and status (COMPLETE or partial). Update success criteria checkboxes.
- [ ] T030 Update `docs/ssids-log.md` with Phase 9.1b changelog entry: what was built, key decisions, performance impact, any issues encountered.
- [ ] T031 [P] Update `CLAUDE.md` in `sparse/` with Phase 9.1b status and any new conventions (workspace reuse pattern, thread-local workspace naming).
- [ ] T032 [P] Run `cargo fmt --check` and `cargo doc --no-deps` to verify formatting and documentation compile without warnings.
- [ ] T033 Run full SuiteSparse validation: `cargo test -- --ignored --test-threads=1` with release profile. Verify all 65 matrices pass backward error < 5e-11. This is the final correctness gate.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: No strict dependency on Phase 1 (baseline collection is informational), but recommended to collect baseline first
- **US1 Tier 1 (Phase 3)**: Depends on Phase 2 (FactorizationWorkspace definition)
- **US1 Tier 2 (Phase 4)**: Depends on Phase 3 (frontal workspace must be working before adding kernel workspace)
- **US2 (Phase 5)**: Can start after Phase 2; independent of Phases 3-4. However, recommended to complete US1 first to measure remaining contribution block impact.
- **Polish (Phase 6)**: Depends on all desired user story phases being complete

### User Story Dependencies

- **US1 (P1)**: Can start after Foundational (Phase 2). Phases 3 and 4 are sequential (Tier 1 before Tier 2).
- **US2 (P2)**: Can start after Foundational (Phase 2). Independent of US1 but recommended to sequence after US1 to measure actual impact.

### Within Each Phase

- Tasks without [P] must run sequentially in listed order
- Tasks with [P] can run in parallel with other [P] tasks in the same phase
- Each phase ends with a checkpoint (test suite pass)

### Parallel Opportunities

- **T014** (AptpKernelWorkspace definition) can run in parallel with T007-T013 (different file: factor.rs vs numeric.rs)
- **T023** (contribution block optimization) can run in parallel with T014-T022 (different functions, no shared state)
- **T031 and T032** (documentation) can run in parallel

---

## Parallel Example: User Story 1

```bash
# Phase 3 and Phase 4 T014 can overlap:
# numeric.rs work (T007-T013) runs in parallel with:
Task: "T014 [P] [US1] Define AptpKernelWorkspace in src/aptp/factor.rs"
# Because T014 modifies factor.rs while T007-T013 modify numeric.rs

# Phase 6 documentation tasks can run in parallel:
Task: "T031 [P] Update CLAUDE.md"
Task: "T032 [P] Run cargo fmt and cargo doc checks"
```

---

## Implementation Strategy

### MVP First (US1 Tier 1 Only — Phases 1-3)

1. Complete Phase 1: Baseline collection
2. Complete Phase 2: Define FactorizationWorkspace
3. Complete Phase 3: Hoist frontal matrix + row indices + g2l into workspace
4. **STOP and VALIDATE**: Run full test suite + baseline comparison
5. If c-71/c-big show significant improvement, assess whether Tier 2 and US2 are needed

### Incremental Delivery

1. Phase 1-2: Setup + Foundational → workspace types ready
2. Phase 3: US1 Tier 1 → frontal matrix workspace reused → measure improvement
3. Phase 4: US1 Tier 2 → APTP kernel temporaries reused → measure incremental improvement
4. Phase 5: US2 → contribution block copy optimized → measure incremental improvement
5. Phase 6: Benchmark validation + documentation → ship

### Key Insight

The MVP (Phase 3 alone) eliminates the single largest allocation — the m x m frontal matrix — which accounts for the majority of allocation overhead per supernode. Phases 4-5 provide incremental gains. If Phase 3 achieves the performance targets, Phases 4-5 may be deferred.

---

## Notes

- All tasks operate on existing files — no new source files are created
- The public API (SparseLDLT) is NOT changed — all modifications are internal
- The constitution requires all existing tests to pass after each task (Principle I: Correctness First)
- The `diagnostic` feature must continue to compile and work correctly (FR-010)
- Commit after each completed phase checkpoint
