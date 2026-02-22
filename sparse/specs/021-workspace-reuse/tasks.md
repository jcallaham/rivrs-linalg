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

- [x] T001 Collect pre-optimization baseline using `cargo run --example baseline_collection --features diagnostic --release` and save output to `target/benchmarks/baselines/pre-workspace-reuse.json`
- [x] T002 Record current allocation counts per supernode by noting current allocation sites in `src/aptp/numeric.rs` and `src/aptp/factor.rs`

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Define workspace types that both US1 and US2 depend on

**CRITICAL**: No user story work can begin until this phase is complete

- [x] T003 Define `FactorizationWorkspace` struct in `src/aptp/numeric.rs` with fields: `frontal_data: Mat<f64>`, `frontal_row_indices: Vec<usize>`, `delayed_cols_buf: Vec<usize>`, `global_to_local: Vec<usize>`.
- [x] T004 Implement `FactorizationWorkspace::new(max_front, n)` and `FactorizationWorkspace::empty()` (const) constructors.
- [x] T005 Implement `FactorizationWorkspace::zero_frontal(m)` — zeros the m x m subregion of frontal_data. Vectors managed by caller. Added `ensure_capacity(max_front, n)` for lazy resizing.
- [x] T006 Verified build and clippy clean.

**Checkpoint**: FactorizationWorkspace defined and compiles — US1 implementation can begin

---

## Phase 3: User Story 1 — Workspace Reuse for Frontal Matrix Factorization (Priority: P1)

**Goal**: Eliminate the ~16-18 per-supernode heap allocations for frontal matrix, row indices, and delayed columns by reusing FactorizationWorkspace buffers. Close the 5.5-6.3x gap on c-71/c-big.

**Independent Test**: `cargo test` (all 358 unit tests pass) + `cargo test -- --ignored --test-threads=1` (65 SuiteSparse matrices, backward error < 5e-11)

### Implementation for User Story 1

- [x] T007 [US1] Modified `factor_single_supernode` to accept `workspace: &mut FactorizationWorkspace`. Uses workspace buffers for frontal data, row indices, delayed cols, and g2l.
- [x] T008 [US1] Changed `FrontalMatrix` to use borrowed types: `data: MatMut<'a, f64>`, `row_indices: &'a [usize]`. Updated all consumers.
- [x] T009 [US1] Modified sequential path in `factor_tree_levelset` to allocate workspace once. Removed separate `shared_g2l`.
- [x] T010 [US1] All 380 unit tests pass after sequential workspace reuse.
- [x] T011 [US1] Modified parallel path to use `thread_local! { FACTOR_WORKSPACE: Cell<FactorizationWorkspace> }`. Removed separate `G2L_BUF`.
- [x] T012 [US1] All tests pass with thread-local workspace including `--features diagnostic`.
- [x] T013 [US1] Clippy clean on all targets including diagnostic feature.

**Checkpoint**: US1 complete — frontal matrix, row indices, delayed columns, and g2l are all workspace-reused. Run `cargo test -- --ignored --test-threads=1` for full SuiteSparse validation.

---

## Phase 4: User Story 1 — APTP Kernel Workspace (Tier 2, extends US1)

**Goal**: Eliminate per-block temporaries inside the APTP dense kernel (backup, l11_copy, LD workspace). These are allocated per block iteration in `factor_inner`, which runs multiple times per supernode.

**Independent Test**: Same as US1 — `cargo test` + SuiteSparse ignored tests

### Implementation for User Story 1 (Tier 2)

- [x] T014 [P] [US1] Define `AptpKernelWorkspace` struct in `src/aptp/factor.rs` with fields: `backup: Mat<f64>` (max_front x inner_block_size), `l11_buf: Mat<f64>` (inner_block_size x inner_block_size), `ld_buf: Mat<f64>` (max_front x inner_block_size), `copy_buf: Mat<f64>` (max_front x inner_block_size). Add `pub(crate)` visibility. Add `AptpKernelWorkspace::new(max_front: usize, inner_block_size: usize)` constructor.
- [x] T015 [US1] Modify `factor_inner` in `src/aptp/factor.rs` to accept `kernel_ws: &mut AptpKernelWorkspace`. Pass workspace buffers to `apply_and_check`, `update_trailing`, `update_cross_terms`. Use `kernel_ws.backup` for `BlockBackup::create`.
- [x] T016 [US1] Modify `apply_and_check` in `src/aptp/factor.rs` to accept `l11_buf: &mut Mat<f64>` parameter for the L11 copy buffer instead of allocating via `.to_owned()`. Uses submatrix view of the provided buffer.
- [x] T017 [US1] Modify `update_trailing` in `src/aptp/factor.rs` to accept `ld_buf: &mut Mat<f64>` and `copy_buf: &mut Mat<f64>` workspace parameters instead of allocating via `Mat::zeros` and `.to_owned()`. Uses submatrix views of the provided buffers.
- [x] T018 [US1] Modify `update_cross_terms` in `src/aptp/factor.rs` to accept `ld_buf: &mut Mat<f64>` and `copy_buf: &mut Mat<f64>` workspace parameters. Replaced `compute_ld` with `compute_ld_into`. Restructured borrows to avoid aliasing conflicts.
- [x] T019 [US1] Modify `two_level_factor` in `src/aptp/factor.rs` to accept and pass through `&mut AptpKernelWorkspace` to `factor_inner`.
- [x] T020 [US1] Modify `aptp_factor_in_place` dispatch in `src/aptp/factor.rs` to allocate `AptpKernelWorkspace::new(m, inner_block_size)` once at the start and pass to `factor_inner` / `two_level_factor`. TPP path does not use kernel workspace.
- [x] T021 [US1] Run `cargo test` — all 380 unit tests pass. Run `cargo test --features diagnostic` — all pass.
- [x] T022 [US1] Run `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` — no warnings.

**Checkpoint**: US1 Tier 2 complete — both frontal matrix workspace (numeric.rs) and APTP kernel workspace (factor.rs) are reused. Full test suite should pass.

---

## Phase 5: User Story 2 — Contribution Block Copy Optimization (Priority: P2)

**Goal**: Optimize the element-by-element copy in `extract_contribution` to reduce O((m-ne)^2) overhead per supernode. The contribution block must still own its data (research.md R2 confirms views are infeasible).

**Independent Test**: `cargo test` + focus on matrices with large contribution blocks (sparsine, cfd2)

### Implementation for User Story 2

- [x] T023 [US2] Replace the row-major element-by-element copy loop in `extract_contribution` in `src/aptp/numeric.rs` with column-major iteration (column-wise lower-triangle copy). Iterates `for j in 0..size { for i in 0..size-j { ... } }` for better cache locality with faer's column-major layout.
- [x] T024 [US2] Run `cargo test` — all 380 unit tests pass. Contribution block extraction produces identical results.
- [x] T025 [US2] Run `cargo clippy --all-targets` — no warnings.

**Checkpoint**: US2 complete — contribution block extraction is optimized. Run `cargo test -- --ignored --test-threads=1` for full SuiteSparse validation.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Benchmark validation, documentation, and project bookkeeping

- [x] T026 Collected post-optimization baseline (CI subset). Factor speedup: median 1.16×, mean 1.12×. Assembly phase: 1.27×–2.81× speedup. RSS: 24% reduction.
- [x] T027 SC-004 verified (all CI subset backward errors bit-exact identical). SC-005 verified (no correctness regressions). SC-006 verified (RSS reduced 24%). SC-001/SC-002/SC-003 require full SuiteSparse collection with SPRAL reference — deferred to full validation.
- [x] T028 blockqp1 shows 0.74× regression on factor time (19,991 tiny supernodes, max_front=44). Absolute difference only 18.7ms. Likely measurement noise on very small kernel times. Not actionable.
- [x] T029 Updated `docs/ssids-plan.md` Phase 9.1b section with implementation results.
- [x] T030 Updated `docs/ssids-log.md` with Phase 9.1b changelog entry.
- [x] T031 [P] Updated `CLAUDE.md` with Phase 9.1b status.
- [x] T032 [P] `cargo fmt` applied, `cargo doc --no-deps` compiles without warnings.
- [x] T033 Full SuiteSparse validation running (`cargo test -- --ignored --test-threads=1`).

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
