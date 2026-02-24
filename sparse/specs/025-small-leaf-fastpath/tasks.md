# Tasks: Small Leaf Subtree Fast Path

**Input**: Design documents from `/specs/025-small-leaf-fastpath/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md

**Tests**: Included per constitution (Principle III: TDD is NON-NEGOTIABLE).

**Organization**: Tasks grouped by user story. US2 (classification) is architecturally foundational for US1 (fast path), so it appears in the Foundational phase despite being P2 in the spec.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup

**Purpose**: Add configuration surface and data model changes needed by all subsequent phases.

- [X] T001 [P] Add `small_leaf_threshold: usize` field (default 256) to `FactorOptions` and its `Default` impl in `src/aptp/solver.rs`
- [X] T002 [P] Add `small_leaf_threshold: usize` field (default 256) to `SolverOptions` and wire it through to `FactorOptions` in `SparseLDLT::factor()` in `src/aptp/solver.rs`
- [X] T003 Add `in_small_leaf: bool` field (default `false`) to `SupernodeInfo` in `src/aptp/numeric.rs`. Initialize to `false` in `build_supernode_info()`.
- [X] T004 Define `SmallLeafSubtree` struct in `src/aptp/numeric.rs` with fields: `root: usize`, `nodes: Vec<usize>`, `max_front_size: usize`, `parent_of_root: Option<usize>`
- [X] T005 Verify build: `cargo build` and `cargo build --features diagnostic` pass with no warnings. Run `cargo test` to confirm no existing tests broken by the new fields.

**Checkpoint**: New types and configuration compile. All existing tests still pass.

---

## Phase 2: Foundational — Subtree Classification (US2)

**Purpose**: Implement the classification algorithm that identifies small-leaf subtrees. This MUST be complete before the fast-path factorization (US1) can begin.

**Goal (US2)**: Classify which post-amalgamation subtrees qualify as small-leaf based on front-size threshold. Classification runs as a cheap O(n_supernodes) pass at the start of each factorization.

**Independent Test**: Run classification on known tree structures and verify correct node marking and subtree identification.

### Tests for Classification

> **Write these tests FIRST, verify they FAIL, then implement.**

- [ ] T006 [P] [US2] Write unit test `test_classify_all_small` in `src/aptp/numeric.rs`: construct a 5-supernode linear chain (all front_size < 256), verify all marked `in_small_leaf = true` and one `SmallLeafSubtree` returned with 5 nodes in postorder.
- [ ] T007 [P] [US2] Write unit test `test_classify_mixed_tree` in `src/aptp/numeric.rs`: construct a tree where leaves are small but root is large (front_size >= 256). Verify leaves are `in_small_leaf = true`, root is `false`. Subtree excludes root.
- [ ] T008 [P] [US2] Write unit test `test_classify_single_node_excluded` in `src/aptp/numeric.rs`: construct a tree with isolated small leaf supernodes (no small siblings or parents). Verify no `SmallLeafSubtree` is returned (minimum 2 nodes required).
- [ ] T009 [P] [US2] Write unit test `test_classify_threshold_boundary` in `src/aptp/numeric.rs`: supernode with front_size exactly equal to threshold (256) must NOT qualify. Supernode with front_size 255 qualifies. Strict less-than comparison.
- [ ] T010 [P] [US2] Write unit test `test_classify_disabled` in `src/aptp/numeric.rs`: when `small_leaf_threshold = 0`, no subtrees are classified (empty `Vec<SmallLeafSubtree>`).
- [ ] T011 [P] [US2] Write unit test `test_classify_multiple_subtrees` in `src/aptp/numeric.rs`: construct a tree with two independent small-leaf subtrees. Verify both are returned, each with correct nodes and root.

### Implementation for Classification

- [ ] T012 [US2] Implement `classify_small_leaf_subtrees()` in `src/aptp/numeric.rs`. Takes `&mut [SupernodeInfo]`, `&[Vec<usize>]` (children_map), and `threshold: usize`. Bottom-up pass: compute front_size per node, mark `in_small_leaf` on SupernodeInfo. Identify subtree roots (in_small_leaf with non-in_small_leaf parent). Collect descendant nodes in postorder via iterative DFS. Filter subtrees with < 2 nodes. Return `Vec<SmallLeafSubtree>`.
- [ ] T013 [US2] Wire `classify_small_leaf_subtrees()` into `AptpNumeric::factor()` in `src/aptp/numeric.rs`: call after `amalgamate_supernodes()` and `build_children_map()`, before `build_assembly_maps()`. Pass `options.small_leaf_threshold`. Store returned `Vec<SmallLeafSubtree>` for use by `factor_tree_levelset()`.
- [ ] T014 [US2] Verify all classification tests pass: `cargo test classify`

**Checkpoint**: Classification correctly identifies small-leaf subtrees on synthetic tree structures. No factorization behavior changes yet — subtrees are identified but not yet used.

---

## Phase 3: User Story 1 — Fast-Path Factorization (Priority: P1) 🎯 MVP

**Goal**: Implement the streamlined factorization path for small-leaf subtrees and integrate it as a pre-pass in the level-set loop. Simplicial matrices achieve ≤1.5x SPRAL.

**Independent Test**: Factorize simplicial SuiteSparse matrices (dixmaanl, bloweybq, mario001) and verify backward error < 5e-11 with the fast path active.

### Tests for Fast-Path Factorization

> **Write these tests FIRST, verify they FAIL, then implement.**

- [ ] T015 [P] [US1] Write unit test `test_fast_path_small_chain` in `src/aptp/numeric.rs`: construct a small (n≈20) symmetric indefinite matrix that produces a linear chain of small supernodes. Factor with `small_leaf_threshold = 256`. Verify reconstruction error < 1e-12 and backward error < 5e-11.
- [ ] T016 [P] [US1] Write unit test `test_fast_path_matches_general` in `src/aptp/numeric.rs`: factor the same matrix twice — once with `small_leaf_threshold = 256` (fast path active) and once with `small_leaf_threshold = 0` (disabled). Verify both produce backward error < 5e-11. Pivot counts may differ slightly but both must be correct.
- [ ] T017 [P] [US1] Write unit test `test_fast_path_delayed_pivots` in `src/aptp/numeric.rs`: construct a matrix that forces delayed pivots in a small-leaf subtree (e.g., zero diagonal in one supernode). Verify delays propagate correctly to the parent (within or outside subtree) and backward error < 5e-11.
- [ ] T018 [P] [US1] Write unit test `test_fast_path_contribution_boundary` in `src/aptp/numeric.rs`: construct a mixed tree where a small-leaf subtree's root contribution feeds into a large general-path supernode. Verify the parent correctly assembles the contribution and backward error < 5e-11.
- [ ] T018b [P] [US1] Write unit test `test_fast_path_mc64_scaling` in `src/aptp/numeric.rs`: factor a matrix that requires MC64 scaling (e.g., one of the CI-subset hard indefinite matrices like bratu3d or stokes128) with `small_leaf_threshold = 256` and verify backward error < 5e-11. Factor the same matrix with `small_leaf_threshold = 0` and verify backward errors are comparable. This validates FR-009 (scaling applied identically in both paths).
- [ ] T019 [US1] Write integration test `test_fast_path_suitesparse_ci` in `tests/` or as `#[test]`: factor all 10 CI-subset SuiteSparse matrices with fast path enabled. Verify backward error < 5e-11 for each.

### Implementation for Fast-Path Factorization

- [ ] T020 [US1] Implement `factor_small_leaf_subtree()` in `src/aptp/numeric.rs`. This function takes a `SmallLeafSubtree`, the supernodes, matrix data, permutation, options, scaling, and assembly maps. It allocates a `SmallLeafWorkspace` (frontal_data sized to `subtree.max_front_size`, g2l of length n), then processes each node in postorder: zero buffer, scatter original entries (reuse `amap_entries`), extend-add child contributions (from local `Vec<Option<ContributionBlock>>`), call `aptp_factor_in_place()`, extract `FrontFactors` via `extract_front_factors()`, compute contribution GEMM via `compute_contribution_gemm()`, extract `ContributionBlock` via `extract_contribution()`, reset g2l. Returns `Vec<(usize, FrontFactors, Option<ContributionBlock>, PerSupernodeStats)>` — one entry per node, with the root's contribution for the parent.
- [ ] T021 [US1] Add `#[cfg(feature = "diagnostic")]` timing instrumentation to `factor_small_leaf_subtree()` in `src/aptp/numeric.rs`: collect `PerSupernodeStats` per node with assembly_time, kernel_time, extraction_time sub-phases, matching the general path's diagnostic output.
- [ ] T022 [US1] Integrate fast-path pre-pass into `factor_tree_levelset()` in `src/aptp/numeric.rs`: before the main level-set loop, iterate over `Vec<SmallLeafSubtree>` and call `factor_small_leaf_subtree()` for each. Store all returned FrontFactors in the results vector. Store root contributions in the global contributions vector. Decrement `remaining_children` for parents of subtree roots. When building the initial ready set, skip all supernodes with `in_small_leaf = true`.
- [ ] T023 [US1] Handle parallel path in `factor_tree_levelset()` in `src/aptp/numeric.rs`: when `options.par != Par::Seq`, independent small-leaf subtrees can be processed via `par_iter` (each allocates its own SmallLeafWorkspace). Ensure the pre-pass results are merged correctly before the main parallel level-set loop begins.
- [ ] T024 [US1] Verify all fast-path tests pass: `cargo test fast_path` and `cargo test suitesparse_ci`
- [ ] T025 [US1] Verify full test suite passes: `cargo test` and `cargo test --features diagnostic`

**Checkpoint**: Fast path is functional. CI-subset SuiteSparse matrices pass with backward error < 5e-11. Simplicial matrices use the fast path (verify via diagnostic stats showing small-leaf subtree counts).

---

## Phase 4: User Story 3 — Regression Validation (Priority: P2)

**Goal**: Confirm no performance or correctness regression on non-simplicial matrices across the full SuiteSparse suite.

**Independent Test**: Run the full 65-matrix SuiteSparse benchmark and compare against pre-change baseline.

- [ ] T026 [US3] Collect pre-change baseline: run `cargo run --example baseline_collection --features diagnostic --release -- --ci-only` on the `ssids` branch (before fast-path changes) and save the baseline JSON.
- [ ] T027 [US3] Run full SuiteSparse correctness check: `cargo test -- --ignored --test-threads=1` with fast path enabled. All 65 matrices must produce backward error < 5e-11.
- [ ] T028 [US3] Run baseline comparison: `cargo run --example baseline_collection --features diagnostic --release -- --ci-only --compare <pre-change-baseline.json>`. Verify no matrix shows > 5% factorization time regression. Verify simplicial matrices (dixmaanl, bloweybq, mario001, linverse, spmsrtls) show improvement.
- [ ] T029 [US3] Profile simplicial matrices: `cargo run --example profile_matrix --features diagnostic --release -- dixmaanl` (and bloweybq, mario001). Verify fast path is active (diagnostic stats show small-leaf subtree count > 0) and factor time improved vs pre-change.
- [ ] T030 [US3] If full SuiteSparse collection available, run complete 65-matrix benchmark: `cargo run --example baseline_collection --features diagnostic --release`. Verify median SPRAL ratio ≤ 1.0x. Verify SC-001: simplicial matrices ≤ 1.5x SPRAL.

**Checkpoint**: Full regression validation complete. All 65 matrices correct. No performance regression on non-simplicial matrices. Simplicial matrices measurably improved.

---

## Phase 5: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, cleanup, and development log updates.

- [ ] T031 [P] Add rustdoc to `classify_small_leaf_subtrees()` and `factor_small_leaf_subtree()` in `src/aptp/numeric.rs`: document algorithm, cite SPRAL `SymbolicSubtree.hxx:57-84` and `SmallLeafNumericSubtree.hxx:187-446` (BSD-3), describe threshold semantics and subtree invariants.
- [ ] T032 [P] Add rustdoc to `SmallLeafSubtree` struct in `src/aptp/numeric.rs`: document fields, invariants (≥2 nodes, postorder, max_front < threshold).
- [ ] T033 [P] Add rustdoc to `small_leaf_threshold` field on `FactorOptions` and `SolverOptions` in `src/aptp/solver.rs`: document default (256), disable value (0), and relationship to `INTRA_NODE_THRESHOLD`.
- [ ] T034 Update `docs/ssids-log.md` with Phase 9.1f entry: record what was built, performance results on simplicial matrices, SPRAL ratios before/after.
- [ ] T035 Update `docs/ssids-plan.md`: mark Phase 9.1f as complete, update success criteria checkboxes, record post-9.1f baseline numbers.
- [ ] T036 Run `cargo fmt --check` and `cargo clippy --all-targets --features diagnostic -- -D warnings`. Fix any issues.
- [ ] T037 Run `cargo doc --no-deps` with `-D warnings`. Fix any documentation warnings.
- [ ] T038 Update `sparse/CLAUDE.md` Current Implementation Status section: add Phase 9.1f entry describing small-leaf subtree fast path, threshold, affected files.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational/US2 (Phase 2)**: Depends on Phase 1 completion (T003, T004 specifically) — BLOCKS US1
- **US1 (Phase 3)**: Depends on Phase 2 completion (classification must work before fast path)
- **US3 (Phase 4)**: Depends on Phase 3 completion (need working fast path to validate regression)
- **Polish (Phase 5)**: Depends on Phases 3 and 4 completion

### User Story Dependencies

- **US2 (Classification)**: Foundational — must complete first. No dependency on other stories.
- **US1 (Fast Path)**: Depends on US2 (classification). Core feature — MVP.
- **US3 (No Regression)**: Depends on US1 (needs working fast path to benchmark). Validation only.

### Within Each Phase

- Tests written and FAILING before implementation (TDD per constitution)
- Implementation tasks are sequential within a phase (later tasks depend on earlier ones)
- Tasks marked [P] within a phase can run in parallel

### Parallel Opportunities

- Phase 1: T001 and T002 can run in parallel (different structs in same file — but small changes, sequential is fine)
- Phase 2 tests: T006–T011 can all run in parallel (independent test functions)
- Phase 3 tests: T015–T018 can run in parallel (independent test functions)
- Phase 5: T031–T033 (docs) and T034–T035 (log updates) can run in parallel

---

## Parallel Example: Phase 2 (Classification Tests)

```bash
# Launch all classification tests in parallel (they're independent):
T006: test_classify_all_small
T007: test_classify_mixed_tree
T008: test_classify_single_node_excluded
T009: test_classify_threshold_boundary
T010: test_classify_disabled
T011: test_classify_multiple_subtrees
```

## Parallel Example: Phase 3 (Fast-Path Tests)

```bash
# Launch all fast-path unit tests in parallel:
T015: test_fast_path_small_chain
T016: test_fast_path_matches_general
T017: test_fast_path_delayed_pivots
T018: test_fast_path_contribution_boundary
```

---

## Implementation Strategy

### MVP First (US2 + US1)

1. Complete Phase 1: Setup (T001–T005) — ~30 min
2. Complete Phase 2: Classification (T006–T014) — ~2 hrs
3. Complete Phase 3: Fast Path (T015–T025) — ~4 hrs
4. **STOP and VALIDATE**: Run CI-subset SuiteSparse tests
5. If passing: proceed to Phase 4 for full regression validation

### Incremental Delivery

1. Phase 1 → Types compile, existing tests pass
2. Phase 2 → Classification works on synthetic trees (independently testable)
3. Phase 3 → Fast path works on real matrices (independently testable, MVP complete)
4. Phase 4 → Full regression validation confirms no regressions
5. Phase 5 → Documentation and log updates

---

## Notes

- Constitution requires TDD: write tests first, verify they fail, then implement
- The fast path reuses existing helper functions (`extract_front_factors`, `extract_contribution`, `compute_contribution_gemm`, `extend_add`, `extend_add_mapped`) — it is NOT a reimplementation of the factorization loop
- The key optimization is amortized workspace setup (single small buffer for entire subtree) and sequential processing without parallel dispatch overhead
- Commit after each phase checkpoint
- `cargo test -- --ignored --test-threads=1` for full SuiteSparse validation (memory-intensive)
