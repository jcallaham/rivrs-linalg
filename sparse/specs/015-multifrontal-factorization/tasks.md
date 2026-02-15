# Tasks: Multifrontal Numeric Factorization

**Input**: Design documents from `/specs/015-multifrontal-factorization/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, contracts/

**Tests**: TDD approach — tests are included per CLAUDE.md and constitution (Principle III: TDD). Tests written before implementation, must fail first.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup

**Purpose**: Module scaffolding and type definitions

- [ ] T001 Add `pub mod numeric;` to `src/aptp/mod.rs` and re-export `AptpNumeric`, `FrontFactors`, `FactorizationStats`
- [ ] T002 Create `src/aptp/numeric.rs` with module-level doc comment, imports, and stub type definitions for `SupernodeInfo`, `FrontalMatrix`, `ContributionBlock`, `FrontFactors`, `FactorizationStats`, `AptpNumeric` per data-model.md. Include accessor methods on `FrontFactors` (`l11()`, `d11()`, `l21()`, `local_perm()`, `num_eliminated()`, `col_indices()`, `row_indices()`) and `AptpNumeric` (`stats()`, `front_factors()`, `n()`). Add `AptpNumeric::factor()` as a stub returning `todo!()`
- [ ] T003 Create empty test file `tests/multifrontal.rs` with imports for `rivrs_sparse::aptp::*`, `faer::sparse::*`, and test utility helpers

**Checkpoint**: `cargo build` succeeds with stub types. `cargo test --test multifrontal` compiles (no tests yet).

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core internal functions that ALL user stories depend on. These are the building blocks used by the factorization loop.

**CRITICAL**: No user story work can begin until this phase is complete.

- [ ] T004 Implement `build_supernode_info(symbolic: &AptpSymbolic) -> Vec<SupernodeInfo>` in `src/aptp/numeric.rs` — supernodal path: extract from `supernode_begin/end/pattern/parent`; simplicial path: derive from `etree()` and CSC structure (`col_ptr()`/`row_idx()`). Postcondition: `info[s].parent.map_or(true, |p| p > s)` for all s
- [ ] T005 Add helper `build_children_map(infos: &[SupernodeInfo]) -> Vec<Vec<usize>>` in `src/aptp/numeric.rs` — iterate all supernodes, collect children for each parent
- [ ] T006 Implement `extend_add(parent_data: &mut Mat<f64>, child: &ContributionBlock, global_to_local: &[usize])` in `src/aptp/numeric.rs` — for each lower-triangle entry in child contribution, map global indices to parent local indices via `global_to_local`, add to parent matrix
- [ ] T007 Implement `extract_front_factors(frontal_data: &Mat<f64>, frontal_row_indices: &[usize], num_fully_summed: usize, result: &AptpFactorResult) -> FrontFactors` in `src/aptp/numeric.rs` — extract L11 (`[0..ne, 0..ne]`), D11 (build truncated MixedDiagonal from first `ne` pivots of `result.d`), L21 (`[k..m, 0..ne]`), local_perm (`result.perm[0..k]`), col_indices (map through `frontal_row_indices`), row_indices (`frontal_row_indices[k..m]`)
- [ ] T008 Implement `extract_contribution(frontal_data: &Mat<f64>, frontal_row_indices: &[usize], num_fully_summed: usize, result: &AptpFactorResult) -> ContributionBlock` in `src/aptp/numeric.rs` — copy trailing `(m - ne) x (m - ne)` submatrix, build row_indices from `result.delayed_cols` (mapped through `result.perm` and `frontal_row_indices`) + `frontal_row_indices[k..m]`, set `num_delayed = k - ne`
- [ ] T009 Implement `reassemble_global_factors(numeric: &AptpNumeric, symbolic: &AptpSymbolic) -> (Mat<f64>, MixedDiagonal)` as a `#[cfg(test)]` utility in `src/aptp/numeric.rs` — allocate dense L (n x n, identity), allocate global MixedDiagonal (n), scatter each supernode's L11/D11/L21 entries to global positions using `col_indices`/`row_indices`, handling local_perm. This enables reconstruction tests.

**Checkpoint**: All internal functions compile. `cargo build` succeeds. No tests yet — those come in the user story phases.

---

## Phase 3: User Story 3 - Assembly of Frontal Matrices (Priority: P1)

**Goal**: Correctly assemble dense frontal matrices from sparse entries and child contribution blocks.

**Independent Test**: Assemble frontal matrices for small matrices with known elimination trees and verify the dense block matches the expected submatrix of the permuted sparse matrix.

**Why first**: Assembly is tested independently of the full factorization loop. Getting assembly right before building the loop ensures correct data flow.

### Tests for User Story 3

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T010 [P] [US3] Test `build_supernode_info` for supernodal symbolic analysis in `tests/multifrontal.rs` — use a hand-constructed matrix, call `AptpSymbolic::analyze`, verify `col_begin/col_end/pattern/parent` for each supernode match expected values
- [ ] T011 [P] [US3] Test `build_supernode_info` for simplicial symbolic analysis in `tests/multifrontal.rs` — use a small matrix where faer produces simplicial, verify each column produces a 1-column SupernodeInfo with correct pattern and parent
- [ ] T012 [P] [US3] Test scatter of original entries for a leaf supernode (no children) in `tests/multifrontal.rs` — hand-constructed 4x4 or 5x5 sparse matrix, assemble one leaf front, verify F11 and F21 entries match expected dense submatrix values, all other entries zero
- [ ] T013 [P] [US3] Test `extend_add` in `tests/multifrontal.rs` — create a known ContributionBlock with specific row_indices and data, call `extend_add` into a parent frontal matrix, verify entries land at correct parent-local positions
- [ ] T014 [US3] Test assembly with children in `tests/multifrontal.rs` — 6x6 sparse matrix with 2 supernodes (child + parent), assemble child, factor child (using `aptp_factor_in_place`), extract contribution, assemble parent with contribution, verify parent frontal matrix matches expected dense submatrix

### Implementation for User Story 3

- [ ] T015 [US3] Implement scatter logic (private function) in `src/aptp/numeric.rs` — for each column j in `[col_begin..col_end)`, map through fill-reducing permutation, iterate nonzeros in that column of the original sparse matrix, map row indices through inverse permutation, scatter into frontal matrix via `global_to_local`. Handle symmetric placement (lower triangle of original → both lower-tri positions in frontal)

**Checkpoint**: Assembly tests pass. `cargo test --test multifrontal` shows all US3 tests green.

---

## Phase 4: User Story 2 - Correctly Handle Delayed Pivots (Priority: P1)

**Goal**: Delayed columns propagate from child to parent supernodes as additional fully-summed columns and are eventually eliminated.

**Independent Test**: Construct matrices where specific pivots fail at a child but succeed at the parent. Verify factorization completes with correct delay count.

**Why before US1**: The factorization loop (US1) needs delayed pivot handling to work correctly. Testing delayed pivots on small controlled examples before building the full loop avoids debugging cascading issues.

### Tests for User Story 2

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T016 [P] [US2] Test `extract_contribution` in `tests/multifrontal.rs` — factor a 4x4 frontal matrix with `num_fully_summed=2`, verify contribution block size, row_indices, and num_delayed match expected values
- [ ] T017 [P] [US2] Test `extract_front_factors` in `tests/multifrontal.rs` — factor a 5x5 frontal matrix with `num_fully_summed=3`, verify L11, D11 (truncated MixedDiagonal), L21, col_indices, row_indices, local_perm
- [ ] T018 [US2] Test delayed pivot propagation in `tests/multifrontal.rs` — construct a 2-supernode matrix where child delays 1 column, assemble parent with delayed column as additional fully-summed, factor parent, verify all columns eventually eliminated and reconstruction error < 1e-12
- [ ] T019 [US2] Test all-delayed front in `tests/multifrontal.rs` — construct a frontal matrix where `aptp_factor_in_place` eliminates 0 columns (all delayed). Verify `FrontFactors` has `num_eliminated == 0` and contribution block contains all original columns

### Implementation for User Story 2

- [ ] T020 [US2] Implement delayed-column collection logic in `src/aptp/numeric.rs` — given children's contribution blocks, collect delayed global indices (`contribution.row_indices[0..num_delayed]`) and merge into parent's fully-summed column set. Compute parent frontal matrix size: `k_parent = sn.ncols() + total_delayed`, `m = k_parent + pattern.len()`

**Checkpoint**: Extraction and delay propagation tests pass. `cargo test --test multifrontal` shows all US2 + US3 tests green.

---

## Phase 5: User Story 1 - Factor a Sparse Symmetric Indefinite Matrix (Priority: P1)

**Goal**: End-to-end multifrontal factorization via `AptpNumeric::factor()`.

**Independent Test**: Factor hand-constructed and CI SuiteSparse matrices, verify reconstruction error < 10^-12.

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T021 [P] [US1] Test single-supernode-matches-dense in `tests/multifrontal.rs` — small dense matrix (n < 10) that produces a single supernode, verify multifrontal result matches Phase 5's `aptp_factor` exactly (reconstruction error identical)
- [ ] T022 [P] [US1] Test simplicial path in `tests/multifrontal.rs` — use a small sparse matrix (n ~ 5-8) where faer produces simplicial analysis, factor via `AptpNumeric::factor()`, verify reconstruction error < 1e-12
- [ ] T023 [US1] Test hand-constructed matrices in `tests/multifrontal.rs` — factor all 15 hand-constructed test matrices via `AptpNumeric::factor()`, verify reconstruction error < 1e-12 for each. Use `reassemble_global_factors` + existing `reconstruction_error` from `validate.rs` (or equivalent)
- [ ] T024 [US1] Test factorization statistics in `tests/multifrontal.rs` — for a small matrix with known structure, verify `FactorizationStats` fields: `total_1x1_pivots`, `total_2x2_pivots`, `total_delayed`, `max_front_size`
- [ ] T025 [US1] Test dimension mismatch error in `tests/multifrontal.rs` — call `AptpNumeric::factor()` with mismatched symbolic/matrix dimensions, verify `SparseError::DimensionMismatch` is returned
- [ ] T026 [US1] Test inertia validation in `tests/multifrontal.rs` — for hand-constructed matrices with known inertia, compute inertia from combined `MixedDiagonal` factors across all supernodes, verify match

### Implementation for User Story 1

- [ ] T027 [US1] Implement `AptpNumeric::factor()` in `src/aptp/numeric.rs` — the full multifrontal factorization loop per plan.md algorithm overview: (1) `build_supernode_info`, (2) `build_children_map`, (3) allocate `global_to_local` array (size n, sentinel `usize::MAX`), (4) allocate `contributions: Vec<Option<ContributionBlock>>`, (5) postorder loop `0..n_supernodes`: collect delayed from children → compute frontal size → set up `global_to_local` → allocate + scatter + extend-add → `aptp_factor_in_place(frontal, k, options)` → `extract_front_factors` → `extract_contribution` (if non-root and ne < m) → cleanup `global_to_local` → accumulate stats. (6) Return `AptpNumeric { front_factors, stats, n }`. Validate dimension match at entry. Handle root delayed columns (return `SparseError::NumericalSingularity` if root has unresolved delays)
- [ ] T028 [US1] Add CI SuiteSparse integration tests in `tests/suitesparse_ci.rs` — add test function that factors all 10 CI-subset matrices via `AptpNumeric::factor()`, verifying reconstruction error < 1e-12 for each
- [ ] T029 [US1] Add full SuiteSparse integration test (ignored) in `tests/suitesparse_ci.rs` or `tests/multifrontal.rs` — `#[ignore]` test that factors all 67 SuiteSparse matrices via `AptpNumeric::factor()`, verifying reconstruction error < 1e-12

**Checkpoint**: `cargo test --test multifrontal` and `cargo test --test suitesparse_ci` pass. Full SuiteSparse test passes when run with `--ignored --test-threads=1`.

---

## Phase 6: User Story 4 - Schur Complement and Contribution Block Computation (Priority: P2)

**Goal**: Verify that the Schur complement and contribution block are computed correctly. (Note: the actual computation is handled implicitly by the Phase 5 kernel per research Decision 1. This user story focuses on validating that the implicit approach produces correct results.)

**Independent Test**: Factor a multi-level tree, verify contribution blocks flow correctly from leaves to root.

### Tests for User Story 4

- [ ] T030 [P] [US4] Test dense equivalence in `tests/multifrontal.rs` — for matrices with n < 100, factor both via `AptpNumeric::factor()` (sparse multifrontal) and `aptp_factor()` (dense Phase 5), verify both produce reconstruction errors within 10x of each other
- [ ] T031 [US4] Test multi-level contribution flow in `tests/multifrontal.rs` — construct a 3-level tree (grandchild → child → root), factor via `AptpNumeric::factor()`, verify reconstruction error < 1e-12 and that contribution blocks correctly accumulated at each level

**Checkpoint**: All user story tests pass. `cargo test --test multifrontal` shows all tests green.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, cleanup, and final validation

- [ ] T032 [P] Add rustdoc comments to all `pub` types and methods in `src/aptp/numeric.rs` — `AptpNumeric`, `FrontFactors`, `FactorizationStats`, `AptpNumeric::factor()`, all accessor methods. Include algorithm references (Duff & Reid 1983, Liu 1992, Hogg et al. 2020)
- [ ] T033 [P] Add rustdoc to `pub(crate)` types in `src/aptp/numeric.rs` — `SupernodeInfo`, `FrontalMatrix`, `ContributionBlock`, `build_supernode_info()`, `extend_add()`, `extract_front_factors()`, `extract_contribution()`
- [ ] T034 Run `cargo clippy` and `cargo fmt --check` on `src/aptp/numeric.rs` and `tests/multifrontal.rs`, fix any warnings
- [ ] T035 Run full SuiteSparse test `cargo test -- --ignored --test-threads=1` and verify all 67 matrices pass with reconstruction error < 1e-12
- [ ] T036 Update `docs/ssids-log.md` with Phase 6 completion entry — what was built, key decisions (pass-entire-front strategy, unified supernode abstraction), statistics, and any lessons learned
- [ ] T037 Update `docs/ssids-plan.md` Phase 6 section to reflect completion status

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 — BLOCKS all user stories
- **US3: Assembly (Phase 3)**: Depends on Phase 2. Tested independently of the factorization loop
- **US2: Delayed Pivots (Phase 4)**: Depends on Phase 2. Can run in parallel with US3 (different test concerns)
- **US1: Full Factorization (Phase 5)**: Depends on Phases 2, 3, and 4 (uses assembly + delay handling)
- **US4: Schur Complement Validation (Phase 6)**: Depends on Phase 5 (needs working `AptpNumeric::factor()`)
- **Polish (Phase 7)**: Depends on all user stories being complete

### User Story Dependencies

- **US3 (Assembly)** and **US2 (Delayed Pivots)**: Can proceed in parallel after Foundational phase
- **US1 (Full Factorization)**: Depends on US3 + US2 completion (the loop calls assembly + delay logic)
- **US4 (Schur Complement)**: Depends on US1 (validates end-to-end results)

### Within Each User Story

- Tests MUST be written and FAIL before implementation
- Internal functions before integration
- Unit tests before integration tests
- Commit after each task or logical group

### Parallel Opportunities

Within Phase 2: T004 and T005 can run in parallel (different functions). T006, T007, T008 can run in parallel (different functions, no dependencies). T009 depends on T007 (uses FrontFactors).

Within Phase 3 tests: T010, T011, T012, T013 can all run in parallel (independent test functions). T014 depends on scatter + extend-add working.

Within Phase 5 tests: T021, T022 can run in parallel. T025 is independent.

Within Phase 7: T032, T033 can run in parallel. T036, T037 can run in parallel.

---

## Parallel Example: Phase 2 (Foundational)

```
# Group 1 — can run in parallel:
T004: build_supernode_info
T005: build_children_map

# Group 2 — can run in parallel (after Group 1):
T006: extend_add
T007: extract_front_factors
T008: extract_contribution

# Group 3 — depends on T007:
T009: reassemble_global_factors
```

## Parallel Example: Phase 3 Tests (US3)

```
# All independent — can run in parallel:
T010: test build_supernode_info (supernodal)
T011: test build_supernode_info (simplicial)
T012: test scatter for leaf supernode
T013: test extend_add
```

---

## Implementation Strategy

### MVP First (US3 + US2 + US1)

1. Complete Phase 1: Setup (stubs compile)
2. Complete Phase 2: Foundational (internal functions)
3. Complete Phase 3: US3 Assembly (tests + scatter implementation)
4. Complete Phase 4: US2 Delayed Pivots (tests + delay logic)
5. Complete Phase 5: US1 Full Factorization (tests + `AptpNumeric::factor()`)
6. **STOP and VALIDATE**: `cargo test` passes, `cargo test -- --ignored` passes
7. Phase 6: US4 Schur Complement Validation
8. Phase 7: Polish

### Incremental Delivery

1. Phase 1 + 2 → All internal functions compile
2. Add US3 → Assembly tested independently → Commit
3. Add US2 → Delayed pivots tested independently → Commit
4. Add US1 → End-to-end factorization working → Commit (MVP!)
5. Add US4 → Schur complement validated → Commit
6. Polish → Documentation + final validation → Commit

---

## Notes

- [P] tasks = different files or functions, no dependencies
- [Story] label maps task to specific user story for traceability
- All source goes in single file `src/aptp/numeric.rs` (~1000-1500 lines)
- All dedicated tests go in `tests/multifrontal.rs`
- CI SuiteSparse integration tests go in existing `tests/suitesparse_ci.rs`
- `reassemble_global_factors` is `#[cfg(test)]` only — O(n^2) memory, not for solve path
- Constitution Principle III (TDD) requires tests before implementation
- Reconstruction error < 10^-12 is the primary correctness oracle (Constitution v1.1.0)
