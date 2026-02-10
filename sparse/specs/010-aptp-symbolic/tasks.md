# Tasks: APTP Symbolic Analysis

**Input**: Design documents from `/specs/010-aptp-symbolic/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, contracts/

**Tests**: Required (constitution Principle III: TDD is NON-NEGOTIABLE). Tests written first, verified failing, then implementation.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3, US4)
- Include exact file paths in descriptions

## Path Conventions

- **Single project**: `src/`, `tests/` at repository root (`/workspace/rivrs-linalg/sparse/`)

---

## Phase 1: Setup

**Purpose**: Module scaffolding and export wiring

- [ ] T001 Create `src/aptp/symbolic.rs` with module-level doc comment citing algorithm references (Liu 1990, Gilbert et al. 1992/1994, Hogg et al. 2016), empty `AptpSymbolic` struct stub, and empty `SymbolicStatistics` struct stub
- [ ] T002 Add `pub mod symbolic;` to `src/aptp/mod.rs` and add re-exports for `AptpSymbolic` and `SymbolicStatistics`
- [ ] T003 Verify `cargo check` passes with the new empty stubs

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Core types and error handling that ALL user stories depend on

**CRITICAL**: No user story work can begin until this phase is complete

- [ ] T004 Define `SymbolicStatistics` struct in `src/aptp/symbolic.rs` with fields: `dimension: usize`, `predicted_nnz: usize`, `average_col_count: f64`, `is_supernodal: bool`, `n_supernodes: Option<usize>`, `total_pivot_buffer: usize`. Derive `Debug, Clone`. Implement `Display` for human-readable diagnostics output.
- [ ] T005 Define `AptpSymbolic` struct in `src/aptp/symbolic.rs` with fields: `inner: SymbolicCholesky<usize>`, `etree: Vec<isize>`, `col_counts: Vec<usize>`, `pivot_buffer: Vec<usize>`. Add necessary faer imports (`factorize_symbolic_cholesky`, `prefactorize_symbolic_cholesky`, `SymbolicCholesky`, `SymbolicCholeskyRaw`, `SymmetricOrdering`, `CholeskySymbolicParams`, `Side`, `SymbolicSparseColMatRef`, `PermRef`, `EliminationTreeRef`). Do NOT implement `analyze()` yet — just the struct definition.
- [ ] T006 Review `src/error.rs` and confirm `SparseError` has adequate variants for symbolic analysis failures (`NotSquare`, `DimensionMismatch`, `AnalysisFailure`). If `AnalysisFailure` needs a `context: String` field or faer error wrapping, add it now. Ensure `FaerError` can be mapped to `SparseError::AnalysisFailure` with a descriptive message.
- [ ] T007 Verify `cargo check` passes with the full struct definitions and imports

**Checkpoint**: Foundation ready — AptpSymbolic struct defined, errors adequate, user story implementation can begin

---

## Phase 3: User Story 1 — Symbolic Analysis with Default Ordering (Priority: P1) MVP

**Goal**: A solver developer provides a sparse symmetric matrix and receives a symbolic analysis result with elimination tree, fill-in prediction, and permutation using AMD ordering.

**Independent Test**: Load any test matrix, run `AptpSymbolic::analyze` with `SymmetricOrdering::Amd`, verify valid statistics.

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation (TDD — Constitution Principle III)**

- [ ] T008 [P] [US1] Write unit test `test_analyze_1x1_matrix` in `src/aptp/symbolic.rs` (inline `#[cfg(test)]` module): construct a 1×1 sparse matrix, call `analyze` with AMD ordering, assert result has dimension 1, predicted_nnz >= 1, valid statistics
- [ ] T009 [P] [US1] Write unit test `test_analyze_diagonal_matrix` in `src/aptp/symbolic.rs`: construct a 5×5 diagonal sparse matrix, call `analyze` with AMD ordering, assert predicted_nnz == 5 (no fill-in for diagonal), etree has all roots (-1)
- [ ] T010 [P] [US1] Write unit test `test_analyze_deterministic` in `src/aptp/symbolic.rs`: analyze the same matrix twice with AMD ordering, assert `predicted_nnz`, `statistics`, `etree`, and `col_counts` are identical across both calls
- [ ] T011 [P] [US1] Write unit test `test_analyze_predicted_nnz_lower_bound` in `src/aptp/symbolic.rs`: construct a matrix with known lower-triangle nnz, verify predicted_nnz >= structural nnz of lower triangle
- [ ] T012 [P] [US1] Write unit test `test_analyze_statistics_consistency` in `src/aptp/symbolic.rs`: analyze a small matrix, verify `statistics()` fields are consistent (dimension matches nrows, average_col_count = predicted_nnz / dimension, total_pivot_buffer = sum of pivot_buffer_estimates)
- [ ] T013 [US1] Write integration test `test_analyze_all_hand_constructed` in `tests/symbolic_analysis.rs`: load all 15 hand-constructed test matrices via `load_test_matrix`, call `AptpSymbolic::analyze` with AMD ordering on each, assert all succeed with valid statistics (dimension matches, predicted_nnz > 0, average_col_count >= 1.0)
- [ ] T014 [US1] Write integration test `test_analyze_suitesparse_ci_subset` in `tests/symbolic_analysis.rs`: load the 10 SuiteSparse CI subset matrices, call `AptpSymbolic::analyze` with AMD ordering on each, assert all succeed with valid statistics
- [ ] T015 [US1] Verify all US1 tests FAIL (compile errors from unimplemented `analyze` method). Document the failure in a commit message.

### Implementation for User Story 1

- [ ] T016 [US1] Implement the private helper for permuted symbolic structure in `src/aptp/symbolic.rs`: given a `SymbolicSparseColMatRef` and a permutation (fwd/inv arrays from `PermRef`), construct the permuted CSC symbolic structure `P^T A P`. Algorithm: for each column `j` in the permuted matrix, gather row indices from column `fwd[j]` of the original matrix, remap each row index `i` to `inv[i]`, collect into new column pointer / row index arrays, and sort row indices within each column. Return a `SparseColMat<usize, ()>` (symbolic only — no values needed) whose `.symbolic()` can be passed to `prefactorize_symbolic_cholesky`. Only the upper triangle needs to be kept (matching `Side::Upper`). Use faer's `MemStack` for workspace allocation. Reference: Gilbert et al. (1992) Section 3.3 for permuted structure computation. See also `faer-0.22.6/src/sparse/linalg/cholesky.rs:4192` where `factorize_symbolic_cholesky` performs this internally.
- [ ] T017 [US1] Implement `AptpSymbolic::analyze()` in `src/aptp/symbolic.rs`: (1) validate input — return `SparseError::NotSquare` if not square; (2) call `factorize_symbolic_cholesky` with `Side::Upper`, the given ordering, and `CholeskySymbolicParams::default()`, mapping `FaerError` to `SparseError::AnalysisFailure`; (3) extract permutation from result and compute permuted symbolic structure; (4) call `prefactorize_symbolic_cholesky` on the permuted structure to get etree and col_counts; (5) compute pivot_buffer using the 10% heuristic (per-supernode for supernodal, per-column for simplicial); (6) return `AptpSymbolic { inner, etree, col_counts, pivot_buffer }`.
- [ ] T018 [US1] Implement faer-delegating accessor methods on `AptpSymbolic` in `src/aptp/symbolic.rs`: `perm()` → `self.inner.perm()`, `predicted_nnz()` → `self.inner.len_val()`, `nrows()` → `self.inner.nrows()`, `ncols()` → `self.inner.ncols()`, `raw()` → `self.inner.raw()`. Add full rustdoc with `# Examples` for `analyze()` and `perm()`.
- [ ] T019 [US1] Implement APTP-specific accessor methods on `AptpSymbolic` in `src/aptp/symbolic.rs`: `etree()` → `&self.etree`, `col_counts()` → `&self.col_counts`, `pivot_buffer_estimates()` → `&self.pivot_buffer`, `total_pivot_buffer()` → sum of pivot_buffer. Add `is_supernodal()` method that pattern-matches on `self.inner.raw()`.
- [ ] T020 [US1] Implement `statistics()` method on `AptpSymbolic` in `src/aptp/symbolic.rs`: construct `SymbolicStatistics` from struct fields. Handle 0-dimension edge case for average_col_count (avoid division by zero).
- [ ] T021 [US1] Run `cargo test` — all US1 tests should now pass. Run `cargo clippy` and `cargo fmt --check`. Fix any issues. Run `cargo doc` and verify all public items are documented.

**Checkpoint**: User Story 1 complete — `AptpSymbolic::analyze` works on all test matrices with AMD ordering, returns valid statistics, is deterministic.

---

## Phase 4: User Story 2 — Custom Ordering Support (Priority: P2)

**Goal**: Accept a pre-computed fill-reducing permutation (e.g., identity, or future MC64/METIS output) and use it instead of AMD.

**Independent Test**: Construct a permutation, pass as `SymmetricOrdering::Custom`, verify analysis completes and permutation is reflected in result.

### Tests for User Story 2

- [ ] T022 [P] [US2] Write unit test `test_analyze_custom_identity_ordering` in `src/aptp/symbolic.rs`: construct an identity permutation via `perm_from_forward`, pass as `SymmetricOrdering::Custom`, verify analysis succeeds and `perm()` returns the identity permutation
- [ ] T023 [P] [US2] Write unit test `test_analyze_custom_reverse_ordering` in `src/aptp/symbolic.rs`: construct a reverse permutation, pass as `SymmetricOrdering::Custom`, verify analysis succeeds, predicted_nnz may differ from AMD result
- [ ] T024 [P] [US2] Write unit test `test_analyze_invalid_perm_dimension` in `src/aptp/symbolic.rs`: construct a permutation with wrong dimension (e.g., n-1 for an n×n matrix), pass as custom ordering, verify `SparseError::DimensionMismatch` is returned (not a panic)
- [ ] T025 [US2] Write integration test `test_custom_ordering_on_hand_constructed` in `tests/symbolic_analysis.rs`: for a representative hand-constructed matrix, run analysis with AMD and with identity ordering, verify both succeed and produce valid (potentially different) statistics

### Implementation for User Story 2

- [ ] T026 [US2] Add input validation for custom ordering in `AptpSymbolic::analyze()` in `src/aptp/symbolic.rs`: when `SymmetricOrdering::Custom(perm)` is provided, check that `perm.len() == matrix.nrows()`, returning `SparseError::DimensionMismatch` with descriptive message if not. (Note: the core `analyze` from US1 already passes ordering through to faer — this task adds the pre-validation.)
- [ ] T027 [US2] Run `cargo test` — all US1 + US2 tests should pass. Run `cargo clippy` and `cargo fmt --check`.

**Checkpoint**: User Story 2 complete — custom orderings accepted, validated, and propagated correctly. Identity and reverse orderings work. Invalid permutations produce descriptive errors.

---

## Phase 5: User Story 3 — Supernodal Structure Access (Priority: P3)

**Goal**: Expose supernodal decomposition (supernode column ranges, assembly tree, row structure per supernode) for use by Phase 6 multifrontal factorization.

**Independent Test**: Run symbolic analysis on a matrix with multiple supernodes, verify supernodal accessors return structurally valid results.

### Tests for User Story 3

- [ ] T028 [P] [US3] Write unit test `test_n_supernodes_simplicial` in `src/aptp/symbolic.rs`: analyze a tiny matrix (e.g., 3×3) that faer classifies as simplicial, verify `n_supernodes()` returns `None` and `is_supernodal()` returns `false`
- [ ] T029 [P] [US3] Write unit test `test_supernode_column_ranges_partition` in `src/aptp/symbolic.rs`: analyze a matrix classified as supernodal, verify that supernode_begin/supernode_end cover all columns without gaps or overlaps (i.e., begin[0] == 0, end[last] == dimension, begin[i+1] == end[i])
- [ ] T030 [P] [US3] Write unit test `test_assembly_tree_valid` in `src/aptp/symbolic.rs`: analyze a SuiteSparse matrix classified as supernodal, call `supernode_parent(s)` for each supernode, verify: (a) each non-root supernode has a parent with index > s (postorder property), (b) exactly one root exists with parent == None, (c) the tree covers all supernodes. This tests US3 acceptance scenario 2.
- [ ] T031 [P] [US3] Write unit test `test_supernode_row_pattern_nonempty` in `src/aptp/symbolic.rs`: analyze a SuiteSparse matrix classified as supernodal, verify `supernode_pattern(s)` returns `Some(pattern)` with `pattern.len() >= 0` for each supernode, and at least one supernode has a non-empty pattern (off-diagonal rows)
- [ ] T032 [US3] Write integration test `test_supernodal_structure_suitesparse` in `tests/symbolic_analysis.rs`: for SuiteSparse CI matrices (which are likely supernodal), verify `n_supernodes().is_some()`, supernode ranges are valid, each supernode's column count is >= 1, assembly tree is valid, and row patterns are accessible

### Implementation for User Story 3

- [ ] T033 [US3] Implement `n_supernodes()` method on `AptpSymbolic` in `src/aptp/symbolic.rs`: pattern-match on `self.inner.raw()` — return `Some(sn.n_supernodes())` for supernodal, `None` for simplicial
- [ ] T034 [US3] Implement `supernode_begin()` and `supernode_end()` methods on `AptpSymbolic` in `src/aptp/symbolic.rs`: pattern-match on `self.inner.raw()` — return `Some(sn.supernode_begin())` / `Some(sn.supernode_end())` for supernodal, `None` for simplicial. Add rustdoc explaining the supernode column range semantics.
- [ ] T035 [US3] Implement `supernode_pattern(s: usize)` method on `AptpSymbolic` in `src/aptp/symbolic.rs`: delegate to faer's `sn.supernode(s).pattern()` for supernodal, return `None` for simplicial. This exposes the row structure per supernode (FR-007). Add rustdoc explaining that the pattern contains row indices of off-diagonal entries for the supernode's columns.
- [ ] T036 [US3] Implement `supernode_parent(s: usize)` method on `AptpSymbolic` in `src/aptp/symbolic.rs`: derive supernodal parent pointers from the column-level etree. Algorithm: for supernode `s` with column range `[begin, end)`, look up `etree[end - 1]` (the etree parent of the last column in the supernode). If it is -1 (root), return `None`. Otherwise, find which supernode contains that parent column using binary search on `supernode_begin`. Cache the parent array on first call (or compute eagerly in `analyze()`). Cite Liu (1992) for assembly tree derivation from elimination tree. This satisfies FR-007 (assembly tree) and US3 acceptance scenario 2.
- [ ] T037 [US3] Run `cargo test` — all US1 + US2 + US3 tests should pass. Run `cargo clippy` and `cargo fmt --check`.

**Checkpoint**: User Story 3 complete — supernodal structure accessible via public API. Column range partition verified. Assembly tree parent pointers derived and validated. Row structure per supernode exposed. Simplicial fallback returns None.

---

## Phase 6: User Story 4 — APTP Pivot Buffer Estimation (Priority: P3)

**Goal**: Provide per-unit (supernode or column) workspace estimates for delayed pivots, enabling the numeric factorization phase to pre-allocate memory.

**Independent Test**: Run symbolic analysis, verify buffer estimates are non-negative, proportional, and reproducible.

### Tests for User Story 4

- [ ] T038 [P] [US4] Write unit test `test_pivot_buffer_nonnegative` in `src/aptp/symbolic.rs`: analyze several test matrices, verify all entries in `pivot_buffer_estimates()` are non-negative (>= 0, which is trivially true for usize, but verify the vec is non-empty and total > 0 for non-trivial matrices)
- [ ] T039 [P] [US4] Write unit test `test_pivot_buffer_proportional` in `src/aptp/symbolic.rs`: analyze a matrix, verify `total_pivot_buffer()` > 0 and is roughly proportional to `predicted_nnz()` (within a reasonable range, e.g., 1%–50% of predicted_nnz)
- [ ] T040 [P] [US4] Write unit test `test_pivot_buffer_reproducible` in `src/aptp/symbolic.rs`: analyze the same matrix twice, verify `pivot_buffer_estimates()` are identical element-wise
- [ ] T041 [US4] Write unit test `test_pivot_buffer_length_matches_structure` in `src/aptp/symbolic.rs`: for supernodal results, verify `pivot_buffer_estimates().len() == n_supernodes().unwrap()`; for simplicial results, verify `pivot_buffer_estimates().len() == nrows()`

### Implementation for User Story 4

- [ ] T042 [US4] Verify that the 10% heuristic implemented in T017 produces correct buffer estimates — the pivot buffer computation is part of `analyze()` (already implemented in US1). This task validates the implementation: review the buffer computation logic, add a doc comment to the private helper explaining the heuristic formula and citing Hogg et al. (2016) Section 2.4 on delayed pivot propagation.
- [ ] T043 [US4] Run `cargo test` — all US1 + US2 + US3 + US4 tests should pass. Run `cargo clippy` and `cargo fmt --check`.

**Checkpoint**: User Story 4 complete — pivot buffer estimates validated as non-negative, proportional, reproducible, and correctly sized.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Edge cases, benchmarks, documentation, and final validation

- [ ] T044 [P] Write unit test `test_analyze_empty_matrix` in `src/aptp/symbolic.rs`: construct a 0×0 sparse matrix, verify `analyze` either returns a trivial result or a descriptive error (not a panic)
- [ ] T045 [P] Write unit test `test_analyze_block_diagonal` in `src/aptp/symbolic.rs`: construct a block-diagonal matrix (two disconnected blocks), verify analysis succeeds and etree reflects the forest structure (multiple roots with parent == -1)
- [ ] T046 [P] Write unit test `test_analyze_not_square_error` in `src/aptp/symbolic.rs`: construct a non-square symbolic matrix (if faer allows), verify `SparseError::NotSquare` is returned
- [ ] T047 Add symbolic analysis benchmark group in `benches/solver_benchmarks.rs`: benchmark `AptpSymbolic::analyze` on representative SuiteSparse matrices with AMD ordering, recording baseline timing for SC-006 (< 5% of future factor time)
- [ ] T048 Ensure all public items in `src/aptp/symbolic.rs` have complete rustdoc: `AptpSymbolic` struct, `SymbolicStatistics` struct, all public methods. Include `# Examples`, `# Errors`, `# Panics` sections where appropriate. Cite algorithm references in module-level doc comment.
- [ ] T049 Run full validation: `cargo test`, `cargo clippy`, `cargo fmt --check`, `cargo doc --no-deps`. All must pass with zero warnings.
- [ ] T050 Update `docs/ssids-log.md` with Phase 3 completion entry: what was built, key design decisions (two-call strategy, dual-variant handling, 10% buffer heuristic), test results, and next steps.
- [ ] T051 Update `docs/ssids-plan.md` Phase 3 success criteria checkboxes to reflect completion status.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Phase 1 (Setup)**: No dependencies — start immediately
- **Phase 2 (Foundational)**: Depends on Phase 1 — BLOCKS all user stories
- **Phase 3 (US1 — P1 MVP)**: Depends on Phase 2 — core analyze() implementation
- **Phase 4 (US2 — P2)**: Depends on Phase 3 (needs working `analyze()` to test custom ordering)
- **Phase 5 (US3 — P3)**: Depends on Phase 3 (needs working `analyze()` to test supernodal access)
- **Phase 6 (US4 — P3)**: Depends on Phase 3 (buffer estimates are computed in `analyze()`)
- **Phase 7 (Polish)**: Depends on all user stories being complete

### User Story Dependencies

- **US1 (P1)**: Foundational dependency only. All other stories depend on US1's `analyze()`.
- **US2 (P2)**: Depends on US1 (needs working analyze to add validation). Can proceed after US1 checkpoint.
- **US3 (P3)**: Can proceed after US1 checkpoint. Independent of US2 and US4.
- **US4 (P3)**: Can proceed after US1 checkpoint. Independent of US2 and US3.

### Within Each User Story

1. Tests MUST be written and FAIL before implementation (TDD)
2. Implementation tasks in dependency order
3. Checkpoint: all tests pass before moving to next story

### Parallel Opportunities

- **Phase 1**: T001 and T002 are sequential (T002 depends on T001)
- **Phase 2**: T004 and T005 can run in parallel (different struct definitions); T006 is independent
- **Phase 3 tests**: T008–T012 can all run in parallel (different test functions, same file)
- **Phase 4 tests**: T022–T024 can all run in parallel
- **Phase 5 tests**: T028–T031 can all run in parallel
- **Phase 6 tests**: T038–T040 can all run in parallel
- **US3 and US4**: Can run in parallel after US1 is complete (independent accessors)
- **Phase 7**: T044–T047 can all run in parallel (different test/bench files)

---

## Parallel Example: User Story 1 Tests

```
# Launch all US1 unit tests together (they are in different test functions, same file):
T008: test_analyze_1x1_matrix in src/aptp/symbolic.rs
T009: test_analyze_diagonal_matrix in src/aptp/symbolic.rs
T010: test_analyze_deterministic in src/aptp/symbolic.rs
T011: test_analyze_predicted_nnz_lower_bound in src/aptp/symbolic.rs
T012: test_analyze_statistics_consistency in src/aptp/symbolic.rs

# Then integration tests (depend on unit test patterns but different file):
T013: test_analyze_all_hand_constructed in tests/symbolic_analysis.rs
T014: test_analyze_suitesparse_ci_subset in tests/symbolic_analysis.rs
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001–T003)
2. Complete Phase 2: Foundational (T004–T007)
3. Complete Phase 3: User Story 1 (T008–T021)
4. **STOP and VALIDATE**: `AptpSymbolic::analyze` works on all 25 test matrices with AMD ordering
5. This is a functional MVP — symbolic analysis is usable by downstream phases

### Incremental Delivery

1. Setup + Foundational → Struct defined, errors ready
2. US1 → Core analyze() works with AMD (MVP!)
3. US2 → Custom ordering validated (Phase 4 integration point ready)
4. US3 + US4 (parallel) → Supernodal access + buffer estimates (Phase 6 integration point ready)
5. Polish → Edge cases, benchmarks, docs, plan updates

---

## Notes

- [P] tasks = different files or independent functions, no dependencies
- [Story] label maps task to specific user story for traceability
- Constitution Principle III (TDD) requires tests first — this is NON-NEGOTIABLE
- Constitution Principle IV requires all academic references cited in rustdoc
- Constitution Principle II (Clean Room) — only faer (MIT) and academic papers consulted
- Commit after each phase checkpoint (Constitution Principle VI)
- The 10% buffer heuristic is a starting point; refine in Phase 5-6 when numeric factorization validates it empirically
