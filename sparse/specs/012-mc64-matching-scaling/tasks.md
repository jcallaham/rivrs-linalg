# Tasks: MC64 Matching & Scaling

**Input**: Design documents from `specs/012-mc64-matching-scaling/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, contracts/api.md

**Tests**: Required by Constitution Principle III (TDD). Tests encode SPRAL's scaling property criteria and are written before implementation per red-green-refactor cycle.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

## Phase 1: Setup

**Purpose**: Module skeleton and type definitions

- [ ] T001 Create `src/aptp/matching.rs` with public types (`Mc64Result`, `Mc64Job`) and stub `mc64_matching` function returning `Err(AnalysisFailure)`. Include internal type stubs (`CostGraph`, `MatchingState`). Add module-level doc comment citing Duff & Koster (2001) and Duff & Pralet (2005).
- [ ] T002 Register module in `src/aptp/mod.rs`: add `pub mod matching;` and `pub use matching::{mc64_matching, Mc64Result, Mc64Job};`
- [ ] T003 Verify `cargo build` compiles cleanly with new module

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Cost graph construction — the shared data structure that all subsequent algorithm steps depend on

- [ ] T004 Write unit tests for `build_cost_graph` in `src/aptp/matching.rs` (`#[cfg(test)]` module). Use a known 3x3 symmetric matrix with hand-computed logarithmic costs. Verify: (a) upper-triangle input is expanded to full bipartite CSC, (b) column maxima are correct, (c) costs are non-negative and match `log(col_max_j) - log|a[i,j]|`, (d) diagonal entries are included. Verify tests FAIL.
- [ ] T005 Implement `build_cost_graph` in `src/aptp/matching.rs`. Accept `&SparseColMat<usize, f64>` (upper-triangular CSC), expand to full symmetric CSC with logarithmic edge costs. Compute column maxima in log domain. Return `CostGraph` struct with `col_ptr`, `row_idx`, `cost`, `col_max` fields. Verify T004 tests PASS.

**Checkpoint**: Cost graph construction works. All subsequent phases depend on this.

---

## Phase 3: User Story 1 — Preprocess an Indefinite Matrix for Better Pivoting (Priority: P1) MVP

**Goal**: MC64 produces a valid maximum-product matching and symmetric scaling factors for sparse symmetric indefinite matrices. Scaled matrix has unit diagonal and off-diagonal entries <= 1.

**Independent Test**: Apply MC64 to hand-constructed indefinite matrices. Verify SPRAL scaling properties: `|s_i * a_ij * s_j| <= 1.0`, row max >= 0.75, unit matched diagonal, matching validity (entries exist, no duplicates), cycle structure (singletons + 2-cycles only).

### Tests for User Story 1

> **NOTE: Write these tests FIRST, ensure they FAIL before implementation**

- [ ] T006 [P] [US1] Write unit tests for `greedy_initial_matching` in `src/aptp/matching.rs`. Use a 4x4 matrix where the greedy heuristic can find at least 3 of 4 matches. Verify: (a) dual variables satisfy feasibility (`u[i] + v[j] <= c[i,j]`), (b) matched edges satisfy extremality (`u[i] + v[j] == c[i,j]`), (c) cardinality > 0. Verify tests FAIL.
- [ ] T007 [P] [US1] Write unit tests for `dijkstra_augment` in `src/aptp/matching.rs`. Use a 3x3 bipartite graph with one unmatched column and a known shortest augmenting path. Verify: (a) augmenting path is found, (b) matching size increases by 1, (c) dual variables remain feasible after update. Verify tests FAIL.
- [ ] T008 [P] [US1] Write unit tests for `symmetrize_scaling` in `src/aptp/matching.rs`. Use known dual variable vectors u, v with hand-computed expected scaling `exp(-(u[i]+v[i])/2)`. Verify: (a) scaling factors are positive, (b) values match expected within 1e-12. Verify tests FAIL.
- [ ] T009 [US1] Write integration tests in `tests/mc64_matching.rs` for end-to-end matching + scaling validation. Create hand-constructed test matrices: (a) 3x3 diagonal (identity matching, unit scaling), (b) 4x4 tridiagonal indefinite, (c) 5x5 arrow indefinite, (d) 1x1 and 2x2 trivial cases. For each matrix verify SPRAL scaling properties: `|s_i * a_ij * s_j| <= 1.0` for all entries, `max_j |s_i * a_ij * s_j| >= 0.75` per row, `matched == n`, matching entries exist in matrix, no duplicate matches, cycle structure is singletons + 2-cycles only. Verify tests FAIL.

### Implementation for User Story 1

- [ ] T010 [US1] Implement `greedy_initial_matching` in `src/aptp/matching.rs`. Compute initial dual variables via column-minimum and row-minimum passes over the cost graph. Build greedy matching from zero-reduced-cost edges. Implement secondary matching pass (2-pass rearrangement from Duff & Koster 2001 Section 4). Return `MatchingState` with `row_match`, `col_match`, `u`, `v` arrays. Verify T006 tests PASS.
- [ ] T011 [US1] Implement `dijkstra_augment` in `src/aptp/matching.rs`. Given an unmatched root column, run Dijkstra on the reduced-cost bipartite graph using `BinaryHeap`. Track distances `d[i]`, parent pointers `pr[j]`, shortest augmenting path length `lsap`. Stop when `lsap <= lsp` or queue empty. If path found: augment matching along alternating path and update dual variables (`u'[i] = u[i] + d[i] - lsap`, `v'[j] = c[i,j] - u'[i]`). Return `true` if augmented, `false` if no path. Verify T007 tests PASS.
- [ ] T012 [US1] Implement `symmetrize_scaling` in `src/aptp/matching.rs`. Compute `scaling[i] = exp(-(u[i] + v[i]) / 2)` from the final dual variables. Handle edge cases: if `u[i] + v[i]` would cause overflow/underflow in exp, clamp to reasonable range. Verify T008 tests PASS.
- [ ] T013 [US1] Wire the main `mc64_matching` function in `src/aptp/matching.rs`. Add input validation: return `NotSquare` error for non-square matrices, `InvalidInput` error for zero dimension or non-finite entries. Handle trivial case (n=1 return identity matching with unit scaling). Main loop: `build_cost_graph` → `greedy_initial_matching` → for each unmatched column call `dijkstra_augment` → `symmetrize_scaling` → construct `Mc64Result` with `Perm<usize>` matching, `Vec<f64>` scaling, and matched count. Verify T009 tests PASS.
- [ ] T014 [US1] Add hand-constructed test matrices that exercise edge cases in `tests/mc64_matching.rs`: (a) badly-scaled matrix with entries spanning 10 orders of magnitude, (b) matrix with zero diagonal entries, (c) well-conditioned PD matrix (should produce near-identity matching and near-unit scaling). Verify SPRAL scaling properties hold for all.

**Checkpoint**: MC64 matching and scaling works on hand-constructed matrices. All SPRAL scaling properties validated. Cycle structure correct. This is the MVP.

---

## Phase 4: User Story 2 — Use MC64 Scaling with Independent METIS Ordering (Priority: P2)

**Goal**: MC64 scaling and METIS ordering compose as independent preprocessing steps. Scaling factors remain valid regardless of ordering. 2-cycle metadata is available.

**Independent Test**: Compute MC64 on a test matrix, compute METIS on the same matrix independently, feed METIS ordering into AptpSymbolic::analyze, verify both succeed and scaling factors are valid.

### Tests for User Story 2

- [ ] T015 [US2] Write integration tests in `tests/mc64_matching.rs` for independent MC64 + METIS composition. Use 2-3 hand-constructed indefinite matrices. For each: (a) compute `mc64_matching`, (b) compute `metis_ordering` on the same matrix, (c) feed METIS ordering into `AptpSymbolic::analyze` via `SymmetricOrdering::Custom`, (d) verify symbolic analysis succeeds with valid fill estimates, (e) verify MC64 scaling factors are still valid (structure-independent). Also verify 2-cycle count from matching is available and non-negative. Verify tests FAIL.

### Implementation for User Story 2

- [ ] T016 [US2] Add a method or helper function to count singletons and 2-cycles from the matching permutation in `src/aptp/matching.rs`. Walk the matching: if `σ(i) == i` → singleton; if `σ(i) == j` and `σ(j) == i` and `i != j` → 2-cycle (count pair once); assert no longer cycles. Return `(singleton_count, two_cycle_count)`. Add unit test verifying correctness on known matching.
- [ ] T017 [US2] Implement the integration test scenarios from T015. Verify all pass: MC64 and METIS compose independently, scaling is structure-independent, symbolic analysis succeeds with METIS ordering.

**Checkpoint**: MC64 and METIS work as independent preprocessing steps. Scaling is structure-independent. 2-cycle metadata available.

---

## Phase 5: User Story 3 — Validate Matching on SuiteSparse Test Suite (Priority: P3)

**Goal**: MC64 produces correct matchings and scaling across the full range of SuiteSparse test matrices, including structurally singular cases.

**Independent Test**: Run MC64 on all SuiteSparse matrices. Verify matching cardinality, scaling properties, and cycle structure. For singular matrices, verify Duff-Pralet correction produces valid scaling.

### Tests for User Story 3

- [ ] T018 [P] [US3] Write unit tests for `duff_pralet_correction` in `src/aptp/matching.rs`. Construct a 4x4 structurally singular matrix (one row/column has no matching). Verify: (a) unmatched scaling is computed as `1.0 / max_k |a[i,k] * s_k|` over matched k, (b) isolated rows get scaling 1.0, (c) all scaling factors are positive and finite. Verify tests FAIL.
- [ ] T019 [P] [US3] Write integration test for structurally singular matrices in `tests/mc64_matching.rs`. Construct a hand-built singular matrix. Verify: (a) `matched < n`, (b) unmatched indices are at end of matching permutation, (c) SPRAL scaling properties hold for matched block, (d) Duff-Pralet correction applied to unmatched indices. Verify tests FAIL.
- [ ] T020 [US3] Write SuiteSparse CI subset integration tests in `tests/mc64_matching.rs` (regular tests, not `#[ignore]`). Load each of the 10 CI matrices. For each: verify (a) `matched == n` for nonsingular matrices, (b) all scaling factors positive and finite, (c) `|s_i * a_ij * s_j| <= 1.0` for all entries, (d) row max >= 0.75, (e) cycle structure is singletons + 2-cycles only. Verify tests FAIL.
- [ ] T021 [US3] Write full SuiteSparse validation tests in `tests/mc64_matching.rs` (`#[ignore]`). Load all 67 matrices. Same validation as T020 plus: (a) compute diagonal dominance metric (min ratio `|d_ii| / sum_{j!=i} |a_ij|`) before and after scaling for each indefinite matrix; assert improvement on at least 80% of indefinite matrices (SC-003), (b) log 2-cycle count per matrix (useful metadata for future condensation assessment). Verify tests FAIL.

### Implementation for User Story 3

- [ ] T022 [US3] Implement structural singularity handling in `src/aptp/matching.rs`. When `matched < n` after the main matching loop in `mc64_matching`: (a) identify matched index set I, (b) extract structurally nonsingular submatrix A[I,I] as a new `SparseColMat`, (c) re-run `mc64_matching` on the submatrix to get proper weighted matching + scaling, (d) map submatrix scaling back to original indices, (e) reorder matching permutation to place unmatched indices at end. Then call `duff_pralet_correction` for the final step: for each unmatched index i compute `scaling[i] = 1.0 / max_k |a[i,k] * scaling[k]|` over matched k (convention: `1/0 = 1.0`). The `duff_pralet_correction` helper handles only the unmatched scaling correction (step e of the API contract); steps (a)-(d) live in the `mc64_matching` orchestration. Verify T018 and T019 tests PASS.
- [ ] T023 [US3] Verify SuiteSparse CI subset tests pass (T020). If any matrices fail, debug and fix the matching algorithm. All 10 CI matrices must produce valid matching and scaling.
- [ ] T024 [US3] Run full SuiteSparse validation (T021) with `cargo test -- --ignored --test-threads=1`. All 67 matrices must pass. Log cycle structure and diagonal dominance statistics. Fix any failures.

**Checkpoint**: MC64 validated across 67 real-world SuiteSparse matrices. Structural singularity handled. All scaling properties verified.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, code quality, and project bookkeeping

- [ ] T025 [P] Add comprehensive rustdoc to `src/aptp/matching.rs`. Module-level doc with algorithm description, academic references (Duff & Koster 2001, Duff & Pralet 2005), SPRAL equivalents, complexity analysis. Function-level docs for `mc64_matching` with `# Arguments`, `# Returns`, `# Errors`, `# Algorithm`, `# References` sections. Document `Mc64Result` fields and `Mc64Job` variants.
- [ ] T026 [P] Run `cargo fmt --check` and `cargo clippy` on the full crate. Fix any warnings or formatting issues in `src/aptp/matching.rs` and `tests/mc64_matching.rs`.
- [ ] T027 Update `docs/ssids-log.md` with Phase 4.2 completion entry. Include: what was built, key decisions (independent MC64+METIS, Duff-Pralet singularity handling), test results summary (CI subset pass count, cycle structure statistics), and notes for future condensation work.
- [ ] T028 Run full test suite (`cargo test` + `cargo test -- --ignored --test-threads=1`) and verify all tests pass. Final compilation check with no warnings.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Phase 1 — BLOCKS all user stories
- **US1 (Phase 3)**: Depends on Phase 2. Core matching algorithm. **MVP target.**
- **US2 (Phase 4)**: Depends on Phase 3 (needs working `mc64_matching`). Tests METIS composition.
- **US3 (Phase 5)**: Depends on Phase 3 (needs working `mc64_matching`). Adds singularity handling and broad validation. Can start T018-T019 (singularity tests) in parallel with US2.
- **Polish (Phase 6)**: Depends on all user stories being complete

### User Story Dependencies

- **US1 (P1)**: Can start after Phase 2 — no dependencies on other stories. MVP.
- **US2 (P2)**: Depends on US1 completion (needs `mc64_matching` to work). Adds independent METIS composition.
- **US3 (P3)**: Depends on US1 completion. Adds structural singularity + broad validation. T018-T021 (test writing) can start in parallel with US2 since they test different concerns.

### Within Each User Story

- Tests MUST be written and FAIL before implementation (Constitution Principle III)
- Internal helpers before main function
- Unit tests before integration tests
- Core algorithm before edge cases

### Parallel Opportunities

- **Phase 3 tests**: T006, T007, T008 can all run in parallel (different internal functions)
- **Phase 5 tests**: T018, T019 can run in parallel (different test scenarios)
- **Phase 6**: T025, T026 can run in parallel (docs vs lint, different concerns)
- **Cross-phase**: T018-T019 (US3 singularity tests) can be written in parallel with T015-T017 (US2 implementation)

---

## Parallel Example: User Story 1

```bash
# Launch all US1 unit tests in parallel (T006, T007, T008 — different functions):
T006: Unit test for greedy_initial_matching in src/aptp/matching.rs
T007: Unit test for dijkstra_augment in src/aptp/matching.rs
T008: Unit test for symmetrize_scaling in src/aptp/matching.rs

# Then sequentially implement each function (T010 → T011 → T012 → T013)
# since each builds on the previous

# T014 (edge case tests) can start after T013 completes
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T003)
2. Complete Phase 2: Foundational — cost graph (T004-T005)
3. Complete Phase 3: User Story 1 — core matching + scaling (T006-T014)
4. **STOP and VALIDATE**: Run `cargo test matching` — all hand-constructed matrices pass SPRAL scaling properties
5. This is a working MC64 implementation suitable for non-singular matrices

### Incremental Delivery

1. Setup + Foundational → Cost graph ready
2. User Story 1 → Core matching + scaling validated on hand-constructed matrices (MVP)
3. User Story 2 → Independent MC64 + METIS composition verified
4. User Story 3 → Singularity handling + SuiteSparse validation (production-ready)
5. Polish → Documentation, cleanup, ssids-log update

---

## Notes

- [P] tasks = different files or different functions, no dependencies
- [Story] label maps task to specific user story for traceability
- Constitution Principle III (TDD) requires tests written and failing before implementation
- SPRAL scaling properties (entries <= 1, row max >= 0.75, unit matched diagonal) are the primary correctness oracle
- No new crate dependencies — pure Rust implementation using std::collections::BinaryHeap
- Commit after each task or logical group (Constitution Principle VI)
- Reference files: Duff & Koster (2001) at `/workspace/rivrs-linalg/references/ssids/duff2001.md`, Duff & Pralet (2005) at `/workspace/rivrs-linalg/references/ssids/duff2005.md`, SPRAL at `/opt/references/spral/src/scaling.f90`
