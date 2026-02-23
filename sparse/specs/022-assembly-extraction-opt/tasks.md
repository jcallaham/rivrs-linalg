# Tasks: Assembly & Extraction Optimization (Phase 9.1c)

**Input**: Design documents from `/specs/022-assembly-extraction-opt/`
**Prerequisites**: plan.md (required), spec.md (required), research.md, data-model.md, quickstart.md

**Tests**: Not explicitly requested. Existing test suite (380+ unit tests + 65-matrix SuiteSparse suite) provides regression coverage. No new test tasks generated.

**Organization**: Tasks are grouped by user story. Implementation order follows plan.md phases (A→B→C→D), which reorders spec priorities for risk management: US2 (bulk extraction, lowest risk) before US3 (bulk zeroing) before US1 (scatter maps, most complex).

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup

**Purpose**: Establish baseline measurements before any optimization changes

- [X] T001 Collect pre-optimization baseline via `cargo run --example baseline_collection --features diagnostic --release` and save output
- [X] T002 Run `profile_matrix` on c-71 and c-big to record pre-optimization assembly/extraction/kernel percentages

**Checkpoint**: Baseline measurements recorded for comparison after optimization

---

## Phase 2: User Story 2 — Bulk Extraction via Column-Slice Copies (Priority: P2, Plan Phase A)

**Goal**: Replace element-by-element extraction loops with per-column `copy_from_slice` using faer's `col_as_slice` API in `extract_contribution` and `extract_front_factors`.

**Independent Test**: `cargo test` (all unit tests pass) + `cargo test -- --ignored --test-threads=1` (65 SuiteSparse matrices, backward error < 5e-11). Backward errors must be identical to pre-optimization baseline.

### Implementation for User Story 2

- [X] T003 [US2] Optimize `extract_contribution` in `src/aptp/numeric.rs` — replace element-by-element lower-triangle column copy with per-column `col_as_slice()`/`col_as_slice_mut()` + `copy_from_slice` for each column of the contribution block
- [X] T004 [US2] Optimize `extract_front_factors` L21 extraction in `src/aptp/numeric.rs` — replace element-by-element copy of L21 columns (rows `ne..m`) with per-column `col_as_slice()` + `copy_from_slice`
- [X] T005 [US2] Optimize `extract_front_factors` L11 extraction in `src/aptp/numeric.rs` — replace element-by-element copy with per-column slice copies where pivot type allows contiguous segments (1x1 pivots); preserve element-by-element for 2x2 pivot diagonal entries
- [X] T006 [US2] Run `cargo test` and `cargo test --features diagnostic` to verify all unit tests pass with bulk extraction
- [X] T007 [US2] Run `cargo test -- --ignored --test-threads=1` to verify all 65 SuiteSparse matrices pass with identical backward errors
- [X] T008 [US2] Run `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` to verify no warnings

**Checkpoint**: Extraction functions use bulk copies. All tests pass. Backward errors identical to baseline.

---

## Phase 3: User Story 3 — Optimized Frontal Matrix Zeroing (Priority: P3, Plan Phase B)

**Goal**: Replace element-by-element zeroing in `zero_frontal` with per-column `fill(0.0)` on contiguous column slices.

**Independent Test**: `cargo test` + SuiteSparse ignored tests. All backward errors unchanged.

### Implementation for User Story 3

- [X] T009 [US3] Optimize `zero_frontal` in `src/aptp/numeric.rs` — replace nested loop with per-column `col_as_slice_mut(j)[j..m].fill(0.0)` for the lower-triangle zeroing
- [X] T010 [US3] Run `cargo test` and `cargo test -- --ignored --test-threads=1` to verify correctness unchanged
- [X] T011 [US3] Run `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` to verify no warnings

**Checkpoint**: Zeroing uses bulk fill operations. All tests pass. Backward errors identical to baseline.

---

## Phase 4: User Story 1 — Precomputed Scatter Maps (Priority: P1, Plan Phase C)

**Goal**: Precompute assembly index mappings during symbolic analysis so that `scatter_original_entries_multi` and `extend_add` avoid per-entry index arithmetic at factorization time.

**Independent Test**: `cargo test` + SuiteSparse ignored tests. All backward errors < 5e-11. `profile_matrix` on c-71 shows assembly percentage decrease.

### Implementation for User Story 1

- [X] T012 [US1] Define `AssemblyMaps` struct in `src/aptp/numeric.rs` (moved from symbolic.rs — amalgamation happens in factor, not analyze) with fields: `amap_entries: Vec<u32>` (4 u32 per entry), `amap_offsets: Vec<usize>`, `ea_map: Vec<u32>`, `ea_offsets: Vec<usize>`, `ea_child_snode: Vec<usize>`, `ea_snode_child_begin: Vec<usize>`
- [X] T013 [US1] Implement amap construction in `build_assembly_maps()` in `src/aptp/numeric.rs` — iterates owned column ranges, maps CSC entries to frontal linear indices with 4-tuple format (src, dest, scale_row, scale_col), handles non-contiguous owned_ranges and upper-triangle dedup
- [X] T014 [US1] Implement extend-add map construction in `build_assembly_maps()` in `src/aptp/numeric.rs` — computes child pattern row → parent local row mappings (zero-delay assumption)
- [X] T015 [US1] Integrate scatter map computation into `AptpNumeric::factor()` in `src/aptp/numeric.rs` — called after amalgamation, passed through `factor_tree_levelset` and `factor_single_supernode`
- [X] T016 [US1] Run `cargo test` to verify factorization still works correctly with added map computation
- [X] T017 [US1] Implement amap-based scatter in `factor_single_supernode` in `src/aptp/numeric.rs` — when `ndelay_in == 0`, uses precomputed entries for direct indexed scatter with scaling support; falls back to `scatter_original_entries_multi` when delays present
- [X] T018 [US1] Implement `extend_add_mapped` in `src/aptp/numeric.rs` — uses precomputed row mapping instead of `global_to_local` lookups; `factor_single_supernode` dispatches to it when child has zero delays
- [X] T019 [US1] Update `factor_single_supernode` and `factor_tree_levelset` call sites in `src/aptp/numeric.rs` — `assembly_maps: &AssemblyMaps` parameter threaded through sequential and parallel paths
- [X] T020 [US1] Run `cargo test` and `cargo test --features diagnostic` to verify all unit tests pass with precomputed maps
- [X] T021 [US1] Run `cargo test -- --ignored --test-threads=1` to verify all 65 SuiteSparse matrices pass with backward error < 5e-11
- [X] T022 [US1] Run `cargo clippy --all-targets` and `cargo clippy --all-targets --features diagnostic` to verify no warnings

**Checkpoint**: Scatter maps precomputed during analysis and used during factorization. All tests pass. Assembly overhead reduced.

---

## Phase 5: Polish & Benchmark Validation (Plan Phase D)

**Purpose**: Verify performance targets and document results

- [X] T023 Run `profile_matrix` on c-71, c-big, ncvxqp7 and document phase breakdown changes vs pre-optimization baseline
- [X] T024 Run full SPRAL comparison benchmark (`spral_benchmark --threads 1 --rivrs`) and verify: c-71 ratio < 3.0x, c-big ratio < 3.0x, median ratio <= 1.01x, all backward errors < 5e-11, no matrix regression > 5%
- [X] T025 Update `docs/ssids-log.md` with Phase 9.1c results (profiling breakdown, benchmark ratios, scatter map memory usage)
- [X] T026 Update `docs/ssids-plan.md` with Phase 9.1c completion status
- [X] T027 Update `CLAUDE.md` Phase 9.1c status from "Next" to "Completed"

**Checkpoint**: All success criteria (SC-001 through SC-007) verified. Documentation updated.

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — collect baseline first
- **US2 — Bulk Extraction (Phase 2)**: Depends on Phase 1 baseline — lowest risk optimization, do first
- **US3 — Bulk Zeroing (Phase 3)**: Depends on Phase 2 completion — incremental low-risk optimization
- **US1 — Scatter Maps (Phase 4)**: Depends on Phase 3 completion — most complex, builds on stable extraction/zeroing
- **Polish (Phase 5)**: Depends on all user stories complete — final validation and documentation

### User Story Dependencies

- **User Story 2 (P2)**: Modified only in `src/aptp/numeric.rs` (extraction functions). No dependencies on other stories.
- **User Story 3 (P3)**: Modified only in `src/aptp/numeric.rs` (zero_frontal). No dependencies on other stories.
- **User Story 1 (P1)**: Modified in both `src/aptp/symbolic.rs` (map construction) and `src/aptp/numeric.rs` (scatter/extend_add). No dependencies on US2/US3 but sequenced last due to complexity.

### Within Each User Story

- Implementation before verification
- All tests must pass before moving to next phase
- Clippy clean before moving to next phase

### Parallel Opportunities

Within Phase 2 (US2):
- T003, T004, T005 modify different functions in the same file — execute sequentially to avoid conflicts

Within Phase 4 (US1):
- T012, T013, T014 are in `symbolic.rs` — execute sequentially
- T017, T018 modify different functions in `numeric.rs` — could potentially parallel but safer sequential due to shared file

---

## Parallel Example: User Story 1

```bash
# Phase 4 tasks are in two files but have logical dependencies:
# symbolic.rs tasks (T012-T015) must complete before numeric.rs tasks (T017-T019)
# because numeric.rs needs the AssemblyMaps struct defined in symbolic.rs

# Sequential within symbolic.rs:
Task: "Define AssemblyMaps struct in src/aptp/symbolic.rs"
Task: "Implement amap construction in src/aptp/symbolic.rs"
Task: "Implement extend-add map construction in src/aptp/symbolic.rs"
Task: "Integrate into AptpSymbolic::analyze() in src/aptp/symbolic.rs"

# Then sequential within numeric.rs:
Task: "Modify scatter_original_entries_multi in src/aptp/numeric.rs"
Task: "Modify extend_add in src/aptp/numeric.rs"
Task: "Update factor_single_supernode call sites in src/aptp/numeric.rs"
```

---

## Implementation Strategy

### MVP First (User Story 2 — Bulk Extraction)

1. Complete Phase 1: Baseline collection
2. Complete Phase 2: User Story 2 (bulk extraction)
3. **STOP and VALIDATE**: Verify backward errors unchanged, measure extraction percentage improvement
4. This is the safest starting point — pure internal optimization with no new data structures

### Incremental Delivery

1. Baseline → record pre-optimization measurements
2. US2 (bulk extraction) → measure extraction improvement → validate
3. US3 (bulk zeroing) → measure zeroing improvement → validate
4. US1 (scatter maps) → measure assembly improvement → validate
5. Full benchmark validation → verify all success criteria met

### Risk Ordering Rationale

The spec priorities (P1=scatter maps, P2=bulk extraction, P3=zeroing) reflect impact ordering. The implementation order (US2→US3→US1) follows the plan's risk ordering:
- **US2/US3 first**: Pure internal loop optimizations, no new data structures, lowest risk of regression
- **US1 last**: Adds new struct to AptpSymbolic, modifies function signatures, introduces fallback logic for delayed columns — most complex but highest reward

---

## Notes

- All changes are internal to `src/aptp/numeric.rs` and `src/aptp/symbolic.rs` — no public API changes
- The `diagnostic` feature timing instrumentation must continue working (FR-009)
- Both sequential and parallel factorization paths must work correctly (FR-010)
- Backward errors must be identical within floating-point tolerance (FR-007)
- Scatter map memory is unbounded, proportional to nnz (per spec clarification)
