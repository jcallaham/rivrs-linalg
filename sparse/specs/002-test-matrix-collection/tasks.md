# Tasks: Test Matrix Collection Assembly

**Input**: Design documents from `/specs/002-test-matrix-collection/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, data-model.md, quickstart.md

**Tests**: Not explicitly requested in the feature specification. Tests are omitted. The validation script (T019) serves as the automated quality gate.

**Organization**: Tasks are grouped by user story to enable independent implementation and testing of each story.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3)
- Include exact file paths in descriptions

---

## Phase 1: Setup (Shared Infrastructure)

**Purpose**: Directory structure, Git LFS, Python environment, metadata schema initialization

- [x] T001 Create test-data directory structure: `sparse/test-data/`, `sparse/test-data/hand-constructed/`, `sparse/test-data/suitesparse/easy-indefinite/`, `sparse/test-data/suitesparse/hard-indefinite/`, `sparse/test-data/suitesparse/positive-definite/`, `sparse/test-data/scripts/`
- [x] T002 Add Git LFS tracking rule `sparse/test-data/**/*.mtx filter=lfs diff=lfs merge=lfs -text` to `/workspace/rivrs-linalg/.gitattributes` and run `git lfs install` (if available; if git-lfs is not installed, document the prerequisite in a `sparse/test-data/README.md` and proceed without LFS for now)
- [x] T003 Create `sparse/test-data/scripts/requirements.txt` listing Python dependencies: `ssgetpy`, `numpy`, `scipy`
- [x] T004 Initialize empty metadata index at `sparse/test-data/metadata.json` with schema skeleton: `{"schema_version": "1.0", "generated": "", "total_count": 0, "matrices": []}`

---

## Phase 2: Foundational (Blocking Prerequisites)

**Purpose**: Shared utility code used by both the hand-constructed generator (US1) and the download script (US2)

**CRITICAL**: No user story work can begin until this phase is complete

- [x] T005 Create Matrix Market writer utility function in `sparse/test-data/scripts/mtx_utils.py` — function `write_symmetric_mtx(filepath, n, rows, cols, vals, comment)` that outputs valid `%%MatrixMarket matrix coordinate real symmetric` files with 1-indexed entries (lower triangle only). Also include `read_mtx_header(filepath)` that returns `(n, nnz)` by parsing the header line.
- [x] T006 Create metadata builder utility in `sparse/test-data/scripts/metadata_utils.py` — functions: `create_matrix_entry(name, source, category, path, size, nnz, properties, **kwargs)` returning a dict matching the MetadataIndex schema from data-model.md; `merge_metadata(existing_index, new_entries)` that merges entries into the index and updates `total_count` and `generated` fields; `write_metadata(index, filepath)` that writes the JSON file with indent=2.
- [x] T007 Create factorization JSON writer utility in `sparse/test-data/scripts/factorization_utils.py` — function `write_factorization(filepath, matrix_name, permutation, l_entries, d_blocks, inertia, notes)` that writes the exact factorization companion `.json` file matching the ExactFactorization schema in data-model.md. Include `verify_factorization(A_dense, P, L, D)` that checks `||P^T A P - L D L^T|| < epsilon` using numpy.

**Checkpoint**: Utility library ready — hand-constructed and SuiteSparse scripts can now be built

---

## Phase 3: User Story 1 — Run Solver Tests Against Known Matrices (Priority: P1) MVP

**Goal**: Produce 15+ hand-constructed matrices (5x5 to 20x20) in Matrix Market format, each with a companion JSON file containing the exact LDL^T factorization (L, D, P, inertia). These are the foundation for all solver correctness testing.

**Independent Test**: Load any `.mtx` file from `sparse/test-data/hand-constructed/`, verify it parses correctly, load the companion `.json`, and confirm `P^T A P = L D L^T` holds numerically.

### Implementation for User Story 1

- [x] T008 [P] [US1] Implement arrow matrix generator in `sparse/test-data/scripts/generate_hand_constructed.py` — function `generate_arrow(n, positive_definite=True)` that creates an n×n arrow matrix (dense first row/column + diagonal). Compute LDL^T factorization analytically via Schur complement. Output `arrow-5-pd.mtx` + `.json` and `arrow-10-indef.mtx` + `.json` and `arrow-15-indef.mtx` + `.json` to `sparse/test-data/hand-constructed/`.
- [x] T009 [P] [US1] Implement tridiagonal matrix generator in `sparse/test-data/scripts/generate_hand_constructed.py` — function `generate_tridiagonal(n, positive_definite=True)` that creates an n×n tridiagonal matrix. Compute LDL^T via the recurrence `d_i = a_i - b_{i-1}^2/d_{i-1}`. Output `tridiag-5-pd.mtx` + `.json`, `tridiag-10-indef.mtx` + `.json`, `tridiag-20-indef.mtx` + `.json` to `sparse/test-data/hand-constructed/`.
- [x] T010 [P] [US1] Implement block diagonal matrix generator in `sparse/test-data/scripts/generate_hand_constructed.py` — function `generate_block_diagonal(block_sizes, block_types)` that creates a block-diagonal matrix with independently factored blocks (one PD, one indefinite, one with 2×2 pivots). Output `block-diag-15.mtx` + `.json` to `sparse/test-data/hand-constructed/`.
- [x] T011 [P] [US1] Implement bordered block diagonal generator in `sparse/test-data/scripts/generate_hand_constructed.py` — function `generate_bordered_block(block_sizes, border_width)` adapted from SPRAL `gen_bordered_block_diag` pattern (BSD-3). Creates blocks with a dense border causing fill-in. Output `bordered-block-20.mtx` + `.json` to `sparse/test-data/hand-constructed/`.
- [x] T012 [P] [US1] Implement stress-test matrix generators in `sparse/test-data/scripts/generate_hand_constructed.py` — three functions: (a) `generate_delayed_pivot_stress(n)` — zero diagonal with large off-diagonal, forcing maximum 2×2 pivots; (b) `generate_fill_in_stress(n)` — matrix with a dense column causing worst-case fill-in; (c) `generate_ill_conditioned(n, cond)` — near-singular matrix with known exact factorization. Output `stress-delayed-pivots.mtx` + `.json`, `stress-fill-in.mtx` + `.json`, `stress-ill-cond.mtx` + `.json` to `sparse/test-data/hand-constructed/`.
- [x] T013 [P] [US1] Implement degenerate case generators in `sparse/test-data/scripts/generate_hand_constructed.py` — four functions adapted from SPRAL test patterns (BSD-3): (a) `generate_zero_diagonal()` — 3×3 from SPRAL `simple_mat_zero_diag`; (b) `generate_singular()` — 3×3 from SPRAL `simple_sing_mat`; (c) `generate_trivial_1x1()` — single element; (d) `generate_trivial_2x2()` — 2×2 requiring 2×2 pivot. Output `zero-diagonal-3.mtx` + `.json`, `singular-3.mtx` + `.json`, `trivial-1x1.mtx` + `.json`, `trivial-2x2.mtx` + `.json` to `sparse/test-data/hand-constructed/`.
- [x] T014 [US1] Add `main()` entry point to `sparse/test-data/scripts/generate_hand_constructed.py` that calls all generators, verifies every factorization numerically (`verify_factorization` from T007), collects metadata entries for all generated matrices, and merges them into `sparse/test-data/metadata.json` using `metadata_utils.py`. Print summary: matrix count, any verification failures.
- [x] T015 [US1] Run `generate_hand_constructed.py` and verify: (a) all 15+ `.mtx` files created in `sparse/test-data/hand-constructed/`; (b) all companion `.json` files created; (c) all factorizations pass numerical verification; (d) `metadata.json` updated with hand-constructed entries.

**Checkpoint**: 15+ hand-constructed matrices with verified exact factorizations. US1 is independently testable: load any `.mtx`, parse it, verify against companion `.json`.

---

## Phase 4: User Story 2 — Validate Solver on Real-World Sparse Problems (Priority: P2)

**Goal**: Download ~67 curated SuiteSparse matrices (30 easy indefinite, 18 hard indefinite, 19 positive definite) organized by difficulty category, with complete metadata for each.

**Independent Test**: Load any downloaded SuiteSparse `.mtx` file, verify it parses correctly, and confirm its metadata entry (size, nnz, category, properties) matches the actual file contents.

### Implementation for User Story 2

- [x] T016 [US2] Create `sparse/test-data/scripts/download_suitesparse.py` with curated matrix lists — define three Python lists: `EASY_INDEFINITE` (~30 entries), `HARD_INDEFINITE` (~18 entries), `POSITIVE_DEFINITE` (~19 entries) containing `"Group/Name"` strings matching the curated lists in plan.md. Include `LARGE_THRESHOLD_ROWS = 200000` for matrices that are metadata-only (download-on-demand). Each list entry should also carry paper references and properties (difficulty, killer_case flag, expected_delayed_pivots) as inline annotations.
- [x] T017 [US2] Implement download logic in `sparse/test-data/scripts/download_suitesparse.py` — function `download_matrix(group_name, category, dest_base, max_rows=None)` that uses `ssgetpy.search()` to find the matrix, calls `mat.download(format='MM', destpath=..., extract=True)` to download, and returns a metadata dict with SuiteSparse properties (id, group, kind, psym, nsym, rows, cols, nnz). Handle download failures gracefully (log warning, skip matrix, continue).
- [x] T018 [US2] Implement CLI interface in `sparse/test-data/scripts/download_suitesparse.py` — `argparse` with options: `--category` (filter by easy-indefinite/hard-indefinite/positive-definite), `--max-rows N` (size filter), `--name "Group/Name"` (specific matrix), `--verify-only` (check existing files without downloading), `--dry-run` (list what would be downloaded). Main function iterates over curated lists, downloads to `sparse/test-data/suitesparse/{category}/`, generates metadata entries, and merges into `sparse/test-data/metadata.json`.
- [x] T019 [US2] Run `download_suitesparse.py --dry-run` to verify curated lists resolve correctly via ssgetpy, then run actual download for a small subset (`--max-rows 50000`) to validate the download pipeline end-to-end. Verify: (a) `.mtx` files appear in correct category subdirectories; (b) metadata.json is updated with SuiteSparse entries; (c) downloaded files parse as valid Matrix Market.
- [x] T020 [US2] Run full download: `download_suitesparse.py` (all categories, default size limits). For matrices exceeding 50MB or `LARGE_THRESHOLD_ROWS`, verify they appear in `metadata.json` with `"in_repo": false` and a valid `"download_command"` field.

**Checkpoint**: ~67 SuiteSparse matrices downloaded and categorized. US2 is independently testable: load any SuiteSparse matrix, verify it parses and its metadata matches the file.

---

## Phase 5: User Story 3 — Browse and Select Matrices by Properties (Priority: P2)

**Goal**: Ensure the metadata index is complete, queryable, and consistent — every matrix has all required fields, and property-based filtering returns correct results.

**Independent Test**: Load `metadata.json`, filter by various property combinations (size, category, difficulty, killer_case), verify results are non-empty and correct.

### Implementation for User Story 3

- [x] T021 [US3] Create collection validation script at `sparse/test-data/scripts/validate_collection.py` — functions: (a) `validate_metadata(index)` — check every entry has all required fields per data-model.md, no placeholder values, names are unique, symmetric=true for all; (b) `validate_files(index, base_dir)` — for entries with `in_repo: true`, verify `.mtx` file exists and `(n, nnz)` from header matches metadata; for hand-constructed entries, verify `.json` factorization file exists; (c) `validate_properties(index)` — check `positive_definite` and `indefinite` are mutually exclusive, `killer_case` entries have `difficulty: "hard"`, hand-constructed entries have `factorization_path`.
- [x] T022 [US3] Add summary statistics to `validate_collection.py` — function `print_summary(index)` that reports: total matrix count, count by category, count by source, size range (min/max), count of killer cases, count of matrices with exact factorizations, count of in-repo vs download-on-demand matrices. Verify SC-001 (>=70 total), SC-004 (4 orders of magnitude size range), SC-005 (>=5 killer cases).
- [x] T023 [US3] Add property query examples to `validate_collection.py` — function `test_queries(index)` that runs and verifies: (a) filter for indefinite matrices < 1000×1000 → returns all hand-constructed indefinite matrices; (b) filter for `killer_case: true` → returns >=5 matrices; (c) filter for `category: "positive-definite"` → returns >=19 matrices; (d) filter for matrices with `factorization_path` → returns >=15 matrices. Print results for each query.
- [x] T024 [US3] Run `validate_collection.py` and fix any issues found — iterate until all validations pass and all summary statistics meet success criteria SC-001 through SC-007.

**Checkpoint**: Complete, consistent metadata index. US3 is independently testable: load `metadata.json` and run property queries.

---

## Phase 6: User Story 4 — Reproduce Academic Paper Results (Priority: P3)

**Goal**: Ensure every matrix referenced in the core APTP papers is present in the collection with proper citation and classification.

**Independent Test**: For each paper (Duff 2020, Hogg 2016), verify that all matrices from their benchmark tables appear in the collection with correct `paper_references` entries.

### Implementation for User Story 4

- [x] T025 [US4] Add paper reference annotations to metadata — update `download_suitesparse.py` curated lists and the `generate_hand_constructed.py` metadata entries to include `paper_references` arrays citing the specific paper and table for each matrix. References from research.md: Duff, Hogg, Lopez (2020) Tables 1-2; Hogg, Ovtchinnikov, Scott (2016) Table III; Duff & Pralet (2005) Tables 3.1-3.2; Schenk & Gärtner (2006).
- [x] T026 [US4] Create paper coverage report in `validate_collection.py` — function `check_paper_coverage(index)` that verifies: (a) all matrices from Duff et al. (2020) Tables 1-2 are present; (b) all matrices from Hogg et al. (2016) Table III are present; (c) all matrices have at least one `paper_references` entry or are documented as "supplementary" (added beyond paper sets for coverage). Print coverage percentage per paper.
- [x] T027 [US4] Run paper coverage check and ensure 100% coverage for the two core APTP papers. Add any missing matrices to the curated download lists and re-run the download.

**Checkpoint**: All paper-referenced matrices present with proper citations. US4 is independently testable: check paper coverage report.

---

## Phase 7: Polish & Cross-Cutting Concerns

**Purpose**: Documentation, plan updates, final validation

- [x] T028 [P] Update `sparse/dev/ssids-plan.md` Phase 0.2 success criteria checkboxes to reflect completion status
- [x] T029 [P] Add Phase 0.2 entry to `sparse/dev/ssids-log.md` documenting: what was built (directory structure, scripts, matrix counts by category), key decisions (curated lists from papers, LFS strategy, 50MB threshold), and any issues encountered
- [x] T030 [P] Create `sparse/test-data/README.md` with setup instructions adapted from quickstart.md: prerequisites (Python 3.8+, ssgetpy, numpy, scipy, git-lfs), first-time setup commands, directory layout, how to add new matrices
- [x] T031 Run final `validate_collection.py` pass and confirm all success criteria: SC-001 (>=70 matrices), SC-002 (all parseable), SC-003 (all metadata complete), SC-004 (4+ orders of magnitude), SC-005 (>=5 killers), SC-006 (verified factorizations), SC-007 (queryable index)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — can start immediately
- **Foundational (Phase 2)**: Depends on Setup completion — BLOCKS all user stories
- **US1 (Phase 3)**: Depends on Foundational (T005, T006, T007) — can start immediately after Phase 2
- **US2 (Phase 4)**: Depends on Foundational (T005, T006) — can run in PARALLEL with US1
- **US3 (Phase 5)**: Depends on US1 (T015) and US2 (T020) — needs complete metadata from both
- **US4 (Phase 6)**: Depends on US2 (T020) — needs SuiteSparse matrices downloaded
- **Polish (Phase 7)**: Depends on US3 (T024) and US4 (T027)

### User Story Dependencies

```
Phase 1 (Setup)
    │
Phase 2 (Foundational: T005, T006, T007)
    │
    ├── Phase 3 (US1: Hand-constructed) ──┐
    │                                      ├── Phase 5 (US3: Metadata validation)
    ├── Phase 4 (US2: SuiteSparse) ───────┤
    │                                      └── Phase 6 (US4: Paper coverage)
    │
    └── Phase 7 (Polish) ← depends on US3 + US4
```

### Parallel Opportunities

**Within Phase 2**: T005, T006, T007 are on different files — all [P]
**Within Phase 3 (US1)**: T008–T013 are independent generator functions in the same file — all [P] (they write different matrix files, but share a source file; can be written as independent functions then assembled)
**Across US1 and US2**: Phase 3 and Phase 4 can run in parallel after Phase 2 completes
**Within Phase 7**: T028, T029, T030 are on different files — all [P]

---

## Parallel Example: Phases 3 and 4

```bash
# After Phase 2 completes, launch US1 and US2 in parallel:

# Agent A: Hand-constructed matrices (US1)
Task: T008 "Implement arrow matrix generator in generate_hand_constructed.py"
Task: T009 "Implement tridiagonal matrix generator in generate_hand_constructed.py"
Task: T010 "Implement block diagonal generator in generate_hand_constructed.py"
Task: T011 "Implement bordered block diagonal generator in generate_hand_constructed.py"
Task: T012 "Implement stress-test generators in generate_hand_constructed.py"
Task: T013 "Implement degenerate case generators in generate_hand_constructed.py"
# Then: T014 (main entry point), T015 (run and verify)

# Agent B: SuiteSparse download (US2)
Task: T016 "Create download script with curated lists"
Task: T017 "Implement download logic"
Task: T018 "Implement CLI interface"
# Then: T019 (test download), T020 (full download)
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001–T004)
2. Complete Phase 2: Foundational (T005–T007)
3. Complete Phase 3: User Story 1 (T008–T015)
4. **STOP and VALIDATE**: Run `generate_hand_constructed.py`, verify 15+ matrices with exact factorizations
5. This is a viable MVP: hand-constructed matrices are sufficient for initial solver development

### Incremental Delivery

1. Setup + Foundational → Infrastructure ready
2. US1 (hand-constructed) → Foundation for solver testing (**MVP**)
3. US2 (SuiteSparse) → Real-world validation data
4. US3 (metadata validation) → Quality gate on complete collection
5. US4 (paper coverage) → Academic reproducibility
6. Polish → Documentation and plan updates

### Task Summary

| Phase | Tasks | Parallel Tasks |
|-------|-------|---------------|
| 1 Setup | T001–T004 (4) | 0 |
| 2 Foundational | T005–T007 (3) | 3 |
| 3 US1 Hand-constructed | T008–T015 (8) | 6 |
| 4 US2 SuiteSparse | T016–T020 (5) | 0 |
| 5 US3 Metadata | T021–T024 (4) | 0 |
| 6 US4 Paper coverage | T025–T027 (3) | 0 |
| 7 Polish | T028–T031 (4) | 3 |
| **Total** | **31 tasks** | **12 parallelizable** |

---

## Notes

- [P] tasks = different files, no dependencies
- [Story] label maps task to specific user story for traceability
- US1 and US2 are independently testable and can run in parallel after Phase 2
- US3 requires both US1 and US2 to complete (needs full metadata)
- US4 can start as soon as US2 completes (only needs SuiteSparse matrices)
- Constitution principle III (TDD) is satisfied: this feature IS the test data, not solver code
- Git LFS may not be available in this container — T002 handles this gracefully with a fallback
- The 50MB threshold for in-repo storage (Decision 9 from research.md) is enforced in T018/T020
