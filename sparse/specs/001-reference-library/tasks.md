# Tasks: Complete Phase 0.1 Reference Library

**Input**: Design documents from `sparse/specs/001-reference-library/`
**Prerequisites**: plan.md (required), spec.md (required for user stories), research.md, quickstart.md

**Tests**: Not applicable — this is a documentation-only feature. Validation is via manual review against spec acceptance criteria.

**Organization**: Tasks are grouped by user story to enable independent completion and verification.

## Format: `[ID] [P?] [Story] Description`

- **[P]**: Can run in parallel (different files, no dependencies)
- **[Story]**: Which user story this task belongs to (e.g., US1, US2, US3, US4)
- Include exact file paths in descriptions

## Path Conventions

- Paper index: `references/ssids/INDEX.md` (relative to monorepo root)
- Synthesis documents: `sparse/dev/references/` (relative to monorepo root)
- SPRAL source: `references/spral/src/ssids/` (read-only reference)
- faer source: `references/faer-rs/faer/src/sparse/` (read-only reference)
- Academic papers: `references/ssids/*.md` and `references/ssids/*.pdf` (read-only reference)

---

## Phase 1: Setup

**Purpose**: Create directory structure for new deliverables

- [X] T001 Create `sparse/dev/references/` directory for synthesis documents
- [X] T002 Verify all 14 papers listed in Phase 0.1 of `sparse/dev/ssids-plan.md` are present in `references/ssids/` with both `.md` and `.pdf` versions; record any gaps

**Checkpoint**: Directory structure ready, paper inventory confirmed

---

## Phase 2: User Story 1 - Reference All Required Papers (Priority: P1)

**Goal**: Create a paper index so any developer can find and identify every reference paper by category, citation, and relevance.

**Independent Test**: Open `references/ssids/INDEX.md` and verify: (1) every paper in `references/ssids/` is listed with full citation, (2) papers are categorized into Core APTP / Multifrontal / Pivoting / Ordering / Supernodal / Infrastructure, (3) the RAL Technical Reports status is documented, (4) markdown quality flags are noted.

### Implementation for User Story 1

- [X] T003 [US1] Read the header/title section of each `.md` file in `references/ssids/` and record: filename, full title, authors, year, journal/venue
- [X] T004 [US1] Write paper index in `references/ssids/INDEX.md` containing: title row, table mapping each filename to full bibliographic citation, category assignment (Core APTP, Multifrontal Foundations, Pivoting Strategies, Ordering & Analysis, Supernodal, Infrastructure), and a one-sentence relevance summary per paper
- [X] T005 [US1] Add a "RAL Technical Reports" entry to `references/ssids/INDEX.md` documenting that the Duff/Hogg/Lopez 2020 SIAM J. Sci. Comput. publication supersedes the RAL Technical Reports on SSIDS v2.0
- [X] T006 [US1] Add a "Markdown Quality" section to `references/ssids/INDEX.md` noting any papers whose `.md` conversion appears incomplete or garbled (flag for re-conversion; PDF is authoritative source)
- [X] T007 [US1] Add a "faer Reference" entry to `references/ssids/INDEX.md` noting that the faer JOSS paper lives at `references/faer-rs/paper.md` (not in the ssids directory)

**Checkpoint**: INDEX.md is complete and self-contained. A developer can find any paper by category or author.

---

## Phase 3: User Story 2 - SPRAL Source Code Review Notes (Priority: P2)

**Goal**: Produce an annotated architecture review of SPRAL's SSIDS implementation so developers can understand the reference solver without reading Fortran source.

**Independent Test**: A developer unfamiliar with SPRAL can read `sparse/dev/references/SPRAL-CODE-REVIEW.md` and describe the three-phase architecture (analyze, factor, solve) and where APTP logic lives within 30 minutes.

### Implementation for User Story 2

- [X] T008 [US2] Write executive summary section (§1) of `sparse/dev/references/SPRAL-CODE-REVIEW.md`: three-phase architecture overview (analyse → factor → solve) with one paragraph per phase. Source: `references/spral/dev/Fortran/ssids.rst` and `references/spral/examples/Fortran/ssids.f90`
- [X] T009 [US2] Write module-by-module overview (§2) of `sparse/dev/references/SPRAL-CODE-REVIEW.md`: table of all key SSIDS source files with file path, line count estimate, responsibility, and key types. Cover: `ssids.f90`, `datatypes.f90`, `akeep.f90`, `fkeep.F90`, `inform.f90`, `contrib.f90`, `subtree.f90`, `anal.F90`, `core_analyse.f90`. Source: `references/spral/src/ssids/` headers
- [X] T010 [US2] Write CPU factorization stack section (§3) of `sparse/dev/references/SPRAL-CODE-REVIEW.md`: table covering CPU-specific files: `cpu/subtree.f90`, `cpu/NumericSubtree.hxx`, `cpu/NumericNode.hxx`, `cpu/SymbolicSubtree.hxx`, and all files in `cpu/kernels/`. Source: `references/spral/src/ssids/cpu/`
- [X] T011 [US2] Write data flow diagram section (§4) of `sparse/dev/references/SPRAL-CODE-REVIEW.md`: ASCII/text flow diagram showing data flow through analyse → factor → solve phases, including key data structures passed between phases (`akeep`, `fkeep`, contribution blocks). Source: research.md §R1
- [X] T012 [US2] Write APTP implementation details section (§5) of `sparse/dev/references/SPRAL-CODE-REVIEW.md`: identify `ldlt_app.cxx` as core APTP file, describe `Column<T>` class, `check_threshold()`, `apply_pivot()`, three pivot methods (APP_AGGRESSIVE, APP_BLOCK, TPP), and failed pivot handling. Source: `references/spral/src/ssids/cpu/kernels/ldlt_app.cxx`
- [X] T013 [US2] Write key data structures section (§6) of `sparse/dev/references/SPRAL-CODE-REVIEW.md`: document `ssids_akeep` fields (sptr, sparent, rptr, rlist, invp), `ssids_fkeep` fields (subtree array, scaling), `node_type` (nelim, ndelay, lcol, perm), `contrib_type` (val, rlist, ndelay, delay_perm, delay_val), and diagonal D storage convention (1x1, 2x2, zero). Source: `references/spral/src/ssids/akeep.f90`, `datatypes.f90`, `contrib.f90`
- [X] T014 [US2] Write external dependencies and configuration section (§7) of `sparse/dev/references/SPRAL-CODE-REVIEW.md`: list BLAS/LAPACK routines used (DGEMM, DSYRK, DPOTRF, DSYTRF, DTRSM, DTRSV), optional deps (METIS, hwloc, CUDA), CPU vs GPU path separation, and key `ssids_options` fields (pivot_method, u, small, cpu_block_size). Source: `references/spral/src/ssids/cpu/cpu_iface.f90`, `datatypes.f90`
- [X] T015 [US2] Add source citations to every section of `sparse/dev/references/SPRAL-CODE-REVIEW.md`: each section must reference the specific SPRAL files consulted with relative paths from `references/spral/`

**Checkpoint**: SPRAL-CODE-REVIEW.md covers all 7 planned sections. All facts are traceable to specific SPRAL source files.

---

## Phase 4: User Story 3 - faer Integration Points Documentation (Priority: P2)

**Goal**: Produce a component-by-component map of faer's sparse infrastructure with reuse classifications for the SSIDS project.

**Independent Test**: The integration notes document identifies at least 5 "direct use" components with specific faer module paths and version info, and explains the cholesky.rs Bunch-Kaufman vs APTP difference.

### Implementation for User Story 3

- [X] T016 [P] [US3] Write version and overview section (§1) of `sparse/dev/references/FAER-INTEGRATION-NOTES.md`: record faer version 0.24.0, commit `8dfccee`, MIT license, and summarize faer's sparse module structure. Source: `references/faer-rs/faer/src/sparse/` directory listing
- [X] T017 [P] [US3] Write component reuse classification table (§2) of `sparse/dev/references/FAER-INTEGRATION-NOTES.md`: table with columns Component | File Path | Size | Classification (Direct use / Adapt / Reference) | APTP Notes. Cover all 12 components from research.md §R2. Source: research.md §R2 component table
- [X] T018 [US3] Write per-component deep dive (§3) of `sparse/dev/references/FAER-INTEGRATION-NOTES.md`: for each "direct use" component (CSC, AMD, COLAMD, elimination tree, triangular solve, permutations, workspace management), provide: key types, public API entry points, and integration notes. Source: `references/faer-rs/faer/src/sparse/` module files
- [X] T019 [US3] Write cholesky.rs analysis section (§4) of `sparse/dev/references/FAER-INTEGRATION-NOTES.md`: document what faer's cholesky.rs implements (LLT + LDLT + LBLT, simplicial + supernodal, regularization), explain the Bunch-Kaufman vs APTP difference (pre-decided vs post-facto pivoting), identify reusable patterns (column-by-column update, ereach algorithm), and list what must be built new for APTP (a posteriori stability check, column delay, hybrid 1x1/2x2). Source: `references/faer-rs/faer/src/sparse/linalg/cholesky.rs`
- [X] T020 [US3] Write integration strategy and workspace patterns section (§5) of `sparse/dev/references/FAER-INTEGRATION-NOTES.md`: recommended approach for using faer (direct dependency vs pattern reuse), MemStack/StackReq pattern description, and suggested build order for SSIDS components. Source: research.md §R2 workspace pattern
- [X] T021 [US3] Add source citations to every section of `sparse/dev/references/FAER-INTEGRATION-NOTES.md`: each entry must reference specific faer file paths relative to `references/faer-rs/`

**Checkpoint**: FAER-INTEGRATION-NOTES.md has 5+ "direct use" components with file paths. The cholesky.rs analysis clearly explains what to reuse vs build new.

---

## Phase 5: User Story 4 - Algorithm Pseudocode Extraction (Priority: P3)

**Goal**: Produce a unified pseudocode reference document that a developer can implement from without flipping between multiple papers.

**Independent Test**: The algorithm document covers the complete simplicial APTP factorization loop (column processing, 1x1 acceptance, 2x2 fallback, column delay), plus symbolic analysis pseudocode, with citations to source papers.

### Implementation for User Story 4

- [X] T022 [US4] Write APTP overview section (§1) of `sparse/dev/references/APTP-ALGORITHM.md`: explain what APTP is, how it differs from TPP and SBK (PARDISO), and why it was chosen for SSIDS. Include comparison table. Source: duff2020.md §2, hogg2016.md §2.3
- [X] T023 [US4] Write main APTP factorization section (§2) of `sparse/dev/references/APTP-ALGORITHM.md`: reproduce Algorithm 3.1 from Duff et al. 2020 in clean pseudocode — 2D block-partitioned loop with Factor, ApplyN, ApplyT, Adjust, UpdateNN, UpdateNT, UpdateTN kernels. Define parameters nb, ib, u, nelim_j. Source: duff2020.md §3
- [X] T024 [US4] Write dense block kernels section (§3) of `sparse/dev/references/APTP-ALGORITHM.md`: pseudocode for Factor kernel (two-level recursion), ApplyN (subdiagonal), ApplyT (superdiagonal), Adjust (synchronize nelim_j, handle 2x2 boundary), and all three Update variants. Source: duff2020.md §3
- [X] T025 [US4] Write complete pivoting algorithm section (§4) of `sparse/dev/references/APTP-ALGORITHM.md`: reproduce Algorithm 4.1 from Duff et al. 2020 — 1x1/2x2 pivot selection on dense blocks with stability bounds and the 2x2 test |Δ| ≥ (1/2)|a_mt|^2. Source: duff2020.md §4
- [X] T026 [US4] Write pivot acceptance and column delay section (§5) of `sparse/dev/references/APTP-ALGORITHM.md`: threshold test |l_ij| > u^{-1}, fail-in-place mechanism (backup/restore), column delay to parent node, and threshold parameter guide (u values: 0.001, 0.01, 0.1 with trade-offs from Table 5). Source: duff2020.md §3, §7
- [X] T027 [US4] Write symbolic analysis section (§6) of `sparse/dev/references/APTP-ALGORITHM.md`: elimination tree computation (Liu's algorithm), column count algorithm from Gilbert et al. 1994 (row/column counts via postorder traversal and LCA), O(mα(m,n)) complexity. Source: gilbert1994.md §2, liu1992.md
- [X] T028 [US4] Write multifrontal structure section (§7) of `sparse/dev/references/APTP-ALGORITHM.md`: frontal matrix assembly, contribution block formation, three-phase solver structure (analyse/factor/solve), and how APTP integrates into the multifrontal framework. Source: liu1992.md, hogg2016.md §2.3, duff2020.md §5
- [X] T029 [US4] Add source citations to every section of `sparse/dev/references/APTP-ALGORITHM.md`: each pseudocode block must reference the paper, section number, and algorithm/equation number it derives from

**Checkpoint**: APTP-ALGORITHM.md is a standalone reference. A developer can implement the simplicial APTP factorization from this document alone.

---

## Phase 6: Polish & Cross-Cutting Concerns

**Purpose**: Final validation and development log update

- [X] T030 Cross-reference all four deliverables to ensure consistent terminology: verify that diagonal D storage convention, threshold parameter notation (u vs α), and pivot method names match across INDEX.md, SPRAL-CODE-REVIEW.md, FAER-INTEGRATION-NOTES.md, and APTP-ALGORITHM.md
- [X] T031 Validate all source citations in all four documents: every claim must trace to a specific paper (by filename in `references/ssids/`), SPRAL file (by path in `references/spral/`), or faer file (by path in `references/faer-rs/`)
- [X] T032 Update `sparse/dev/ssids-log.md` with a Phase 0.1 completion entry documenting: what was built, key findings (faer Bunch-Kaufman vs APTP difference, SPRAL APTP kernel location), and readiness for Phase 0.2
- [X] T033 Verify Phase 0.1 exit criteria from `sparse/dev/ssids-plan.md`: all key papers obtained and organized (INDEX.md), algorithm pseudocode extracted (APTP-ALGORITHM.md), SPRAL source code reviewed (SPRAL-CODE-REVIEW.md), faer integration points identified (FAER-INTEGRATION-NOTES.md)

---

## Dependencies & Execution Order

### Phase Dependencies

- **Setup (Phase 1)**: No dependencies — start immediately
- **User Story 1 (Phase 2)**: Depends on T002 (paper inventory) from Setup
- **User Story 2 (Phase 3)**: No dependency on US1; can start after Setup
- **User Story 3 (Phase 4)**: No dependency on US1 or US2; can start after Setup
- **User Story 4 (Phase 5)**: No dependency on US1, US2, or US3; can start after Setup
- **Polish (Phase 6)**: Depends on ALL user stories being complete

### User Story Dependencies

- **User Story 1 (P1)**: Independent — only reads existing paper files
- **User Story 2 (P2)**: Independent — only reads SPRAL source
- **User Story 3 (P2)**: Independent — only reads faer source
- **User Story 4 (P3)**: Independent — only reads academic papers

All four user stories produce separate files with no cross-dependencies. They CAN be executed in parallel.

### Within Each User Story

- Earlier sections (overview, tables) before later sections (deep dives, citations)
- Citation pass (final task in each story) depends on all content being written
- Each story is complete when all its tasks are done

### Parallel Opportunities

- **US2 + US3 + US4 can all start immediately** after Setup (Phase 1)
- Within US3: T016 and T017 are marked [P] (independent overview + table)
- US1 is low-effort and can be completed quickly as an MVP

---

## Parallel Example: All User Stories

```
# After Setup (Phase 1) completes, launch all stories in parallel:

Stream A (US1 - Paper Index):    T003 → T004 → T005 → T006 → T007
Stream B (US2 - SPRAL Review):   T008 → T009 → T010 → T011 → T012 → T013 → T014 → T015
Stream C (US3 - faer Notes):     T016+T017 → T018 → T019 → T020 → T021
Stream D (US4 - APTP Algorithm): T022 → T023 → T024 → T025 → T026 → T027 → T028 → T029

# After all streams complete:
Polish: T030 → T031 → T032 → T033
```

---

## Implementation Strategy

### MVP First (User Story 1 Only)

1. Complete Phase 1: Setup (T001-T002)
2. Complete Phase 2: User Story 1 - Paper Index (T003-T007)
3. **STOP and VALIDATE**: INDEX.md is usable — developers can find any paper
4. Commit and continue to remaining stories

### Incremental Delivery

1. US1 (Paper Index) — immediate value, low effort
2. US2 (SPRAL Review) + US3 (faer Notes) — parallel, both P2 priority
3. US4 (APTP Algorithm) — depends on understanding from papers but not on other deliverables
4. Polish — consistency check after all documents exist

### Single Developer Strategy

Recommended order for sequential execution:

1. Setup → US1 (quick win, builds familiarity with papers)
2. US4 (algorithm extraction — deepens understanding needed for US2/US3)
3. US2 (SPRAL review — now informed by algorithm understanding)
4. US3 (faer notes — now informed by SPRAL architecture knowledge)
5. Polish

---

## Notes

- All tasks produce Markdown files — no code compilation or testing needed
- Source materials are read-only references in `references/` — never modify them
- Clean room compliance: only consult BSD-3 (SPRAL) and MIT (faer) licensed code; cite all sources
- Commit after completing each user story (natural commit points)
- Each deliverable should be 2000-5000 words based on plan.md scope estimate
