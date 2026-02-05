# Implementation Plan: Complete Phase 0.1 Reference Library

**Branch**: `001-reference-library` | **Date**: 2026-02-05 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `sparse/specs/001-reference-library/spec.md`

## Summary

Complete the Phase 0.1 reference library by producing four synthesis documents from existing materials: a paper index (INDEX.md), SPRAL code review (SPRAL-CODE-REVIEW.md), faer integration notes (FAER-INTEGRATION-NOTES.md), and APTP algorithm pseudocode extraction (APTP-ALGORITHM.md). All 14 academic papers are already compiled; no external acquisition needed. This is a documentation-only feature — no Rust code is written or modified.

## Technical Context

**Language/Version**: N/A (documentation-only feature; no code written)
**Primary Dependencies**: N/A
**Storage**: Markdown files in `references/ssids/` and `sparse/docs/references/`
**Testing**: Manual review against spec acceptance criteria
**Target Platform**: N/A
**Project Type**: Documentation / research synthesis
**Performance Goals**: N/A
**Constraints**: Must maintain clean room audit trail; all sources must be cited
**Scale/Scope**: 4 documents, ~2000-5000 words each

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Documentation accuracy verified against source papers and code |
| II. Clean Room Implementation | PASS | Only consults BSD-3 (SPRAL, LAPACK) and MIT (faer) licensed code; all sources cited |
| III. Test-Driven Development | N/A | No code written; documentation validated against acceptance criteria |
| IV. Algorithm Documentation & Academic Attribution | PASS | Core deliverable — every document cites sources |
| V. Numerical Stability & Robustness | N/A | No numerical code; algorithm descriptions note stability properties |
| VI. Structured Development Discipline | PASS | Phase 0.1 must complete before Phase 0.2 or any solver code |
| VII. Code Quality & Rust Best Practices | N/A | No Rust code written |

**Gate result**: PASS — no violations. This feature is pure documentation aligned with Phase 0.1 exit criteria.

## Project Structure

### Documentation (this feature)

```text
sparse/specs/001-reference-library/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Phase 0 research findings
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Deliverables (repository)

```text
references/ssids/
├── INDEX.md                          # NEW: Paper index with citations and categories

sparse/docs/references/
├── SPRAL-CODE-REVIEW.md              # NEW: Annotated SPRAL SSIDS architecture review
├── FAER-INTEGRATION-NOTES.md         # NEW: faer component reuse map
└── APTP-ALGORITHM.md                 # NEW: Extracted pseudocode from papers
```

**Structure Decision**: Documentation files only. INDEX.md goes in `references/ssids/` alongside the papers it indexes. The three synthesis documents go in `sparse/docs/references/` as they are specific to the sparse solver project.

## Research Findings

All research is consolidated in [research.md](research.md). Key findings:

### SPRAL SSIDS Architecture (from code review)

- **Three-phase design**: `ssids_analyse()` → `ssids_factor()` → `ssids_solve()`
- **Core APTP kernel**: `src/ssids/cpu/kernels/ldlt_app.cxx` (~2600 lines)
- **Subtree decomposition**: Assembly tree partitioned for NUMA-aware parallelism
- **Three pivot methods**: APP_AGGRESSIVE (single-pass), APP_BLOCK (default, block-rollback), TPP (serial fallback)
- **Threshold test**: `|l_ij| > 1/u` where default `u = 0.01`
- **Diagonal D storage**: 2-element encoding per pivot (1x1: `[a, 0.0]`, 2x2: `[a, b, INF, c]`, zero: `[0, 0]`)
- **External deps**: BLAS/LAPACK, optional METIS ordering, optional CUDA

### faer Sparse Infrastructure (from code review)

| Component | File | Reuse | Notes |
|-----------|------|-------|-------|
| CSC storage | `sparse/csc/` | Direct use | Full API, symmetric matrix support |
| AMD ordering | `sparse/linalg/amd.rs` (24 KB) | Direct use | Directly callable |
| COLAMD ordering | `sparse/linalg/colamd.rs` (19 KB) | Direct use | Column AMD variant |
| Elimination trees | `sparse/linalg/cholesky.rs` | Direct use | `ghost_prefactorize_symbolic_cholesky()` |
| Triangular solve | `sparse/linalg/triangular_solve.rs` | Direct use | Forward/backward substitution |
| Permutation utilities | `perm/` | Direct use | Forward/inverse permutations |
| Workspace management | `MemStack`/`StackReq` | Direct use | Stack-based temp allocation |
| Sparse LDL^T (symbolic) | `sparse/linalg/cholesky.rs` | Adapt | Reuse ereach, etree; adapt numeric phase |
| Sparse LDL^T (numeric) | `sparse/linalg/cholesky.rs` (~164 KB) | Adapt | Column-by-column pattern; replace Bunch-Kaufman with APTP |
| Dense Bunch-Kaufman | dense module | Reference only | Understand 2x2 block handling |

**Critical finding**: faer's sparse LDL^T uses **pre-decided Bunch-Kaufman pivoting**, not post-facto APTP. The column-by-column update pattern can be adapted, but the pivot decision logic must be replaced entirely with:
- A posteriori stability check after elimination
- Column delay mechanism for problematic pivots
- Hybrid 1x1 + 2x2 decision based on post-factor stability

**faer version reviewed**: 0.24.0 (commit `8dfccee`)

### APTP Algorithm (from paper extraction)

Key algorithms extracted from Duff/Hogg/Lopez (2020) and Hogg et al. (2016):

1. **Main APTP factorization** (Algorithm 3.1): 2D block-partitioned loop with Factor, ApplyN, ApplyT, Adjust, and Update kernels
2. **Complete pivoting** (Algorithm 4.1): Dense block-level 1x1/2x2 pivot selection
3. **Pivot acceptance test**: `|l_ij| > u^{-1}` checked after applying pivot
4. **Column delay**: Fail-in-place with backup/restore; pass to parent node if still failing
5. **Threshold parameter**: `u = 0.01` default; smaller u → sparser factor, fewer delays
6. **2x2 pivot test**: `|Δ| ≥ (1/2)|a_mt|^2` where `Δ = a_mm*a_tt - a_mt^2`

From Gilbert et al. (1994):
7. **Symbolic column counts**: O(mα(m,n)) algorithm for predicting nonzero structure

## Deliverable Specifications

### D1: INDEX.md (Paper Index)

**Location**: `references/ssids/INDEX.md`
**Satisfies**: FR-001, FR-002, FR-003

**Structure**:
- Table mapping filename → full citation → category → relevance summary
- Categories: Core APTP, Multifrontal Foundations, Pivoting Strategies, Ordering & Analysis, Supernodal, Infrastructure
- Note on RAL Technical Reports superseded by Duff/Hogg/Lopez 2020 journal publication
- Quality flags for any markdown conversions that need re-review

**Source data**: Existing files in `references/ssids/`, paper headers already verified during spec phase.

### D2: SPRAL-CODE-REVIEW.md

**Location**: `sparse/docs/references/SPRAL-CODE-REVIEW.md`
**Satisfies**: FR-004

**Structure**:
1. Executive summary: three-phase architecture
2. Module-by-module overview (file path, responsibility, key types)
3. Data flow diagram: analyse → factor → solve
4. APTP implementation details: which files, pivot decision flow, delayed column handling
5. Key data structures: `ssids_akeep`, `ssids_fkeep`, `node_type`, `contrib_type`
6. Diagonal D storage format
7. External dependencies (BLAS/LAPACK, METIS)
8. CPU vs GPU path separation
9. Configuration options (`ssids_options` key fields)

**Sources**: SPRAL BSD-3 source at `references/spral/src/ssids/`, SPRAL documentation at `references/spral/docs/Fortran/ssids.rst`, SPRAL examples.

### D3: FAER-INTEGRATION-NOTES.md

**Location**: `sparse/docs/references/FAER-INTEGRATION-NOTES.md`
**Satisfies**: FR-005

**Structure**:
1. faer version and repository commit reviewed
2. Component-by-component map with reuse classification table
3. For each component: file path, API summary, reuse assessment, APTP adaptation needed
4. Deep dive: cholesky.rs — Bunch-Kaufman vs APTP differences
5. Recommended integration strategy: which faer APIs to call, which patterns to adapt, what to build from scratch
6. Workspace management patterns to follow

**Sources**: faer source at `references/faer-rs/faer/src/sparse/`, faer JOSS paper at `references/faer-rs/paper.md`.

### D4: APTP-ALGORITHM.md

**Location**: `sparse/docs/references/APTP-ALGORITHM.md`
**Satisfies**: FR-006

**Structure**:
1. Overview: what APTP is and why it differs from traditional TPP
2. Main APTP factorization (Algorithm 3.1 from Duff et al. 2020) — full pseudocode
3. Dense block kernels: Factor, ApplyN, ApplyT, Adjust, Update variants
4. Complete pivoting algorithm (Algorithm 4.1) — 1x1 and 2x2 selection
5. Pivot acceptance test — threshold formula and stability bounds
6. Column delay mechanism — fail-in-place, backup/restore, pass to parent
7. Threshold parameter guide — values, trade-offs, empirical results
8. Symbolic analysis: elimination tree traversal, column count algorithm (Gilbert et al. 1994)
9. Multifrontal structure: frontal matrix assembly, contribution blocks
10. Comparison table: TPP vs APTP vs SBK (PARDISO)

**Sources**: Duff/Hogg/Lopez (2020) §3-4-7, Hogg et al. (2016) §2-3, Gilbert et al. (1994) §2, Liu (1992).

## Complexity Tracking

No constitution violations — no complexity tracking needed.
