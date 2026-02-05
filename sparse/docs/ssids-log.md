# SSIDS Development Log

## Phase 0.1: Literature Review & Reference Library

**Status**: Complete
**Branch**: `001-reference-library`
**Date**: 2026-02-05

### What Was Built

Four synthesis documents and a paper library completing the Phase 0.1 reference
library, organized following the structure from `ssids-plan.md`:

```
sparse/docs/references/
├── INDEX.md                          # Paper index with citations and categories
├── papers/                           # 14 academic papers (.md format)
│   ├── davis2016.md ... schenk2006.md
├── algorithms/
│   └── APTP-ALGORITHM.md            # Extracted pseudocode from papers
└── notes/
    ├── SPRAL-CODE-REVIEW.md          # Annotated SPRAL SSIDS architecture review
    └── FAER-INTEGRATION-NOTES.md     # faer component reuse map
```

1. **`docs/references/INDEX.md`** (~1400 words) -- Paper index cataloging all 14
   academic papers by category (Core APTP, Multifrontal Foundations, Pivoting
   Strategies, Ordering & Analysis, Infrastructure) with full bibliographic
   citations, relevance summaries, markdown quality flags, and RAL Technical
   Reports status.

2. **`docs/references/notes/SPRAL-CODE-REVIEW.md`** (~5100 words) -- Annotated
   architecture review of SPRAL's SSIDS implementation covering: executive
   summary of three-phase design, module-by-module overview (9 top-level
   modules), CPU factorization stack (24 files), ASCII data flow diagrams, APTP
   kernel internals (`ldlt_app.cxx`), key data structures (`akeep`, `fkeep`,
   `node_type`, `contrib_type`, diagonal D storage), and external dependencies
   with `ssids_options` configuration fields.

3. **`docs/references/notes/FAER-INTEGRATION-NOTES.md`** (~3500 words) --
   Component-by-component map of faer 0.24.0 sparse infrastructure with reuse
   classifications: 8 "direct use" components (CSC, AMD, COLAMD, elimination
   tree, triangular solve, permutations, workspace management, CSR), 2 "adapt"
   components (symbolic/numeric Cholesky), 2 "reference only" (dense
   Bunch-Kaufman, sparse LU). Includes deep dive on cholesky.rs
   Bunch-Kaufman vs APTP difference and integration strategy with build order.

4. **`docs/references/algorithms/APTP-ALGORITHM.md`** (~5000 words) -- Unified
   pseudocode reference covering: APTP overview with TPP/SBK comparison table,
   Algorithm 3.1 (main APTP loop), all 7 dense block kernels, Algorithm 4.1
   (complete pivoting), pivot acceptance/column delay mechanism, symbolic
   analysis (elimination tree + Gilbert-Ng-Peyton column counts), and
   multifrontal structure integration.

### Key Findings

1. **faer Bunch-Kaufman vs APTP**: faer's sparse LDLT uses pre-decided
   Bunch-Kaufman pivoting. The column-by-column update pattern and `ereach`
   traversal are reusable, but the pivot decision logic must be replaced entirely
   with a posteriori stability checks, column delay mechanism, and hybrid
   1x1/2x2 handling.

2. **SPRAL APTP kernel location**: The core APTP implementation lives in
   `src/ssids/cpu/kernels/ldlt_app.cxx` (~2600 lines). Key abstractions:
   `Column<T>` class tracking per-block-column elimination state,
   `check_threshold()` for a posteriori stability, and three pivot methods
   (APP_AGGRESSIVE, APP_BLOCK, TPP).

3. **Diagonal D storage convention**: SPRAL stores D^{-1} as a flat `2n` array
   with an `Inf` sentinel for 2x2 pivots (`d[2k+2] = Inf`). faer uses a
   different approach (`subdiag` + permutation arrays). The SPRAL convention is
   more compact for sparse solvers.

4. **Threshold parameter**: Default `u = 0.01` provides good stability/fill-in
   trade-off. Smaller `u` (0.001) reduces fill-in but may compromise accuracy;
   larger `u` (0.1) improves stability at the cost of more delayed pivots.

5. **faer reuse estimate**: ~70% of needed sparse infrastructure is available
   from faer as direct-use components. The remaining ~30% (APTP kernel, column
   delay, block-diagonal D handling) must be built from academic paper
   references.

### Readiness for Phase 0.2

All Phase 0.1 exit criteria are met:
- All key papers obtained and organized (INDEX.md)
- Algorithm pseudocode extracted (APTP-ALGORITHM.md)
- SPRAL source code reviewed (SPRAL-CODE-REVIEW.md)
- faer integration points identified (FAER-INTEGRATION-NOTES.md)

Phase 0.2 (Test Matrix Collection) can proceed. The reference library provides
sufficient algorithmic understanding to select appropriate test matrices and
define expected behaviors.
