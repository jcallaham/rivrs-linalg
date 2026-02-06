# SSIDS Development Log

## Phase 0.2: Test Matrix Collection Assembly

**Status**: Complete
**Branch**: `002-test-matrix-collection`
**Date**: 2026-02-05

### What Was Built

A comprehensive test matrix collection with 82 matrices (15 hand-constructed +
67 SuiteSparse), organized under `sparse/test-data/`:

```
sparse/test-data/
├── metadata.json                    # Complete index (82 matrices)
├── README.md                        # Setup instructions
├── hand-constructed/                # 15 matrices with exact LDL^T factorizations
│   ├── {name}.mtx                   # Matrix Market format
│   └── {name}.json                  # Exact L, D, P, inertia
├── suitesparse/
│   ├── easy-indefinite/             # 30 matrices from APTP paper benchmarks
│   ├── hard-indefinite/             # 18 matrices (includes killer cases)
│   └── positive-definite/           # 19 matrices for fast-path validation
├── suitesparse-ci/                  # 10 representative matrices (plain git, ~73MB)
│   ├── easy-indefinite/{name}.mtx
│   ├── hard-indefinite/{name}.mtx
│   └── positive-definite/{name}.mtx
└── scripts/
    ├── requirements.txt             # Python dependencies (ssgetpy, numpy, scipy)
    ├── mtx_utils.py                 # Matrix Market I/O
    ├── metadata_utils.py            # Metadata index builder
    ├── factorization_utils.py       # Factorization writing and verification
    ├── generate_hand_constructed.py  # Hand-constructed matrix generator
    ├── download_suitesparse.py      # SuiteSparse download with curated lists
    └── validate_collection.py       # Collection validation and reporting
```

### Hand-Constructed Matrices (15)

Six categories of matrices with analytically verified LDL^T factorizations:

1. **Arrow matrices** (3): 5x5 PD, 10x10 indefinite, 15x15 indefinite
2. **Tridiagonal matrices** (3): 5x5 PD, 10x10 indefinite, 20x20 indefinite
3. **Block diagonal** (1): 15x15 with PD, indefinite, and 2x2-pivot blocks
4. **Bordered block diagonal** (1): 20x20 adapted from SPRAL pattern (BSD-3)
5. **Stress tests** (3): zero diagonal (forces 2x2 pivots), dense column (fill-in), ill-conditioned
6. **Degenerate cases** (4): 3x3 zero-diagonal, 3x3 singular, 1x1 trivial, 2x2 requiring 2x2 pivot

All 15 factorizations verified: `||P^T A P - L D L^T|| / ||A|| < 1e-10`.

### SuiteSparse Matrices (67)

Curated from APTP benchmark papers:
- **Easy indefinite (30)**: From Duff, Hogg, Lopez (2020) Table 1 and Hogg et al. (2016)
- **Hard indefinite (18)**: From Duff et al. (2020) Table 2 (KKT/saddle-point systems)
- **Positive definite (19)**: From Hogg et al. (2016) Table III (dagger-marked)

Large matrices (>200K rows) stored as metadata-only with download commands.

### Key Decisions

1. **Curated lists from papers** (not broad queries): Provides gold-standard difficulty
   classification from actual APTP benchmark results.
2. **Tiered storage (no LFS)**: Hand-constructed + CI subset committed to plain git;
   full SuiteSparse collection gitignored and extracted from archive at container build.
   Originally planned Git LFS, replaced with three-tier approach to avoid LFS costs
   and complexity for 4GB of immutable reference data.
3. **CI subset (10 matrices, ~73MB)**: Representative set in `suitesparse-ci/` committed
   to plain git (~19MB in pack). Avoids 40KB/s SuiteSparse API throttle in CI. Covers
   all categories including 2 killer cases.
4. **Full collection via archive**: `references/ssids/suitesparse.tar.gz` (1.3GB) stored
   outside git as static reference data. Docker entrypoint extracts automatically.
5. **Python scripts with uv**: Virtual environment at `scripts/.venv/`, dependencies via `uv pip`.
6. **ssgetpy API**: `ssgetpy.search(name)` with exact group/name matching; download with
   manual tar.gz extraction for reliable file placement.
7. **vibrobox group correction**: SuiteSparse has `Cote/vibrobox`, not `GHS_indef/vibrobox`.
8. **cont-201 killer case**: Added `killer_case: True` to `cont-201` (large optimization problem
   with many zero diagonals forcing 2×2 pivots) to reach the SC-005 threshold of ≥5 killer cases.
   Total killer cases: stokes128, cont-201, ncvxqp3, c-big, c-71.

### Issues Encountered

- **ssgetpy API change**: The `matname` keyword argument no longer works; must use positional
  `name_or_id` parameter.
- **ssgetpy extract=True**: Doesn't always place the correct `.mtx` file; rewrote to use
  `extract=False` + manual `tarfile` extraction with name-based file selection.
- **Slow download speeds**: SuiteSparse downloads throttled to ~40KB/s from this environment.

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
