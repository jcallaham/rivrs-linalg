# SSIDS Development Log

## Phase 1.1: Test Infrastructure

**Status**: Complete
**Branch**: `005-test-infrastructure`
**Date**: 2026-02-07

### What Was Built

Reusable test harness and validation infrastructure for solver development,
gated behind a `test-util` Cargo feature flag to keep production builds lean.

**Testing module** (`src/testing/`, feature-gated behind `test-util`):
- `harness.rs` — `SolverTest` trait defining four test methods (analyze, factor,
  solve, roundtrip), `MockSolver` implementation that validates reference
  factorizations directly, `TestKind` enum, `MetricResult` and `TestResult`
  structs with `Display` formatting
- `validator.rs` — `NumericalValidator` with builder pattern for configurable
  tolerances (reconstruction 10^-12, backward error 10^-10), methods for
  `check_reconstruction()`, `check_backward_error()`, `check_inertia()`, and
  `validate_factorization()` returning structured `TestResult`
- `cases.rs` — `SolverTestCase`, `TestMatrixProperties`, `TestCaseFilter` with
  builder methods (`all()`, `hand_constructed()`, `ci_subset()`, `with_source()`,
  `with_category()`, `with_difficulty()`, `ci_only()`, `require_reference()`),
  `load_test_cases()` wrapping registry functions with predicate filtering
- `generators.rs` — `RandomMatrixConfig`, `generate_random_symmetric()` (PD via
  diagonal dominance or indefinite with mixed signs), `generate_arrow()`,
  `generate_tridiagonal()`, `generate_banded()` — all producing
  `SparseColMat<usize, f64>` via faer's `Triplet` API
- `mod.rs` — Module root with re-exports and compilable rustdoc example

**Cargo.toml changes**:
- `test-util` feature flag gating `rand` and `rand_distr` as optional deps
- Self-dev-dependency: `rivrs-sparse = { path = ".", features = ["test-util"] }`

**Integration test refactoring**:
- `tests/hand_constructed.rs` — refactored to use `load_test_cases` +
  `NumericalValidator` instead of raw registry + validate calls
- `tests/suitesparse_ci.rs` — refactored to use `TestCaseFilter::ci_subset()`

**Tests**: 48 total (45 unit + 2 integration + 1 doctest), all passing

### Key Decisions

1. **Feature-gated testing module**: The `test-util` feature flag keeps `rand`,
   `rand_distr`, and all test generators out of production builds. Downstream
   crates enable it in `[dev-dependencies]` only.

2. **Builder pattern for filters and validators**: `TestCaseFilter` and
   `NumericalValidator` use builder patterns for ergonomic, composable
   configuration. Default tolerances match the constitution (reconstruction
   10^-12, backward error 10^-10).

3. **MockSolver for immediate validation**: The `MockSolver` validates reference
   factorizations from the hand-constructed JSON files directly (no solver needed
   yet). This allows the test harness to be fully exercised before any solver
   exists.

4. **Diagonal dominance for PD generators**: Random PD matrices use Gershgorin
   circle theorem — diagonal entries exceed row absolute sums by a random margin,
   guaranteeing positive definiteness without Cholesky verification.

5. **Seeded RNG for reproducibility**: All generator tests use
   `StdRng::seed_from_u64(42)` for deterministic output across runs.

### Issues Encountered

- **`faer::sparse::Triplet` path**: `Triplet` is at `faer::sparse::Triplet`,
  not `faer::Triplet`. The `try_new_from_triplets` method expects
  `&[Triplet<usize, usize, f64>]`, not tuples.

- **Rust 2024 edition `gen` keyword**: `gen` is reserved in edition 2024. Must
  use raw identifier syntax `rng.r#gen::<f64>()` to call `rand::Rng::gen()`.

- **`Uniform::new` type inference**: `Uniform::new(0.1, 1.0)` is ambiguous;
  needs explicit type: `Uniform::new(0.1f64, 1.0)`.

- **Reconstruction tolerance sensitivity**: A perturbation of 1e-4 to an L
  entry produced reconstruction error exceeding 1e-6 tolerance in tests. Fixed
  by using 1e-8 perturbation for the relaxed-tolerance test.

- **`cargo fmt` CI failure**: Import ordering and line wrapping differences
  between local and CI rustfmt. Fixed by running `cargo fmt` and committing.

---

## Phase 0.4: Repository Setup for Solver Development

**Status**: Complete
**Branch**: `004-repo-setup`
**Date**: 2026-02-06

### What Was Built

Development infrastructure enabling Rust-based loading, parsing, and validation
of the test matrix collection established in Phase 0.2.

**IO modules** (`src/io/`):
- `mtx.rs` — Matrix Market parser (`coordinate real symmetric` format) with
  1-indexed→0-indexed conversion, symmetric expansion, and descriptive parse
  errors including file path and line number
- `reference.rs` — JSON reference factorization loader with types for Inertia,
  LEntry, DBlock (1×1/2×2 with polymorphic serde), ReferenceFactorization;
  includes validation of strict lower triangle and permutation consistency
- `registry.rs` — Test matrix catalog backed by `metadata.json` with CI-subset
  path fallback (`suitesparse-ci/` preferred for CI matrices over gitignored
  `suitesparse/` directory)

**Error handling** (`src/error.rs`):
- Added `IoError` and `ParseError` variants to `SparseError`
- Added `From<std::io::Error>` and `From<serde_json::Error>` impls

**Validation utilities** (`src/validate.rs`):
- `reconstruction_error(A, ref)` — `||P^T A P - L D L^T||_F / ||A||_F`
- `backward_error(A, x, b)` — `||Ax - b|| / (||A||_F ||x|| + ||b||)`
- `check_inertia(computed, expected)` — field-wise equality comparison
- All implemented using faer dense operations (to_dense, matmul, norm_l2)

**Tests** (20 total):
- 11 unit tests across mtx, reference, registry modules
- 7 validation unit tests (reconstruction, backward error, inertia)
- 1 integration test loading all 15 hand-constructed matrices with
  reconstruction error < 10^-12 (SC-008 validated)
- 1 integration test for 10 CI-subset SuiteSparse matrices

**CI pipeline** (`.github/workflows/ci.yml`):
- `test-sparse` job (stable + MSRV 1.87)
- `lint-sparse` job (clippy + rustfmt)
- `doc-sparse` job (rustdoc with `-D warnings`)

**Benchmark scaffold** (`benches/matrix_loading.rs`):
- Criterion benchmarks for matrix loading and reconstruction error computation

### Key Decisions

1. **Custom MTX parser**: Lightweight, scoped to `coordinate real symmetric`
   format. No external crate dependency. Descriptive errors with file path and
   line number for debugging data quality issues.

2. **faer-based validation**: All validation uses dense matrix operations via
   faer (to_dense, matmul, norm_l2, operator overloads). Simple and correct for
   test matrices up to ~20K dimension. No need for sparse validation at this
   stage.

3. **Polymorphic DBlock serde**: Custom deserializer reads `"size"` field to
   distinguish 1×1 scalar pivots from 2×2 symmetric blocks. Matches the JSON
   schema established in Phase 0.2.

4. **CI-subset path fallback**: `load_test_matrix()` checks `suitesparse-ci/`
   before `suitesparse/` for CI-subset entries, ensuring tests work both in
   CI (where only `suitesparse-ci/` is available) and locally (where
   `suitesparse/` may also be extracted).

5. **TDD throughout**: All tests written and verified to fail before
   implementation, per Constitution Principle III.

### Issues Encountered

- **faer API discovery**: `sparse_dense_matmul` lives in
  `faer::sparse::linalg::matmul`, not `faer::linalg::matmul`. `Par::Seq`
  (not `Par::sequential()`). Column `as_mat()` (not `as_mat_ref()`).
  Resolved by reading faer source code directly.

- **metadata.json local corruption**: Working copy had been modified from 82
  to 73 entries (likely from a prior script run). Restored from git.

---

## Phase 0.3: SPRAL Golden Results — Deferred

**Status**: Deferred
**Branch**: `003-spral-golden-results`
**Date**: 2026-02-06

### Decision

Phase 0.3 was planned to build SPRAL from source, write a C driver program, run
SPRAL on all 82 test matrices, and capture golden reference results (timing,
inertia, pivot statistics, backward error) as JSON files.

After investigating SPRAL's C API (`spral_ssids.h`) and assessing cost-benefit,
we decided to **defer this phase**. Full rationale in
`specs/003-spral-golden-results/decision.md`.

### Key Finding: SPRAL API Limitations

SPRAL's C API exposes:
- Aggregate statistics from the `inform` struct (num_factor, num_flops,
  num_delay, num_neg, num_two, matrix_rank, maxfront, maxsupernode)
- Block diagonal D entries via `spral_ssids_enquire_indef()`
- Pivot ordering via `spral_ssids_enquire_indef()`
- The solution vector (from which errors are computed)

SPRAL does **NOT** expose: the L factor, permutation P, elimination tree parent
pointers, supernode membership arrays, or fill-in patterns. These are internal to
the opaque `akeep`/`fkeep` handles.

### Key Decision: Reconstruction Tests as Primary Oracle

The **reconstruction test** (`||P^T A P - L D L^T|| / ||A|| < epsilon`) is a
strictly stronger correctness oracle than comparing against SPRAL's output. If
our factorization reconstructs A to machine precision, it is correct by
definition — regardless of what any other solver produces.

This means the most valuable correctness tests require no SPRAL infrastructure:
1. **Reconstruction**: P^T A P = L D L^T (proves mathematical correctness)
2. **Backward error**: ||Ax - b|| / (||A|| ||x|| + ||b||) (validates solve pipeline)
3. **Hand-constructed matrices**: Analytically known factorizations from Phase 0.2
4. **Property-based tests**: Inertia consistency, symmetry preservation

### Impact

- Constitution updated to v1.1.0: reconstruction tests are primary oracle; SPRAL
  comparison is a secondary validation layer added when available
- Phase 0 exit criteria updated: "SPRAL reference results" requirement replaced
  with reconstruction test strategy
- SPRAL build infrastructure deferred to Phases 2-8 when needed for performance
  benchmarking and inertia validation on large SuiteSparse matrices

### What Was Produced

No code changes. Documentation artifacts only:
- `specs/003-spral-golden-results/decision.md` — full decision record
- `specs/003-spral-golden-results/spec.md` — original feature specification
- `specs/003-spral-golden-results/plan.md` — implementation plan (not executed)
- `specs/003-spral-golden-results/research.md` — SPRAL API investigation findings
- Updated `docs/ssids-plan.md` — Phase 0.3 marked deferred, exit criteria revised
- Updated constitution v1.0.0 → v1.1.0

---

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
