# Implementation Plan: SPRAL Golden Results Generation

**Branch**: `003-spral-golden-results` | **Date**: 2026-02-06 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/003-spral-golden-results/spec.md`

## Summary

Build SPRAL from source (CPU-only, Meson) in the development container, write a C driver program that runs the SSIDS analyze-factorize-solve pipeline on each of the 82 test matrices from Phase 0.2, and capture comprehensive golden reference results (inertia, pivot statistics, errors, timing) as JSON files. Results follow the three-tier storage strategy: hand-constructed and CI subset committed to git, full SuiteSparse set gitignored.

## Technical Context

**Language/Version**: C (driver program), Shell (orchestration scripts); building against SPRAL Fortran library
**Primary Dependencies**: SPRAL (BSD-3, built from source at `/opt/references/spral/`), libopenblas, liblapack, libmetis, libhwloc
**Storage**: JSON files in `test-data/spral-reference/` (three-tier: git, git-CI, gitignored)
**Testing**: Shell-based validation (jq for JSON, diff for reproducibility)
**Target Platform**: Linux arm64 (development container, Debian 12)
**Project Type**: Test infrastructure (not part of the Rust solver)
**Performance Goals**: N/A (correctness-focused; timing data is captured but not optimized)
**Constraints**: No Python in container (validation via jq/shell). No GPU. Container must have gfortran, gcc already installed.
**Scale/Scope**: 82 matrices, ~82 JSON result files, 4 shell scripts, 1 C driver program (~500 LOC)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Golden results provide the correctness oracle for all future phases. The driver validates its own output (checks SPRAL return codes, computes backward error). |
| II. Clean Room Implementation | PASS | Driver links against SPRAL (BSD-3) via its public C API. No restricted code consulted. SPRAL source is freely consultable. |
| III. Test-Driven Development | PASS | This feature *creates* the test data that TDD depends on. The driver itself is validated by cross-checking hand-constructed matrices with known factorizations. |
| IV. Algorithm Documentation | PASS | Research.md documents all decisions. Driver source will cite SPRAL API docs and examples. |
| V. Numerical Stability | PASS | Uses SPRAL's own APTP implementation. Backward error computation follows SPRAL's scaled residual formula. |
| VI. Structured Development | PASS | This is Phase 0.3, executed after Phase 0.2 completion. Exit criteria from ssids-plan.md are addressed. |
| VII. Code Quality | N/A | Driver is C infrastructure code, not Rust solver code. Standard C quality applies. |

**No violations. All gates pass.**

## Project Structure

### Documentation (this feature)

```text
specs/003-spral-golden-results/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0 research findings
├── data-model.md        # Entity definitions and schemas
├── quickstart.md        # Getting started guide
├── contracts/
│   └── driver-cli.md    # CLI interface contracts
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # (Phase 2 - created by /speckit.tasks)
```

### Source Code (repository root)

```text
test-data/
├── scripts/
│   ├── spral_driver.c           # C driver program (new)
│   ├── mmio.c                   # Matrix Market I/O (new, public domain)
│   ├── mmio.h                   # Matrix Market I/O header (new, public domain)
│   ├── run-spral.sh             # Orchestration script (new)
│   ├── generate-summary.sh      # Summary generator (new)
│   └── verify-reproducibility.sh # Reproducibility checker (new)
├── spral-reference/
│   ├── config.json              # SPRAL build configuration record (new)
│   ├── summary.json             # Aggregated summary report (new)
│   ├── hand-constructed/        # Tier 1: 15 result files (new, git-tracked)
│   ├── suitesparse-ci/          # Tier 2: 10 result files (new, git-tracked)
│   └── suitesparse/             # Tier 3: 67 result files (new, gitignored)
├── hand-constructed/            # (existing from Phase 0.2)
├── suitesparse-ci/              # (existing from Phase 0.2)
├── suitesparse/                 # (existing from Phase 0.2, gitignored)
└── metadata.json                # (existing from Phase 0.2)

docker/
└── Dockerfile                   # Updated to install SPRAL build deps and build SPRAL
```

**Structure Decision**: All new code lives under `test-data/scripts/` alongside existing Phase 0.2 scripts. Results go in `test-data/spral-reference/` parallel to the matrix storage directories. This keeps test infrastructure co-located with test data.

## Complexity Tracking

No violations to justify. All constitution gates pass without exception.

## Implementation Phases

### Phase A: SPRAL Build Infrastructure (Dockerfile Build Stage)

**Goal**: Get SPRAL building and passing its own tests as a Dockerfile build stage, so the library is pre-installed in the container image.

**Steps**:
1. Add a new build stage to `docker/Dockerfile` that installs build dependencies (python3, meson, ninja, libopenblas-dev, libmetis-dev, libhwloc-dev)
2. In that stage, run meson setup (CPU-only, release, OpenMP), compile, and install SPRAL to `/usr/local`
3. Copy the installed SPRAL artifacts (libs, headers) into the final development stage
4. Verify SPRAL's own `meson test` passes during the build
5. Record the SPRAL build configuration (commit hash, compiler version, options) in `config.json`

**Key Decisions**:
- CPU-only build (`-Dgpu=false`) — no CUDA dependency
- OpenMP enabled for realistic threading behavior
- Release build for representative timing data
- METIS ordering (default, `ordering=1`)
- No scaling (`scaling=0`) for simplest reproducible baseline

### Phase B: Driver Program

**Goal**: C program that runs SSIDS on a single .mtx file and outputs JSON.

**Steps**:
1. Obtain/write `mmio.c`/`mmio.h` for Matrix Market reading (public domain NIST implementation or minimal custom reader)
2. Write `spral_driver.c`:
   - Read Matrix Market file (symmetric coordinate format)
   - Convert to CSC lower triangle with 1-based indexing (SPRAL's expected format)
   - Initialize SPRAL options (default + command-line overrides)
   - Run analyze phase, capture inform fields and timing
   - Run factorize phase, capture inform fields and timing
   - Generate RHS as `b = A * ones(n)` (known solution `x = ones(n)`)
   - Run solve phase, capture timing
   - Compute forward error: `max |x_i - 1.0|`
   - Compute backward error: `||b - Ax|| / (||A|| ||x|| + ||b||)` (infinity norm, scaled)
   - Compute inertia: positive = rank - num_neg, negative = num_neg, zero = n - rank
   - Output JSON result to stdout or `--output` file
   - Handle failures at any phase (write partial result with status)
3. Compile and test on `arrow-5-pd.mtx` (known positive definite, trivial case)
4. Test on `arrow-10-indef.mtx` (known indefinite case with 2x2 pivots expected)
5. Test on `singular-3.mtx` (expected singularity detection)

**Key Design Choices**:
- C rather than Fortran for easier JSON output and string handling
- Uses SPRAL's C API (32-bit pointer variant `spral_ssids_analyse_ptr32` if Matrix Market uses 32-bit indices, or 64-bit `spral_ssids_analyse` for large matrices)
- Environment variables `OMP_CANCELLATION=TRUE` and `OMP_PROC_BIND=TRUE` set by run scripts

### Phase C: Orchestration and Batch Execution

**Goal**: Shell scripts to run the driver across all matrices and produce summary.

**Steps**:
1. Write `run-spral.sh`:
   - Accepts tier argument (hand-constructed, suitesparse-ci, suitesparse)
   - Reads `metadata.json` to get list of matrices for the tier
   - Runs `spral_driver` on each matrix, writing results to `spral-reference/<tier>/`
   - Reports progress and captures any driver crashes
2. Write `generate-summary.sh`:
   - Reads all result JSON files
   - Cross-references with `metadata.json` for matrix dimensions
   - Classifies difficulty based on backward error and pivot statistics
   - Outputs `summary.json`
3. Write `verify-reproducibility.sh`:
   - Runs driver twice per matrix (on specified tier)
   - Compares structural fields (exact match)
   - Compares numerical fields (within 1e-14 tolerance)
   - Reports pass/fail per matrix

### Phase D: Execution and Validation

**Goal**: Generate all golden results, validate, and commit.

**Steps**:
1. Run on hand-constructed matrices (15) — all should succeed, cross-check inertia with known factorizations from `.json` files
2. Run on CI subset (10) — expect all to succeed
3. Run on full SuiteSparse set (67) — expect ≥95% success
4. Verify reproducibility on CI subset
5. Generate summary report
6. Update `.gitignore` for `test-data/spral-reference/suitesparse/`
7. Commit Tier 1 and Tier 2 results to git
8. Document any failures or surprising results in development log

### Phase E: Documentation and Cleanup

**Goal**: Update project documentation to reflect Phase 0.3 completion.

**Steps**:
1. Update `docs/ssids-plan.md` success criteria checkboxes for Phase 0.3
2. Add Phase 0.3 entry to `docs/ssids-log.md`
3. Update `CLAUDE.md` with new test infrastructure details
4. Update `metadata.json` `reference_results` fields for matrices that now have SPRAL results

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| SPRAL fails to build on arm64 | Low | High | SPRAL CI tests on multiple platforms; fallback to cross-compilation or x86 container |
| MeTiS not available for arm64 | Low | Medium | Can build without MeTiS (uses inferior ordering); or build MeTiS from source |
| Some SuiteSparse matrices fail | Medium | Low | Expected — singular/extreme matrices may fail; FR-010 requires capturing failures |
| Large matrices exceed container memory | Medium | Low | Record OOM as known limitation; most matrices are moderate-sized |
| Reproducibility fails due to OpenMP non-determinism | Low | Medium | Pin thread count to 1 if needed; structural results are always deterministic |
