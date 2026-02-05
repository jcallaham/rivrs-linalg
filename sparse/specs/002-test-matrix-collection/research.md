# Research: Test Matrix Collection Assembly

**Date**: 2026-02-05
**Feature**: `002-test-matrix-collection`

## Decision 1: SuiteSparse Matrix Retrieval Tool

**Decision**: Use `ssgetpy` (Python) to query and download matrices from the SuiteSparse Matrix Collection in Matrix Market format.

**Rationale**: ssgetpy provides a programmatic interface to the SuiteSparse collection with filtering by properties (`dtype`, `isspd`, `rowbounds`, `nzbounds`, `kind`). It supports Matrix Market download format (`format='MM'`). The user has already successfully used it to find 412 real symmetric indefinite matrices.

**Alternatives considered**:
- Manual download from sparse.tamu.edu: Not scriptable, not reproducible
- suitesparse-matrix-collection-website API: Lower-level, ssgetpy wraps it conveniently
- Direct HTTP download of .tar.gz: Possible but ssgetpy handles caching and extraction

**Key finding — ssgetpy limitations**:
- No explicit "symmetric" boolean filter. Must use `isspd=False` + post-filter `nsym >= 0.999` for symmetric indefinite matrices.
- The `kind` parameter uses SQL LIKE substring matching (e.g., `kind="optimization"` matches "optimization problem").
- `limit` parameter defaults to 10 — must be set explicitly for larger result sets.

## Decision 2: Matrix Selection Strategy — Curated Lists Over Broad Queries

**Decision**: Use curated matrix lists from the APTP/SSIDS benchmark papers as the primary selection, supplemented by targeted discovery queries for additional coverage.

**Rationale**: The 412 results from a broad "real, symmetric, indefinite" query include many matrices unsuitable for SSIDS testing (graph problems, eigenvalue problems, least squares). The APTP papers (Duff, Hogg, Lopez 2020; Hogg, Ovtchinnikov, Scott 2016) provide gold-standard test sets with known difficulty classifications based on actual delayed pivot counts.

**Alternatives considered**:
- Broad automated query with automated classification: Would require running a solver to classify difficulty — that's Phase 0.3, not 0.2.
- Only paper-referenced matrices: Too narrow, would miss important edge cases and coverage.

**Matrix selection by category**:

### Easy Indefinite (target: ~30)

From Duff et al. (2020) Table 1 and Hogg et al. (2016) Table III:
- Oberwolfach/t2dal (4.3k), GHS_indef/dixmaanl (60k), Oberwolfach/rail_79841 (79.8k)
- GHS_indef/dawson5 (51.5k), Boeing/bcsstk39 (46.8k), Boeing/pct20stif (52.3k)
- GHS_indef/copter2 (55.5k), GHS_indef/helm2d03 (392.3k), Boeing/crystk03 (24.7k)
- Oberwolfach/filter3D (106.4k), Koutsovasilis/F2 (71.5k), McRae/ecology1 (1000k)
- Cunningham/qa8fk (66.1k), Oberwolfach/gas_sensor (66.9k), Oberwolfach/t3dh (79.2k)
- Lin/Lin (256k), PARSEC/H2O (67k), GHS_indef/sparsine (50k)
- From Hogg 2016: Newman/astro-ph (16.7k), Andrianov/mip1 (66.5k), PARSEC/SiNa (5.7k)
- Schmid/thermal2 (1230k), INPRO/msdoor (416k), Boeing/pwtk (218k)
- BenElechi/BenElechi1 (246k), PARSEC/Si10H16 (17.1k), PARSEC/Si5H12 (19.9k)
- From Duff 2005: GHS_indef/helm3d01 (32.2k), GHS_indef/vibrobox (12.3k)

**Characteristics**: Well-conditioned after AMD ordering, few delayed pivots (<1% of columns), no matching-based prescaling needed. Domains: structural mechanics, model reduction, thermal.

### Hard Indefinite (target: ~15)

From Duff et al. (2020) Table 2 and Duff (2005):
- TSOPF/TSOPF_FS_b39_c7 (28.2k), TSOPF/TSOPF_FS_b162_c1 (10.8k)
- GHS_indef/cont-201 (80.6k), GHS_indef/stokes128 (49.7k)
- GHS_indef/ncvxqp1 (12.1k), GHS_indef/darcy003 (389.9k)
- GHS_indef/cont-300 (180.9k), GHS_indef/bratu3d (27.8k)
- GHS_indef/cvxqp3 (17.5k), GHS_indef/d_pretok (182.7k)
- GHS_indef/turon_m (189.9k), GHS_indef/ncvxqp5 (62.5k)
- GHS_indef/ncvxqp3 (75k), GHS_indef/ncvxqp7 (87.5k)
- Schenk_IBMNA/c-big (345.2k)
- From Hogg 2016 APTP-sensitive: GHS_indef/c-71 (76.6k), ND/nd6k (18k), ND/nd12k (36k)

**Characteristics**: Augmented/saddle-point structure, large off-diagonal entries, many delayed pivots. Domains: optimization (KKT systems), power systems (TSOPF), nonconvex QP, Stokes flow.

### Positive Definite (target: ~20)

From Hogg et al. (2016) Table III (dagger-marked):
- GHS_psdef/crankseg_1 (52.8k), Rothberg/cfd2 (123k), DNVS/thread (29.7k)
- DNVS/shipsec1 (141k), DNVS/shipsec8 (115k), Oberwolfach/boneS01 (127k)
- GHS_psdef/crankseg_2 (63.8k), Schenk_AFE/af_shell7 (505k)
- DNVS/shipsec5 (180k), AMD/G3_circuit (1590k), GHS_psdef/bmwcra_1 (149k)
- Schenk_AFE/af_0_k101 (504k), GHS_psdef/ldoor (952k), DNVS/ship_003 (122k)
- Um/offshore (260k), ND/nd6k (18k), GHS_psdef/inline_1 (504k)
- GHS_psdef/apache2 (715k), ND/nd12k (36k)

### "Killer" Cases (matrices that defeat static pivoting)

From Hogg 2016 (high delayed pivots) and Duff 2005:
- GHS_indef/ncvxqp3: Nonconvex QP, APTP-sensitive
- GHS_indef/c-71: Nonlinear optimization, APTP-sensitive
- ND/nd6k and ND/nd12k: 3D mesh, APTP-sensitive
- GHS_indef/stokes128: Stokes flow, saddle-point structure
- Schenk_IBMNA/c-big: Nonlinear optimization, large

## Decision 3: SPRAL Test Matrices — Procedural Generation, Not Files

**Decision**: Port SPRAL's hand-constructed matrices and random generation patterns into the test suite as Rust code, rather than extracting static matrix files.

**Rationale**: SPRAL contains no pre-built .mtx files. All test matrices are generated procedurally at runtime via:
1. **Hand-constructed small matrices** (4x4 SPD, 3x3 zero-diagonal, 3x3 singular) in `tests/ssids/ssids.f90`
2. **Structured generators**: `gen_bordered_block_diag()` for block diagonal with dense border
3. **Random generators**: `gen_random_posdef()`, `gen_random_indef()`, `gen_random_indef_fred()` (with scattered zero diagonals)
4. Test coverage: 100+ random problems per suite, sizes 1 to 2000, 50/50 SPD/indefinite split

**What to extract for our hand-constructed set**:
- 4x4 SPD matrix (`simple_mat_lower`)
- 3x3 zero-diagonal matrix (`simple_mat_zero_diag`)
- 3x3 singular matrix (`simple_sing_mat`)
- 3x3 null-column matrix (`simple_sing_mat2`)
- Bordered block diagonal pattern (configurable sizes and border width)
- Random indefinite with scattered zero diagonals pattern

**Alternatives considered**:
- Run SPRAL to produce .mtx files: Adds SPRAL build dependency; SPRAL doesn't naturally output .mtx format.
- Skip SPRAL patterns entirely: Loses well-tested structural patterns that exercise specific solver behaviors.

## Decision 4: Hand-Constructed Matrix Design

**Decision**: Build 15+ hand-constructed matrices (5x5 to 20x20) covering 5 structural categories, each with analytically computed exact LDL^T factorizations.

**Categories**:

1. **Arrow matrices** (known supernode structure):
   - Dense first row/column, diagonal otherwise
   - Sizes: 5x5, 10x10, 15x15
   - Both positive definite and indefinite variants

2. **Block diagonal** (independent subsystems):
   - 2-4 dense blocks of varying sizes
   - With and without dense border (from SPRAL pattern)
   - Tests solver's ability to detect independent subproblems

3. **Tridiagonal** (simple elimination tree):
   - Sizes: 5x5, 10x10, 20x20
   - Known exact LDL^T factorization computable by recurrence
   - Both positive definite and indefinite

4. **Stress-test structures**:
   - Matrix requiring maximum delayed pivots (zero diagonal, large off-diagonal)
   - Matrix with worst-case fill-in (dense column pattern)
   - Ill-conditioned matrix (large condition number, known exact solution)

5. **Degenerate cases** (from SPRAL):
   - Zero diagonal matrix
   - Singular matrix (zero rows/columns)
   - 1x1 and 2x2 trivial cases

**Factorization documentation**: Each matrix will include a companion `.json` file with the exact L, D, and permutation P such that P^T A P = L D L^T, computed symbolically or by hand.

## Decision 5: Git LFS for Matrix Files

**Decision**: Set up Git LFS to track `.mtx` files before adding any matrix data to the repository.

**Rationale**: SuiteSparse matrices range from a few KB (small problems) to hundreds of MB (large problems like G3_circuit with 1.59M rows). Storing these directly in Git would bloat the repository history. The root `.gitattributes` already tracks `*.pdf` via LFS.

**Current state**:
- Root `.gitattributes` exists with `*.pdf filter=lfs diff=lfs merge=lfs -text`
- Git LFS is NOT installed in this Docker container
- No `.mtx` files exist yet
- No `test-data/` directory exists yet

**Plan**:
1. Install `git-lfs` in the Docker container (or document as a setup prerequisite)
2. Add `*.mtx filter=lfs diff=lfs merge=lfs -text` to root `.gitattributes`
3. Run `git lfs install` in the repo
4. Commit `.gitattributes` changes before adding any `.mtx` files

**Size management strategy**:
- Hand-constructed matrices (5x5–20x20): Tiny, a few KB each — LFS is overkill but harmless
- Small/medium SuiteSparse matrices (up to ~100K rows): Typically 1-50 MB — LFS is appropriate
- Large SuiteSparse matrices (>100K rows): Can be 100MB+ — store metadata only with download script, not in repo
- Threshold: Matrices with .mtx files > 50MB should be download-on-demand via the Python script

**Alternatives considered**:
- No LFS, just commit directly: Unacceptable for 70+ matrices
- External hosting only (S3, etc.): Adds infrastructure dependency, harder for contributors
- Compressed .mtx.gz: Matrix Market is already reasonably compact; adds complexity for marginal savings

## Decision 6: Directory Structure

**Decision**: Use the directory structure from ssids-plan.md Phase 0.2 with minor adjustments.

```
sparse/test-data/
├── metadata.json                # Complete index of all matrices
├── hand-constructed/
│   ├── arrow-5-pd.mtx          # 5x5 arrow, positive definite
│   ├── arrow-5-pd.json         # Exact factorization for arrow-5-pd
│   ├── arrow-10-indef.mtx
│   ├── arrow-10-indef.json
│   ├── tridiag-10-pd.mtx
│   ├── tridiag-10-pd.json
│   ├── tridiag-10-indef.mtx
│   ├── tridiag-10-indef.json
│   ├── block-diag-15.mtx
│   ├── block-diag-15.json
│   ├── bordered-block-20.mtx
│   ├── bordered-block-20.json
│   ├── stress-delayed-pivots.mtx
│   ├── stress-delayed-pivots.json
│   ├── stress-fill-in.mtx
│   ├── stress-fill-in.json
│   ├── stress-ill-conditioned.mtx
│   ├── stress-ill-conditioned.json
│   ├── zero-diagonal-3.mtx
│   ├── zero-diagonal-3.json
│   ├── singular-3.mtx
│   ├── singular-3.json
│   ├── trivial-1x1.mtx
│   ├── trivial-1x1.json
│   └── trivial-2x2.mtx
├── suitesparse/
│   ├── easy-indefinite/         # ~30 matrices
│   ├── hard-indefinite/         # ~18 matrices
│   └── positive-definite/       # ~19 matrices
├── spral-inspired/              # Matrices adapted from SPRAL test patterns
│   └── (generated at test time, not stored)
└── scripts/
    ├── download_suitesparse.py  # ssgetpy download script
    └── generate_hand_constructed.py  # Matrix generation script
```

**Rationale**: Separates hand-constructed (committed) from SuiteSparse (LFS-tracked) from generated (runtime). Scripts directory keeps tooling co-located.

## Decision 7: Metadata Schema

**Decision**: Use a single `metadata.json` at the test-data root, following the schema from ssids-plan.md with extensions.

```json
{
  "matrices": [
    {
      "name": "arrow-10-indef",
      "source": "hand-constructed",
      "category": "hand-constructed",
      "path": "hand-constructed/arrow-10-indef.mtx",
      "size": 10,
      "nnz": 28,
      "properties": {
        "symmetric": true,
        "positive_definite": false,
        "indefinite": true,
        "structure": "arrow",
        "difficulty": "trivial"
      },
      "factorization_path": "hand-constructed/arrow-10-indef.json",
      "paper_references": []
    },
    {
      "name": "GHS_indef/ncvxqp3",
      "source": "suitesparse",
      "category": "hard-indefinite",
      "path": "suitesparse/hard-indefinite/ncvxqp3/ncvxqp3.mtx",
      "size": 75000,
      "nnz": 275000,
      "properties": {
        "symmetric": true,
        "positive_definite": false,
        "indefinite": true,
        "structure": "augmented",
        "difficulty": "hard",
        "kind": "optimization problem",
        "expected_delayed_pivots": "high",
        "killer_case": true
      },
      "paper_references": [
        "Hogg, Ovtchinnikov, Scott (2016) Table III #7",
        "Duff, Hogg, Lopez (2020)"
      ],
      "suitesparse_id": null,
      "reference_results": {}
    }
  ],
  "schema_version": "1.0",
  "generated": "2026-02-05",
  "total_count": 0
}
```

`reference_results` is intentionally empty — that will be populated in Phase 0.3 (SPRAL golden results).

## Decision 8: Interior Point Test Problems

**Decision**: Include a small set of KKT/augmented system matrices from interior point methods, sourced from the SuiteSparse collection.

**Rationale**: Interior point methods are a primary consumer of sparse indefinite solvers. The APTP algorithm was specifically designed for these saddle-point systems. Several are already in the hard-indefinite set (ncvxqp*, cont-201, cont-300, cvxqp3).

**Additional candidates from Duff (2005)**:
- GHS_indef/aug3dcqp (35.5k): Expanded system from 3D PDE
- GHS_indef/blockqp1 (60k): QP block structure
- GHS_indef/mario001 (38.4k): Stokes equation
- GHS_indef/qpband (20k): Banded QP

These will be classified under `hard-indefinite` with a `kind: "optimization problem"` or `kind: "fluid mechanics"` tag in metadata.

## Decision 9: Size Cutoff for In-Repo Storage

**Decision**: Store matrices up to 50MB in the repository via LFS. For matrices exceeding 50MB, store only metadata with a download URL/command.

**Rationale**: Several matrices in the curated lists are very large (G3_circuit: 1.59M rows, ecology1: 1M rows). Their .mtx files can exceed 100MB. Including all of them via LFS would be costly for cloning and CI. The download script already supports selective downloading.

**Large matrices (metadata-only, download-on-demand)**:
- AMD/G3_circuit (1590k rows)
- McRae/ecology1 (1000k rows)
- GHS_psdef/ldoor (952k rows)
- GHS_psdef/apache2 (715k rows)
- Schmid/thermal2 (1230k rows)
- Schenk_AFE/af_shell7 (505k rows)
- Schenk_AFE/af_0_k101 (504k rows)
- GHS_psdef/inline_1 (504k rows)
- GHS_indef/helm2d03 (392k rows)
- GHS_indef/darcy003 (390k rows)

These will have `"in_repo": false` in metadata and a `"download_command"` field.
