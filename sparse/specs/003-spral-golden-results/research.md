# Research: SPRAL Golden Results Generation

**Feature**: 003-spral-golden-results
**Date**: 2026-02-06

## R1: Building SPRAL from Source in the Container

**Decision**: Build SPRAL using Meson with CPU-only configuration (no GPU) in the existing Docker container.

**Rationale**:
- SPRAL uses Meson as its primary build system, well-maintained with CI across Linux/macOS/Windows
- GPU support is unnecessary (we only need reference results, not GPU performance data)
- Container already has `gfortran` (12.2.0) and `gcc` — the core compilers
- Missing: python3, meson, ninja, libopenblas-dev, liblapack-dev, libmetis-dev, libhwloc-dev

**Alternatives Considered**:
- Autotools build: Deprecated by SPRAL maintainers, less reliable
- Pre-built binaries: Not available for arm64 container architecture
- Building without MeTiS: Possible but would use inferior default ordering

**Required Installation Steps**:
```bash
apt-get update && apt-get install -y \
    python3 python3-pip ninja-build \
    libopenblas-dev liblapack-dev \
    libmetis-dev libhwloc-dev

pip3 install meson
```

**Build Configuration**:
```bash
cd /opt/references/spral
meson setup builddir --prefix=/usr/local \
    --buildtype=release \
    -Dgpu=false \
    -Dopenmp=true \
    -Dexamples=true \
    -Dtests=true

cd builddir && meson compile && meson install
```

**Runtime Environment**:
```bash
export OMP_CANCELLATION=TRUE
export OMP_PROC_BIND=TRUE
```

## R2: Driver Program Language and Approach

**Decision**: Write the driver program in C, linking against the SPRAL C API.

**Rationale**:
- SPRAL provides a clean, well-documented C API in `spral_ssids.h`
- The C example (`examples/C/ssids.c`) demonstrates the complete workflow
- C is simpler than Fortran for JSON output generation and string handling
- C compiles alongside the SPRAL Meson build (can be added as a custom target or compiled separately)
- The existing SPRAL Fortran driver (`driver/spral_ssids.F90`) reads Rutherford-Boeing format, not Matrix Market — we need Matrix Market support

**Alternatives Considered**:
- Fortran driver (adapting existing `spral_ssids.F90`): Would work but Fortran string/JSON handling is cumbersome; the existing driver uses Rutherford-Boeing format, not Matrix Market
- Rust driver with FFI: Adds unnecessary complexity; C is the natural bridge language for SPRAL
- Shell script wrapping the SPRAL driver: The existing driver outputs unstructured text, would require fragile parsing

## R3: Matrix Market File Reading

**Decision**: Use the standard `mmio.c`/`mmio.h` reference implementation from NIST for reading Matrix Market files, or write a minimal CSC reader in C.

**Rationale**:
- SPRAL does NOT include a Matrix Market reader — it uses Rutherford-Boeing format
- Our test matrices from Phase 0.2 are all in Matrix Market format (`.mtx` files)
- The Matrix Market I/O library (`mmio.c`) is public domain and widely used
- Only need to read symmetric matrices in coordinate format and convert to CSC lower triangle (SPRAL's expected input)

**Alternatives Considered**:
- Convert all matrices to Rutherford-Boeing: Would allow using SPRAL's built-in reader, but adds a conversion step and changes the canonical format from Phase 0.2
- Embed matrices in C source: Only viable for hand-constructed matrices, not 67 SuiteSparse matrices

## R4: Output Format for Golden Results

**Decision**: JSON output, one file per matrix, with a machine-readable schema.

**Rationale**:
- JSON is human-readable and easily parsed by Rust (serde_json), Python, jq, and shell scripts
- One file per matrix allows incremental generation and selective re-runs
- The existing metadata.json from Phase 0.2 already uses JSON — consistency with established patterns
- jq is available in the container for validation (no Python needed)

**Alternatives Considered**:
- CSV summary only: Insufficient detail (can't capture variable-length data like pivot orders)
- Binary format: Not human-inspectable, harder to debug
- SQLite database: Over-engineered for ~82 result files

**Schema** (per-matrix result file):
```json
{
  "matrix_name": "string",
  "matrix_path": "string",
  "spral_config": {
    "version": "string (git commit hash)",
    "build_type": "release",
    "compiler": "gfortran 12.2.0",
    "options": {
      "ordering": 1,
      "scaling": 0,
      "pivot_method": 1,
      "u": 0.01,
      "nemin": 32,
      "action": true,
      "posdef": false
    }
  },
  "status": "success | factorization_failed | analysis_failed | error",
  "error_info": {
    "flag": 0,
    "message": "string or null"
  },
  "analysis": {
    "time_seconds": 0.0,
    "predicted_num_factor": 0,
    "predicted_num_flops": 0,
    "nparts": 0
  },
  "factorization": {
    "time_seconds": 0.0,
    "num_factor": 0,
    "num_flops": 0,
    "num_delay": 0,
    "num_neg": 0,
    "num_two": 0,
    "matrix_rank": 0,
    "maxfront": 0,
    "maxsupernode": 0,
    "num_sup": 0,
    "maxdepth": 0
  },
  "solve": {
    "time_seconds": 0.0,
    "rhs_method": "b = A * ones(n)",
    "rhs_seed": null,
    "forward_error": 0.0,
    "backward_error": 0.0
  },
  "inertia": {
    "positive": 0,
    "negative": 0,
    "zero": 0
  },
  "matrix_diagnostics": {
    "duplicate_entries": 0,
    "missing_diagonal": 0,
    "out_of_range": 0
  }
}
```

## R5: SPRAL API Surface for Data Extraction

**Decision**: Use the C API functions and inform structure to extract all required data.

**Findings from SPRAL source investigation**:

**After `spral_ssids_analyse()`**:
- `inform.flag` — status (0=success, <0=error, >0=warning)
- `inform.num_factor` — predicted number of factor entries
- `inform.num_flops` — predicted number of flops
- `inform.maxfront` — maximum front size
- `inform.maxdepth` — maximum elimination tree depth
- `inform.num_sup` — number of supernodes
- `inform.matrix_dup` — duplicate entries detected
- `inform.matrix_missing_diag` — missing diagonal entries
- `inform.matrix_outrange` — out-of-range entries

**After `spral_ssids_factor()`**:
- `inform.num_factor` — actual factor entries
- `inform.num_flops` — actual flops
- `inform.num_delay` — delayed columns
- `inform.num_neg` — negative eigenvalue count
- `inform.num_two` — 2x2 pivots used
- `inform.matrix_rank` — numerical rank
- `inform.maxfront` — maximum front size (may update)
- `inform.maxsupernode` — maximum supernode size

**After `spral_ssids_solve1()`**:
- Solution vector x (compute errors from this)
- `inform.flag` — solve status

**Via `spral_ssids_enquire_indef()`**:
- `piv_order` — pivot permutation (optional, not stored in golden results to save space)

**Inertia Computation** (from driver source):
- `num_neg` = `inform.num_neg`
- `num_zero` = `n - inform.matrix_rank`
- `num_positive` = `inform.matrix_rank - inform.num_neg`

**Backward Error Computation** (from driver source):
- RHS generated as `b = A * ones(n)`, so exact solution is `x = ones(n)`
- Forward error: `max |x_computed - 1.0|`
- Backward error: scaled residual `||b - Ax|| / (||A|| ||x|| + ||b||)`

## R6: Three-Tier Storage Strategy for Results

**Decision**: Follow the same three-tier pattern established in Phase 0.2.

**Mapping**:
- **Tier 1 (git-tracked)**: `test-data/spral-reference/hand-constructed/` — 15 result files (~small JSON, <1KB each)
- **Tier 2 (git-tracked CI)**: `test-data/spral-reference/suitesparse-ci/` — 10 result files for CI subset
- **Tier 3 (gitignored)**: `test-data/spral-reference/suitesparse/` — full 67 result files (generated at container build or on-demand)

**Rationale**: Consistent with Phase 0.2 storage decisions. Golden result JSON files are small (a few KB each), so even the full set could fit in git, but keeping the three-tier pattern maintains consistency and keeps the repo lean.

## R7: Reproducibility Strategy

**Decision**: Use deterministic RHS generation (`b = A * ones(n)`) and fixed SPRAL options. Verify reproducibility by running twice on CI subset.

**Rationale**:
- `b = A * ones(n)` gives a known exact solution `x = ones(n)` without requiring random number generation
- This is exactly what SPRAL's own driver does (see `driver/spral_ssids.F90` lines 74-85)
- Fixed options (no scaling, default METIS ordering, default pivot threshold) ensure deterministic behavior
- SPRAL's analysis phase is deterministic given the same ordering
- SPRAL's factorization is deterministic on same hardware with same thread count (when `OMP_PROC_BIND=TRUE`)

**Verification**: Run driver twice on CI subset matrices, diff the result files. Structural fields must match exactly; numerical fields (forward/backward error) must agree within 1e-14.
