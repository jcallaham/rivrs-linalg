# Implementation Plan: Repository Setup for Solver Development

**Branch**: `004-repo-setup` | **Date**: 2026-02-06 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/004-repo-setup/spec.md`

## Summary

Establish the development infrastructure for the SSIDS sparse solver: Matrix Market parser, JSON factorization loader, test matrix registry, numerical validation utilities (reconstruction error, backward error, inertia), CI pipeline for sparse domain, and benchmark scaffold. This completes Phase 0.4 of the SSIDS development plan, making the project ready for Phase 1 (Symbolic Analysis).

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (sparse/dense LA), serde + serde_json (JSON parsing)
**Storage**: Filesystem (test-data/ directory with .mtx and .json files)
**Testing**: cargo test (integration tests), criterion 0.5 (benchmarks), approx 0.5 (float comparison)
**Target Platform**: Linux (CI: ubuntu-latest), development in Docker container
**Project Type**: Single Rust library crate (rivrs-sparse)
**Performance Goals**: Load all 15 hand-constructed matrices in < 1 second; reconstruction error < 10^-12
**Constraints**: No external Matrix Market parsing crate; no Python; clean room implementation
**Scale/Scope**: 82 test matrices (15 hand-constructed, 10 CI-subset, 67 full SuiteSparse)

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Reconstruction error validation implemented and tested against known factorizations |
| II. Clean Room Implementation | PASS | No algorithm implementation in this phase — infrastructure only. All code is original. |
| III. Test-Driven Development | PASS | Integration tests written for all modules; validation functions tested against known results |
| IV. Algorithm Documentation | PASS | No algorithms implemented; module-level rustdoc for all public functions |
| V. Numerical Stability | PASS | Validation utilities use numerically sound formulas (relative norms, scaled errors) |
| VI. Structured Development | PASS | Phase 0.4 follows directly from completed Phase 0.2; all Phase 0 exit criteria addressed |
| VII. Code Quality | PASS | CI enforces clippy, fmt, doc; Result types for all fallible operations |

**Post-design re-check**: All gates still pass. No algorithm implementation means Principles I, II, V are satisfied by infrastructure correctness (loaders, validators). Phase 1 will be the first test of algorithm-level compliance.

## Project Structure

### Documentation (this feature)

```text
specs/004-repo-setup/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0 research decisions
├── data-model.md        # Entity definitions
├── quickstart.md        # Build/test/usage guide
├── contracts/
│   └── modules.md       # Module API contracts
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # (created by /speckit.tasks)
```

### Source Code (repository root)

```text
sparse/
├── Cargo.toml              # Add serde, serde_json dependencies
├── src/
│   ├── lib.rs              # Add pub mod io, validate
│   ├── error.rs            # Extend with IO error variants
│   ├── io.rs               # Module root: pub mod mtx, reference, registry
│   ├── io/
│   │   ├── mtx.rs          # Matrix Market parser (coordinate real symmetric)
│   │   ├── reference.rs    # JSON factorization loader (serde-based)
│   │   └── registry.rs     # Test matrix registry (metadata.json loader)
│   └── validate.rs         # Reconstruction error, backward error, inertia check
├── tests/
│   ├── common/
│   │   └── mod.rs          # Shared test helpers (existing, extend)
│   ├── hand_constructed.rs # Load all 15 matrices, verify properties, validate factorizations
│   └── suitesparse_ci.rs   # Load all 10 CI-subset matrices, verify properties
├── benches/
│   └── matrix_loading.rs   # Criterion benchmark: matrix load timing
└── .github/workflows/
    └── ci.yml              # Extended with sparse domain jobs
```

**Structure Decision**: Single library crate, flat module layout under `src/`. Test infrastructure (IO, validation) lives in the main crate because it will be reused by all future solver phases. No separate test crate — the modules are small enough to coexist with future solver code. If the crate grows large, IO and validation can be extracted to a `rivrs-sparse-testkit` crate in the workspace migration (Phase 2 of the monorepo plan).

## Implementation Approach

### Phase A: Core IO (FR-001, FR-002, FR-003)

1. **Extend `error.rs`** with IO-specific variants:
   - `IoError { source, path }` — wraps `std::io::Error`
   - `ParseError { reason, path, line }` — MTX/JSON parse failures
   - Add `From<std::io::Error>` and `From<serde_json::Error>` impls

2. **Implement `src/io/mtx.rs`** — Matrix Market parser:
   - Parse header: validate `%%MatrixMarket matrix coordinate real symmetric`
   - Skip comment lines (`%`)
   - Parse size line: `nrows ncols nnz`
   - Parse triplets: `row col value` (1-indexed → 0-indexed)
   - Symmetrize: for each off-diagonal (i,j,v), also store (j,i,v)
   - Construct `SparseColMat<usize, f64>` via `try_new_from_triplets`
   - Error on: unsupported format, malformed lines, out-of-bounds indices

3. **Implement `src/io/reference.rs`** — JSON factorization loader:
   - Define `ReferenceFactorization`, `LEntry`, `DBlock`, `Inertia` structs with `#[derive(Deserialize)]`
   - `DBlock` as enum: `OneByOne { value: f64 }` | `TwoByTwo { values: [[f64; 2]; 2] }`
   - Custom serde for `d_blocks` (size field determines variant)
   - Validate: l_entries have col < row, permutation is valid permutation of 0..n

4. **Implement `src/io/registry.rs`** — Test matrix registry:
   - Define `MatrixMetadata`, `MatrixProperties` structs with `#[derive(Deserialize)]`
   - `load_registry()` → parse `metadata.json` relative to `CARGO_MANIFEST_DIR`
   - `load_test_matrix(name)` → find in registry, load .mtx, optionally load .json
   - Return `Ok(None)` when .mtx file doesn't exist on disk (gitignored matrices)

5. **Wire up `src/io.rs`** and update `src/lib.rs`

### Phase B: Numerical Validation (FR-004)

1. **Implement `src/validate.rs`**:
   - `reconstruction_error(a, reference)`:
     - Build dense L (n×n unit lower triangular from l_entries)
     - Build dense D (n×n block diagonal from d_blocks)
     - Build permutation from reference.permutation
     - Convert sparse A to dense
     - Apply permutation: `P^T * A_dense * P` via column/row reindexing
     - Compute `L * D * L^T` via dense multiply
     - Return `||PAP - LDL^T||_F / ||A||_F`
   - `backward_error(a, x, b)`:
     - Compute r = A*x - b via `sparse_dense_matmul`
     - Return `||r|| / (||A||_F * ||x|| + ||b||)`
   - `check_inertia(computed, expected)`:
     - Exact field-wise comparison

### Phase C: Integration Tests (FR-006, FR-008)

1. **`tests/hand_constructed.rs`**:
   - Load all 15 hand-constructed matrices via registry
   - For each: verify dimensions, nnz match metadata
   - For each with reference: compute reconstruction error, assert < 1e-12
   - For each with reference: verify inertia matches

2. **`tests/suitesparse_ci.rs`**:
   - Load all 10 CI-subset matrices
   - Verify dimensions and nnz match metadata
   - Verify matrix structure (symmetric, correct index ranges)

3. **Graceful skip for full SuiteSparse**:
   - In any test that attempts to load a gitignored matrix, check path existence first
   - Print skip message if not present

### Phase D: CI Pipeline (FR-005)

1. **Extend `.github/workflows/ci.yml`**:
   - Add `test-sparse` job (matrix: stable + 1.87 MSRV)
   - Add `lint-sparse` job (clippy + fmt)
   - Add `doc-sparse` job
   - Mirror existing control/ job structure

### Phase E: Benchmark Scaffold (FR-007)

1. **`benches/matrix_loading.rs`**:
   - Criterion benchmark group
   - Benchmark: load `arrow-10-indef.mtx` from disk
   - Benchmark: load + reconstruction error computation
   - Add `[[bench]]` entry to Cargo.toml if needed

## Complexity Tracking

No constitution violations to justify. All modules are straightforward infrastructure with no algorithm complexity.
