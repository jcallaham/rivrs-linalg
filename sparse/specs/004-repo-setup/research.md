# Research: 004-repo-setup

**Date**: 2026-02-06
**Feature**: Repository Setup for Solver Development

## R1: Matrix Market Parsing Strategy

**Decision**: Custom lightweight parser scoped to `coordinate real symmetric` format.

**Rationale**: All test matrices (hand-constructed and SuiteSparse) use the same Matrix Market coordinate format with header `%%MatrixMarket matrix coordinate real symmetric`. The format is simple (header line, optional comments, size line, triplet lines). A custom parser gives full control over error messages and direct construction of faer `SparseColMat` via `try_new_from_triplets`. The `matrixmarket-rs` crate adds a dependency for minimal gain and may not produce faer types directly.

**Alternatives considered**:
- `matrixmarket-rs` crate — adds dependency, requires conversion to faer types, may not handle all comment styles
- Build script code generation — over-engineered for this use case, complicates debugging

**Format spec** (from inspection of actual files):
```
%%MatrixMarket matrix coordinate real symmetric
% optional comment lines
<nrows> <ncols> <nnz>
<row> <col> <value>    # 1-indexed, lower triangle only for symmetric
...
```

## R2: JSON Factorization Schema

**Decision**: Use `serde_json` for deserialization into typed Rust structs. Schema is stable and well-understood from Phase 0.2 output.

**Rationale**: The JSON files are authored by Phase 0.2 tooling and follow a consistent schema. `serde` + `serde_json` are standard Rust ecosystem choices with zero controversy. The `d_blocks` field has a polymorphic structure (size 1 → scalar values, size 2 → 2x2 matrix values) that needs a custom deserializer or enum representation.

**Schema summary** (from inspection of all 15 JSON files):
```json
{
  "matrix_name": "string",
  "permutation": [int, ...],           // 0-indexed
  "l_entries": [{"row": int, "col": int, "value": float}, ...],  // strict lower triangle
  "d_blocks": [
    {"size": 1, "values": [float]},           // 1x1 pivot
    {"size": 2, "values": [[f,f],[f,f]]}      // 2x2 pivot (row-major)
  ],
  "inertia": {"positive": int, "negative": int, "zero": int},
  "notes": "string"
}
```

**New dependencies needed**: `serde`, `serde_json` (dev-dependencies for test infrastructure).

## R3: faer API for Validation Arithmetic

**Decision**: Leverage faer's built-in dense arithmetic for validation. No custom matrix math needed.

**Rationale**: faer 0.22 provides all required operations:
- `SparseColMat::try_new_from_triplets()` — construct sparse matrices from MTX data
- `.to_dense()` — convert to dense `Mat<f64>` for small matrices
- `&A * &B`, `&A - &B` — dense matrix multiply and subtract via operator overloads
- `.transpose()` — zero-cost transpose view
- `norm_l2()` — Frobenius norm (what faer calls `norm_l2` for matrices)
- `sparse_dense_matmul()` — sparse matrix × dense vector for backward error
- `Perm` / `PermRef` — permutation types with `permute_rows` / `permute_cols`

For hand-constructed matrices (max 20×20), dense conversion is trivially efficient. The validation pattern:
1. Build L as dense lower-triangular from `l_entries` + identity diagonal
2. Build D as dense block-diagonal from `d_blocks`
3. Build P as `Perm` from `permutation` array
4. Compute `L * D * L^T` via dense multiply
5. Convert sparse A to dense, apply permutation: `P^T * A_dense * P`
6. Compute `||P^T A P - L D L^T||_F / ||A||_F`

**Alternatives considered**:
- `permute_self_adjoint` on sparse A — works but requires workspace allocation setup; for tiny matrices, dense permutation via row/col reordering is simpler
- `nalgebra` for dense ops — unnecessary second dependency when faer already provides everything

## R4: Test Data Registry Design

**Decision**: Parse `metadata.json` at test time using `serde_json`. Provide helper functions that resolve matrix names to absolute file paths relative to `CARGO_MANIFEST_DIR`.

**Rationale**: The metadata.json file (82 entries) is the single source of truth established in Phase 0.2. Parsing it into a `Vec<MatrixMetadata>` struct gives typed access to all matrix properties. Using `CARGO_MANIFEST_DIR` env var (set by cargo at compile time) ensures paths resolve correctly in both `cargo test` and IDE test runners.

**Three-tier loading strategy**:
- **Hand-constructed** (`test-data/hand-constructed/`): Always available, always tested. Have companion `.json` files.
- **CI-subset** (`test-data/suitesparse-ci/`): Always available (committed to git), tested in CI. No `.json` files — load as sparse matrix only, verify dimensions/nnz from metadata.
- **Full SuiteSparse** (`test-data/suitesparse/`): Gitignored, may not be present. Tests skip gracefully with `println!("Skipping: ...")` when directory absent.

## R5: CI Pipeline Extension

**Decision**: Add `test-sparse`, `lint-sparse`, and `doc-sparse` jobs to existing `.github/workflows/ci.yml`, mirroring the structure of the existing `control/` jobs.

**Rationale**: The existing CI already has a clean pattern (test matrix with stable + MSRV, separate lint job, separate doc job). Extending it with parallel sparse jobs keeps the pattern consistent and the CI fast (sparse and control jobs run in parallel).

**Key consideration**: The CI-subset SuiteSparse matrices (~73MB) are committed to git. GitHub Actions checkout will include them. No special caching or artifact download needed.

## R6: Benchmark Scaffold

**Decision**: Single criterion benchmark group in `benches/` that loads a hand-constructed matrix and times the load operation, with a placeholder for future solver timing.

**Rationale**: Criterion is already a dev-dependency. A minimal benchmark that loads `arrow-10-indef.mtx` establishes the infrastructure without pretending to benchmark nonexistent solver operations. Future phases will add benchmark groups for analyze, factor, and solve.
