# Data Model: 004-repo-setup

**Date**: 2026-02-06
**Feature**: Repository Setup for Solver Development

## Entities

### MatrixMetadata

Represents one entry from `test-data/metadata.json`. Read-only after Phase 0.2.

| Field | Type | Description |
|-------|------|-------------|
| name | String | Matrix identifier (e.g., "arrow-10-indef") |
| source | String | Origin ("hand-constructed" or "SuiteSparse") |
| category | String | Classification ("hand-constructed", "easy-indefinite", "hard-indefinite", "positive-definite") |
| path | String | Relative path from test-data/ to .mtx file |
| size | usize | Matrix dimension (n for n×n) |
| nnz | usize | Number of stored nonzeros (lower triangle for symmetric) |
| in_repo | bool | Whether the .mtx file is committed to git |
| ci_subset | bool (optional) | Whether this is in the CI-subset (from properties or inferred from path) |
| properties | MatrixProperties | Structural/numerical properties |
| factorization_path | Option\<String\> | Relative path to companion .json (hand-constructed only) |

### MatrixProperties

Embedded within MatrixMetadata.

| Field | Type | Description |
|-------|------|-------------|
| symmetric | bool | Always true for this project |
| positive_definite | bool | Whether matrix is positive definite |
| indefinite | bool | Whether matrix is indefinite |
| difficulty | String | "trivial", "easy", "hard" |
| structure | Option\<String\> | "arrow", "tridiagonal", "block-diagonal", etc. |

### ReferenceFactorization

Loaded from companion `.json` files for hand-constructed matrices. Represents the analytically known LDL^T factorization.

| Field | Type | Description |
|-------|------|-------------|
| matrix_name | String | Must match the MatrixMetadata name |
| permutation | Vec\<usize\> | Pivot permutation (0-indexed) |
| l_entries | Vec\<LEntry\> | Strict lower-triangular entries of L |
| d_blocks | Vec\<DBlock\> | Block diagonal D (1×1 or 2×2 blocks) |
| inertia | Inertia | Eigenvalue sign counts |
| notes | String | Human-readable description |

### LEntry

Single entry in the strict lower triangle of L.

| Field | Type | Description |
|-------|------|-------------|
| row | usize | Row index (0-indexed) |
| col | usize | Column index (0-indexed), col < row |
| value | f64 | Entry value |

### DBlock

One block of the block diagonal D. Polymorphic: either 1×1 scalar or 2×2 matrix.

| Variant | Fields | Description |
|---------|--------|-------------|
| OneByOne | value: f64 | Scalar pivot |
| TwoByTwo | values: [[f64; 2]; 2] | 2×2 symmetric pivot block (row-major) |

### Inertia

Eigenvalue sign classification of a symmetric matrix.

| Field | Type | Description |
|-------|------|-------------|
| positive | usize | Count of positive eigenvalues |
| negative | usize | Count of negative eigenvalues |
| zero | usize | Count of zero eigenvalues |

**Invariant**: `positive + negative + zero == matrix dimension`

### TestMatrix

Runtime composite combining a loaded sparse matrix with its metadata and optional reference factorization.

| Field | Type | Description |
|-------|------|-------------|
| metadata | MatrixMetadata | From metadata.json |
| matrix | SparseColMat\<usize, f64\> | Loaded from .mtx file |
| reference | Option\<ReferenceFactorization\> | From .json file (hand-constructed only) |

## Relationships

```
metadata.json (82 entries)
    │
    ├── hand-constructed/ (15 matrices)
    │   ├── *.mtx → SparseColMat<usize, f64>
    │   └── *.json → ReferenceFactorization
    │
    ├── suitesparse-ci/ (10 matrices, in-repo)
    │   └── *.mtx → SparseColMat<usize, f64>  (no .json)
    │
    └── suitesparse/ (67 matrices, gitignored)
        └── *.mtx → SparseColMat<usize, f64>  (no .json)
```

## Validation Functions (Stateless)

These operate on faer dense/sparse types and return scalar error values. Not entities, but core to the feature.

| Function | Input | Output | Formula |
|----------|-------|--------|---------|
| reconstruction_error | A (sparse), L (dense), D (dense), P (perm) | f64 | \|\|P^T A P - L D L^T\|\|_F / \|\|A\|\|_F |
| backward_error | A (sparse), x (dense vec), b (dense vec) | f64 | \|\|Ax - b\|\| / (\|\|A\|\| \|\|x\|\| + \|\|b\|\|) |
| check_inertia | computed: Inertia, expected: Inertia | bool | Exact match of all three counts |
