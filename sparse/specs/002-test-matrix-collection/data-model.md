# Data Model: Test Matrix Collection

**Feature**: `002-test-matrix-collection`

## Entities

### TestMatrix

A sparse matrix stored in Matrix Market format with associated metadata.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| name | string | yes | Unique identifier (e.g., "arrow-10-indef" or "GHS_indef/ncvxqp3") |
| source | enum | yes | Origin: "hand-constructed", "suitesparse", "spral-inspired", "paper-referenced" |
| category | enum | yes | Classification: "hand-constructed", "easy-indefinite", "hard-indefinite", "positive-definite" |
| path | string | yes | Relative path from test-data/ to the .mtx file |
| size | integer | yes | Matrix dimension n (square: n x n) |
| nnz | integer | yes | Number of stored nonzeros (lower triangle for symmetric) |
| in_repo | boolean | yes | Whether the .mtx file is committed to the repository |
| download_command | string | no | ssgetpy command to fetch this matrix (if in_repo=false) |

### MatrixProperties

Structural and numerical properties of a test matrix.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| symmetric | boolean | yes | Always true for this collection |
| positive_definite | boolean | yes | True for PD validation set |
| indefinite | boolean | yes | True for easy/hard indefinite sets |
| structure | enum | no | Known structure: "arrow", "block-diagonal", "bordered-block", "tridiagonal", "augmented", "general" |
| difficulty | enum | yes | "trivial" (hand-constructed), "easy", "hard" |
| kind | string | no | SuiteSparse problem domain (e.g., "optimization problem", "structural problem") |
| expected_delayed_pivots | enum | no | Expected APTP behavior: "none", "low", "medium", "high" |
| killer_case | boolean | no | True if matrix defeats static pivoting strategies |
| singular | boolean | no | True if matrix is structurally or numerically singular |
| condition_info | string | no | Notes on conditioning (e.g., "ill-conditioned, cond ~ 1e12") |

### ExactFactorization

Pre-computed exact LDL^T factorization for hand-constructed matrices.

| Field | Type | Required | Description |
|-------|------|----------|-------------|
| matrix_name | string | yes | Reference to parent TestMatrix.name |
| permutation | array[int] | yes | Permutation vector P such that P^T A P = L D L^T |
| l_entries | array[{row, col, value}] | yes | Lower triangular factor L (unit diagonal) entries |
| d_blocks | array[{size, values}] | yes | Block diagonal D: 1x1 blocks [d] or 2x2 blocks [[d11, d21], [d21, d22]] |
| inertia | {positive, negative, zero} | yes | Eigenvalue sign counts |
| notes | string | no | Derivation method (e.g., "computed by hand", "verified with sympy") |

### MetadataIndex

The root-level index file containing all matrix entries.

| Field | Type | Description |
|-------|------|-------------|
| schema_version | string | Metadata format version (start at "1.0") |
| generated | string | ISO date of last generation |
| total_count | integer | Total number of matrices in collection |
| matrices | array[MatrixEntry] | All matrix records |

Where MatrixEntry combines TestMatrix fields, MatrixProperties (nested under `properties`), optional `factorization_path`, `paper_references` (array of citation strings), and optional `reference_results` (empty object, populated in Phase 0.3).

### SuiteSparseMeta

Additional metadata specific to SuiteSparse-sourced matrices.

| Field | Type | Description |
|-------|------|-------------|
| suitesparse_id | integer | SuiteSparse collection ID |
| group | string | SuiteSparse group (e.g., "GHS_indef") |
| psym | float | Pattern symmetry (0.0-1.0) |
| nsym | float | Numerical symmetry (0.0-1.0) |

## Relationships

```
MetadataIndex 1──* MatrixEntry
MatrixEntry 1──1 MatrixProperties
MatrixEntry 1──0..1 ExactFactorization  (hand-constructed only)
MatrixEntry 1──0..1 SuiteSparseMeta     (suitesparse only)
MatrixEntry *──* PaperReference         (citation strings)
```

## Validation Rules

1. Every TestMatrix MUST have `symmetric: true`
2. `positive_definite: true` implies `indefinite: false` and vice versa
3. If `source == "hand-constructed"`, then `factorization_path` MUST be provided
4. If `in_repo == false`, then `download_command` MUST be provided
5. `name` MUST be unique across the entire collection
6. `path` MUST point to an existing file (if `in_repo == true`) or be a valid expected path
7. `size` and `nnz` MUST match the actual matrix file contents
8. `inertia.positive + inertia.negative + inertia.zero` MUST equal `size`

## State Transitions

Not applicable — matrices are static data artifacts, not stateful objects.

## File Format Reference

### Matrix Market (.mtx)

```
%%MatrixMarket matrix coordinate real symmetric
% Comment lines
n n nnz
row col value
row col value
...
```

- 1-indexed row/column indices
- For symmetric matrices, only lower triangle entries are stored
- Real values in scientific notation (e.g., `1.234567e+02`)

### Factorization JSON (.json)

```json
{
  "matrix_name": "arrow-10-indef",
  "permutation": [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
  "l_entries": [
    {"row": 1, "col": 0, "value": 0.5},
    ...
  ],
  "d_blocks": [
    {"size": 1, "values": [2.0]},
    {"size": 2, "values": [[1.0, 0.5], [0.5, -1.0]]},
    ...
  ],
  "inertia": {"positive": 6, "negative": 4, "zero": 0},
  "notes": "Computed by hand. Verified with sympy."
}
```
