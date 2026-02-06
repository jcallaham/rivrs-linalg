# Data Model: Core Test Infrastructure

**Branch**: `005-test-infrastructure` | **Date**: 2026-02-06

## Entities

### SolverTestCase

A complete test scenario for validating a solver component.

| Field | Type | Description | Source |
|-------|------|-------------|--------|
| name | String | Matrix identifier (e.g., "arrow-5-pd") | metadata.json |
| matrix | SparseColMat<usize, f64> | The sparse symmetric matrix | .mtx file |
| properties | TestMatrixProperties | Structural/numerical properties | metadata.json |
| reference | Option<ReferenceFactorization> | Known-correct factorization (hand-constructed only) | .json file |

**Lifecycle**: Created (load from registry or generate) → Used (passed to SolverTest methods) → Discarded (no mutation).

**Identity**: Uniquely identified by `name`. Registry matrices are singletons; generated matrices use caller-provided names.

### TestMatrixProperties

Structural and numerical properties of a test matrix. Extends/wraps existing `registry::MatrixProperties`.

| Field | Type | Description |
|-------|------|-------------|
| size | usize | Matrix dimension (n for n×n) |
| nnz | usize | Number of stored nonzeros |
| symmetric | bool | Always true for this project |
| positive_definite | bool | Whether matrix is PD |
| indefinite | bool | Whether matrix has mixed eigenvalue signs |
| difficulty | String | "trivial", "easy", "hard" |
| structure | Option<String> | "arrow", "tridiagonal", "block-diagonal", etc. |
| source | String | "hand-constructed", "suitesparse", or "generated" |
| category | String | Classification for filtering |

### TestResult

Structured outcome of a solver test.

| Field | Type | Description |
|-------|------|-------------|
| passed | bool | Overall pass/fail |
| metrics | Vec<MetricResult> | Individual metric outcomes |
| diagnostics | Vec<String> | Human-readable diagnostic messages |
| matrix_name | String | Which matrix was tested |
| test_kind | TestKind | Which test method produced this result |

### MetricResult

A single numerical metric check.

| Field | Type | Description |
|-------|------|-------------|
| name | String | e.g., "reconstruction_error", "backward_error" |
| value | f64 | Computed metric value |
| threshold | f64 | Pass/fail threshold |
| passed | bool | value < threshold |

### TestKind (enum)

| Variant | Description |
|---------|-------------|
| Analyze | Symbolic analysis phase validation |
| Factor | Numeric factorization phase validation |
| Solve | Triangular solve phase validation |
| Roundtrip | Full pipeline end-to-end validation |

### NumericalValidator

Configurable validation engine.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| reconstruction_tol | f64 | 1e-12 | Max allowed reconstruction error |
| backward_error_tol | f64 | 1e-10 | Max allowed backward error |

**Methods**:
- `default()` → NumericalValidator with constitution-mandated defaults
- `with_reconstruction_tol(tol)` → Builder method
- `with_backward_error_tol(tol)` → Builder method
- `check_reconstruction(matrix, reference) → MetricResult`
- `check_backward_error(matrix, x, b) → MetricResult`
- `check_inertia(computed, expected) → MetricResult`
- `validate_factorization(test_case, reference) → TestResult` — runs all applicable checks

### TestCaseFilter

Criteria for loading subsets of test matrices.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| source | Option<String> | None | Filter by "hand-constructed" or "suitesparse" |
| category | Option<String> | None | Filter by category |
| difficulty | Option<String> | None | Filter by difficulty |
| ci_only | bool | false | Only CI-subset matrices |
| require_reference | bool | false | Only matrices with reference factorizations |

## Relationships

```
TestCaseFilter  ──loads──►  Vec<SolverTestCase>
                                   │
                                   ├── matrix: SparseColMat
                                   ├── properties: TestMatrixProperties
                                   └── reference: Option<ReferenceFactorization>
                                          │
                                          ├── permutation: Vec<usize>
                                          ├── l_entries: Vec<LEntry>
                                          ├── d_blocks: Vec<DBlock>
                                          └── inertia: Inertia

SolverTest trait
    ├── test_analyze(case) → TestResult
    ├── test_factor(case)  → TestResult
    ├── test_solve(case)   → TestResult
    └── test_roundtrip(case) → TestResult
                                   │
                                   ├── passed: bool
                                   ├── metrics: Vec<MetricResult>
                                   └── diagnostics: Vec<String>

NumericalValidator
    ├── check_reconstruction(matrix, reference) → MetricResult
    ├── check_backward_error(matrix, x, b)      → MetricResult
    └── check_inertia(computed, expected)        → MetricResult
```

## Existing Types (Unchanged)

These types from Phase 0.4 remain as-is. The new harness wraps them but does not modify them.

- `registry::MatrixMetadata` — Metadata from metadata.json
- `registry::MatrixProperties` — Structural properties (wrapped by TestMatrixProperties)
- `registry::TestMatrix` — Loaded matrix with reference (wrapped by SolverTestCase)
- `reference::ReferenceFactorization` — Known-correct L, D, P, inertia
- `reference::Inertia` — Eigenvalue sign counts
- `reference::LEntry` — Single L factor entry
- `reference::DBlock` — 1x1 or 2x2 diagonal block
- `error::SparseError` — Error type for I/O and parse failures
