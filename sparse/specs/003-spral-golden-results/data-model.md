# Data Model: SPRAL Golden Results Generation

**Feature**: 003-spral-golden-results
**Date**: 2026-02-06

## Entities

### 1. SPRALConfig

Captures the exact SPRAL build and runtime configuration used to generate results. Recorded once and referenced by all result files.

**Attributes**:
- `spral_commit`: string — Git commit hash of SPRAL source used
- `build_type`: string — "release" or "debug"
- `compiler`: string — Compiler identification (e.g., "gfortran 12.2.0")
- `gpu_enabled`: boolean — Always false for this feature
- `openmp_enabled`: boolean — Whether OpenMP parallelism is active
- `ordering`: integer — Ordering strategy (1=METIS default)
- `scaling`: integer — Scaling strategy (0=none)
- `pivot_method`: integer — Pivot method (1=APP aggressive)
- `u`: float — Pivot tolerance (0.01 default)
- `nemin`: integer — Supernode amalgamation parameter (32 default)
- `action`: boolean — Continue on singularity (true)

**Validation Rules**:
- `spral_commit` must be a valid 40-character hex string
- `u` must be in range [0.0, 0.5]
- `nemin` must be positive

---

### 2. SPRALResult

Per-matrix output containing all phases of SPRAL execution. One instance per test matrix.

**Attributes**:
- `matrix_name`: string — Unique identifier matching metadata.json
- `matrix_path`: string — Relative path to source .mtx file
- `status`: enum — {success, analysis_failed, factorization_failed, solve_failed, error}
- `error_flag`: integer — SPRAL inform.flag value
- `error_message`: string (nullable) — Human-readable error description

**Analysis Phase Attributes**:
- `analyze_time_seconds`: float
- `predicted_num_factor`: integer — Predicted factor entries
- `predicted_num_flops`: integer — Predicted flops
- `nparts`: integer — Number of tree partitions

**Factorization Phase Attributes**:
- `factor_time_seconds`: float
- `num_factor`: integer — Actual factor entries
- `num_flops`: integer — Actual flops
- `num_delay`: integer — Delayed columns
- `num_neg`: integer — Negative eigenvalue count
- `num_two`: integer — 2x2 pivots used
- `matrix_rank`: integer — Numerical rank
- `maxfront`: integer — Maximum front size
- `maxsupernode`: integer — Maximum supernode size
- `num_sup`: integer — Number of supernodes
- `maxdepth`: integer — Maximum elimination tree depth

**Solve Phase Attributes**:
- `solve_time_seconds`: float
- `rhs_method`: string — Always "b = A * ones(n)"
- `forward_error`: float — max |x_computed - 1.0|
- `backward_error`: float — ||b - Ax|| / (||A|| ||x|| + ||b||)

**Inertia Attributes**:
- `inertia_positive`: integer — matrix_rank - num_neg
- `inertia_negative`: integer — num_neg
- `inertia_zero`: integer — n - matrix_rank

**Matrix Diagnostics**:
- `duplicate_entries`: integer — from inform.matrix_dup
- `missing_diagonal`: integer — from inform.matrix_missing_diag
- `out_of_range`: integer — from inform.matrix_outrange

**Validation Rules**:
- `inertia_positive + inertia_negative + inertia_zero` must equal matrix dimension `n`
- `backward_error` should be < 1e-10 for well-conditioned successful factorizations
- `num_factor >= 0` and `num_flops >= 0`
- If `status == success`, all phase attributes must be populated
- If `status == analysis_failed`, only error fields are populated
- If `status == factorization_failed`, analysis fields are populated, factorization fields may be partial

**Relationships**:
- References one `SPRALConfig` (shared across all results)
- Corresponds 1:1 to a matrix entry in `metadata.json`

---

### 3. SummaryReport

Aggregated view across all SPRALResults for quick reference.

**Attributes**:
- `generated_date`: string — ISO 8601 date
- `total_matrices`: integer
- `successful`: integer
- `failed`: integer
- `entries`: list of SummaryEntry

**SummaryEntry Attributes**:
- `matrix_name`: string
- `n`: integer — Matrix dimension
- `nnz`: integer — Number of nonzeros
- `status`: string
- `category`: string — {hand-constructed, suitesparse-ci, suitesparse}
- `difficulty`: string — {trivial, easy, moderate, hard, failed}
- `analyze_time`: float
- `factor_time`: float
- `solve_time`: float
- `backward_error`: float (nullable)
- `num_delay`: integer
- `num_neg`: integer
- `num_two`: integer

**Difficulty Classification Rules**:
- `trivial`: n < 100 AND backward_error < 1e-14
- `easy`: backward_error < 1e-12 AND num_delay == 0
- `moderate`: backward_error < 1e-10 OR num_delay > 0
- `hard`: backward_error >= 1e-10 OR num_two > 0 OR factorization warnings
- `failed`: status != success

---

## State Transitions

### Result Generation Pipeline

```
Matrix (.mtx file)
    │
    ▼
[Load Matrix Market] ──── error ──→ status=error
    │
    ▼
[Convert to CSC Lower Triangle]
    │
    ▼
[spral_ssids_analyse] ──── flag<0 ──→ status=analysis_failed
    │
    ▼
[spral_ssids_factor] ──── flag<0 ──→ status=factorization_failed
    │
    ▼
[Generate RHS: b = A*ones]
    │
    ▼
[spral_ssids_solve1] ──── flag<0 ──→ status=solve_failed
    │
    ▼
[Compute Errors]
    │
    ▼
[Write JSON Result] ──→ status=success
```

## File Layout

```
test-data/spral-reference/
├── config.json                    # SPRALConfig (shared)
├── summary.json                   # SummaryReport
├── hand-constructed/              # Tier 1: git-tracked
│   ├── arrow-5-pd.json
│   ├── arrow-10-indef.json
│   └── ... (15 files)
├── suitesparse-ci/                # Tier 2: git-tracked CI subset
│   ├── LFAT5.json
│   ├── bcsstk14.json
│   └── ... (10 files)
└── suitesparse/                   # Tier 3: gitignored
    ├── aug3dcqp.json
    └── ... (67 files)
```
