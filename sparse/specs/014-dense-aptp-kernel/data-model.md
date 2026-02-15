# Data Model: Dense APTP Factorization Kernel

**Feature Branch**: `014-dense-aptp-kernel`
**Date**: 2026-02-14

## Entities

### AptpOptions

Configuration for the APTP factorization kernel.

| Field       | Type          | Description                                                                 |
|-------------|---------------|-----------------------------------------------------------------------------|
| threshold   | f64           | Stability threshold u: entries must satisfy \|l_ij\| < 1/u. Default: 0.01  |
| small       | f64           | Singularity detection: pivots with \|d\| < small treated as zero. Default: 1e-20 |
| fallback    | AptpFallback  | Strategy when 1x1 pivot fails stability check                              |

**Validation Rules**:
- `threshold` must be in (0.0, 1.0] — 0.0 would mean infinite growth allowed (use separate flag for "no checking")
- `small` must be >= 0.0
- Provides `Default` implementation: threshold=0.01, small=1e-20, fallback=BunchKaufman

### AptpFallback (enum)

Strategy for handling columns that fail the 1x1 stability check.

| Variant      | Description                                                                     |
|--------------|---------------------------------------------------------------------------------|
| BunchKaufman | Attempt 2x2 pivot with best partner column; delay if 2x2 also fails            |
| Delay        | Immediately mark column as delayed without attempting 2x2 pivot                 |

**State Transitions**:
```
Column k arrives for pivoting
  → Try 1x1 pivot (divide by diagonal)
    → PASS: accept as OneByOne
    → FAIL (BunchKaufman fallback):
        → Try 2x2 with partner column
          → PASS: accept as TwoByTwo { partner }
          → FAIL: mark as Delayed
    → FAIL (Delay fallback):
        → mark as Delayed
```

### AptpFactorResult (in-place API result)

Returned by the in-place factorization. The L factor is stored in the mutated input matrix.

| Field         | Type              | Description                                                          |
|---------------|-------------------|----------------------------------------------------------------------|
| d             | MixedDiagonal     | Block diagonal factor with mixed 1x1/2x2 blocks (from Phase 2)      |
| perm          | Vec\<usize\>      | Column permutation: perm[i] = original column index at position i    |
| num_eliminated| usize             | Number of successfully eliminated columns (q <= num_fully_summed)    |
| delayed_cols  | Vec\<usize\>      | Original column indices that were not eliminated                     |
| stats         | AptpStatistics    | Summary statistics for diagnostics                                   |
| pivot_log     | Vec\<AptpPivotRecord\> | Per-column diagnostic log                                       |

**Invariants**:
- `num_eliminated + delayed_cols.len() == num_fully_summed` (total fully-summed columns)
- `stats.num_1x1 + 2 * stats.num_2x2 + stats.num_delayed == num_fully_summed`
- `perm` is a valid permutation of `0..matrix_dimension`
- L factor (in lower triangle of mutated input) has unit diagonal for eliminated columns

### AptpFactorization (convenience API result)

Returned by the convenience wrapper that copies input and extracts L.

| Field         | Type              | Description                                                          |
|---------------|-------------------|----------------------------------------------------------------------|
| l             | Mat\<f64\>        | Unit lower triangular factor (extracted from in-place result)        |
| d             | MixedDiagonal     | Block diagonal factor (same as AptpFactorResult.d)                   |
| perm          | Perm\<usize\>     | Column permutation as faer Perm type                                 |
| delayed_cols  | Vec\<usize\>      | Original column indices not eliminated                               |
| stats         | AptpStatistics    | Summary statistics                                                   |
| pivot_log     | Vec\<AptpPivotRecord\> | Per-column diagnostic log                                       |

**Relationship to AptpFactorResult**: Constructed from AptpFactorResult by extracting the L factor from the mutated matrix and converting `Vec<usize>` permutation to `Perm<usize>`.

### AptpStatistics

Summary of the factorization outcome.

| Field        | Type  | Description                                                    |
|--------------|-------|----------------------------------------------------------------|
| num_1x1      | usize | Count of 1x1 pivots accepted                                  |
| num_2x2      | usize | Count of 2x2 pivot pairs (each pair counts as 1)              |
| num_delayed  | usize | Count of delayed columns                                       |
| max_l_entry  | f64   | Maximum absolute value across all L entries (stability metric) |

**Invariant**: `num_1x1 + 2 * num_2x2 + num_delayed == total_fully_summed_columns`

### AptpPivotRecord

Per-column diagnostic for a single pivot decision.

| Field        | Type      | Description                                              |
|--------------|-----------|----------------------------------------------------------|
| col          | usize     | Original column index                                    |
| pivot_type   | PivotType | Classification from Phase 2 (OneByOne, TwoByTwo, Delayed)|
| max_l_entry  | f64       | Worst stability metric for this column's L entries       |
| was_fallback | bool      | True if 2x2 fallback was attempted (regardless of outcome)|

## Entity Relationships

```
AptpOptions ──controls──> aptp_factor_in_place()
                               │
                               ▼
                         AptpFactorResult
                               │
              ┌────────────────┼────────────────┐
              ▼                ▼                ▼
        MixedDiagonal    AptpStatistics   Vec<AptpPivotRecord>
        (Phase 2)              │                │
              │                │                │
              ▼                ▼                ▼
         PivotType        num_1x1/2x2      PivotType
         Block2x2         num_delayed       (Phase 2)
         (Phase 2)        max_l_entry
              │
              ▼
          Inertia
          (Phase 2)
```

## Existing Types Used (Phase 2, unchanged)

| Type           | Module          | Role in Phase 5                                      |
|----------------|-----------------|------------------------------------------------------|
| MixedDiagonal  | aptp/diagonal   | D factor storage with mixed 1x1/2x2 blocks           |
| PivotType      | aptp/pivot      | Column classification (OneByOne, TwoByTwo, Delayed)   |
| Block2x2       | aptp/pivot      | 2x2 symmetric block {first_col, a, b, c}             |
| Inertia        | aptp/inertia    | Eigenvalue sign counts (positive, negative, zero)     |
| SparseError    | error           | Error type for factorization failures                 |
