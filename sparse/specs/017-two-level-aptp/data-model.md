# Data Model: Two-Level APTP Factorization

**Feature**: 017-two-level-aptp
**Date**: 2026-02-16

## Type Changes Summary

### Modified Types

#### `AptpOptions` (src/aptp/factor.rs)

**Current** (lines 64-83):
```
AptpOptions {
    threshold: f64,        // 0.01
    small: f64,            // 1e-20
    fallback: AptpFallback // BunchKaufman
}
```

**Proposed** — add block size fields:
```
AptpOptions {
    threshold: f64,           // 0.01 (unchanged)
    small: f64,               // 1e-20 (unchanged)
    fallback: AptpFallback,   // BunchKaufman (unchanged)
    outer_block_size: usize,  // NEW: 256 (default)
    inner_block_size: usize,  // NEW: 32 (default)
}
```

**Validation rules**:
- `outer_block_size > 0`
- `inner_block_size > 0`
- `inner_block_size <= outer_block_size`
- `outer_block_size` should be a multiple of `inner_block_size` for clean blocking (not enforced, just recommended)

**Impact**: AptpOptions is passed to `aptp_factor_in_place` (factor.rs), constructed in `AptpNumeric::factor` (numeric.rs), and built from `FactorOptions` (solver.rs). All three sites need updates.

#### `FactorOptions` (src/aptp/solver.rs)

**Current** (lines 62-76):
```
FactorOptions {
    threshold: f64,
    fallback: AptpFallback,
}
```

**Proposed** — add block size fields:
```
FactorOptions {
    threshold: f64,
    fallback: AptpFallback,
    outer_block_size: usize,  // NEW: 256 (default)
    inner_block_size: usize,  // NEW: 32 (default)
}
```

### Unchanged Types (verified no changes needed)

| Type | Location | Why unchanged |
|------|----------|---------------|
| `AptpFactorResult` | factor.rs:107-120 | Output interface — two-level produces the same result structure |
| `AptpFactorization` | factor.rs:124-137 | Convenience wrapper — delegates to AptpFactorResult |
| `AptpStatistics` | factor.rs:143-152 | Aggregation unchanged — stats.num_1x1/num_2x2/num_delayed still valid |
| `AptpPivotRecord` | factor.rs:156-165 | Per-column log — still one record per eliminated column |
| `AptpFallback` | factor.rs:91-96 | BunchKaufman/Delay — applies to middle level, no new variants |
| `MixedDiagonal` | diagonal.rs | D storage — same 1x1/2x2 blocks regardless of blocking strategy |
| `Block2x2` | pivot.rs | 2x2 pivot data — unchanged |
| `PivotType` | pivot.rs | OneByOne/TwoByTwo/Delayed — unchanged |
| `AptpNumeric` | numeric.rs | Multifrontal result — kernel call interface unchanged |
| `FrontFactors` | numeric.rs | Per-supernode L/D — unchanged (produced from same AptpFactorResult) |
| `SparseLDLT` | solver.rs | User-facing struct — unchanged |

### New Types

#### `BlockBackup` (src/aptp/factor.rs — internal)

**Purpose**: Store backup of matrix entries for one outer block column, enabling restore on pivot failure.

```
BlockBackup {
    col_start: usize,       // First column of the outer block
    block_cols: usize,       // Number of columns in this block
    data: Mat<f64>,          // Copy of A[col_start:, col_start:col_start+block_cols]
                             // (block_cols columns, m-col_start rows)
}
```

**Operations**:
- `create(a: MatRef, col_start, block_cols, m) -> Self`: Copy relevant columns
- `restore_failed(a: MatMut, col_start, nelim, block_cols, m)`: Restore columns [col_start+nelim, col_start+block_cols) from backup

**Visibility**: Private to factor.rs (not exported).

**Memory**: O(nb × m) per block. Since blocks are processed sequentially, only one BlockBackup is live at a time. Total backup memory is O(nb × m), not O(m²).

## Data Flow

```
FactorOptions (user)
    ↓ (construct)
AptpOptions (with block sizes)
    ↓ (pass to)
aptp_factor_in_place
    ↓ (dispatch based on p vs outer_block_size)
    ├─ p ≤ nb: factor_inner (single block)
    └─ p > nb: two_level_factor (block loop)
         ↓ (per outer block)
         ├─ BlockBackup::create
         ├─ factor_inner (nb×nb diagonal)
         │    ↓ (per ib sub-block)
         │    ├─ complete_pivoting_factor (ib×ib leaf)
         │    └─ update_schur (inner Schur complement)
         ├─ apply_and_check (TRSM + threshold check → nelim)
         ├─ BlockBackup::restore_failed (if nelim < block_cols)
         ├─ update_trailing (GEMM)
         └─ update_delayed (UpdateNT/UpdateTN)
    ↓ (accumulate)
AptpFactorResult (unchanged output)
```

## Entity Relationship

```
AptpOptions ──configures──→ aptp_factor_in_place
    │                            │
    ├── threshold, small ────→ factor_inner (inner APTP)
    ├── fallback ────────────→ factor_inner (middle-level fallback)
    ├── outer_block_size ────→ dispatch + block loop stride
    └── inner_block_size ────→ factor_inner sub-block stride
                                 │
                                 └─→ complete_pivoting_factor (leaf)
                                      │
                                      └─→ AptpFactorResult (per-block)
                                           │
                                           ├── d: MixedDiagonal
                                           ├── perm: Vec<usize>
                                           ├── num_eliminated: usize
                                           └── delayed_cols: Vec<usize>
```
