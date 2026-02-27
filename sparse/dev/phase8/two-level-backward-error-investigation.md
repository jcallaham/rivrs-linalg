# Two-Level APTP Backward Error Investigation

## Problem Statement

Phase 8.1 introduced a two-level APTP kernel (`factor_inner` + `factor_block_diagonal`)
that regressed backward error on hard indefinite matrices compared to Phase 7's
single-level approach:

| Matrix | Phase 7 (single-level) | Phase 8.1 (two-level, ib=32) |
|--------|----------------------|------------------------------|
| bratu3d | 1.00e-9 | 8.32e-4 |
| stokes128 | 1.38e-9 | 6.10e-6 |
| bloweybq | 4.27e-11 | 4.27e-11 |

## Root Cause

The regression was caused by **restricted pivot search scope** in the two-level code.

- **Single-level** (`aptp_factor_in_place`): searches ALL remaining rows/columns of the
  frontal matrix for each pivot. For a front of size 1496 (bratu3d's max), the search
  covers 1496 rows.

- **Two-level** (`factor_inner` → `factor_block_diagonal`): inner blocking with
  `inner_block_size = 32` meant each call to `factor_block_diagonal` searched only
  32 rows for pivots. This produced suboptimal pivots, leading to massive threshold
  failures (53K delays on bratu3d with plain METIS, still poor with MatchOrderMetis).

In SPRAL's equivalent code (`block_ldlt.hxx`), the inner kernel processes the **entire
BLOCK_SIZE tile** (256 rows) at once with per-column rank-1/rank-2 updates, searching
all tile rows for each pivot. The BLAS-3 benefit comes from inter-tile operations
(TRSM for panel, GEMM for trailing), not from inner blocking within the tile.

## Approaches Tried

### 1. Cross-Update Fix (UpdateNT from Algorithm 3.1)

**Hypothesis**: Missing Schur complement update between factored columns and delayed
columns from earlier blocks (lines 10-13 of Duff et al. 2020, Algorithm 3.1).

**Implementation**: Added ~50 lines after `update_trailing` to apply the Schur complement
from the current block's factorization to previously-delayed columns.

**Result**: No improvement. bratu3d BE unchanged at 8.32e-4. The cross-update was
mathematically correct but not the dominant issue — only 3 delays triggered it.

### 2. Extended Search to Full Matrix Height

**Hypothesis**: Extend `factor_block_diagonal` to search ALL rows (m, not just
block_size) for pivots, compute L entries for all rows, and apply Schur complement
for block columns × all rows.

**Implementation**: Added `m: usize` parameter to `factor_block_diagonal`. Extended
search, L computation, and Schur complement. Extended `local_perm` from `block_size`
to `m - col_start` entries.

**Result**: 5 test failures:
- 3 tests: "permutation has duplicate index" — `local_perm` extended to full matrix
  size caused invalid permutations when mapped through `col_order` in `factor_inner`
- 2 tests: reconstruction error ~0.19 — factorization numerically incorrect
- bratu3d BE improved slightly (3.14e-4) but stokes128 and bloweybq REGRESSED

**Root cause of failure**: Two fundamental issues:
1. **Perm handling**: `local_perm` of size `m - col_start` includes panel positions.
   When applied to `col_order[k..m]`, the perm must be a valid permutation of all
   those positions. But `factor_block_diagonal` only swaps positions reachable from
   the search, so positions beyond the search scope are identity-mapped. When the
   outer loop applies the perm across iterations, the identity-mapped tail positions
   can produce duplicates.

2. **BlockBackup incompatibility**: `swap_symmetric(a, block_pos, panel_pos)` modifies
   entries in BOTH the block column (backed up) and the trailing column (NOT backed
   up). On restore, the block columns are correct but trailing columns retain
   corrupted data from the swap. This makes the matrix inconsistent.

   Specifically: `swap_symmetric` applied as a similarity transform P^T A P modifies
   entries in rows/columns `block_pos` and `panel_pos`. For entries at
   `(panel_pos, trailing_col)` where `trailing_col > block_end`, these are outside
   BlockBackup's scope. Restoring only the block columns leaves the trailing columns
   with data from the wrong row.

### 3. Full-Tile Processing (Final Fix)

**Hypothesis**: Match SPRAL's architecture — process the entire tile as one block
instead of breaking it into ib-sized sub-blocks. The search scope equals the tile
size (up to `outer_block_size`, typically 256).

**Implementation**: One-line change in `factor_inner`:
```rust
// Before:
let block_size = (end_pos - k).min(ib);
// After:
let block_size = end_pos - k;
```

**Result**: All backward errors match single-level exactly:
- bratu3d: 1.00e-9 (was 8.32e-4)
- stokes128: 1.38e-9 (was 6.10e-6)
- bloweybq: 4.27e-11 (unchanged)

**Why it works**:
1. Search scope = tile size (up to 256 rows) — same as SPRAL's block_ldlt
2. BlockBackup covers the ENTIRE tile: `a[k..m, k..k+block_size]` where
   `block_size = end_pos - k`. All swaps are within [k, end_pos), so all
   swap-affected columns are in the backup.
3. apply_and_check handles panel rows beyond the tile (TRSM + threshold)
4. update_trailing provides BLAS-3 GEMM between tiles
5. The while loop typically runs once per tile (all pivots succeed), giving
   the correct SPRAL architecture: Factor → Apply → Update per tile.

## SPRAL Architecture Comparison

The fix aligns our code with SPRAL's two-level structure:

| Component | SPRAL | Our Code (after fix) |
|-----------|-------|---------------------|
| Outer loop | `run_elim_pivoted` (BLOCK_SIZE tiles) | `two_level_factor` (outer_block_size tiles) |
| Tile factor | `block_ldlt::factor(m, n, ...)` | `factor_block_diagonal(a, k, tile_size, ...)` |
| Panel solve | `apply_pivot` (TRSM) | `apply_and_check` (TRSM + threshold) |
| Trailing update | `Block::update` (GEMM) | `update_trailing` (GEMM) |
| Inner blocking | INNER_BLOCK_SIZE=32 (optimization only) | Not used (tile processed at once) |

Key difference: SPRAL's `block_ldlt` uses INNER_BLOCK_SIZE=32 for some cache blocking
within the tile, but the search scope still covers all tile rows. Our code currently
processes the tile without inner blocking (pure BLAS-2 within tile). Inner blocking
could be re-added later as a performance optimization, but the search scope must
always cover the full tile.

## Remaining Considerations

### Inner Blocking for Performance

The current fix removes inner blocking entirely within each tile. For large tiles
(outer_block_size = 256), this means 256 rank-1/rank-2 updates, each touching up to
(m - k) rows — pure BLAS-2. The BLAS-3 payoff comes only from inter-tile operations.

Re-adding inner blocking would require:
- `factor_block_diagonal` processes ib columns but searches all tile rows
- L computation extended to all tile rows (not just ib)
- Schur complement limited to ib block columns × all tile rows
- `apply_and_check` does TRSM for remaining tile rows (tile_start + ib .. end_pos)
  plus panel rows (end_pos .. m)
- Requires correct perm handling for extended search with ib-sized blocks

This is captured in `dev/blas3-refactor.md` and deferred to a future optimization
phase. The current approach is correct and matches SPRAL's performance profile.

### `inner_block_size` Parameter

The `inner_block_size` field in `AptpOptions` / `FactorOptions` is now effectively
unused. It could be:
- Removed (breaking change for API)
- Kept for future inner blocking re-addition
- Repurposed for the complete-pivoting leaf block size within `factor_block_diagonal`

Current decision: keep as `_ib` (unused) for forward compatibility.

### SPRAL's Threshold Check

Investigation confirmed that SPRAL uses a unified threshold bound `1.0 / u` for
ALL pivot types (both 1x1 and 2x2). The `blas3-refactor.md` document's claim of
`sqrt(2)` bound for 2x2 pivots was incorrect. Our `apply_and_check` uses the correct
unified bound.

## Academic References

- Duff, Hogg & Lopez (2020): "A new sparse LDL^T solver using a posteriori
  threshold pivoting." Algorithm 3.1 defines the two-level blocked APTP.
- Hogg & Scott (2013): Algorithm 4.1 defines complete pivoting within inner blocks.
- SPRAL source: `block_ldlt.hxx` (inner kernel), `ldlt_app.cxx` (outer loop)
