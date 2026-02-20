# Phase 8.1 BLAS-3 factor_inner Investigation

**Date**: 2026-02-18
**Branch**: `017-two-level-aptp`
**Context**: Debugging backward error after refactoring `factor_inner` to BLAS-3
Factor/Apply/Update architecture matching SPRAL's `run_elim_pivoted_notasks`.

## Background

Commit `0bdf9d8` ("Refactor factor_inner to BLAS-3 Factor/Apply/Update architecture")
introduced a backward error regression on bratu3d: BE=1e-9 (Phase 7) to BE=7.76e-4
(Phase 8.1). The plan identified two bugs:

1. **BUG #1 (PRIMARY)**: Failure handling did full restore + retry instead of SPRAL-style
   partial restore + advance (keeping passed columns' L/D/Schur contributions).
2. **BUG #2 (SECONDARY)**: Missing cross-term Schur updates for failed columns
   (`update_trailing` only covered trailing x trailing, missing failed x failed and
   trailing x failed regions).

## Implementation

Five changes were made to `factor_inner` in `src/aptp/factor.rs`:

1. **`BlockBackup::restore_diagonal_with_perm`** (new method): Restores only failed columns
   from pre-factor backup with permutation mapping. Two regions: diagonal block (symmetric
   perm lookup) and panel below (column perm only). Matches SPRAL's
   `restore_part_with_sym_perm` (ldlt_app.cxx:562-574).

2. **`compute_ld`** (new function): Computes W = L * D for both 1x1 and 2x2 pivots.
   Matches SPRAL's `calcLD<OP_N>` (calc_ld.hxx:41+).

3. **`update_cross_terms`** (replaced `update_delayed`): Applies Schur complement from
   passed columns to failed and trailing regions. Two regions:
   - Failed x failed diagonal: `A[k+e..k+bs, k+e..k+bs] -= W_blk * L_blk^T`
   - Trailing x failed cross-term: `A[ts..m, k+e..k+bs] -= W_panel * L_blk^T`

4. **Row perm propagation moved before `apply_and_check`**: Steps 5-6 (row perm to
   columns 0..k + col_order update) now happen unconditionally before TRSM, matching
   SPRAL where `apply_rperm_and_backup` happens before `apply_pivot_app`.

5. **Unified success/failure path**: On threshold failure, partial restore + advance
   (no retry). k += effective_nelim where effective_nelim may be < block_size.

## Investigation Results (4 Parallel Agents)

### Agent A: Delay-Triggering Tests

Built tests using SPRAL's `cause_delays` pattern (multiply n/8 random rows by 1000x).
Initial test failures turned out to be caused by Agent B's finding (upper-triangle
test artifact), not actual factorization bugs.

**Tests added**:
- `test_factor_inner_with_delays`: 6 configurations, various (m, p, ib, seed, scale)
- `test_factor_inner_with_delays_targeted`: 12 hand-picked configurations
- `test_factor_inner_with_delays_aggressive`: 20 configurations with scale=1e6

All pass after Agent B's test fix.

### Agent B: `restore_diagonal_with_perm` Verification

**Finding**: The restore function is CORRECT. The test helper
`check_partial_factorization_in_place` had a bug — it was checking the UPPER triangle
of the factored matrix, which `swap_symmetric` doesn't maintain. Production code only
reads the lower triangle.

**Fix**: Modified `check_partial_factorization_in_place` to only check lower triangle
(`for j in 0..=i`) with weight factor 2.0 for off-diagonal entries in the Frobenius
norm. All 342 lib tests pass.

### Agent C: `update_cross_terms` Math Verification

**Finding**: The math in `compute_ld` and `update_cross_terms` is correct. The debug
test (`debug_partial_factor_config1`) showed:
- Eliminated columns: error = 2.22e-16 (machine precision)
- Cross terms: error = 2.84e-14 (fine)
- Schur "errors" up to 5.67e2 — BUT these are exclusively in the **upper triangle**
  of the Schur submatrix, which is never updated by `update_trailing`'s
  `BlockStructure::TriangularLower` GEMM.

The `check_partial_factorization_in_place` function (lower-triangle only) reports
error = 6.75e-20, confirming correctness.

### Agent D: Swap-to-End Timing Verification

**Finding**: No bug in swap-to-end timing. The delayed column swaps happen AFTER all
Schur updates are complete, ensuring the data at the new position includes the correct
Schur contributions. The loop `for i in (0..n_delayed).rev()` processes delays in
reverse order, which correctly handles cascading position adjustments.

## Key Findings

### 1. Lower-Triangle Consistency

All operations are consistently lower-triangle only:
- `factor_block_diagonal`: searches `for i in j..search_end` (lower triangle)
- Within-block Schur: `for i in j..block_end` (lower triangle)
- `update_trailing`: `BlockStructure::TriangularLower`
- `update_cross_terms`: `for i in j..n_failed` (lower triangle for diagonal region)

The upper triangle contains stale values after updates, but this is harmless because
no production code reads it.

### 2. sparsine Backward Error is Pre-Existing

The `GHS_indef/sparsine` backward error of 3.07e-3 is **NOT a regression from Phase
8.1**. Evidence:

- Phase 7 (`ssids-log.md`, Phase 7.5): sparsine listed with BE = 6.86e-4
- Phase 8.1 (`backward-error-investigation.md`): sparsine BE = 3.03e-3 with both
  METIS and MatchOrderMetis
- Classified as "Large-Front Accuracy" failure mode — large fronts (up to 916 rows)
  with ~1e-3 backward error regardless of ordering
- The 4/8 CI subset pass rate was established in Phase 7

### 3. Root Cause of Large-Front Accuracy Issues

Per `backward-error-investigation.md`, the accuracy gap vs SPRAL is likely due to:
- SPRAL uses inner block size of 32 with BLAS-3 (GEMM/TRSM) between inner blocks
  within each tile, giving better numerical error accumulation
- Our current `factor_inner` uses ib=32 inner blocks (matching SPRAL) after the
  Phase 8.1 refactoring, but the numerical properties may differ in the
  complete-pivoting implementation details
- SPRAL stores D^{-1} (inverse) and multiplies; we store D and divide
- These are systematic accuracy differences, not bugs

## SPRAL Architecture Comparison

### swap-to-end vs In-Place

| Aspect | Our swap-to-end | SPRAL's in-place |
|--------|----------------|------------------|
| During factorization | O(m) per delay (swap row+col) | No swaps, rfrom/cfrom skip |
| After factorization | Nothing — already arranged | ~70 lines of rearrangement |
| LEFT apply needed? | No (no in-place failed cols) | Yes (threshold on left blocks) |
| Cross-term handling | Swap AFTER updates -> correct | In-place during update loop |

### Block Processing

| Aspect | Our Code | SPRAL |
|--------|----------|-------|
| Inner block size | `options.inner_block_size` (default 32) | `INNER_BLOCK_SIZE = 32` (hardcoded) |
| Inner factoring | `factor_block_diagonal` (complete pivoting, BLAS-2) | `block_ldlt` (complete pivoting, BLAS-2) |
| Panel solve | `apply_and_check` (TRSM + threshold) | `apply_pivot_app` (TRSM + threshold) |
| Trailing update | `update_trailing` (GEMM, lower triangle) | `Block::update` (GEMM with rfrom/cfrom) |
| Failure handling | Partial restore + advance | Partial restore + advance |
| Cross-term updates | `update_cross_terms` (separate function) | Within `Block::update` loop |

## Test Results

- **342 lib tests**: All pass
- **22 integration tests**: All pass
- **SuiteSparse CI**: sparsine fails at 3.07e-3 (pre-existing, same as Phase 7)
- **Delay-triggering tests**: All pass with reconstruction error < 1e-10

## Remaining Work

1. **Verify bratu3d improvement**: The plan's motivation was bratu3d regressing from
   BE=1e-9 to 7.76e-4. Need to run bratu3d with current code to see if partial restore
   + cross-terms fix the regression.

2. **Large-front accuracy**: sparsine, helm3d01, dawson5, copter2, astro-ph all fail
   at ~1e-3 regardless of ordering. Potential fixes:
   - Iterative refinement (one step of `x += A^{-1}(b - Ax)`)
   - Store D^{-1} instead of D (match SPRAL)
   - Profile-guided inner block size tuning

3. **MatchOrderMetis ordering regression**: linverse degrades from 6.47e-18 (Metis)
   to 1.11e-5 (MatchOrderMetis) due to condensed graph quality. Options: heuristic
   default, try-both strategy, or expose as user choice.
