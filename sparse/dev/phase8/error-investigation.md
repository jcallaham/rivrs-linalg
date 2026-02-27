# Systematic Backward Error Investigation: SPRAL Comparison

## Context

bratu3d (n=27792) shows persistent backward error ~5e-3 with ~80K delays regardless of blocking strategy (single-level, two-level with nb=128..1024). stokes128 shows 3.65e-5 with 93K delays. bloweybq is clean at 3.26e-11 with 65 delays. The consistency across blocking strategies points to systematic bugs rather than numerical drift.

## Primary Findings

### Finding 1: `bratu3d_debug.rs` Uses Wrong Ordering

**Status: CONFIRMED via code reading**

`AnalyzeOptions::default()` returns `OrderingStrategy::Metis` (plain METIS), NOT
`MatchOrderMetis`. The CLAUDE.md memory claims "Default ordering changed to MatchOrderMetis"
but this was **never implemented** — both `AnalyzeOptions::default()` (solver.rs:67) and
`SolverOptions::default()` (solver.rs:110) return `Metis`.

The `bratu3d_debug.rs` example (line 16) uses `AnalyzeOptions::default()` = plain METIS.
Phase 7 investigation already showed bratu3d with plain METIS gives ~53K delays / 5e-3 BE,
while MatchOrderMetis gives 1 delay / 1e-9 BE.

**However**, the SuiteSparse full test (`tests/solve.rs:961`) correctly routes hard-indefinite matrices (including bratu3d) to `MatchOrderMetis`. So the integration test results should be correct. The benchmark data may be from the example using the wrong ordering.

### Finding 2: Missing Schur Cross-Term Update for Delayed Columns (PROBABLE BUG)

**Status: HIGH-CONFIDENCE HYPOTHESIS from code analysis — needs verification**

In `factor_inner`, when `block_nelim < block_size` (some columns delayed within a block):

1. `factor_block_diagonal` applies the Schur update WITHIN the block `[k..k+block_size]` only
2. `apply_and_check` computes L21 for panel rows `[k+block_size..m]` for eliminated columns
3. `update_trailing` applies `A22 -= L21 * D * L21^T` to `a[k+block_size..m, k+block_size..m]`

**BUT**: The cross-terms `a[k+block_size..m, k+block_nelim..k+block_size]` (connecting panel
rows to delayed columns within the block) are **NEVER updated**. They should receive:

```
a[trailing_i, delayed_j] -= sum_e L21[trailing_i, e] * D[e] * L_block[delayed_j, e]
```

In SPRAL, `run_elim_pivoted` updates ALL trailing blocks including cross-terms. Our
`update_trailing` only updates the square trailing block.

The placeholder `update_delayed` function (factor.rs:1446, currently a no-op) was intended
for this but never implemented.

**Impact pattern matches observations**:
- bloweybq (65 delays, max_front=13): Tiny cross-terms → negligible error (3e-11)
- stokes128 (93K delays, max_front=597): Moderate cross-terms → moderate error (4e-5)
- bratu3d (80K delays, max_front=7240): Large cross-terms → large error (5e-3)

## Investigation Plan

### Phase 0: Confirm the Symptom [30 min]

**Goal**: Determine if the 5e-3 BE comes from wrong ordering or a real regression.

1. Run bratu3d with `MatchOrderMetis` explicitly and record BE + delays
2. Run bratu3d with plain `Metis` and confirm it matches the 80K delays / 5e-3 data
3. If MatchOrderMetis gives ~1e-9: the example was using wrong ordering; fix the default
4. If MatchOrderMetis ALSO gives ~5e-3: there's a real regression (investigate further)

**Files**: `examples/bratu3d_debug.rs`, `src/aptp/solver.rs:64-69`

### Phase 1: Fix Ordering Defaults [15 min]

**Goal**: Ensure MatchOrderMetis is the default (matching the documented intent).

- `src/aptp/solver.rs:67`: Change `AnalyzeOptions::default()` to `MatchOrderMetis`
- `src/aptp/solver.rs:110`: Change `SolverOptions::default()` to `MatchOrderMetis`
- Note: `SparseLDLT::analyze()` (symbolic-only) returns an error for MatchOrderMetis since it
  needs numeric values. Default for symbolic-only API may need to stay as `Metis`.
- Update `bratu3d_debug.rs` to use `MatchOrderMetis` explicitly

### Phase 2: Implement Cross-Term Update [2-3 hours]

**Goal**: Fix the missing Schur complement update for delayed columns.

#### 2.1 Add `update_delayed_cross_terms` function

New function in `src/aptp/factor.rs` (replace the no-op `update_delayed`):

```rust
/// Apply Schur complement update from eliminated columns to the cross-region
/// connecting panel rows with delayed columns within the current block.
///
/// Computes: a[trailing..m, delayed_cols] -= (L21 * D) * L_block^T
/// where L21 = a[trailing..m, col_start..col_start+nelim] (panel L entries)
/// and L_block = a[delayed_cols, col_start..col_start+nelim] (within-block L entries)
fn update_delayed_cross_terms(
    mut a: MatMut<'_, f64>,
    col_start: usize,       // k
    nelim: usize,            // block_nelim
    block_cols: usize,       // block_size
    m: usize,
    d: &MixedDiagonal,
)
```

Algorithm:
1. Compute W = L21 * D (same as in `update_trailing`, trailing_size × nelim)
2. Copy L_block from `a[col_start+nelim..col_start+block_cols, col_start..col_start+nelim]`
   (n_delayed × nelim matrix of L entries for delayed columns)
3. GEMM: `a[trailing..m, col_start+nelim..col_start+block_cols] -= W * L_block^T`

This is a rectangular GEMM (trailing_size × n_delayed -= trailing_size × nelim × nelim × n_delayed).

#### 2.2 Insert call in `factor_inner` loop body

After `update_trailing` (step 9 in current code), add:

```rust
// Cross-term update: panel rows × delayed columns within block
if block_nelim > 0 && block_nelim < block_size {
    update_delayed_cross_terms(a.rb_mut(), k, block_nelim, block_size, m, &block_d);
}
```

This only triggers when there are both eliminated AND delayed columns in the same block.

#### 2.3 Operation order in `factor_inner` (updated)

```
1.  block_size = min(end_pos - k, ib)
2.  BlockBackup::create(a, k, block_size, m)
3.  factor_block_diagonal(a, k, block_size, small, block_end)  // block-scoped swap
4.  Permute panel column entries by block_perm
5.  Zero D off-diagonals
6.  apply_and_check → effective_nelim
7.  adjust_for_2x2_boundary
8.  IF fail: restore from backup, delay, continue
9.  Propagate row perm to columns 0..k
10. Apply block_perm to col_order
11. update_trailing (GEMM on square trailing block)
12. update_delayed_cross_terms (GEMM on rectangular cross-region)   ← NEW
13. Accumulate D, stats, pivot log
14. Handle singular columns (swap to end)
15. k += block_nelim
```

### Phase 3: Verify with Targeted Tests [1 hour]

1. **Unit test**: 8×8 matrix with block_size=4 that forces 2 eliminations + 2 delays.
   Verify reconstruction error with and without cross-term update.

2. **Regression test**: Run `cargo test --lib -- factor` (49 tests) and
   `cargo test --test solve` (22 tests).

3. **Accuracy comparison**: Run bratu3d, stokes128, bloweybq with both orderings.
   Record BE before and after fix.

4. **SuiteSparse CI subset**: `cargo test --test solve -- test_solve_suitesparse_ci`

### Phase 4: Secondary Comparisons with SPRAL [if needed]

Only pursue these if Phase 2 fix doesn't resolve the accuracy issues:

#### 4.1 D Storage Convention
- SPRAL stores D^{-1} (inverse); we store D itself
- Numerically equivalent for 1x1; slight rounding difference for 2x2
- **Files**: SPRAL `block_ldlt.hxx:362-406` vs our `diagonal.rs`

#### 4.2 Scatter Skip Logic for Delayed Columns
- Our skip at `numeric.rs:575` (`local_row >= sn_ncols && local_row < k`)
- Verify entries connecting parent columns to delayed columns are correctly
  covered by child contribution blocks
- **Files**: `numeric.rs:529-592` vs SPRAL `assemble.hxx:50-79`

#### 4.3 Two-Level Row Permutation Propagation
- `two_level_factor` propagates only `block_cols` rows (lines 1915-1936)
- Full `swap_symmetric` for delayed columns (line 1899) permutes ALL rows
- Verify pattern rows beyond block_cols are correctly handled
- **Files**: `factor.rs:1869-1986`

#### 4.4 Solve Forward/Backward Index Mapping
- Compare gather/scatter with SPRAL's `ldlt_app_solve_fwd/bwd`
- **Files**: `solve.rs:120-269` vs SPRAL `ldlt_app.cxx:2556-2589`

## Key Files

| File | What to examine |
|------|----------------|
| `src/aptp/factor.rs:1543-1856` | `factor_inner` — main target for cross-term fix |
| `src/aptp/factor.rs:1371-1435` | `update_trailing` — reference for W=L21*D computation |
| `src/aptp/factor.rs:1446-1507` | `update_delayed` — placeholder to replace |
| `src/aptp/solver.rs:64-69` | `AnalyzeOptions::default()` — ordering fix |
| `src/aptp/solver.rs:107-114` | `SolverOptions::default()` — ordering fix |
| `src/aptp/numeric.rs:529-592` | `scatter_original_entries` — secondary check |
| `src/aptp/numeric.rs:762-799` | `extract_contribution` — verify delayed data |

## Verification

1. `cargo test --lib -- factor` — all factor tests pass
2. `cargo test --lib` — all 336+ tests pass
3. `cargo test --test solve` — all 22 integration tests pass
4. `cargo clippy --lib --tests` — no new warnings
5. Run bratu3d/stokes128/bloweybq with MatchOrderMetis — compare BE before/after
6. Full SuiteSparse run if Phase 0-3 look promising
