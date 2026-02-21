# Accuracy Investigation: Post-Ordering Pipeline

Tracking the investigation into backward error discrepancies between our solver and
SPRAL on "easy-indefinite" and "hard-indefinite" SuiteSparse matrices.

## Status

**Phase**: Diagnostic (comparing our pipeline against SPRAL stage-by-stage)

**Failing matrices** (13 of 65 SuiteSparse, all with MatchOrderMetis ordering):

| Matrix | n | Backward Error | Category |
|--------|---|---------------|----------|
| sparsine | 50,000 | 2.0e-3 | easy-indefinite |
| copter2 | 55,476 | 6.4e-4 | easy-indefinite |
| dawson5 | 51,537 | 6.8e-4 | easy-indefinite |
| helm3d01 | 32,226 | 9.8e-4 | easy-indefinite |
| astro-ph | 16,706 | 1.5e-3 | easy-indefinite (singular) |
| TSOPF_FS_b162_c1 | 10,798 | 2.2e-7 | hard-indefinite |
| cont-300 | 180,895 | 2.4e-5 | hard-indefinite |
| cvxqp3 | 17,500 | 1.2e-7 | hard-indefinite |
| d_pretok | 182,730 | 2.5e-6 | hard-indefinite |
| ncvxqp3 | 75,000 | 4.5e-8 | hard-indefinite |
| ncvxqp5 | 62,500 | 1.1e-7 | hard-indefinite |
| ncvxqp7 | 87,500 | 1.7e-8 | hard-indefinite |
| thread | 29,736 | 1.6e-4 | positive-definite |

Target: backward error < 5e-11 (SPRAL's threshold).

## Key Findings

### 1. SPRAL Uses TPP Fallback After APP Failure

**Source**: `spral/src/ssids/cpu/factor.hxx`, `factor_node_indef()`

SPRAL's per-supernode factorization has a two-pass strategy:

1. **First pass**: APP (A Posteriori Threshold Pivoting) via `ldlt_app_factor`
2. **If APP doesn't eliminate all columns** (`nelim < ncol`): TPP (Threshold Partial
   Pivoting) via `ldlt_tpp_factor` on the remaining columns

The `inform` struct tracks both passes:
- `not_first_pass`: pivots that failed APP
- `not_second_pass`: pivots that also failed TPP (truly delayed to parent)

**Our code only does APP.** When APP can't eliminate a pivot, we delay it unconditionally
to the parent supernode. SPRAL salvages many of these with TPP first. This means:

- We delay more pivots than SPRAL does
- Delayed pivots accumulate error through extend-add assembly
- This is a plausible root cause for the accuracy gap, especially on easy-indefinite
  matrices where plain METIS (no pre-pairing) produces more initial delays

**Hypothesis**: Adding a TPP fallback pass would significantly reduce delays and improve
backward error on the failing matrices.

### 2. SPRAL's Default Ordering is Plain METIS

**Source**: `spral/src/ssids/datatypes.f90` line 209, Duff et al. (2020) paper

SPRAL defaults to `ordering=1` (plain METIS, no matching/scaling). The Duff 2020 paper
is explicit:

> "In the case of easy-indefinite matrices, no scaling is performed prior to the
> factorization whereas hard-indefinite matrices are scaled and ordered using a
> matching-based ordering and scaling."

This means SPRAL achieves good backward error on easy-indefinite matrices **with plain
METIS and no scaling** — relying on the APTP kernel (with TPP fallback) to handle delays.

Our solver changed the default to MatchOrderMetis specifically because plain METIS
caused unacceptable backward error. This is consistent with Finding #1: our APP-only
kernel can't handle the higher delay counts that plain METIS produces.

### 3. MC64 Scaling: Unmatched Index Handling (astro-ph)

**Source**: `spral/src/scaling.f90` lines 688-801

SPRAL comparison on astro-ph (structurally singular, matched=15677/16706) shows scaling
vectors diverge by ~1e217 relative diff on unmatched indices. Fill quality still matches
(ratio=0.990).

Root cause: SPRAL has dead code in the singular path of `hungarian_wrapper`. The guard
`match(i) < 0` on line 699 is never true because `hungarian_match` uses 0 (not negative)
for unmatched rows. Effectively:

- SPRAL re-matches the **full** graph (not the matched subgraph)
- The Duff-Pralet correction for unmatched indices **never fires**
- All indices get `scaling[i] = exp((u[i] + v[i] - cmax[i]) / 2)` using duals from the
  full matching, where unmatched rows have non-zero duals from Dijkstra iterations

Our `DualsDirect` strategy explicitly zeros `u[i] = 0` for unmatched rows before
computing scaling. For matched indices the duals and scaling are identical. For unmatched
indices, SPRAL's duals retain values from Dijkstra iterations while ours are zeroed.

This is a cosmetic difference — it only affects ~6% of indices on astro-ph and the fill
quality is unaffected. The accuracy issue on astro-ph is from its singularity and the
general kernel accuracy problem (Finding #1), not from scaling differences.

**Note**: We use `SingularScalingStrategy::DualsDirect`, not `RematchDuffPralet`. SPRAL's
dead code means SPRAL also effectively uses the DualsDirect approach (no subgraph
extraction, no Duff-Pralet correction), just with different zero-initialization for
unmatched duals.

## SPRAL Full Solve Comparison Results

**Tool**: `examples/spral_solve_comparison.rs` + `tools/spral_full_solve.f90`

Both solvers receive identical matrix + ordering (from our `match_order_metis`) + scaling + RHS.
SPRAL uses `ordering=0` (user-supplied) and `scaling=0` (user-supplied via `scale` argument).

### Failing Matrices

| Matrix | n | our_BE | spral_BE | ratio | delays | 2x2 | !1st | !2nd | maxfrt |
|--------|---|--------|----------|-------|--------|-----|------|------|--------|
| sparsine | 50,000 | 7.3e-2 | 1.3e-12 | 5.7e10 | 23 | 9,722 | 669 | 23 | 11,124 |
| copter2 | 55,476 | 5.5e-1 | 6.6e-13 | 8.4e11 | 19 | 20,229 | 519 | 14 | 1,675 |
| dawson5 | 51,537 | 1.5e0 | 3.5e-13 | 4.2e12 | 12 | 19,954 | 66 | 4 | 976 |
| helm3d01 | 32,226 | 3.5e-1 | 2.8e-13 | 1.2e12 | 7 | 12,202 | 258 | 6 | 1,493 |
| astro-ph | 16,706 | 3.1e-2 | 6.6e-17 | 4.7e14 | 6,350 | 5,850 | 6,586 | 6,184 | 2,282 |
| TSOPF_FS_b162_c1 | 10,798 | 6.0e-6 | 7.0e-15 | 8.5e8 | 299 | 3,607 | 376 | 299 | 326 |
| cont-300 | 180,895 | 6.5e-1 | 1.3e-13 | 4.9e12 | 0 | 89,624 | 4 | 0 | 906 |
| cvxqp3 | 17,500 | 4.2e-6 | 2.7e-11 | 1.6e5 | 64 | 6,464 | 104 | 64 | 1,827 |
| d_pretok | 182,730 | 2.9e-16 | 9.8e-16 | 0.3 | 27 | 77,468 | 2 | 0 | 1,180 |
| ncvxqp3 | 75,000 | 9.4e-7 | 9.2e-10 | 1,026 | 84 | 27,263 | 529 | 64 | 3,837 |
| ncvxqp5 | 62,500 | 7.0e-6 | 7.9e-13 | 8.9e6 | 55 | 24,167 | 765 | 37 | 2,585 |
| ncvxqp7 | 87,500 | 2.0e-6 | 4.4e-10 | 4,493 | 115 | 32,974 | 469 | 106 | 4,632 |
| thread | 29,736 | 7.4e-3 | 1.6e-15 | 4.6e12 | 2,890 | 2,203 | 9,893 | 2,890 | 3,504 |

### Passing Matrices (Baseline)

| Matrix | n | our_BE | spral_BE | ratio | delays | 2x2 | !1st | !2nd | maxfrt |
|--------|---|--------|----------|-------|--------|-----|------|------|--------|
| bloweybq | 10,001 | 9.1e-17 | 8.3e-17 | 1.1 | 5 | 2,523 | 5 | 0 | 59 |
| t2dal | 4,257 | 3.1e-16 | 4.6e-16 | 0.7 | 0 | 1,993 | 0 | 0 | 134 |
| bratu3d | 27,792 | 4.0e-15 | 4.2e-13 | 0.01 | 29 | 11,506 | 296 | 2 | 1,507 |
| stokes128 | 49,666 | 2.6e-15 | 2.4e-15 | 1.1 | 1 | 20,144 | 0 | 0 | 555 |
| ncvxqp1 | 12,111 | 1.0e-13 | 8.4e-15 | 12.4 | 93 | 4,733 | 115 | 92 | 1,454 |
| cfd2 | 123,440 | 8.0e-16 | 3.0e-15 | 0.3 | 0 | 27,227 | 0 | 0 | 2,146 |

### Analysis

**Key finding: SPRAL achieves dramatically better backward error on all failing matrices
(except d_pretok where we're already better).** The gap ranges from 10^3 to 10^14.

This **confirms the kernel hypothesis**: given identical ordering+scaling, SPRAL's
factorization kernel (APP + TPP fallback) produces much more accurate results than our
APP-only kernel.

**Patterns in the data:**

1. **`not_first_pass` correlates with our failures**: Matrices where SPRAL reports many
   `not_first_pass` events (APP failures salvaged by TPP) are exactly where our solver
   fails. Our code delays these pivots to the parent; SPRAL eliminates them with TPP.

2. **`not_second_pass` = SPRAL's true delays**: The `not_second_pass` count represents
   pivots that also fail TPP and are genuinely delayed. Compare:
   - thread: our delays=2890, SPRAL not_second_pass=2890 (same — all our delays also fail TPP)
   - cvxqp3: our delays=64, SPRAL not_second_pass=64 (same — these are genuine delays)
   - sparsine: our delays=23, SPRAL not_second_pass=23 (same)
   - **But** SPRAL salvaged 646/669 in sparsine, 505/519 in copter2, 62/66 in dawson5

3. **cont-300 is revealing**: SPRAL has `not_first_pass=4, not_second_pass=0` meaning
   only 4 pivots failed APP and all were salvaged by TPP. But cont-300 has **zero delays**
   in our code and SPRAL! Yet our BE is 6.5e-1 vs SPRAL's 1.3e-13. This suggests the
   accuracy issue is NOT just about delay handling — there may be a numerical precision
   issue in the BLAS-3 factorization path itself (TRSM/GEMM accumulation).

4. **d_pretok is the outlier**: We beat SPRAL (0.3x ratio). This has only 2 APP failures
   and 0 TPP failures — trivially easy for both solvers.

5. **Passing matrices have low `not_first_pass`**: bloweybq (5), t2dal (0), stokes128 (0),
   cfd2 (0). These matrices have few or no APP failures, so the TPP fallback is irrelevant.
   bratu3d is interesting: 296 APP failures, but our code still achieves better BE (0.01x).

**Revised hypothesis**: Two independent accuracy issues:

- **Issue A**: Missing TPP fallback. This affects matrices with many `not_first_pass`
  events (thread, astro-ph, ncvxqp*, cvxqp3, etc.). Adding TPP would reduce delays.

- **Issue B**: Numerical precision in BLAS-3 path. cont-300 has zero delays and only 4
  APP failures, yet our BE is 10^12 times worse. This suggests a bug or precision issue in
  our factor_inner TRSM/GEMM pipeline, possibly in how we handle the lower triangle or
  in partial restore paths. The large maxfront matrices (sparsine: 11,124) may also suffer.

### 4. Inner Block Size Directly Controls Accuracy (CONFIRMED)

**Tool**: `examples/block_size_experiment.rs`

Testing cont-300 and dawson5 with different `(outer_block_size, inner_block_size)`:

**cont-300** (n=180,895, maxfront=906, zero delays):

| Configuration | Backward Error |
|---|---|
| default (ob=256, ib=32) | **6.48e-1** |
| single-level (ob=huge, ib=32) | 6.88e-5 |
| **single-level unblocked (ob=huge, ib=huge)** | **1.21e-14** |
| two-level (ob=256, ib=256) | 3.06e-1 |
| ob=128, ib=32 | 6.46e-1 |
| ob=512, ib=32 | 1.44e-1 |
| ob=256, ib=8 | 6.48e-1 |
| ob=256, ib=64 | 6.63e-1 |

**dawson5** (n=51,537, maxfront=976, 12 delays):

| Configuration | Backward Error |
|---|---|
| default (ob=256, ib=32) | **1.45e0** |
| single-level (ob=huge, ib=32) | 3.47e-1 |
| **single-level unblocked (ob=huge, ib=huge)** | **5.79e-14** |
| **two-level (ob=256, ib=256)** | **7.92e-14** |
| ob=128, ib=32 | 1.30e0 |
| ob=512, ib=32 | 1.08e0 |
| ob=256, ib=8 | 9.80e-1 |
| **ob=256, ib=64** | **2.24e-13** |

**Critical conclusions:**

1. **The unblocked path achieves SPRAL-level accuracy** (~1e-14). The old
   `complete_pivoting_factor` path (column-by-column BLAS-2) is correct.

2. **The inner block size is the critical variable.** On dawson5:
   - ib=32 → 1.45e0 (terrible)
   - ib=64 → 2.24e-13 (excellent!)
   - ib=256 → 7.92e-14 (excellent!)

3. **The two-level outer loop also contributes.** On cont-300:
   - ob=huge, ib=32 → 6.88e-5 (better than ob=256, ib=32 → 6.48e-1)
   - ob=huge, ib=huge → 1.21e-14 (perfect)

   So the outer loop multiplies the inner-block error by ~10^4x on cont-300.

4. **The bug is NOT a pure numerical precision issue** — the error difference
   between ib=32 and ib=unblocked is too large (10^12x) for normal floating
   point accumulation. There is likely a correctness bug in `factor_inner`'s
   inner blocking loop that only manifests when `ib < num_fully_summed`.

### 5. Two Distinct Bugs + Missing Feature (CONFIRMED)

Full ib=ob and unblocked tests on all 13 failing matrices reveal three categories:

| Matrix | ib=32 | ib=256 (=ob) | unblocked | Root Cause |
|--------|-------|-------------|-----------|------------|
| dawson5 | 1.4e0 | **7.9e-14** | 5.8e-14 | Inner blocking bug |
| TSOPF_FS_b162_c1 | 6.0e-6 | **4.6e-15** | 5.0e-15 | Inner blocking bug |
| copter2 | 5.5e-1 | 3.7e-1 | **1.7e-13** | Two-level outer loop bug |
| helm3d01 | 3.5e-1 | 3.1e-1 | **5.3e-14** | Two-level outer loop bug |
| astro-ph | 3.1e-2 | 2.5e-2 | **3.4e-16** | Two-level outer loop bug |
| ncvxqp5 | 7.0e-6 | 6.9e-6 | **7.0e-14** | Two-level outer loop bug |
| thread | 7.4e-3 | 5.1e-4 | **2.7e-15** | Two-level outer loop bug |
| cont-300 | 6.5e-1 | 3.1e-1 | **1.2e-14** | Two-level outer loop bug |
| cvxqp3 | 4.2e-6 | 4.0e-6 | **1.2e-11** | Two-level + missing TPP |
| ncvxqp3 | 9.4e-7 | 1.1e-6 | **4.8e-11** | Two-level + missing TPP |
| ncvxqp7 | 2.0e-6 | 5.0e-7 | **7.0e-12** | Two-level + missing TPP |

**Bug 1: Inner blocking loop in `factor_inner`** (affects dawson5, TSOPF_FS_b162_c1):
Fixed by setting ib=ob (single block per factor_inner call). The multi-block
inner loop (ib < p) has a correctness issue when iterating through multiple
ib-sized blocks within a single factor_inner call.

**Bug 2: Two-level outer loop in `two_level_factor`** (affects copter2, helm3d01,
astro-ph, ncvxqp5, thread, cont-300): NOT fixed by ib=ob, but fixed by
unblocked (ob=huge, ib=huge). The outer loop that calls factor_inner on
successive 256-column blocks has its own precision issue, likely in how it
propagates row permutations to prior blocks (lines 2043-2056) or handles
delayed columns from inner calls (lines 2058-2071).

**Missing feature: TPP fallback** (affects cvxqp3, ncvxqp3, ncvxqp7): Even the
unblocked path achieves ~1e-11, not the ~1e-15 that SPRAL gets. These need TPP
to salvage pivots that fail the APP threshold check (SPRAL's `not_first_pass`).

### 6. Detailed SPRAL Code Comparison (factor_inner + two_level_factor)

Two agents performed line-by-line comparison of our `factor_inner` and `two_level_factor`
against SPRAL's `run_elim_pivoted_notasks` (ldlt_app.cxx) and `factor_node_indef`
(factor.hxx). Key architectural differences and suspect deviations documented below.

#### SPRAL's Architecture (Reference)

SPRAL uses a **tiled block layout**: the matrix is divided into `block_size × block_size`
rectangular tiles, each with its own metadata (`CBlockData` with per-block `nelim`).
The inner loop iterates over diagonal tiles:

1. **Factor diagonal tile** → produces `lperm`, `nelim`, D
2. **ApplyT** (left tiles, jblk < blk): `apply_rperm_and_backup` + `apply_pivot<OP_T>`
   - Row permutation via `lperm`, then TRSM + D-scaling
   - `from = cdata[jblk].nelim` skips already-processed columns
3. **ApplyN** (below tiles, iblk > blk): `apply_cperm_and_backup` + `apply_pivot<OP_N>`
   - Column permutation via `lperm`, then TRSM + D-scaling
4. **Update**: Schur complement for all remaining tiles, using `rfrom`/`cfrom` skip

Our code uses a **monolithic matrix** with column ranges instead of tiles.

#### Deviation A: Missing rfrom/cfrom Skip Logic (two_level_factor — CRITICAL)

**SPRAL** (`Block::update`, ldlt_app.cxx:1082-1153):
```cpp
int rfrom = (i_ == elim_col) ? cdata_[i_].nelim : 0;
int cfrom = (j_ == elim_col) ? cdata_[j_].nelim : 0;
// GEMM only on rows [rfrom..nrow) and cols [cfrom..ncol)
```

When computing Schur updates, SPRAL skips rows/columns that were already eliminated
in the current column's diagonal block. This avoids redundant computation AND prevents
applying the Schur complement to already-committed L entries.

**Our code** (`update_trailing`, factor.rs:1447-1511):
```rust
let trailing_start = col_start + block_cols;
// GEMM covers ALL rows [trailing_start..m) — no skip logic
```

Our Schur update always starts from `trailing_start = col_start + block_cols` (the full
block boundary), with no awareness of which rows/columns within the block were eliminated
vs delayed.

**Impact**: In the two-level outer loop, when `factor_inner` returns with some delayed
columns, the subsequent outer-loop Schur update applies to the entire trailing region
without distinguishing between rows that belong to eliminated vs delayed columns within
the current outer block. This could corrupt L entries that are already committed.

**Severity**: HIGH — this is the most likely cause of Bug 2 (copter2, helm3d01, cont-300,
etc.). The error disappears when ob=huge because the two-level outer loop never runs.

#### Deviation B: No Per-Block Elimination Tracking (two_level_factor — STRUCTURAL)

**SPRAL**: Each tile has `CBlockData` with a `nelim` field tracking how many columns
in that tile have been eliminated. This metadata is essential for Deviation A's skip logic.

**Our code**: No per-block metadata. We only track the global `k` position (next column
to eliminate) and `end_pos` (boundary of delayed columns). Without per-block nelim, we
cannot implement rfrom/cfrom skip logic.

**Impact**: This is the structural prerequisite for fixing Deviation A. We would need to
either:
1. Add per-outer-block metadata tracking eliminated count, OR
2. Restructure to a tiled layout matching SPRAL

#### Deviation C: Row Permutation Propagation Timing (two_level_factor — SUSPECT)

**SPRAL**: Row permutations from the diagonal block's `factor()` are propagated to
left tiles via `apply_rperm_and_backup` DURING the ApplyT phase — one left tile at a
time, interleaved with the TRSM.

**Our code** (`two_level_factor`, factor.rs:2043-2056): Row permutations are propagated
to ALL prior outer blocks at once, AFTER `factor_inner` returns. The propagation loops
over all prior columns [0..ob_start) in a single pass.

**Impact**: MEDIUM — the mathematical effect should be the same (permutation is the same
regardless of when it's applied), but the interleaved approach in SPRAL may interact with
the backup/restore mechanism differently.

#### Deviation D: Panel Permutation Scope in factor_inner (INVESTIGATED — NOT A BUG)

Agent a6af333 flagged that Step 3 (panel permutation, lines 1752-1760) operates on all
`block_size` columns while only `block_nelim` were factored. After investigation, this
matches SPRAL's `apply_cperm_and_backup` which also permutes all columns of below-
diagonal blocks. The TRSM subsequently only operates on the `block_nelim` passed columns.

Similarly, Step 5 (row perm to prior columns, lines 1783-1793) was flagged for reading
from `a[(k + block_perm[i], c)]` after Step 3 physically permuted panel entries. After
investigation, Steps 3 and 5 operate on **disjoint matrix regions** (panel rows vs left
columns), so there is no interaction. This matches SPRAL's separate `apply_rperm` and
`apply_cperm` operations.

#### Summary of Deviations

| ID | Location | Severity | Description |
|----|----------|----------|-------------|
| A | two_level_factor | **CRITICAL** | Missing rfrom/cfrom skip in Schur updates |
| B | two_level_factor | STRUCTURAL | No per-block elimination tracking (prerequisite for A) |
| C | two_level_factor | MEDIUM | Row perm propagation timing differs from SPRAL |
| D | factor_inner | NOT A BUG | Panel/row perm scope matches SPRAL after analysis |

**Conclusion**: The most promising fix for Bug 2 is implementing SPRAL's rfrom/cfrom skip
logic in the two-level Schur update. This requires adding per-outer-block elimination
tracking (Deviation B) and using it to restrict the GEMM range in `update_trailing` and
`update_cross_terms` when called from `two_level_factor`.

For Bug 1 (inner blocking in factor_inner), the code comparison did not conclusively
identify the root cause. The multi-block inner loop's correctness should be verified by
constructing minimal reproducing cases with known factorizations, and by comparing
intermediate matrix state (after each inner block) against the unblocked path.

## Diagnostic Plan

### Completed

- [x] `spral_full_solve.f90` — Confirms kernel is the bottleneck, not ordering/scaling
- [x] `spral_comparison.rs` — Confirms ordering+fill quality matches SPRAL
- [x] Block size experiment — Isolates three root causes

### Next Steps (Priority Order)

1. **Fix Bug 2: Implement rfrom/cfrom skip logic in `two_level_factor`**
   (6 matrices). Deviation A is the most likely root cause. Requires:
   a. Add per-outer-block elimination tracking (Deviation B)
   b. Pass rfrom/cfrom parameters to `update_trailing` and `update_cross_terms`
      when called from `two_level_factor`'s outer loop
   c. GEMM range restricted to only un-eliminated rows/columns

2. **Fix Bug 1: Debug `factor_inner` multi-block inner loop**
   (2 matrices). The code comparison found no conclusive deviation from SPRAL.
   Next approach: compare intermediate matrix state after each inner block
   against the unblocked path on a small (n~100) test case. Diff the matrix
   snapshots to pinpoint where they diverge.

3. **Interim workaround**: Set `inner_block_size = outer_block_size` and/or
   increase default `inner_block_size` from 32 to ≥64 in the default options.
   On dawson5, ib=64 already achieves 2.24e-13 (vs 1.45e0 with ib=32).

4. **Add TPP fallback**: After APP pass, if `nelim < ncol`, run TPP on remaining
   columns. This is the only fix for cvxqp3/ncvxqp3/ncvxqp7.

## SPRAL Public API for Comparison

The full Fortran-level pipeline is callable with user-supplied ordering+scaling:

```fortran
options%ordering = 0    ! user-supplied ordering
options%scaling = 0     ! user-supplied scaling via scale(:) arg to ssids_factor

call ssids_analyse(check, n, ptr, row, akeep, options, inform, order=order)
call ssids_factor(posdef, val, akeep, fkeep, options, inform, scale=scale)
call ssids_solve(x, akeep, fkeep, options, inform)
call ssids_enquire_indef(akeep, fkeep, options, inform, piv_order=po, d=d)
```

Key `inform` fields: `num_delay`, `num_two`, `num_neg`, `maxfront`, `matrix_rank`,
`num_factor`, `num_flops`, `not_first_pass`, `not_second_pass`.

D from `ssids_enquire_indef`: `d(1,i)` = diagonal of D^{-1}, `d(2,i)` = off-diagonal.
Entries in pivot order (use `piv_order` to map back to original variables).
