# MC64 Pipeline Investigation: Matching, Scaling, and Condensation

**Date**: 2026-02-18
**Branch**: `017-two-level-aptp`
**Context**: After fixing the extract_front_factors bug, remaining backward error failures
correlate strongly with MC64 matching quality. dawson5 passes with plain METIS but FAILS
with MatchOrderMetis — our MC64 is actively harmful for that matrix.

## Motivation

Post-fix backward error results show two categories of remaining failures:

1. **MC64-quality correlated** (astro-ph, copter2, helm3d01, TSOPF×2, dawson5 with MatchOrder):
   All have imperfect matchings (unmatched rows and/or low min_rmax)
2. **Large-front** (sparsine): Perfect MC64 but max_front=11131

Duff et al. solved all of these to machine precision with SPRAL, which uses HSL_MC64
(the reference implementation). Our MC64 is a clean-room reimplementation.

Key evidence that this is a bug, not a fundamental limitation:
- dawson5: BE=7.45e-17 with METIS, FAILS with MatchOrderMetis (MC64 harmful)
- pwtk: min_rmax=3.68e-5 (terrible MC64 quality) but PASSES after extraction fix
- Duff et al. report machine-precision results on all these matrices

## MC64 Matching Quality on Failing Matrices

| Matrix | BE | matched | min_rmax | MC64 helpful? |
|--------|-----|---------|----------|:---:|
| astro-ph | 2.79e-3 | 15677/16706 (94%) | 0.013 | Harmful? |
| copter2 | 6.87e-4 | 55421/55476 (99.9%) | 0.010 | Harmful? |
| dawson5 | FAIL w/ MO | 51500/51537 (99.9%) | 0.015 | **Actively harmful** |
| helm3d01 | 1.24e-3 | 32225/32226 (99.99%) | 0.124 | Harmful? |
| TSOPF_b162 | 2.34e-7 | 10797/10798 (99.99%) | 0.020 | Harmful? |
| TSOPF_b39 | 5.98e-8 | 28216/28216 (100%) | 0.167 | Harmful? |
| sparsine | 2.18e-3 | 50000/50000 (100%) | 1.0 | N/A (front size) |

## Audit Agents

### Agent A: MC64 Matching Algorithm (`matching.rs` vs SPRAL `scaling.f90`)

**Scope**: Compare our Hungarian algorithm implementation against SPRAL's `hungarian_match`

**Key questions**:
- Does our Dijkstra shortest-path augmentation match SPRAL's implementation?
- Is the initialization heuristic correct? (SPRAL's `hungarian_init_heuristic`)
- Are dual variables (u, v) updated correctly during augmentation?
- Does the algorithm correctly handle the symmetric case?
- Why does astro-ph get only 94% matching? Is this correct or suboptimal?
- What is the matching weight compared to what SPRAL would produce?

**Our files**: `src/aptp/matching.rs` (1672 lines)
**SPRAL files**: `/opt/references/spral/src/scaling.f90` lines 584-1171 (`hungarian_match`,
`hungarian_init_heuristic`, `hungarian_wrapper`)

**Status**: COMPLETE
**Findings**:

#### 1. Initialization Heuristic (greedy_initial_matching vs hungarian_init_heuristic)

Both implementations use a two-phase greedy initialization. The overall structure is equivalent:
Pass 1 finds the minimum-cost column for each row, matches greedily where possible (skipping
dense columns); Pass 2 scans unmatched columns for low reduced-cost rows, with length-2
augmenting path fallback. Two differences were found:

**Difference 1a: Tie-breaking in Pass 2 missing RINF guard**

SPRAL (scaling.f90 lines 883-885):
```fortran
if ((di .eq. vj) .and. (di .ne. RINF)) then
   if ((iperm(i) .ne. 0) .or. (iperm(i0) .eq. 0)) cycle
end if
```
SPRAL checks `di /= RINF` before applying the "prefer unmatched row" tie-breaking rule.
This prevents spurious replacements when both rows have infinite reduced cost (which happens
for rows with no valid edges).

Our code (matching.rs lines 492-494):
```rust
if rc == best_rc && row_match[i] == UNMATCHED && row_match[best_i] != UNMATCHED
```
Missing the RINF guard. However, our cost graph contains only finite costs (zeros are
filtered out in `build_cost_graph`), and `u[i]` is set to 0 for empty rows, so reduced costs
are always finite for existing edges. The guard cannot trigger in practice.

**Classification: CORRECT DIFFERENCE** -- defensive programming in SPRAL; no practical impact.

**Difference 1b: `d_col[j]` initialization**

SPRAL (scaling.f90 line 869): `d(1:n) = 0.0` -- initializes the per-column best reduced-cost
array to zero before Pass 2. Our code (matching.rs line 470): `let mut d_col = vec![0.0_f64; n]`
-- also initializes to zero. Equivalent behavior.

For columns already matched from Pass 1, their `d(j)` / `d_col[j]` remains 0, meaning the
2-augmentation check requires `reduced_cost <= 0` (zero-reduced-cost edges only). Both
implementations use the same condition.

**Classification: CORRECT DIFFERENCE** -- Identical behavior.

#### 2. Dijkstra Shortest-Path Augmentation

The core Dijkstra augmentation follows Duff & Koster's Algorithm MPD in both implementations.
The overall logic is: scan root column, build initial heap/Q2, iterate (extract min from Q2
or heap, scan matched column, relax neighbors), augment along shortest path, update duals.

**Difference 2a: Priority queue -- decrease-key heap vs lazy-deletion heap**

SPRAL uses a custom binary min-heap with `heap_update` (decrease-key/bubble-up),
`heap_delete` (remove-at-position), and `heap_pop` (remove-root). Position tracking via
the `l(i)` array enables O(log n) decrease-key. Additionally:
- `l(i) >= up` means row i is visited/settled -- skip (line 1084).
- `l(i) >= low` means row i is in Q2 or settled -- skip update (line 1098).
- `l(i) == 0` means row i has never been touched.

Our code uses Rust's `BinaryHeap` which does not support decrease-key. Instead:
- When a row's distance decreases, a new `(Reverse(dnew), i)` entry is pushed (line 713).
- Old entries with higher distances remain as stale entries.
- When popping, stale entries are discarded via `!visited[i] && d[i] <= dmin` (line 645)
  and `d[q0] >= csp || visited[q0]` (line 658).

This "lazy deletion" approach is a well-known correct alternative. Each row is processed
at most once (the `visited` flag prevents re-processing), and the first time a row is
popped, its `d[i]` holds the correct shortest distance. The heap may contain O(m) entries
(where m = number of relaxations) vs SPRAL's O(n), but this does not affect correctness.

**Missing Q2-skip optimization:** SPRAL checks `if (l(i) >= low) cycle` (line 1098) to
skip rows already in Q2. Our code may push a row to Q1 even if it is already there.
However, the `visited` check prevents double-processing, so this is a performance-only
difference.

**Classification: CORRECT DIFFERENCE** -- Both produce correct shortest-path augmentations.

**Difference 2b: Q2 processing order (LIFO vs FIFO)**

SPRAL's Q2 set `q(low:up-1)` is processed from right to left: take `q(up-1)`, decrement
`up`. New entries prepended at `low-1`. Result: new entries processed LAST.

Our Q1 Vec is processed via `q1.pop()` (LIFO, takes from end). New entries appended via
`q1.push(i)` (go to end). Result: new entries processed FIRST.

Since all rows in Q2/Q1 have the same distance (dmin), processing order does not affect
which augmenting path length is found (all have equal cost). Different order may lead to
different matchings, but both are equally optimal.

**Classification: CORRECT DIFFERENCE** -- No impact on optimality.

**Difference 2c: `vj` computation -- O(1) vs O(degree)**

SPRAL (line 1081): `vj = dq0 - val(jperm(j)) + dualu(q0)`

SPRAL stores `jperm(j)` = CSC position of the matched entry in column j, giving O(1) access
to the matched entry's cost.

Our code (matching.rs lines 676-685): Scans column j to find row q0's entry, costing
O(degree of column j).

Mathematically identical: `vj = d[q0] - cost[q0,j] + u[q0]` in both cases.

**Classification: CORRECT DIFFERENCE** -- Performance only; O(degree) per settled row vs O(1).

**Difference 2d: State management between augmentations**

SPRAL uses persistent arrays (`d(1:m)`, `l(1:m)`) initialized once, reset only for touched
rows after each augmentation (lines 1144-1153).

Our code allocates fresh `d`, `visited`, `pr`, `out`, `heap`, `q1` arrays in each
`dijkstra_augment` call. Costs O(n) allocation per augmentation.

**Classification: CORRECT DIFFERENCE** -- Performance only; no correctness impact.

**Difference 2e: Augmentation traceback**

Both trace backward from `(isp, jsp)` to root column using `pr[j]` and `out[j]`. The only
difference is representation: SPRAL's `isp` is a CSC position (requires `row(isp)` to get
row index); ours is the row index directly. SPRAL also stores `jperm(j) = isp` during
traceback for future O(1) matched-entry lookups. The path logic is identical.

**Classification: CORRECT DIFFERENCE** -- Representation only.

#### 3. Symmetric Case Handling

**This is the most significant area of divergence between our implementations.**

**SPRAL `hungarian_wrapper` (scaling.f90 lines 679-801):**

For symmetric matrices with FULL matching (`matched == n`):
- Returns `rscaling = dualu`, `cscaling = dualv - cmax`.
- `match_postproc` applies a magnitude balance adjustment that has NO net effect on the
  final symmetric scaling (row/column shift cancels in `(rscaling + cscaling)/2`).
- Final: `scaling = exp((dualu + dualv - cmax) / 2)`.

For symmetric matrices with INCOMPLETE matching (`matched < n`):
1. **Identify matched subset**: Build `old_to_new` / `new_to_old` mappings (lines 699-714).
   Matched rows get new consecutive indices 1..nn. Unmatched rows get negative sentinels.
2. **Extract full-rank submatrix**: Overwrite `ptr2`, `row2`, `val2` with the nn x nn
   submatrix containing only matched rows/columns (lines 716-735). The cost values are
   preserved from the first `hungarian_wrapper` pass (already `cmax - log|a|`).
3. **Re-match**: Call `hungarian_match` on this nn x nn submatrix (line 738). This produces
   **fresh, optimal dual variables** for the matched portion.
4. **Compute scaling from second matching's duals**: For matched index i:
   `rscaling[i] = (dualu_2nd[new_i] + dualv_2nd[new_i] - cmax[old_i]) / 2` (line 751).
   Note: `cmax` uses the ORIGINAL column index because column maxima reflect the full matrix.
5. **Sentinel for unmatched**: `rscaling[i] = -huge(rscaling)` (line 748).
6. **Duff-Pralet correction** in log space (lines 777-797): For unmatched i connected to
   matched k: `rscaling[i] = max(rscaling[i], log|a_ik| + rscaling[k])`. Then negate.
   If still at sentinel, set to 0 (exp(0)=1).
7. **Symmetrize**: `cscaling = rscaling` (line 800). Final: `scaling = exp(rscaling)`.

**Our code (matching.rs lines 196-250):**
- Full matching path (lines 196-207): Compute column duals via complementary slackness,
  symmetrize scaling, enforce scaling bound. Matches SPRAL for the full-matching case.
- Incomplete matching path (lines 209-250):
  1. Zero `u[i]` for unmatched rows.
  2. Compute column duals from the partial matching's u values.
  3. Symmetrize scaling.
  4. Apply Duff-Pralet correction in linear space.
  5. `enforce_scaling_bound` to fix violations.

**DEFINITE BUG: Missing re-matching for structurally singular symmetric matrices.**

When the first matching is incomplete, the dual variables are NOT optimal for the matched
subgraph. SPRAL addresses this by running the Hungarian algorithm a second time on the
extracted full-rank submatrix. This second run produces optimal duals, which directly
determine the scaling quality for the matched portion.

Our code skips this second matching entirely and uses the partial-matching duals directly.
These duals were computed on the full bipartite graph where the augmenting path search was
influenced by edges to/from rows that ended up unmatched.

**Impact assessment:**
- astro-ph (94% matched, 1029 unmatched): Most severe. Large unmatched fraction means the
  first matching's duals are significantly distorted.
- copter2 (99.9% matched, 55 unmatched): Moderate. Few unmatched rows but min_rmax=0.010
  suggests duals are still significantly off.
- dawson5 (99.9% matched, 37 unmatched): Explains why MC64 is actively harmful -- suboptimal
  scaling introduces worse diagonal dominance than the raw matrix.
- helm3d01 (99.99% matched, 1 unmatched): Even 1 unmatched row triggers the suboptimal path.
- TSOPF_b162 (99.99% matched, 1 unmatched): Same.

#### 4. Edge Weight Computation

Both implementations use identical edge weight construction:
```
cost[i,j] = log(max_k |a[k,j]|) - log|a[i,j]|
```

This transforms the maximum-product matching into a minimum-sum assignment. SPRAL computes
`log(abs(val))` first, then `colmax - log_val`. Our code computes `col_max_log - val.ln()`.
Mathematically identical. Both handle symmetry expansion (upper triangle to full storage)
before computing column maxima. Both skip zero entries. Both produce non-negative costs.

**Classification: CORRECT DIFFERENCE** -- No discrepancy.

#### 5. Unmatched Row Handling and Reporting

| Aspect | SPRAL | Ours | Verdict |
|--------|-------|------|---------|
| Unmatched sentinel | `iperm(i) = 0` | `row_match[i] = UNMATCHED` (usize::MAX) | CORRECT DIFFERENCE |
| Dual zeroing | `where(iperm==0) dualu=0` (line 1169) | `u[i] = 0 for unmatched` (line 218) | Equivalent |
| Matching cardinality | `num` output parameter | `matched` field in `Mc64Result` | Equivalent |
| Incomplete-match scaling | Re-match submatrix, then D-P | Direct D-P from first matching | **DEFINITE BUG** (see #3) |
| `enforce_scaling_bound` | Not present in SPRAL | Iterative cap to enforce `|s_i a_ij s_j| <= 1` | CORRECT DIFFERENCE (safety net) |

The `enforce_scaling_bound` function (lines 785-837) has no SPRAL equivalent. SPRAL's dual
variables naturally satisfy dual feasibility (`u[i] + v[j] <= c[i,j]`), which translates
to `|s_i * a_ij * s_j| <= 1`. Our duals may violate this due to the missing re-matching
step, necessitating the safety net. The function itself is correct and convergent.

#### 6. Additional Observations

**`match_postproc` for Full-Rank Symmetric Case:**
SPRAL's `match_postproc` (lines 1608-1615) shifts `rscaling` and `cscaling` to balance
their averages. Since the final symmetric scaling uses `(rscaling + cscaling)/2` and the
shift cancels, this has ZERO effect. Our omission is correct.

**O(1) matched-edge lookup via `jperm`:**
SPRAL maintains `jperm(j)` = CSC position of the matched entry in column j. This is updated
during augmentation traceback (`jperm(jj) = klong`, line 1135) and used for O(1) `vj`
computation. Our code does not maintain this, requiring O(degree) column scans. This is
purely a performance difference. For very wide columns in structurally irregular matrices,
this could add noticeable overhead but does not affect correctness.

#### Summary of Findings

| # | Finding | Classification | Severity |
|---|---------|---------------|----------|
| 1a | Missing RINF guard in tie-breaking | CORRECT DIFFERENCE | None |
| 1b | `d_col` initialization | CORRECT DIFFERENCE | None |
| 2a | Lazy-deletion heap vs decrease-key heap | CORRECT DIFFERENCE | Performance only |
| 2b | Q2 processing order (LIFO vs FIFO) | CORRECT DIFFERENCE | None |
| 2c | `vj` computation O(degree) vs O(1) | CORRECT DIFFERENCE | Performance only |
| 2d | Per-augmentation allocation vs persistent arrays | CORRECT DIFFERENCE | Performance only |
| 2e | Augmentation traceback representation | CORRECT DIFFERENCE | None |
| **3** | **Missing re-matching for structurally singular symmetric case** | **DEFINITE BUG** | **HIGH -- root cause** |
| 4 | Edge weight computation | CORRECT DIFFERENCE | None |
| 5 | Unmatched row handling (except re-matching) | CORRECT DIFFERENCE | None |
| 6 | Missing `match_postproc` magnitude adjustment | CORRECT DIFFERENCE | None |

**Primary finding**: The core Hungarian algorithm (initialization + Dijkstra augmentation +
dual updates + traceback) is correctly implemented and matches SPRAL's logic. The Dijkstra
shortest-path computation produces correct optimal matchings.

**The single significant bug** is the missing re-matching step for structurally singular
symmetric matrices (Finding #3). When `matched < n`, SPRAL extracts the matched submatrix
and runs the Hungarian algorithm a second time to obtain optimal dual variables for the
matched portion. Our code skips this and computes scaling from the first (suboptimal) duals.
This directly explains the poor scaling quality (low min_rmax) on all failing matrices with
incomplete matchings and likely explains why dawson5 degrades with MatchOrderMetis.

**Recommended fix**: In `mc64_matching`, when `matched < n`:
1. Extract the `matched x matched` subgraph from `CostGraph`, remapping row/column indices.
2. Run `greedy_initial_matching` + `dijkstra_augment` loop on this subgraph.
3. Use the second matching's dual variables (u2, v2) with the ORIGINAL `col_max_log` values
   to compute scaling for matched indices.
4. Apply Duff-Pralet correction for unmatched indices using the second-matching scaling.
5. The `enforce_scaling_bound` step may become unnecessary if the second matching produces
   clean dual feasibility, but can be kept as a safety net.

### Agent B: MC64 Scaling Computation (`matching.rs` + `ordering.rs`)

**Scope**: Verify that dual variables are correctly converted to scaling factors

**Key questions**:
- Are the dual variables (u for rows, v for columns) correctly extracted from the matching?
- Is the scaling formula correct? (SPRAL: `scaling(i) = exp(u(i) + v(i))/2` for symmetric)
- For the symmetric case, does SPRAL's `hungarian_wrapper` do post-processing?
- How does SPRAL handle unmatched rows in the scaling? (lines 689-738 of scaling.f90)
- Does SPRAL's "build a full rank submatrix" approach differ from ours?
- Compare our scaling values against SPRAL's on a small test matrix

**Our files**: `src/aptp/matching.rs` (scaling computation), `src/aptp/ordering.rs` (scaling application)
**SPRAL files**: `/opt/references/spral/src/scaling.f90` lines 134-170 (`hungarian_scale_sym_int64`),
597-801 (`hungarian_wrapper`)

**Status**: COMPLETE
**Findings**:

#### 1. Dual Variable to Scaling Factor Conversion

**SPRAL's full-rank symmetric pipeline** (tracing through the call chain):

```
hungarian_scale_sym_int64 (line 134-170):
  calls hungarian_wrapper(sym=.true., ...)
  then: scaling(1:n) = exp( (rscaling(1:n) + cscaling(1:n)) / 2 )

hungarian_wrapper (lines 679-685, full-rank path):
  rscaling(1:m) = dualu(1:m)
  cscaling(1:n) = dualv(1:n) - cmax(1:n)
  calls match_postproc(...)

match_postproc (lines 1608-1615, square case):
  ravg = sum(rscaling(1:m)) / m
  cavg = sum(cscaling(1:n)) / n
  adjust = (ravg - cavg) / 2
  rscaling(1:m) = rscaling(1:m) - adjust
  cscaling(1:n) = cscaling(1:n) + adjust
```

So SPRAL's complete formula for the full-rank symmetric case is:
```
rscaling_raw[i] = dualu[i]
cscaling_raw[j] = dualv[j] - cmax[j]
adjust = (mean(rscaling_raw) - mean(cscaling_raw)) / 2
rscaling[i] = rscaling_raw[i] - adjust
cscaling[j] = cscaling_raw[j] + adjust
final_scaling[i] = exp( (rscaling[i] + cscaling[i]) / 2 )
```

**Our formula** (matching.rs `symmetrize_scaling`, line 766):
```
scaling[i] = exp( (u[i] + v[i] - col_max_log[i]) / 2 )
```

Before the `match_postproc` adjustment, SPRAL computes:
```
rscaling[i] + cscaling[i] = dualu[i] + (dualv[i] - cmax[i])
                           = dualu[i] + dualv[i] - cmax[i]
```
which equals our `u[i] + v[i] - col_max_log[i]`. So the raw formula matches.

**POTENTIAL BUG: Missing `match_postproc` magnitude adjustment.**
SPRAL's `match_postproc` applies a magnitude balancing adjustment for square
matrices: it shifts rscaling/cscaling by `(ravg - cavg)/2` to balance the
average row and column scalings. Because the final symmetric formula is
`exp((rscaling + cscaling)/2)` and the shift cancels out in the sum
(rscaling decreases by adjust, cscaling increases by adjust), the sum
`rscaling[i] + cscaling[i]` is unchanged. Therefore this adjustment has
**NO EFFECT** on the final symmetric scaling. **CORRECT DIFFERENCE** -- our
omission of `match_postproc` does not change the output for the symmetric
full-rank case.

#### 2. Symmetric Case: Incomplete Matching (Structurally Singular)

**DEFINITE BUG: Missing second matching on the full-rank submatrix.**

SPRAL's `hungarian_wrapper` (lines 688-801) handles the structurally singular
symmetric case with a fundamentally different strategy than ours:

**SPRAL's approach:**
1. Run `hungarian_match` on the full (symmetrized) cost graph. Gets a partial
   matching of `nn < n` rows.
2. Build a **full-rank submatrix** keeping only the `nn` matched rows/columns
   (lines 716-735). This submatrix includes the already-computed cost values
   (in log-space).
3. Run `hungarian_match` **a second time** on this `nn x nn` submatrix
   (line 738). This produces **fresh dual variables** (dualu, dualv) on the
   reduced system.
4. Compute scaling from the **second matching's duals**:
   `rscaling[i] = (dualu[new_i] + dualv[new_i] - cmax[old_i]) / 2`
   (line 751). Note the key detail: `cmax` uses the **original** column index
   because the column maxima come from the original graph structure.
5. Set unmatched indices to `-huge(rscaling)` (a sentinel for "log-space
   negative infinity", line 748).
6. Apply Duff-Pralet correction in **log space** (lines 777-797):
   - For unmatched `i` connected to matched `k`:
     `rscaling[i] = max(rscaling[i], log|a_ik| + rscaling[k])`
   - Then negate: `rscaling[i] = -rscaling[i]`
   - Convention: if still at sentinel, set to `0.0` (will become `exp(0)=1`)
7. Copy rscaling to cscaling (symmetric: line 800).
8. Final: `scaling[i] = exp((rscaling[i] + cscaling[i])/2) = exp(rscaling[i])`

**Our approach** (matching.rs lines 209-250):
1. Run matching once. If `matched < n`:
2. Zero `u[i]` for unmatched rows (line 218).
3. Compute column duals via complementary slackness on the **partial** matching
   (line 223).
4. Compute scaling from these duals via `symmetrize_scaling` (line 224).
5. Apply Duff-Pralet correction (line 237) operating on the **linear-domain**
   scaling factors.
6. Apply `enforce_scaling_bound` to fix violations (line 240).

**Critical differences:**

a) **No second matching.** SPRAL re-runs the Hungarian algorithm on the
   full-rank submatrix to get **optimal** dual variables for the matched
   portion. We skip this and use the dual variables from the original partial
   matching. These duals are NOT optimal for the matched subgraph because
   the augmenting path search was influenced by edges to/from unmatched
   rows. The Dijkstra-based dual updates (`u[i] += d[i] - csp`) are only
   applied to visited rows, and the search paths may have been distorted by
   the presence of unmatched rows.

   **Classification: DEFINITE BUG** -- For matrices with incomplete matchings
   (dawson5, astro-ph, copter2, helm3d01, TSOPF_b162), the dual variables
   are suboptimal, producing poor scaling for the matched portion.

b) **Zeroing unmatched u values.** We set `u[i] = 0` for unmatched rows
   (line 218). SPRAL never touches the unmatched rows' duals from the first
   matching; it sets them to `-huge` as a sentinel after the second matching.
   Our zeroing interacts with `compute_column_duals` and `symmetrize_scaling`
   in a way that could produce incorrect values when an unmatched row is
   incident to a matched column's matched edge (since `v[j] = cost[matched_row,j] - u[matched_row]`
   and the matched row is always matched, the zeroing of unmatched rows
   should not affect column duals for matched columns directly). However,
   the fundamental problem is (a), not the zeroing per se.

   **Classification: POTENTIAL BUG** -- Indirect effect; may compound with (a).

#### 3. Weight Transformation

**SPRAL:**
- `hungarian_wrapper` line 647: `val2(k) = log(abs(val(k)))` initially
- Line 656: `val2(k) = colmax - val2(k)` (where colmax = max log|a_kj| in col)
- So cost `c[i,j] = log(col_max_j) - log|a[i,j]|` (non-negative)
- Dual variables `dualu`, `dualv` are in **cost space** (same as our `u`, `v`)
- Final scaling: `rscaling = dualu`, `cscaling = dualv - cmax`
  i.e., SPRAL subtracts cmax to convert from cost-space duals to log-space
  scaling, then `exp()` converts to linear space.

**Ours:**
- `build_cost_graph` line 393: `c = col_max_log[j] - abs_val.ln()`
- Identical cost definition.
- `symmetrize_scaling` line 766: `exp((u[i] + v[i] - col_max_log[i]) / 2)`
- We subtract `col_max_log[i]` from dual `v[i]` inline.

**Classification: CORRECT DIFFERENCE** -- The weight transformation is
mathematically identical. Both use `c[i,j] = log(col_max_j) - log|a_ij|`
as costs, and both subtract `col_max_log` when converting duals to scaling.

#### 4. Scaling for Unmatched Rows

**SPRAL (structurally singular path, lines 767-797):**
1. Matched indices: `rscaling[i] = (dualu_2nd[new_i] + dualv_2nd[new_i] - cmax[old_i]) / 2`
   (computed from the second matching's optimal duals).
2. Unmatched indices initially set to `-huge(rscaling)` (sentinel).
3. Duff-Pralet correction in **log space**:
   `rscaling[i] = max(rscaling[i], log|a_ik| + rscaling[k])` for matched `k`.
4. Then negate: `rscaling[i] = -rscaling[i]`.
5. If still at sentinel: `rscaling[i] = 0.0` (exp(0) = 1.0).
6. Final: `scaling[i] = exp(rscaling[i])` (since rscaling = cscaling for symmetric).

**Ours (lines 843-900):**
1. `duff_pralet_correction` works in **linear space**:
   - `log_max[i] = max(ln|a_ik| + ln(scaling[k]))` for matched `k`.
   - `scaling[i] = exp(-log_max[i])`.
   - If isolated (no matched neighbors): `scaling[i] = 1.0`.

The Duff-Pralet formula itself is mathematically equivalent in log vs linear
form: SPRAL computes `rscaling[i] = -max_k(log|a_ik| + rscaling[k])` where
`rscaling[k]` is in log space, then `exp()`. We compute
`scaling[i] = exp(-max_k(ln|a_ik| + ln(scaling[k])))`. These are the same
provided `rscaling[k]` = `ln(scaling[k])`.

However, the INPUT scaling values differ due to finding (2a) above: SPRAL
uses optimal duals from the second matching, we use suboptimal duals from
the partial first matching. So even though the Duff-Pralet formula is
implemented correctly, the matched-portion scaling it operates on is wrong.

**Classification: CORRECT DIFFERENCE** in formula, but **inputs are wrong**
due to (2a).

#### 5. How Scaling Is Applied in the Solver

**Factorization** (`solver.rs` line 279 + `numeric.rs` lines 579-583):
- Scaling is stored in **elimination order** (`solver.rs` line 244:
  `elim_scaling[i] = result.scaling[perm_fwd[i]]`).
- During frontal matrix assembly, each entry is scaled:
  `val *= s[perm_inv[orig_row]] * s[perm_inv[orig_col]]`
- This is correct: for row `r` and column `c`, the scaled entry is
  `s[elim(r)] * a[r,c] * s[elim(c)]` where `elim(r) = perm_inv[r]`.
  Since `s[elim(r)] = scaling[perm_fwd[elim(r)]] = scaling[r]`, this
  correctly computes `scaling[r] * a[r,c] * scaling[c]`.

**Solve** (`solver.rs` lines 346-372):
1. Permute RHS into elimination order.
2. Scale: `rhs_perm[i] *= scaling[i]` (pre-multiply by S).
3. Core solve: L D L^T x_perm = S * P * b.
4. Unscale: `rhs_perm[i] *= scaling[i]` (post-multiply by S).
5. Unpermute.

The overall solve computes: `x = P^T * S * (LDL^T)^{-1} * S * P * b`.
Since factorization computed `P * (S*A*S) * P^T = LDL^T`, i.e.,
`S*A*S = P^T * LDL^T * P`, so `(S*A*S)^{-1} = P^T * (LDL^T)^{-1} * P`.
Then `x = P^T * S * P^T * (LDL^T)^{-1} * P * S * P * b`... wait, this
needs careful checking.

Actually, let me re-trace. The factorization assembles `S*A*S` into fronts
using the fill-reducing permutation P. So we have `P*(S*A*S)*P^T = LDL^T`.
The solve needs `(S*A*S)^{-1} = P^T*(LDL^T)^{-1}*P`.
Then `A^{-1} = S*(S*A*S)^{-1}*S = S*P^T*(LDL^T)^{-1}*P*S`.
So `x = A^{-1}*b = S * P^T * (LDL^T)^{-1} * P * S * b`.

The solve code does:
1. `rhs_perm = P * b` (permute: `rhs_perm[new] = b[perm_fwd[new]]`)
2. `rhs_perm *= S` (scale in elimination order)
3. `rhs_perm = (LDL^T)^{-1} * rhs_perm` (core solve)
4. `rhs_perm *= S` (unscale)
5. `x = P^T * rhs_perm` (unpermute: `x[perm_fwd[new]] = rhs_perm[new]`)

This computes `x = P^T * S * (LDL^T)^{-1} * S * P * b`, which equals
`S * P^T * (LDL^T)^{-1} * P * S * b` **only if S and P^T commute**.
But S is a diagonal matrix in original index space, and the code applies
S in elimination (permuted) order. Let me re-check.

The scaling vector stored on `SparseLDLT` is `elim_scaling[i] = scaling[perm_fwd[i]]`.
So step 2 applies `rhs_perm[i] *= elim_scaling[i] = scaling[perm_fwd[i]]`.
After step 1, `rhs_perm[i] = b[perm_fwd[i]]`.
After step 2, `rhs_perm[i] = scaling[perm_fwd[i]] * b[perm_fwd[i]]`.
This is `(S*b)[perm_fwd[i]]` in the original ordering, which is `P*(S*b)`.
After step 3, `rhs_perm = (LDL^T)^{-1} * P * S * b`.
After step 4, `rhs_perm[i] *= elim_scaling[i]`, giving
`rhs_perm[i] = elim_scaling[i] * [(LDL^T)^{-1} * P * S * b]_i`.
This is S_elim * (LDL^T)^{-1} * P * S * b where S_elim = diag(elim_scaling).
After step 5, `x[perm_fwd[i]] = rhs_perm[i]`, i.e., `x = P^T * rhs_perm`.
So `x = P^T * S_elim * (LDL^T)^{-1} * P * S * b`.

Now, `S_elim` in elimination order corresponds to `S` in original order:
`(S_elim * v)[i] = elim_scaling[i] * v[i] = scaling[perm_fwd[i]] * v[i]`.
And `(P^T * S_elim * v)[perm_fwd[i]] = elim_scaling[i] * v[i] = scaling[perm_fwd[i]] * v[i]`.
Setting `k = perm_fwd[i]`, `(P^T * S_elim * v)[k] = scaling[k] * v[perm_inv[k]]`.
Meanwhile, `(S * P^T * v)[k] = scaling[k] * (P^T * v)[k] = scaling[k] * v[perm_inv[k]]`.
So `P^T * S_elim = S * P^T`. Good, they commute in this sense.

Therefore `x = S * P^T * (LDL^T)^{-1} * P * S * b = A^{-1} * b`.

**Classification: CORRECT** -- Scaling application in both factorization and
solve is mathematically correct. The symmetric scaling `S*A*S` is correctly
assembled, and the solve correctly applies `S * P^T * (LDL^T)^{-1} * P * S`.

#### Summary of Agent B Findings

| # | Issue | Classification | Impact |
|---|-------|---------------|--------|
| 1 | Raw scaling formula (`exp((u+v-cmax)/2)`) | CORRECT DIFFERENCE | None |
| 2a | **Missing second matching on full-rank submatrix** | **DEFINITE BUG** | Suboptimal duals for all matrices with incomplete matchings (dawson5, astro-ph, copter2, helm3d01, TSOPF_b162). This is the most likely root cause of MC64 being harmful. |
| 2b | Zeroing unmatched row duals (vs SPRAL's sentinel approach) | POTENTIAL BUG | May compound with 2a; benign in isolation for matched-column dual computation |
| 3 | Weight transformation (log-space costs) | CORRECT DIFFERENCE | None |
| 4 | Duff-Pralet correction formula (linear vs log space) | CORRECT DIFFERENCE | None (but operates on wrong inputs due to 2a) |
| 5 | Scaling application in factorize/solve | CORRECT | None |
| 6 | Missing `match_postproc` magnitude adjustment | CORRECT DIFFERENCE | Zero net effect for symmetric case (row/col shift cancels in sum) |
| 7 | `enforce_scaling_bound` (not in SPRAL) | CORRECT DIFFERENCE | Conservative safety net; cannot make things worse |

**Root cause hypothesis**: For structurally singular matrices, SPRAL runs
the Hungarian algorithm **twice** -- once on the full graph to find which
rows/columns participate, then again on the matched subgraph to get optimal
dual variables. We only run it once, using suboptimal duals from a partial
matching. The fix is to implement the second matching pass on the full-rank
submatrix, following SPRAL lines 688-801.

**Explains dawson5**: dawson5 has 37 unmatched rows (51500/51537). With
suboptimal duals, the scaling factors for the 51500 matched rows are wrong,
corrupting the scaled matrix. Plain METIS (no scaling) gives BE=7.45e-17
because the matrix is well-conditioned without scaling.

**Explains astro-ph**: 1029 unmatched rows (15677/16706) -- the largest
fraction of unmatched rows in our test set. Suboptimal duals from the
partial matching would severely degrade scaling quality for the matched
portion, explaining the 2.79e-3 backward error.

### Agent C: Condensation / match_order_metis (`ordering.rs` vs SPRAL `match_order.f90`)

**Scope**: Compare our condensation pipeline against SPRAL's `match_order_metis`

**Key questions**:
- How does SPRAL build the condensed graph from the matching?
- How are unmatched rows handled? (singleton nodes? merged with nearest matched pair?)
- Does SPRAL apply the matching permutation before or after METIS?
- What ordering does SPRAL produce for the condensed graph?
- How is the condensed ordering expanded back to the original dimension?
- Compare our condensation against SPRAL's on dawson5 (37 unmatched rows)

**Our files**: `src/aptp/ordering.rs` (1297 lines) — `match_order_metis` function
**SPRAL files**: `/opt/references/spral/src/match_order.f90` (632 lines)

**Status**: COMPLETE
**Findings**:

#### 1. Condensed Graph Construction

**SPRAL** (`mo_split`, lines 220-396 of `match_order.f90`):

SPRAL's pipeline has two main phases in `mo_split`:

**Phase 1 — Cycle decomposition** (lines 259-290): Walk cycles in the matching
permutation `cperm`. After decomposition, `cperm(i)` is overwritten with:
- `-2`: unmatched by MC64
- `-1`: singleton (self-matched, or leftover from odd-length cycle)
- `> 0`: partner index (2-cycle pair: `cperm(i) = j, cperm(j) = i`)

**Phase 2 — Build condensed numbering** (lines 295-303):
```fortran
k = 1
do i = 1,n
   j = cperm(i)
   if ((j .lt. i) .and. (j .gt. 0)) cycle   ! skip 2-cycle second member
   old_to_new(i) = k
   new_to_old(k) = i
   if (j .gt. 0) old_to_new(j) = k           ! pair shares same index
   k = k + 1
end do
ncomp_matched = k-1
```

**Critical observation**: The skip condition `(j < i) AND (j > 0)` only skips the
second member of 2-cycle pairs. Nodes with `cperm(i) = -1` (singletons) and
`cperm(i) = -2` (unmatched) are **NOT skipped** because `-1 > 0` and `-2 > 0` are
both FALSE. Every such node receives its own condensed index. Therefore
`ncomp_matched` includes **all** nodes: 2-cycle pairs + singletons + unmatched.

**Our code** (`split_matching_cycles`, lines 312-401 of `ordering.rs`):

We first mark unmatched indices in a separate pass:
```rust
for i in 0..n {
    if !is_matched[i] {
        partner[i] = -2;
        visited[i] = true;   // <-- prevents cycle walking
    }
}
```

Then we walk only matched indices, building condensed indices for singletons and
2-cycle pairs. Unmatched indices are **excluded from the condensed graph entirely**.
Our `condensed_dim` = singletons + two_cycle_pairs (matched nodes only).

**DEFINITE BUG**: Unmatched nodes are excluded from the condensed graph and from
METIS ordering. SPRAL includes them as singleton condensed nodes and lets METIS
optimally place them in the fill-reducing ordering.

**Phase 3 — Build condensed adjacency** (SPRAL lines 310-355):

SPRAL iterates over all `n` original nodes (skipping only 2-cycle second members),
scanning all edges and mapping row indices through `old_to_new`. The filter on
line 320 (`if (krow .gt. ncomp_matched) cycle`) is effectively a no-op because all
nodes have indices in `[1, ncomp_matched]`. Duplicate edges within a condensed
column are removed via the marker `iwork(krow) = i`.

For 2-cycle pairs, SPRAL scans BOTH columns (node `i` and partner `j`) to collect
all edges into the single condensed node. This is done at lines 325-334:
```fortran
if (j .gt. 0) then
   do klong = ptr2(j), ptr2(j+1)-1
      ...
   end do
end if
```

SPRAL then extracts the **lower triangular** part of the condensed graph
(lines 341-355) for input to its METIS wrapper, which expects lower-triangular
CSC. Our code builds the full symmetric adjacency, which is correct for
`METIS_NodeND`.

**Our code** (`build_condensed_adjacency`, lines 415-481):

Our implementation correctly handles 2-cycle pair edge merging: it iterates over
all `n` columns, and for each column that belongs to a pair, both original columns'
edges get mapped to the same condensed node via `old_to_new[col]`. Deduplication
is handled via the `marker` array plus a final `sort_unstable() + dedup()` pass.

However, all edges from/to unmatched nodes are skipped:
```rust
if decomp.partner[col] == -2 { continue; }  // line 429
if decomp.partner[row] == -2 { continue; }  // line 439
```

This is a consequence of the condensed graph not including unmatched nodes.

**CORRECT** for matched nodes; **DEFINITE BUG** for unmatched nodes.

#### 2. Unmatched Row Handling (CRITICAL)

**SPRAL**: Unmatched rows become **singleton condensed nodes** that participate
fully in METIS ordering. METIS sees their connectivity to other condensed nodes
and places them at positions that minimize fill. The condensed graph includes all
edges from/to unmatched nodes.

**Our code**: Unmatched rows are excluded from METIS entirely. In `expand_ordering`
(lines 517-522), they are unconditionally appended at the end of the elimination
ordering in natural index order:
```rust
for i in 0..n {
    if decomp.partner[i] == -2 {
        fwd.push(i);
    }
}
```

**Impact on dawson5**: The 37 unmatched rows are placed at elimination positions
`[51500, 51501, ..., 51536]` — the very last positions. If any of these rows have
connections to many other rows (which is likely for a structurally near-singular
matrix), eliminating them last creates catastrophic fill-in. METIS would have
interleaved them via nested dissection to minimize fill.

**Impact on astro-ph**: With 1029 unmatched rows (6% of 16706), the impact is even
more dramatic — over a thousand rows are forced to the worst possible ordering
positions.

**DEFINITE BUG**: This is the most likely primary cause of dawson5 failing with
MatchOrderMetis while passing with plain METIS (BE=7.45e-17).

#### 3. Condensed Ordering Expansion

**SPRAL** (lines 376-395):
```fortran
! Step 1: invert METIS ordering to get iwork[position] = condensed_node
do i = 1, ncomp
   j = order(i)       ! order(i) = position of condensed node i
   iwork(j) = i       ! iwork[position] = condensed_node
end do

! Step 2: walk positions in order, map to original indices
k = 1
do i = 1, ncomp
   j = new_to_old( iwork(i) )   ! first original index for condensed node
   order(j) = k                  ! assign elimination position
   k = k + 1
   if (cperm(j) .gt. 0) then    ! 2-cycle pair has a partner
      j = cperm(j)
      order(j) = k              ! partner gets next consecutive position
      k = k + 1
   end if
end do
```

Key properties of SPRAL's expansion:
- ALL `ncomp` condensed nodes participate (including unmatched-as-singletons)
- 2-cycle pairs (`cperm(j) > 0`) get consecutive positions
- Singletons (`cperm(j) = -1`) and unmatched-as-singletons (`cperm(j) = -2`)
  get single positions
- After the loop, `k = n + 1` (all original indices assigned)
- Output is `order(j) = elimination_position_of_original_index_j` (the inverse perm)

**Our code** (`expand_ordering`, lines 495-541):

```rust
// Walk matched condensed nodes
for &cnode in &inv_order {
    let orig = decomp.new_to_old[cnode];
    fwd.push(orig);
    if decomp.partner[orig] >= 0 {
        fwd.push(decomp.partner[orig] as usize);
    }
}
// Append unmatched at the end
for i in 0..n {
    if decomp.partner[i] == -2 {
        fwd.push(i);
    }
}
```

Our expansion correctly handles matched nodes but unconditionally appends unmatched
nodes at the end rather than letting METIS decide their position.

**DEFINITE BUG** (same root cause as finding #2).

#### 4. Matching Permutation and `is_matched` Interaction

**SPRAL's approach**: The matching permutation returned by `mo_match` uses the
sentinel value `-1` for unmatched rows (`perm(1:n) = -1` initially, then
`perm(new_to_old(i)) = new_to_old(cperm(i))` for matched rows). In `mo_split`,
cycle walking checks `cperm(j) == -1` to detect unmatched nodes and immediately
marks them as `iwork(j) = -2`. The cycle walker never follows the permutation
of an unmatched row.

**Our approach**: `build_singular_permutation` (matching.rs lines 301-334) assigns
unmatched rows to arbitrary free columns to create a valid permutation (every row
maps to some column). This is necessary because `Perm<usize>` requires a valid
permutation. Then `is_matched[i] = is_row_matched[i] || is_col_matched[i]` is
used to flag which indices are truly matched.

**POTENTIAL BUG**: The `is_matched` uses an OR of row-matched and column-matched.
Consider an index `i` where:
- Row `i` is NOT matched (no augmenting path reached it)
- Column `i` IS matched (some other row `r` has `row_match[r] = i`)

Then `is_col_matched[i] = true`, so `is_matched[i] = true`. In
`split_matching_cycles`, index `i` is treated as matched, but its permutation
entry `matching_fwd[i]` was filled in by `build_singular_permutation` with an
arbitrary free column. The cycle walker follows this fake entry, creating a
spurious cycle that includes index `i` and whatever it was arbitrarily assigned to.

For symmetric matrices, `row_match[r] = i` means "row r is matched to column i",
and symmetrically column i appears in the matching. But row i might NOT be matched
to any column. In this case, index i should be treated as unmatched for the
condensation pipeline, but `is_matched[i] = true` prevents this.

The `is_matched` field was originally designed for the Duff-Pralet scaling correction
(where the OR semantics are arguably correct — index i participates in the matching
ecosystem if it is the target column of another row's match). But reusing it for
cycle decomposition introduces the bug described above.

**Severity**: For dawson5 with 37 unmatched rows, it is possible that some of the
37 unmatched row-indices have their column matched by other rows, causing them to
be treated as matched and creating spurious 2-cycles/singletons in the condensed
graph. This would compound the "append at end" bug — not only are truly unmatched
nodes placed at the end, but some notionally "unmatched" nodes are placed into fake
2-cycles, corrupting the condensed graph structure.

**POTENTIAL BUG** — could cause incorrect cycle decomposition for matrices where
the matching is asymmetric (row i unmatched but column i matched).

#### 5. Edge Cases

**Empty matching (all rows unmatched)**:
- SPRAL: All nodes become singleton condensed nodes (cperm = -2 for all). The
  condensed graph is the original graph. METIS orders the full graph. Equivalent
  to plain METIS.
- Our code: `condensed_dim = 0`, METIS is skipped (trivial identity ordering for
  dimension 0), all nodes are appended at the end in natural index order. The
  output is the identity ordering.
- **DEFINITE BUG**: Identity ordering instead of fill-reducing ordering. No
  fill reduction at all. (Low severity: unlikely in practice since MC64 usually
  finds a substantial matching.)

**Perfect matching (all rows matched)**:
- SPRAL: All nodes are singletons or 2-cycle pairs. No unmatched nodes.
- Our code: Same behavior. No unmatched nodes to mishandle.
- **CORRECT**: Both implementations produce equivalent results.

**Self-matches (i matched to i)**:
- SPRAL: After cycle decomposition, `cperm(i) = -1` (singleton). Gets own
  condensed node. During expansion, `cperm(j) > 0` is FALSE, so gets single
  position.
- Our code: `matching_fwd[i] == i` detected as singleton, `partner[i] = -1`.
  Gets own condensed node. During expansion, `partner[orig] >= 0` is FALSE,
  so gets single position.
- **CORRECT DIFFERENCE**: Both handle self-matches identically.

#### 6. METIS Input Format

**SPRAL** (lines 341-355): After building the full condensed adjacency, SPRAL
extracts only the **lower triangular** part for its METIS wrapper:
```fortran
do i = 1, ncomp
   do k = j1, j2-1
      krow = row3(k)
      if (krow .lt. i) cycle   ! skip upper triangle
   end do
end do
```
This is because SPRAL's `metis_order` wrapper (`spral_metis_wrapper`) presumably
expects lower-triangular input and symmetrizes internally.

**Our code**: We pass a full symmetric adjacency to `METIS_NodeND`, which expects
the complete adjacency list for each vertex (both upper and lower triangles).

**CORRECT DIFFERENCE**: Both approaches are valid for their respective METIS
interfaces. `METIS_NodeND` requires full symmetric input; SPRAL's wrapper handles
the conversion.

#### 7. Cycle Decomposition Algorithm Comparison

Both SPRAL and our code use the same cycle-walking approach for decomposition.
The key steps are:

1. Walk each cycle in the matching permutation
2. Pair consecutive members into 2-cycles
3. Odd-length cycles produce one leftover singleton

SPRAL's cycle walker (lines 261-286):
```fortran
j = i
do
   if (cperm(j) .eq. -1) then     ! unmatched
      iwork(j) = -2; exit
   else if (cperm(j) .eq. i) then ! finished cycle
      iwork(j) = -1; exit         ! leftover singleton
   end if
   jj = cperm(j)
   iwork(j) = jj                  ! pair j with jj
   iwork(jj) = j
   j = cperm(jj)
   if (j .eq. i) exit             ! even-length cycle complete
end do
```

Our cycle walker (lines 350-391):
```rust
let mut j = i;
loop {
    let k = matching_fwd[j];
    if visited[k] || k == i {
        // leftover singleton
        if !visited[j] { partner[j] = -1; ... }
        break;
    }
    // pair j and k
    partner[j] = k as isize;
    partner[k] = j as isize;
    let next = matching_fwd[k];
    if next == i { break; }
    j = next;
    if visited[j] { break; }
}
```

The core algorithm is equivalent except for unmatched detection:
- SPRAL checks `cperm(j) == -1` (sentinel for unmatched)
- Our code uses pre-marking via `visited[i] = true` for `!is_matched[i]`

**CORRECT DIFFERENCE** for the cycle walking algorithm itself (aside from the
`is_matched` OR bug discussed in finding #4).

#### 8. Summary of Findings

| # | Finding | Classification | Severity | Impact |
|---|---------|---------------|----------|--------|
| 1 | Unmatched nodes excluded from condensed graph | **DEFINITE BUG** | **HIGH** | Unmatched nodes invisible to METIS; ordering quality degraded |
| 2 | Unmatched nodes appended at end of ordering | **DEFINITE BUG** | **HIGH** | Worst-case fill-in for unmatched rows; primary cause of dawson5 failure |
| 3 | Expansion only covers matched condensed nodes | **DEFINITE BUG** | **HIGH** | Direct consequence of bugs 1-2 |
| 4 | `is_matched` OR logic creates spurious cycle entries | **POTENTIAL BUG** | **MEDIUM** | Corrupts cycle decomposition when row-unmatched but column-matched |
| 5a | Empty matching gives identity ordering | **DEFINITE BUG** | LOW | Edge case; unlikely in practice |
| 5b | Perfect matching | CORRECT | N/A | No issue |
| 5c | Self-match handling | CORRECT | N/A | Identical to SPRAL |
| 6 | METIS input format (full vs lower-tri) | CORRECT DIFFERENCE | N/A | Both valid |
| 7 | Cycle decomposition algorithm | CORRECT | N/A | Aside from is_matched issue |

#### 9. Recommended Fix

The fix requires changes to three functions:

**Fix 1 — `split_matching_cycles`**: Include unmatched nodes as singleton condensed
nodes. Instead of:
```rust
if !is_matched[i] {
    partner[i] = -2;
    visited[i] = true;  // WRONG: prevents condensed index assignment
}
```
Change to still mark them as unmatched (`partner[i] = -2`) but DO assign them a
condensed index and add them to `new_to_old`. Alternatively, use a new marker value
(e.g., `-3`) to distinguish "unmatched-singleton-in-condensed-graph" from
"truly excluded". But the simplest approach is to give them the same treatment as
singletons (`partner[i] = -1`) but still track them separately for diagnostics:
```rust
if !is_matched[i] {
    partner[i] = -2;       // still tagged as unmatched for diagnostics
    old_to_new[i] = condensed_idx;
    new_to_old.push(i);
    condensed_idx += 1;
    visited[i] = true;
}
```

**Fix 2 — `build_condensed_adjacency`**: Remove the skip for unmatched nodes:
```rust
// DELETE: if decomp.partner[col] == -2 { continue; }
// DELETE: if decomp.partner[row] == -2 { continue; }
```
Since unmatched nodes now have condensed indices, their edges will be correctly
mapped and included.

**Fix 3 — `expand_ordering`**: Remove the "append unmatched at end" block:
```rust
// DELETE the block:
// for i in 0..n {
//     if decomp.partner[i] == -2 {
//         fwd.push(i);
//     }
// }
```
Since ALL nodes (including unmatched) now have condensed indices, the main
expansion loop `for &cnode in &inv_order` will already emit them.

**Fix 4 (optional, addresses finding #4)**: Change `is_matched` from OR to
`is_row_matched` only, or better yet, use a sentinel in the matching permutation
(matching_fwd[i] = usize::MAX for unmatched rows) and detect it directly in
`split_matching_cycles` instead of relying on a separate `is_matched` vector.

After these fixes, the condensed graph for dawson5 would have:
- ~25732 condensed nodes for 2-cycle pairs (51500 - 37 unmatched = 51463 matched,
  but some are singletons; exact count depends on cycle structure)
- Plus ~37 singleton condensed nodes for unmatched rows
- METIS would place all of them optimally

### Agent D: Empirical Isolation (matching vs scaling vs condensation)

**Scope**: Run failing matrices with different MC64 configurations to isolate the bug

**Key experiments**:
1. **MatchOrderMetis vs plain METIS** on all failing matrices (dawson5 already known)
2. **MC64 matching perm only** (apply matching permutation, skip scaling, run METIS)
3. **MC64 scaling only** (keep scaling, use identity matching perm, run METIS)
4. **MC64 matching quality**: compute matching weight sum, compare unmatched count
5. **Condensation bypass**: apply MC64 matching, then run METIS on full matrix (not condensed)

**Our files**: `src/aptp/ordering.rs`, `src/aptp/solver.rs`, `src/aptp/matching.rs`
**Test matrices**: dawson5 (smoking gun), astro-ph, copter2, helm3d01, TSOPF_b162, TSOPF_b39

**Status**: COMPLETE
**Findings**:

#### Experiment Setup

Created `examples/mc64_isolation.rs` with four experiment configurations for each failing matrix:

| Config | Ordering | Scaling | Purpose |
|--------|----------|---------|---------|
| MatchOrderMetis (ord+scale) | MC64+condensed METIS | MC64 scaling | Full MatchOrderMetis (baseline) |
| MatchOrder ord, NO scale | MC64+condensed METIS | None | Isolates condensed ordering effect |
| Metis ord + MC64 scale | Plain METIS | MC64 scaling | Isolates scaling effect |
| Metis ord, NO scale (reference) | Plain METIS | None | Reference (no MC64 involvement) |

Also computed MC64 quality diagnostics per matrix: matching cardinality, cycle structure,
scaling factor distribution, diagonal dominance after scaling, and condensation ratio.

#### 1. dawson5: THE SMOKING GUN

**Verdict: MC64 scaling is the primary cause, not the condensed ordering.**

| Configuration | BE | Status | delays | max_front |
|--------------|---:|:------:|-------:|----------:|
| MatchOrderMetis (ord+scale) | 6.91e-4 | FAIL | 91 | 916 |
| MatchOrder ord, NO scale | 4.14e-4 | FAIL | 94 | 916 |
| **Metis ord + MC64 scale** | **7.50e-4** | **FAIL** | 199 | 565 |
| Metis ord, NO scale (reference) | 7.45e-17 | PASS | 198 | 565 |

Critical observation: **Metis ord + MC64 scale** (row 3) fails with BE=7.50e-4 while
the identical ordering WITHOUT scaling (row 4) passes with BE=7.45e-17. The scaling
factors alone degrade backward error by 13 orders of magnitude.

The condensed ordering is a secondary issue: MatchOrder ord WITHOUT scale (row 2, BE=4.14e-4)
still fails but with a different ordering. Both orderings fail when scaling is applied.

MC64 diagnostics for dawson5:
- 51500/51537 matched (99.9%), 0 unmatched rows (all `is_matched` = true)
- **634 longer-cycles** (not just 2-cycles!) — indicates the matching quality is poor
  for the symmetric case
- 4105 singletons, 20721 two-cycles → condensation ratio = 54.5%
- Scaling range: [0.110, 14.3], ratio=130 — moderate, no extreme values
- **Only 4809/51537 (9.3%) diagonal-dominant after scaling, 46728 (90.7%) weak**
- Scaling bound: max |s_i * a_ij * s_j| = 1.0 (correctly enforced)
- No scaling violations

The 634 longer-cycles are a red flag: for a symmetric matrix, the optimal matching
should decompose purely into singletons and 2-cycles. Longer cycles indicate the
matching is operating on an asymmetric cost graph (due to how entries are expanded
from upper triangle), and the resulting dual variables do not produce good symmetric
scaling. The suboptimal duals (Agent A Finding #3: missing re-matching) are the root
cause — the scaling factors computed from these duals make 90.7% of the diagonal weak,
destroying the natural diagonal dominance that plain METIS exploits.

#### 2. vibrobox and crystk02: MC64 Scaling is ESSENTIAL

These two matrices show the opposite pattern — MC64 scaling FIXES backward error:

**vibrobox** (n=12328):

| Configuration | BE | Status | delays | max_front |
|--------------|---:|:------:|-------:|----------:|
| MatchOrderMetis (ord+scale) | 1.23e-18 | PASS | 0 | 937 |
| MatchOrder ord, NO scale | 4.22e-4 | FAIL | 7168 | 1704 |
| **Metis ord + MC64 scale** | **1.23e-18** | **PASS** | 0 | 937 |
| Metis ord, NO scale (reference) | 4.22e-4 | FAIL | 7168 | 1704 |

**crystk02** (n=13965):

| Configuration | BE | Status | delays | max_front |
|--------------|---:|:------:|-------:|----------:|
| MatchOrderMetis (ord+scale) | 3.42e-18 | PASS | 0 | 756 |
| MatchOrder ord, NO scale | 1.40e-5 | FAIL | 6 | 756 |
| **Metis ord + MC64 scale** | **3.42e-18** | **PASS** | 0 | 756 |
| Metis ord, NO scale (reference) | 1.40e-5 | FAIL | 6 | 756 |

For both matrices:
- Perfect matching (100%), 0 longer-cycles, all singletons
- vibrobox: scaling range [2.0e-6, 2.8e-3], ratio=1355 — significant magnitude variation
- crystk02: scaling range [1.0e5, 5.7e5], ratio=5.6 — uniform but large
- Scaling eliminates ALL delayed pivots (7168 → 0 for vibrobox, 6 → 0 for crystk02)
- Scaling also reduces max_front (1704 → 937 for vibrobox)
- The **ordering is irrelevant** — condensed ordering = plain METIS (no 2-cycles to condense)
- **Scaling alone** accounts for the entire improvement

These matrices have well-structured matchings (all singletons, no longer-cycles), so
the dual variables from the first matching are optimal. The scaling correctly normalizes
the matrix entries to improve diagonal dominance.

#### 3. TSOPF matrices: Scaling Helps, But Not Enough

**TSOPF_b162** (n=10798):

| Configuration | BE | Status | delays | max_front |
|--------------|---:|:------:|-------:|----------:|
| MatchOrderMetis (ord+scale) | 2.34e-7 | FAIL | 485 | 297 |
| MatchOrder ord, NO scale | 5.57e-7 | FAIL | 9139 | 1089 |
| Metis ord + MC64 scale | 2.14e-7 | FAIL | 7917 | 322 |
| Metis ord, NO scale (reference) | 1.38e-6 | FAIL | 15481 | 1144 |

**TSOPF_b39** (n=28216):

| Configuration | BE | Status | delays | max_front |
|--------------|---:|:------:|-------:|----------:|
| MatchOrderMetis (ord+scale) | 5.98e-8 | FAIL | 1798 | 322 |
| MatchOrder ord, NO scale | 2.19e-6 | FAIL | 15858 | 1057 |
| **Metis ord + MC64 scale** | **1.19e-10** | **~PASS** | 135661 | 285 |
| Metis ord, NO scale (reference) | 1.48e-6 | FAIL | 325969 | 1012 |

For TSOPF_b39, the best result is **Metis ord + MC64 scale** (BE=1.19e-10), which
nearly passes the 5e-11 threshold. This suggests:
- MC64 scaling helps (1.48e-6 → 1.19e-10 with Metis ordering)
- The condensed ordering HURTS compared to plain METIS when scaling is applied
  (5.98e-8 with condensed vs 1.19e-10 with plain METIS)
- Both TSOPF matrices have ~50% condensation ratio and massive 2-cycle counts
- 14 and 2 longer-cycles respectively — nearly all 2-cycles, so matching is good
- Very wide scaling range (TSOPF_b162: ratio=2.5e5; TSOPF_b39: ratio=6.7e4)

The TSOPF matrices are hard-indefinite with many zero diagonals (7874 and 21098 missing
respectively). MC64 scaling is genuinely needed but the condensed ordering introduces
fill that the large-front factorization cannot handle accurately.

#### 4. Large-Front Failures: Factorization Quality Limit

**copter2** (n=55476), **helm3d01** (n=32226), **sparsine** (n=50000):

All fail with ALL configurations at ~1e-3 backward error. The four-way isolation shows
no significant difference between configurations — BE stays in the [6.87e-4, 9.42e-4]
range for copter2, [8.12e-4, 1.24e-3] for helm3d01, [1.91e-3, 2.53e-3] for sparsine.

Common characteristic: **very large fronts** (copter2: 1218-1801, helm3d01: 1206-1424,
sparsine: 11131-11135). These are factorization quality issues — the BLAS-2 inner
kernel accumulates too much rounding error over hundreds/thousands of rank-1 updates
in a single tile. The BLAS-3 inner blocking (Phase 8.1's unused `ib` parameter) would
help here by limiting rank-1 updates to 32-column blocks with GEMM between them.

**astro-ph** (n=16706): All-zero diagonal, fails with all configurations at ~1e-3.
The all-zero diagonal is a fundamental challenge — every pivot must be 2x2, generating
massive delays (9K-34K) regardless of ordering/scaling.

#### 5. Pattern Classification

The empirical isolation reveals three distinct categories:

**Category A: Scaling-harmful matrices** (dawson5)
- MC64 scaling CAUSES the failure (7.45e-17 → 7.50e-4)
- Root cause: suboptimal dual variables from partial matching (Agent A Finding #3)
  produce bad scaling that destroys natural diagonal dominance
- Fix: implement SPRAL's second-matching pass for structurally singular matrices

**Category B: Scaling-essential matrices** (vibrobox, crystk02)
- MC64 scaling FIXES the failure (4.22e-4 → 1.23e-18)
- Matrices have well-structured matchings (all singletons, no longer-cycles)
- Scaling correctly normalizes entry magnitudes

**Category C: Factorization-limited matrices** (copter2, helm3d01, sparsine, astro-ph)
- Neither ordering nor scaling helps — BE stays at ~1e-3 with all configurations
- Root cause: large fronts (1000-11000) with BLAS-2 inner kernel
- Fix: restore BLAS-3 inner blocking (`ib=32`) in `factor_inner`

**Category D: Scaling-helps-but-insufficient** (TSOPF_b162, TSOPF_b39)
- Scaling improves BE by 1-4 orders of magnitude but not enough to pass
- Condensed ordering may be counterproductive (TSOPF_b39: 5.98e-8 condensed vs 1.19e-10 plain)
- May benefit from both better scaling (second matching) AND BLAS-3 inner blocking

#### 6. Key Finding: Scaling Is More Important Than Ordering

The most surprising result: **for dawson5, the condensed ordering is NOT the primary
problem — the scaling is.** The prior hypothesis (Agent C) was that unmatched nodes
being appended at the end of the ordering was the root cause. But the empirical data
shows:

- "MatchOrder ord, NO scale" (condensed ordering, no scaling): BE=4.14e-4 (FAIL)
- "Metis ord + MC64 scale" (plain ordering, MC64 scaling): BE=7.50e-4 (FAIL)
- "Metis ord, NO scale" (plain ordering, no scaling): BE=7.45e-17 (PASS)

The condensed ordering alone (without scaling) produces a worse ordering than plain
METIS — this is real and consistent with Agent C's findings about unmatched nodes
being excluded from METIS. But the scaling is the MORE DAMAGING component: even with
the optimal METIS ordering, MC64 scaling degrades dawson5 from machine precision to
1e-4.

This reorders the priority of fixes:
1. **First**: Fix MC64 re-matching for structurally singular case (Agent A Finding #3)
2. **Second**: Fix condensed graph to include unmatched nodes (Agent C Findings #1-3)
3. **Third**: Restore BLAS-3 inner blocking for large-front accuracy (separate issue)

#### 7. Diagnostic Tool

The diagnostic example is at `examples/mc64_isolation.rs`. Run with:
```
cargo run --example mc64_isolation --release
```

It tests 9 matrices covering all failure categories and produces per-matrix:
- Baseline backward error (Metis vs MatchOrderMetis)
- MC64 quality diagnostics (matching, cycles, scaling distribution, diagonal dominance)
- Four-way scaling isolation (ordering x scaling)

## Resolution

### Root Cause Summary

Three independent issues contribute to backward error failures:

1. **Missing MC64 re-matching** (Agent A Finding #3, Agent B Finding #2a):
   For structurally singular matrices (matched < n), SPRAL runs the Hungarian algorithm
   twice — once on the full graph, then again on the matched subgraph to get optimal
   duals. Our code skips the second pass. Suboptimal duals produce scaling that actively
   destroys diagonal dominance on matrices like dawson5 (90.7% of diagonal entries become
   weak after scaling, vs machine-precision accuracy without scaling).

2. **Unmatched nodes excluded from condensed graph** (Agent C Findings #1-3):
   SPRAL includes unmatched nodes as singleton condensed nodes in METIS. Our code
   excludes them and appends at the end, forcing worst-case fill-in. This degrades
   the condensed ordering quality.

3. **BLAS-2 inner kernel for large fronts** (separate from MC64):
   Fronts with 1000+ columns accumulate too much rounding error with per-column rank-1
   updates. SPRAL uses BLAS-3 inner blocking (ib=32) with GEMM/TRSM between 32-column
   blocks. Our code has the parameter but does not use it.

### Priority of Fixes

| Priority | Fix | Matrices Affected | Expected Impact |
|:--------:|-----|-------------------|-----------------|
| 1 | MC64 re-matching for structurally singular | dawson5, astro-ph, copter2, helm3d01, TSOPF x2 | dawson5: 7.50e-4 → ~1e-17 |
| 2 | Include unmatched nodes in condensed graph | Same as above | Better ordering quality |
| 3 | BLAS-3 inner blocking (ib=32) | copter2, helm3d01, sparsine, astro-ph | ~1e-3 → ~1e-11 (matching SPRAL) |
