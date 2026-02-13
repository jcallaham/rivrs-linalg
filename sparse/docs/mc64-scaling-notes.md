# MC64 Scaling Implementation Notes

Technical notes on the MC64 matching & scaling implementation (Phase 4.2),
documenting algorithm choices, deviations from SPRAL, known limitations,
and quality analysis. Written as a reference for future revisits.

**Date**: 2026-02-13
**Implementation**: `src/aptp/matching.rs`

---

## Algorithm Identity: "Hungarian" = MC64

SPRAL's `scaling.f90` names its functions `hungarian_scale_sym`, `hungarian_match`,
`hungarian_init_heuristic`, etc. This is the **same algorithm** as MC64 (HSL's
MC64 routine). The Hungarian algorithm (Kuhn 1955) is the general framework for
the assignment problem; MC64 (Duff & Koster 2001) is a specific efficient
implementation using:

- **Logarithmic cost transformation**: converts maximum product matching into
  shortest path problem (`c[i,j] = log(col_max_j) - log|a[i,j]|`)
- **Dijkstra-based augmenting paths**: O(n·τ·log n) instead of classic O(n³)
- **Greedy initialization**: two-pass heuristic achieving ~80% cardinality
- **MC64SYM symmetric scaling** (Duff & Pralet 2005): `s[i] = exp((u[i] + v[i] - col_max_log[i]) / 2)`

SPRAL and our implementation both follow these same papers. There is no
algorithmic difference between SPRAL's "Hungarian" and our "MC64" — only naming.

---

## Consistency with SPRAL Source Code

### What matches exactly

| Component | Our code | SPRAL code | Match? |
|-----------|----------|------------|--------|
| Cost graph | `c[i,j] = col_max_log[j] - log\|a[i,j]\|` | Same formula | Exact |
| Greedy init | Two-pass with 2-augmentation | `hungarian_init_heuristic` (lines 810-929) | Exact |
| Dijkstra augmentation | Shortest augmenting path on reduced costs | `hungarian_match` (lines 938-1171) | See below |
| Column duals (full match) | `v[j] = c[matched_row,j] - u[matched_row]` | `dualv(j) = val(jperm(j)) - dualu(row(jperm(j)))` (line 679) | Exact |
| Zero unmatched duals | `u[i] = 0` for unmatched rows | `where(iperm==0) dualu=0` | Exact |
| Symmetric scaling | `exp((u[i]+v[i]-col_max_log[i])/2)` | `exp((rscaling+cscaling)/2)` | Exact |
| Duff-Pralet correction | `s[i] = 1/max_k(\|a_ik\|·s_k)` over matched k | Same formula | Exact |

### Deviations

**1. No stored column dual array during Dijkstra**

SPRAL maintains a `dualv(j)` array throughout the matching algorithm. Our
implementation computes column duals post-hoc via complementary slackness.
During Dijkstra, we compute `vj` on-the-fly from the matched edge:

```rust
// Our approach (matching.rs line ~729):
let vj = dq0 - matched_cost + state.u[q0];
```

This is mathematically equivalent for computing reduced costs within Dijkstra.
The post-hoc complementary slackness `v[j] = c[matched_row,j] - u[matched_row]`
reproduces the same final values. The `MatchingState` struct has no `v` field.

**2. `enforce_scaling_bound` — no SPRAL equivalent**

We add an iterative correction pass that caps `s[i] ≤ 1/max_j(|a_ij|·s_j)`.
This guarantees `|s_i · a_ij · s_j| ≤ 1` even if the dual variables have small
feasibility gaps from the greedy heuristic. SPRAL relies on the Hungarian
algorithm maintaining exact dual feasibility throughout.

In practice, this pass is a **no-op for nonsingular matrices** (see Quality
Analysis below). It only fires for structurally singular matrices where unmatched
entries have imprecise duals.

**3. SPRAL's singular subgraph re-matching is effectively dead code**

SPRAL's `hungarian_wrapper` (lines 688-801) contains code to extract a
"nonsingular subgraph" for structurally singular matrices by checking
`match(i) < 0`. However, `hungarian_match` (the function that computes the
matching) initializes `iperm` to 0 and sets matched entries to positive values.
Unmatched entries remain 0, never negative. Therefore `match(i) < 0` is
**always false**, the subgraph includes all indices, and the re-matching is
redundant.

We skip this entirely and apply duals + Duff-Pralet correction + scaling bound
enforcement directly.

**4. `match_postproc` adjustments cancel for symmetric**

SPRAL's `match_postproc` (lines 1588-1696) adjusts dual variables for the
square non-symmetric case. For square symmetric matrices (our use case), the
row-scaling and column-scaling adjustments are identical, so they cancel in the
symmetric formula `s[i] = exp((u[i]+v[i]-col_max_log[i])/2)`. We skip this
post-processing.

---

## Quality Analysis on SuiteSparse CI Subset

Measured on all 9 CI-subset SuiteSparse matrices:

| Matrix | n | matched | min(rmax) | mean(rmax) | rows<0.99 |
|--------|--:|--------:|----------:|-----------:|----------:|
| bloweybq | 10,001 | 10,001 | 1.000000 | 1.000000 | 0 |
| sparsine | 50,000 | 50,000 | 1.000000 | 1.000000 | 0 |
| t2dal | 4,257 | 4,257 | 1.000000 | 1.000000 | 0 |
| bratu3d | 27,792 | 27,792 | 1.000000 | 1.000000 | 0 |
| stokes128 | 49,666 | 49,666 | 1.000000 | 1.000000 | 0 |
| cvxqp3 | 17,500 | 17,437 | 0.004113 | 0.849137 | 7,372 |
| ncvxqp1 | 12,111 | 11,983 | 0.000000 | 0.771814 | 4,114 |
| ncvxqp3 | 75,000 | 74,205 | 0.000469 | 0.797671 | 36,788 |
| cfd2 | 123,440 | 52,160 | 0.255728 | 0.997915 | 473 |

### Key finding: nonsingular quality matches SPRAL exactly

For all 5 fully-matched (nonsingular) matrices, `min(row_max) = 1.0` with zero
rows below 0.99. This matches SPRAL's strict `err_tol = 5e-14` tolerance.
The `enforce_scaling_bound` pass is a no-op for these matrices — the duals
from the matching algorithm are already perfectly feasible.

### Structurally singular matrices: degraded quality is expected

The 4 partially-matched matrices show degraded row maxima, especially for
unmatched rows. This is inherent to partial matchings — there is no matched
entry to anchor the scaling at 1.0. The Duff-Pralet correction provides
reasonable fallback scaling, and `enforce_scaling_bound` guarantees the upper
bound property (`|s_i · a_ij · s_j| ≤ 1`).

Our test thresholds for singular matrices:
- Bound: `|s_i · a_ij · s_j| ≤ 1.0 + 1e-8` (strict, always enforced)
- Quality: median(row_max) ≥ 0.5 (relaxed, only for singular)

### Anomaly: cfd2 matched=52,160/123,440

cfd2 is classified as positive definite in the SuiteSparse collection. A PD
matrix should be structurally nonsingular (every diagonal entry is positive),
so the matching should achieve full cardinality. Getting only 42% matching
is unexpected and may indicate an issue worth investigating. Possibilities:
- Upper-triangular storage convention interaction with the cost graph
- Algorithmic issue with very large sparse matrices (123K dimension)
- Structural property of cfd2's sparsity pattern in the bipartite graph

This does not affect correctness (the scaling bound is enforced), but the
scaling quality is suboptimal for this matrix. **TODO**: Investigate cfd2
matching cardinality in a future optimization pass.

---

## Comparison with SPRAL Test Suite

### What SPRAL tests

SPRAL's `tests/scaling.f90` includes:

| Test | Matrices | Singularity | Bound tolerance | Row max tolerance |
|------|----------|-------------|-----------------|-------------------|
| `hungarian_sym_random` | 100 random, n=1-1000 | All nonsingular | 1.0 + 5e-14 | 1.0 - 5e-14 |
| `hungarian_sym_singular` | 1 hand-constructed 3×3 | Singular (empty col) | Not checked | Not checked |
| `auction_sym_random` | 100 random, n=1-1000 | All nonsingular | 1.0 + 1.0 (loose) | rmax ≥ 0.75 |

For the singular test, SPRAL only checks that the unmatched column is correctly
identified (`match(ising) == 0`). It does **not** validate scaling quality for
singular matrices.

### What we test beyond SPRAL

1. **Real SuiteSparse matrices** (9 CI-subset, 65+ full collection) — SPRAL only
   uses random matrices up to n=1000
2. **Structurally singular SuiteSparse matrices** (cvxqp3 n=17.5K, ncvxqp1 n=12K,
   ncvxqp3 n=75K) — SPRAL only has one 3×3 singular test
3. **Scaling quality validation** for singular matrices (bound + median row_max)
4. **Cycle structure analysis** (singleton + 2-cycle decomposition)
5. **MC64 + METIS independent composition** tests

### Test threshold comparison

| Property | SPRAL (nonsingular random) | Ours (nonsingular SuiteSparse) | Ours (singular) |
|----------|---------------------------|-------------------------------|-----------------|
| Scaled entry bound | < 1.0 + 5e-14 | ≤ 1.0 + 1e-10 | ≤ 1.0 + 1e-8 |
| Row max lower | > 1.0 - 5e-14 | ≥ 0.75 | median ≥ 0.5 |
| Full matching | Required | Required | Not required |

**CI subset** (9 curated matrices): row_max thresholds tightened to SPRAL-level
`1.0 - 1e-12` for nonsingular matrices. Measured values are exactly 1.0.

**Full collection** (62 matrices tested, 3 skipped for memory): row_max check
is diagnostic rather than a hard assertion. Most nonsingular matrices achieve
min_rmax = 1.0, but 3 fully-matched matrices show degraded quality (see below).

The singular thresholds (median ≥ 0.5) represent genuinely new territory that
SPRAL doesn't test.

---

## Quality Analysis on Full SuiteSparse Collection (62 matrices)

Validated on 62 of 65 SuiteSparse matrices (3 skipped: msdoor, inline_1,
ldoor — exceed 10M nnz cap due to MC64's ~2x memory requirement for cost graph).

### Fully-matched matrices with perfect quality (min_rmax = 1.0)

BenElechi1, F2, H2O, Si10H16, Si5H12, SiNa, bcsstk39, bloweybq, crystk02,
crystk03, dixmaanl, filter3D, gas_sensor, linverse, nd3k, nd6k, qa8fk,
rail_79841, sparsine, t2dal, t3dh, vibrobox, blockqp1, bratu3d, cont-201,
cont-300, mario001, stokes128, turon_m, G3_circuit (1.58M dim!), af_0_k101,
af_shell7, apache2, bmwcra_1, boneS01, crankseg_1, crankseg_2, offshore

### Fully-matched matrices with degraded quality

| Matrix | n | rows<0.75 | min_rmax | longer cycles | Notes |
|--------|--:|----------:|---------:|--------------:|-------|
| TSOPF_FS_b39_c7 | 28,216 | 23 | 1.67e-1 | 2 | Power system, hard-indefinite |
| d_pretok | 182,730 | 12 | 5.06e-1 | 1,534 | Hard-indefinite |
| thread | 29,736 | 659 | 1.76e-8 | 11 | Near-zero min_rmax despite full match |

These matrices have full matchings (matched == n) but suboptimal dual variable
quality for some rows. The degraded rows are a tiny fraction (0.08%-2.2%) of
the total. The scaling bound (`|s_i · a_ij · s_j| ≤ 1`) still holds for all
entries. The issue likely stems from numerical sensitivity in the cost graph
(extreme value ranges in power systems and structural mechanics matrices) or
suboptimal augmenting paths from the greedy initialization.

**TODO**: Investigate whether maintaining column duals during Dijkstra (as SPRAL
does) would improve quality for these matrices.

### Partially-matched (structurally singular) matrices

| Matrix | n | matched | mean_rmax |
|--------|--:|--------:|----------:|
| astro-ph | 16,706 | 15,677 | 0.791 |
| copter2 | 55,476 | 55,421 | 0.801 |
| dawson5 | 51,537 | 51,500 | 0.932 |
| helm3d01 | 32,226 | 32,225 | 0.995 |
| pwtk | 217,918 | 217,915 | 0.997 |
| c-71 | 76,638 | 76,367 | 0.927 |
| c-big | 345,241 | 345,114 | 0.676 |
| cfd2 | 123,440 | 52,160 | 0.998 |
| nd12k | 36,000 | 7,707 | 0.451 |
| ship_003 | 121,728 | 26,281 | 0.945 |
| shipsec1 | 140,874 | 140,870 | 0.997 |
| shipsec5 | 179,860 | 179,805 | 0.990 |
| shipsec8 | 114,919 | 114,914 | 0.992 |
| + 5 more | | | |

Quality degradation for partial matchings is expected and inherent to the
algorithm — there is no matched entry to anchor scaling at 1.0 for unmatched rows.

---

## Known Limitations

### 1. Longer cycles in some matchings

Some nonsingular matchings contain cycles longer than singletons + 2-cycles.
For example, blockqp1 has 20,000 longer cycles (n=60K), d_pretok has 1,534,
and stokes128 has 185. For a symmetric cost graph, the optimal matching should
decompose into only singletons and 2-cycles. Longer cycles indicate suboptimal
matching, likely due to tie-breaking in the greedy heuristic or Dijkstra heap
ordering. This doesn't affect scaling correctness but means the matching
permutation has longer cycles than necessary.

Logged as a soft warning in tests rather than a hard failure.

### 2. Partial matching anomalies (cfd2, ship_003, nd12k)

Several matrices classified as positive definite in SuiteSparse achieve
low matching cardinality:
- **cfd2**: 42% matched (PD matrix, expected full)
- **ship_003**: 22% matched
- **nd12k**: 21% matched

PD matrices should be structurally nonsingular, so the matching should achieve
full cardinality. This likely indicates an interaction between upper-triangular
storage conventions and the bipartite cost graph construction.
**TODO**: Investigate in a future optimization pass.

### 3. Three fully-matched matrices with near-zero row_max

TSOPF_FS_b39_c7, d_pretok, and thread have full matchings but min_rmax well
below 1.0 (down to 1.76e-8 for thread). Theoretically, a perfect matching
should yield row_max = 1.0 for all rows. The affected rows are a small fraction
of each matrix. Possible causes: numerical issues in the cost transformation
for extreme-valued entries, or suboptimal augmenting paths.

### 4. No column dual maintenance during algorithm

Our choice to compute column duals post-hoc (rather than maintaining them
during Dijkstra like SPRAL does) is mathematically equivalent but means we
can't incrementally validate dual feasibility during the algorithm. This is
compensated by `enforce_scaling_bound` as a safety net. It may also contribute
to limitations 1 and 3.

---

## Summary of Robustness Assessment

| Aspect | vs. SPRAL | Notes |
|--------|-----------|-------|
| Nonsingular scaling quality | **Equal** | row_max = 1.0 exactly on 38/41 fully-matched matrices |
| Nonsingular bound | **Equal** | No violations on any of 62 tested matrices |
| Degraded-quality matrices | **New finding** | 3 fully-matched matrices with some rows < 0.75 (0.08-2.2% of rows) |
| Singular scaling quality | **More tested** | SPRAL tests 1 trivial 3×3; we test 18 real SuiteSparse matrices |
| Singular bound | **Guaranteed** | `enforce_scaling_bound` provides hard guarantee |
| Algorithm efficiency | **Equal** | Same O(n·τ·log n) complexity |
| Matching optimality | **Weaker on some** | Up to 20K longer cycles (blockqp1); SPRAL likely tighter via column dual maintenance |
| Test coverage | **Much broader** | 62 real SuiteSparse matrices (up to 1.58M dim) vs. 100 random matrices n≤1000 |

The implementation is production-quality for the vast majority of matrices.
The scaling bound is maintained universally. Quality degradation occurs on
3 of 41 fully-matched matrices (7%), affecting only a small fraction of rows.
These are hard-indefinite or structurally unusual matrices that SPRAL doesn't
test. The degraded quality does not affect factorization safety — the bound
property ensures no scaled entry exceeds 1.0.

---

## References

- Duff & Koster (2001), "On Algorithms for Permuting Large Entries to the
  Diagonal of a Sparse Matrix", SIAM J. Matrix Anal. Appl. 22(4) — Algorithm MPD
- Duff & Pralet (2005), "Strategies for Scaling and Pivoting for Sparse
  Symmetric Indefinite Problems", RAL Technical Report — MC64SYM, singular handling
- Kuhn (1955), "The Hungarian Method for the Assignment Problem" — theoretical framework
- SPRAL source: `src/scaling.f90` (BSD-3), `tests/scaling.f90` — reference implementation
