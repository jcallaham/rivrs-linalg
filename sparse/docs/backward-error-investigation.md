# Backward Error Investigation: Full SuiteSparse Failures

**Date**: 2026-02-17
**Branch**: `017-two-level-aptp`
**Context**: First full SuiteSparse backward error evaluation (67 matrices)

## Summary

Running the full SuiteSparse collection through `SparseLDLT::solve_full` with default
`MatchOrderMetis` ordering revealed two distinct failure modes:

1. **Ordering regression**: MatchOrderMetis condensed ordering degrades accuracy on
   matrices that don't need MC64 pairing (linverse, and likely others)
2. **Large-front accuracy**: Several matrices fail with ALL orderings at ~1e-3 backward
   error, likely due to accumulated rounding in large frontal factorizations

## Failure Mode 1: MatchOrderMetis Ordering Regression

### Observation

`GHS_indef/linverse` (n=11,999, nnz=95,977) achieves machine-precision backward error
with plain METIS or AMD, but degrades to 1e-5 with MatchOrderMetis:

| Ordering | nb | BE | delays | 2x2 | max_front |
|----------|---:|---:|-------:|----:|----------:|
| MatchOrderMetis | 256 | 1.11e-5 | 42 | 26 | 12 |
| Same ordering, NO scaling | 256 | 1.11e-5 | 43 | 26 | 12 |
| Metis | 256 | 6.47e-18 | 19 | 14 | 11 |
| Amd | 256 | 2.11e-18 | 0 | 0 | 5 |

### Root Cause: Condensed Ordering, Not Scaling

The critical experiment: using the MatchOrderMetis ordering **without** MC64 scaling
gives the exact same backward error (1.11e-5). Scaling is not the cause.

The problem is the condensed METIS ordering itself:
- MC64 finds 1,909 two-cycles on linverse (surprisingly many)
- The condensed graph has dimension 10,090 (vs 11,999 original)
- METIS on the condensed graph produces a worse elimination order
- More delays (42 vs 19) and more 2x2 pivots (26 vs 14) result

### Matrix Properties

linverse has:
- **All positive diagonal** (11,999 positive, 0 negative, 0 zero)
- Well-conditioned diagonal: min=7.37e-2, max=9.16e0, ratio=124
- MC64 scaling range: min=0.33, max=3.68, ratio=11.1
- Small fronts (max_front=12) — rounding error is negligible

Despite having no zero diagonal entries, MC64 finds many off-diagonal entries
that are larger than the diagonal, creating unnecessary 2-cycle pairings that
constrain the ordering.

### SPRAL Comparison

**SPRAL defaults to `ordering=1` (plain METIS), not matching+METIS.**

Key findings from SPRAL source analysis (`ssids.f90`, `match_order.f90`):
- `ordering=1` (default): Plain METIS on full graph
- `ordering=2` (optional): MC64 matching + METIS on condensed graph
- When `ordering=2` is used, SPRAL always runs it unconditionally — no fallback
- No quality comparison between orderings — the user chooses
- No detection of "unnecessary" matchings

SPRAL's documentation describes `ordering=2` as for "numerically difficult problems."
The default plain METIS works well for most matrices.

### Fix Options

1. **Change default to Metis**: Match SPRAL's default. MatchOrderMetis becomes opt-in
   for hard indefinite problems. This fixes linverse but regresses bratu3d (53K delays,
   BE=5e-3 with plain METIS).

2. **Heuristic-based default**: Check diagonal properties:
   - If matrix has zero diagonal entries → MatchOrderMetis
   - If all diagonal entries are nonzero → plain Metis
   This matches the intended use of MC64 (stabilize matrices with zero diagonals).

3. **Try-both strategy**: Run both orderings, compare predicted fill, pick better.
   More expensive but robust.

4. **Expose as user choice**: Like SPRAL, let the user decide and document guidance.

## Failure Mode 2: Large-Front Accuracy

### Observation

Several matrices fail with ALL orderings at ~1e-3 backward error:

| Matrix | n | MatchOrder BE | Metis BE | AMD BE | max_front |
|--------|--:|-------------:|---------:|-------:|----------:|
| helm3d01 | 32,226 | 1.39e-3 | 2.14e-3 | 1.19e-3 | 1,207-1,424 |
| dawson5 | 51,537 | 1.85e-3 | 2.21e-3 | 1.67e-3 | 565-916 |
| copter2 | 55,476 | 1.09e-3 | 1.12e-3 | 1.43e-3 | 1,218-1,799 |
| astro-ph | 16,706 | 1.71e-3 | 1.99e-3 | 3.46e-3 | 2,429-3,240 |

Common characteristics:
- Large fronts (565-3,240)
- Significant 2x2 pivot usage (6-50% of total pivots)
- Backward error ~1e-3 regardless of ordering

### astro-ph: Special Case

Newman/astro-ph has ALL ZERO diagonal entries (pure graph adjacency matrix, no
self-loops). This is an extreme case:
- 16,706 zero diagonals
- 35,707 total delays with MatchOrderMetis
- 969 zero pivots
- Max front = 2,429

This matrix is genuinely difficult for LDL^T — every pivot must be 2x2.

### SPRAL Comparison: Inner Factorization Architecture

**SPRAL's inner kernel is also BLAS-2 (same as ours)**, but with a critical structural
difference in the two-level blocking:

**SPRAL architecture** (`ldlt_app.cxx`):
- Outer block size: `cpu_block_size` (typically 256)
- Inner block size: `INNER_BLOCK_SIZE = 32` (hardcoded)
- Within each outer tile: recurse with inner block size
  - Factor 32 columns (BLAS-2 via `block_ldlt`)
  - Apply pivot: TRSM (BLAS-3)
  - Update trailing: GEMM (BLAS-3)
  - Repeat for next 32 columns
- Between outer tiles: TRSM + GEMM (BLAS-3)

**Our architecture** (after Phase 8.1):
- Outer block size: `outer_block_size` (256)
- `factor_inner` processes entire tile at once (BLAS-2 per-column)
- Between outer tiles: TRSM + GEMM (BLAS-3, just implemented)

The key difference: for a 256-column tile, SPRAL does ~8 rounds of inner BLAS-3
(32-column blocks with GEMM/TRSM between them). We do 256 columns of BLAS-2
updates sequentially. This means SPRAL's trailing updates within a tile benefit
from BLAS-3 numerical properties (block operations with better error accumulation).

**This was the intended Phase 8.1 design** — `inner_block_size` (default 32) was
supposed to control sub-tile BLAS-3 processing. The current code sets
`block_size = end_pos - k` (full tile) instead of `(end_pos - k).min(ib)` (inner
blocks of size ib=32). The `_ib` variable is unused (`factor.rs:1536`).

The tradeoff: larger inner blocks give better pivot search scope (complete pivoting
over more candidate columns), but lose the BLAS-3 structure between inner blocks.
SPRAL chose ib=32 as the fixed inner block size, accepting the narrower search
scope in exchange for better BLAS-3 numerical properties. This is likely why SPRAL
achieves better accuracy on large fronts — the trailing matrix gets GEMM updates
every 32 columns, keeping the numerical state cleaner than 256+ columns of
BLAS-2 rank-1 updates.

### SPRAL: No Iterative Refinement

**SPRAL does NOT perform iterative refinement.** The solve path is:

1. Permute + scale RHS
2. Forward solve (TRSM + GEMM per supernode, leaves to root)
3. Diagonal solve (D^{-1} multiply, handles 1x1 and 2x2)
4. Backward solve (GEMM + TRSM per supernode, root to leaves)
5. Unscale + unpermute

No residual computation, no backward error check, no refinement step.

**This means SPRAL achieves its accuracy purely through the factorization quality.**
The difference in our backward error must come from factorization differences,
not from iterative refinement.

### SPRAL: D^{-1} Storage

One notable SPRAL convention: D stores **D^{-1}** (the inverse), not D itself.
The diagonal solve multiplies by D^{-1} directly. Our code stores D and divides.
This is a minor numerical difference but should not account for 1e-3 backward error.

### Fix Options

1. **Restore inner block processing in `factor_inner`**: Change `block_size` back to
   `(end_pos - k).min(ib)` where `ib=32` (matching SPRAL). This enables BLAS-3
   between inner blocks within each tile. This was the original Phase 8.1 design.

2. **Iterative refinement**: Not strictly needed (SPRAL doesn't have it), but would
   paper over accuracy issues. One step of `x += A^{-1}(b - Ax)` typically improves
   backward error by the condition number ratio. Easy to implement at the SparseLDLT
   level (requires one extra solve + one matvec).

3. **Complete pivoting quality**: SPRAL's `block_ldlt` uses `find_maxloc` (complete
   pivoting) to find the best pivot in the remaining block. Our `factor_inner` uses
   `try_1x1` on the diagonal followed by `try_2x2` with the best off-diagonal partner.
   The pivot search scope may differ.

## Full SuiteSparse Results by Category

Results with `MatchOrderMetis` (current default):

### Easy-Indefinite (28 matrices)

| Matrix | n | BE | Status | kind |
|--------|--:|---:|:------:|------|
| BenElechi1 | 245,874 | 7.33e-19 | PASS | structural |
| F2 | 71,505 | 1.44e-18 | PASS | structural |
| H2O | 67,024 | 4.72e-18 | PASS | quantum chemistry |
| Si10H16 | 17,077 | 1.31e-17 | PASS | quantum chemistry |
| Si5H12 | 19,896 | 8.21e-18 | PASS | quantum chemistry |
| SiNa | 5,743 | 9.26e-18 | PASS | quantum chemistry |
| bcsstk39 | 46,772 | 1.24e-18 | PASS | structural |
| bloweybq | 10,001 | 4.27e-11 | PASS | optimization |
| crystk02 | 13,965 | 3.53e-18 | PASS | structural |
| crystk03 | 24,696 | 1.74e-18 | PASS | structural |
| dixmaanl | 60,000 | 3.39e-18 | PASS | optimization |
| filter3D | 106,437 | 1.13e-18 | PASS | model reduction |
| gas_sensor | 66,917 | 1.32e-18 | PASS | model reduction |
| msdoor | 415,863 | 5.03e-19 | PASS | structural |
| nd3k | 9,000 | 1.26e-17 | PASS | 2D/3D problem |
| qa8fk | 66,127 | 1.98e-18 | PASS | acoustics |
| rail_79841 | 79,841 | 8.21e-19 | PASS | model reduction |
| t2dal | 4,257 | 2.64e-18 | PASS | model reduction |
| t3dh | 79,171 | 1.58e-18 | PASS | model reduction |
| vibrobox | 12,328 | 2.27e-18 | PASS | structural |
| **astro-ph** | 16,706 | 1.71e-3 | **FAIL** | undirected graph |
| **copter2** | 55,476 | 1.09e-3 | **FAIL** | CFD |
| **dawson5** | 51,537 | 1.85e-3 | **FAIL** | structural |
| **helm3d01** | 32,226 | 1.39e-3 | **FAIL** | CFD |
| **linverse** | 11,999 | 1.11e-5 | **FAIL** | inverse problem |
| **pwtk** | 217,918 | 6.63e-5 | **FAIL** | structural |
| **sparsine** | 50,000 | 9.92e-4 | **FAIL** | synthetic |
| **spmsrtls** | 29,995 | 3.23e-6 | **FAIL** | least squares |

**Score: 20/28 pass (71%).**

### Hard-Indefinite (6 matrices so far)

| Matrix | n | BE | Status | kind |
|--------|--:|---:|:------:|------|
| blockqp1 | 60,012 | 1.11e-13 | PASS | optimization |
| aug3dcqp | 35,543 | 3.45e-10 | RELAXED | optimization |
| bratu3d | 27,792 | 1.00e-9 | RELAXED | nonlinear PDE |
| **c-71** | 76,638 | 2.37e-5 | **FAIL** | optimization |
| **TSOPF_FS_b162_c1** | 10,798 | 4.40e-7 | **FAIL** | power network |
| **TSOPF_FS_b39_c7** | 28,216 | 2.94e-6 | **FAIL** | power network |

**Score: 3/6 pass or relaxed (50%).**

### Pattern Analysis

1. **Easy-indefinite failures fall with ALL orderings**: copter2, dawson5, helm3d01
   fail at ~1e-3 with MatchOrderMetis, Metis, AND AMD (see diagnostics above). This
   is a factorization quality issue, not an ordering issue.

2. **linverse is pure ordering regression**: machine precision with Metis, 1e-5 with
   MatchOrderMetis. Others (pwtk, spmsrtls) may be similar — needs verification.

3. **Hard-indefinite matrices benefit from MatchOrderMetis**: bratu3d went from 5e-3
   (plain METIS) to 1e-9 (MatchOrderMetis). blockqp1 achieves 1e-13.

4. **astro-ph is a structural outlier**: all-zero diagonal, graph adjacency matrix.
   Terrible with all orderings.

## SPRAL Source Analysis

### Condensed Graph Implementation: No Bug Found

Detailed comparison of our `build_condensed_adjacency` with SPRAL's `mo_split`:
- **Core condensation logic is correct**: marker-based dedup, self-loop prevention,
  shared-neighbor handling all match SPRAL semantics
- **Minor difference**: our code excludes unmatched nodes from condensed graph and
  appends at end; SPRAL includes them. Only affects structurally singular matrices.
- **Both produce the same METIS input** for structurally nonsingular matrices

### Ordering Strategy: Configuration, Not Bug

The Duff 2020 paper is explicit (Section 7):
> "In the case of easy-indefinite matrices, **no scaling is performed prior to the
> factorization** whereas hard-indefinite matrices are scaled and ordered using a
> matching-based ordering and scaling."

This means:
- **Easy-indefinite**: ordering=1 (plain METIS), no scaling
- **Hard-indefinite**: ordering=2 (MC64+METIS), with scaling

SPRAL's default is ordering=1 (plain METIS). MC64+METIS is the specialist option.

### Domain-Specific Guidance

| Domain | Recommended ordering | Why |
|--------|---------------------|-----|
| Structural FEM (stiffness, mass) | Plain METIS | Naturally dominant diagonal |
| Thermal/acoustic FEM | Plain METIS | Good diagonal structure |
| Model reduction (Oberwolfach) | Plain METIS | Same |
| Quantum chemistry (PARSEC) | Plain METIS | Dense blocks, good structure |
| **Interior point QP/LP (KKT)** | **MC64+METIS** | Saddle-point, zero diagonal blocks |
| **Optimal control** | **MC64+METIS** | KKT structure |
| **Power network TSOPF** | **MC64+METIS** | Optimization over networks |
| **Mixed FEM / Stokes** | **MC64+METIS** | Constrained FE = saddle-point |

The dividing line: **saddle-point structure with zero or near-zero diagonal blocks
needs MC64+METIS. Everything else is better with plain METIS.**

## Recommended Action Plan

### Immediate: Fix ordering default

**Option A — Auto-detection heuristic** (recommended):
Count zero diagonal entries. If any exist, use MatchOrderMetis. Otherwise, use Metis.
This matches the Duff 2020 paper's practice and handles both linverse (no zeros →
plain METIS → machine precision) and bratu3d (many zeros → MC64+METIS → 1e-9).

**Option B — Change default to Metis**:
Match SPRAL's default. Document MatchOrderMetis as opt-in for hard indefinite.
Simpler, but requires users to know their matrix structure.

### Short-term: Restore inner block processing

Change `factor_inner` to use `ib`-sized sub-blocks (matching SPRAL's
`INNER_BLOCK_SIZE=32`). This should improve large-front accuracy for matrices like
copter2 (max_front=1799), helm3d01 (max_front=1424), dawson5 (max_front=916).

### Medium-term: Iterative refinement

SPRAL does NOT use iterative refinement, so accuracy improvements must come from
factorization quality. However, one step of IR is cheap insurance and could be added
as an option at the SparseLDLT level.

## Diagnostic Examples

- `examples/diagnose_failures.rs` — multi-ordering comparison with pivot stats
- `examples/scaling_diagnosis.rs` — isolates scaling vs ordering effects

Run with: `cargo run --example diagnose_failures --release`


## Before audit

Partial result on `test_solve_suitesparse_full` using the Duff et al (2020) protocol (`MatchOrderMetis` for hard-indefinite, otherwise plain METIS):

```
Matrix                                          n           BE  Status
--------------------------------------------------------------------------------
BenElechi/BenElechi1                       245874     6.80e-19  PASS
Koutsovasilis/F2                            71505     1.45e-18  PASS
PARSEC/H2O                                  67024     4.66e-18  PASS
PARSEC/Si10H16                              17077     3.00e-17  PASS
PARSEC/Si5H12                               19896     7.82e-18  PASS
PARSEC/SiNa                                  5743     9.56e-18  PASS
Newman/astro-ph                             16706      1.68e-3  FAIL
Boeing/bcsstk39                             46772     1.15e-18  PASS
GHS_indef/bloweybq                          10001     3.26e-11  PASS
GHS_indef/copter2                           55476      1.47e-3  FAIL
Boeing/crystk02                             13965      1.40e-5  FAIL
Boeing/crystk03                             24696      1.98e-6  FAIL
GHS_indef/dawson5                           51537      2.61e-3  FAIL
GHS_indef/dixmaanl                          60000     4.89e-18  PASS
Oberwolfach/filter3D                       106437     1.05e-18  PASS
Oberwolfach/gas_sensor                      66917     1.12e-18  PASS
GHS_indef/helm3d01                          32226      3.41e-3  FAIL
GHS_indef/linverse                          11999     6.47e-18  PASS
INPRO/msdoor                               415863     4.83e-19  PASS
ND/nd3k                                      9000     1.30e-17  PASS
Boeing/pwtk                                217918      1.38e-4  FAIL
Cunningham/qa8fk                            66127     1.86e-18  PASS
Oberwolfach/rail_79841                      79841     7.34e-19  PASS
GHS_indef/sparsine                          50000      3.03e-3  FAIL
GHS_indef/spmsrtls                          29995      1.72e-6  FAIL
Oberwolfach/t2dal                            4257     2.66e-18  PASS
Oberwolfach/t3dh                            79171     1.52e-18  PASS
Cote/vibrobox                               12328      1.55e-5  FAIL
TSOPF/TSOPF_FS_b162_c1                      10798      3.90e-7  FAIL
TSOPF/TSOPF_FS_b39_c7                       28216      4.55e-6  FAIL
GHS_indef/aug3dcqp                          35543     3.45e-10  RELAXED
GHS_indef/blockqp1                          60012     1.11e-13  PASS
GHS_indef/bratu3d                           27792      7.76e-4  FAIL
GHS_indef/c-71                              76638      1.76e-5  FAIL
```

Result with `MatchOrderMetis` only:

```
Matrix                                          n           BE  Status
--------------------------------------------------------------------------------
BenElechi/BenElechi1                       245874     7.17e-19  PASS
Koutsovasilis/F2                            71505     1.51e-18  PASS
PARSEC/H2O                                  67024     4.72e-18  PASS
PARSEC/Si10H16                              17077     2.62e-17  PASS
PARSEC/Si5H12                               19896     8.10e-18  PASS
PARSEC/SiNa                                  5743     8.82e-18  PASS
Newman/astro-ph                             16706      2.65e-3  FAIL
Boeing/bcsstk39                             46772     1.24e-18  PASS
GHS_indef/bloweybq                          10001     4.27e-11  PASS
GHS_indef/copter2                           55476      1.26e-3  FAIL
Boeing/crystk02                             13965     3.42e-18  PASS
Boeing/crystk03                             24696     1.48e-18  PASS
GHS_indef/dawson5                           51537      1.13e-3  FAIL
GHS_indef/dixmaanl                          60000     3.39e-18  PASS
Oberwolfach/filter3D                       106437     1.09e-18  PASS
Oberwolfach/gas_sensor                      66917     1.24e-18  PASS
GHS_indef/helm3d01                          32226      1.86e-3  FAIL
GHS_indef/linverse                          11999      1.11e-5  FAIL
INPRO/msdoor                               415863     4.99e-19  PASS
ND/nd3k                                      9000     1.26e-17  PASS
Boeing/pwtk                                217918      1.48e-4  FAIL
Cunningham/qa8fk                            66127     1.81e-18  PASS
Oberwolfach/rail_79841                      79841     8.21e-19  PASS
```


Result with plain METIS only:

```
Matrix                                          n           BE  Status
--------------------------------------------------------------------------------
BenElechi/BenElechi1                       245874     6.80e-19  PASS
Koutsovasilis/F2                            71505     1.45e-18  PASS
PARSEC/H2O                                  67024     4.66e-18  PASS
PARSEC/Si10H16                              17077     3.00e-17  PASS
PARSEC/Si5H12                               19896     7.82e-18  PASS
PARSEC/SiNa                                  5743     9.56e-18  PASS
Newman/astro-ph                             16706      1.68e-3  FAIL
Boeing/bcsstk39                             46772     1.15e-18  PASS
GHS_indef/bloweybq                          10001     3.26e-11  PASS
GHS_indef/copter2                           55476      1.47e-3  FAIL
Boeing/crystk02                             13965      1.40e-5  FAIL
Boeing/crystk03                             24696      1.98e-6  FAIL
GHS_indef/dawson5                           51537      2.61e-3  FAIL
GHS_indef/dixmaanl                          60000     4.89e-18  PASS
Oberwolfach/filter3D                       106437     1.05e-18  PASS
Oberwolfach/gas_sensor                      66917     1.12e-18  PASS
GHS_indef/helm3d01                          32226      3.41e-3  FAIL
GHS_indef/linverse                          11999     6.47e-18  PASS
INPRO/msdoor                               415863     4.83e-19  PASS
ND/nd3k                                      9000     1.30e-17  PASS
Boeing/pwtk                                217918      1.38e-4  FAIL
Cunningham/qa8fk                            66127     1.86e-18  PASS
Oberwolfach/rail_79841                      79841     7.34e-19  PASS
GHS_indef/sparsine                          50000      3.03e-3  FAIL
GHS_indef/spmsrtls                          29995      1.72e-6  FAIL
Oberwolfach/t2dal                            4257     2.66e-18  PASS
Oberwolfach/t3dh                            79171     1.52e-18  PASS
Cote/vibrobox                               12328      1.55e-5  FAIL
TSOPF/TSOPF_FS_b162_c1                      10798      9.54e-8  FAIL
TSOPF/TSOPF_FS_b39_c7                       28216      7.78e-6  FAIL
GHS_indef/aug3dcqp                          35543     2.21e-10  RELAXED
GHS_indef/blockqp1                          60012     1.75e-14  PASS
GHS_indef/bratu3d                           27792      6.15e-3  FAIL
GHS_indef/c-71                              76638      1.33e-6  FAIL
```


## After `extract_front_factors` fix


```
Matrix                                          n           BE  Status
--------------------------------------------------------------------------------
BenElechi/BenElechi1                       245874     6.80e-19  PASS
Koutsovasilis/F2                            71505     1.45e-18  PASS
PARSEC/H2O                                  67024     4.66e-18  PASS
PARSEC/Si10H16                              17077     3.00e-17  PASS
PARSEC/Si5H12                               19896     7.82e-18  PASS
PARSEC/SiNa                                  5743     9.56e-18  PASS
Newman/astro-ph                             16706      2.79e-3  FAIL
Boeing/bcsstk39                             46772     1.15e-18  PASS
GHS_indef/bloweybq                          10001     2.68e-18  PASS
GHS_indef/copter2                           55476      6.87e-4  FAIL
Boeing/crystk02                             13965      1.40e-5  FAIL
Boeing/crystk03                             24696      1.98e-6  FAIL
GHS_indef/dawson5                           51537     7.45e-17  PASS
GHS_indef/dixmaanl                          60000     4.89e-18  PASS
Oberwolfach/filter3D                       106437     1.05e-18  PASS
Oberwolfach/gas_sensor                      66917     1.12e-18  PASS
GHS_indef/helm3d01                          32226      1.24e-3  FAIL
GHS_indef/linverse                          11999     6.47e-18  PASS
INPRO/msdoor                               415863     4.83e-19  PASS
ND/nd3k                                      9000     1.30e-17  PASS
Boeing/pwtk                                217918     6.55e-19  PASS
Cunningham/qa8fk                            66127     1.86e-18  PASS
Oberwolfach/rail_79841                      79841     7.34e-19  PASS
GHS_indef/sparsine                          50000      2.18e-3  FAIL
GHS_indef/spmsrtls                          29995     1.30e-17  PASS
Oberwolfach/t2dal                            4257     2.66e-18  PASS
Oberwolfach/t3dh                            79171     1.52e-18  PASS
Cote/vibrobox                               12328      4.22e-4  FAIL
TSOPF/TSOPF_FS_b162_c1                      10798      2.34e-7  FAIL
TSOPF/TSOPF_FS_b39_c7                       28216      5.98e-8  FAIL
GHS_indef/aug3dcqp                          35543     7.00e-19  PASS
GHS_indef/blockqp1                          60012     1.11e-13  PASS
GHS_indef/bratu3d                           27792     2.56e-18  PASS
```

`MatchOrderMetis` only
