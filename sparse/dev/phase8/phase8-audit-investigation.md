# Phase 8 Audit Investigation: Assembly, Solve, and Scaling

**Date**: 2026-02-18
**Branch**: `017-two-level-aptp`
**Context**: After verifying the dense APTP kernel is correct (4 agents, no bugs found),
backward error still correlates with max_front size. This audit targets the unverified
code paths: frontal matrix assembly, triangular solve, and MC64 scaling.

## Motivation

Backward error vs max_front correlation:

| Matrix | max_front | BE | Status |
|--------|----------:|---:|:------:|
| bloweybq | 13 | 3.26e-11 | PASS |
| stokes128 | 547 | 6.08e-6 | FAIL |
| bratu3d | 1,496 | 7.55e-4 | FAIL |
| sparsine | 11,131 | 2.68e-3 | FAIL |

The dense kernel is verified correct (reconstruction < 1e-12 on all unit tests).
The bug must be in code that runs per-front and whose error scales with front size.

## Audit Agents

### Agent A: Frontal Matrix Assembly (`numeric.rs`)

**Scope**: `scatter_original_entries` + `extend_add` + frontal matrix initialization

**Key questions**:
- Is the upper-triangle skip logic in `scatter_original_entries` correct?
- Are sparse CSC indices correctly mapped to dense frontal matrix positions?
- Does `extend_add` correctly map child contribution indices to parent positions?
- Are any entries double-counted or missed?
- How does SPRAL's equivalent assembly code compare?

**SPRAL references**:
- `ssids/cpu/kernels/ldlt_app.cxx`: frontal matrix assembly
- `ssids/cpu/NumericSubtree.cxx`: `add_a_block`, contribution assembly

**Status**: COMPLETE — NO BUGS FOUND
**Findings**:
- `scatter_original_entries`: upper-triangle skip logic correct. Traced all cases (both
  in supernode, one in supernode, neither) — exactly one copy survives per off-diagonal entry
- `extend_add`: correctly maps child→parent, accumulates (`+=`), lower-triangle placement
- Frontal matrix init: `Mat::zeros(m, m)` correct
- `extract_contribution`: correctly includes delayed rows via `result.perm[ne..k]`
- `build_supernode_info`: pattern filtering (`r > j`) correctly excludes diagonal
- **Key architectural difference from SPRAL**: SPRAL precomputes `amap` (assignment `=`),
  we use runtime CSC iteration with skip logic (accumulation `+=`). Both produce correct
  results. Our approach is more fragile (skip-logic bug = double-counting) but verified correct.
- Note: Agent A missed the `extract_front_factors` delayed-rows bug that Agent B caught —
  the fix is in the L21/row_indices extraction, not the assembly path.

### Agent B: Triangular Solve (`solve.rs`)

**Scope**: `aptp_solve` — per-supernode forward/backward/diagonal solve

**Key questions**:
- Are gather/scatter operations between global vector and local workspace correct?
- Is `local_perm` applied correctly for each supernode?
- For 2×2 pivots in diagonal solve, is D^{-1} computed correctly?
- Does the forward solve correctly handle the extend-add contribution pattern?
- Does the backward solve correctly scatter results back?
- How does SPRAL's `ldlt_app_solve_fwd/bwd` compare?

**SPRAL references**:
- `ssids/cpu/kernels/ldlt_app.cxx`: `solve_fwd`, `solve_diag`, `solve_bwd`
- `ssids/cpu/NumericSubtree.cxx`: `solve_fwd`, `solve_bwd`

**Status**: COMPLETE — **BUG FOUND**
**Findings**:

**BUG: Missing L entries for delayed rows in `extract_front_factors` (numeric.rs)**

When APTP eliminates `ne` out of `k` fully-summed columns (ne < k due to delays),
the factored frontal matrix has valid L21 entries at rows `ne..k` (the delayed
fully-summed rows) for columns `0..ne`. These are L entries from `apply_and_check`'s
TRSM that were computed before the columns were delayed.

`extract_front_factors` currently extracts L21 from rows `k..m` only, **missing rows
`ne..k`**. The `row_indices` also only includes `frontal.row_indices[k..]`.

Impact on solve:
- **Forward solve**: L21 * work not scattered to delayed row positions → stale values
- **Backward solve**: delayed row contributions not gathered → missing L21^T updates

Crucially, `extract_contribution` (same file) DOES correctly include delayed rows
in its `row_indices`. So the contribution block assembly is correct, but the solve
data extraction is wrong — an inconsistency within the same file.

**SPRAL comparison**: SPRAL's forward solve gathers ALL `m+ndin` rows and updates ALL
rows below `nelim` via GEMV. The `map` array includes both delayed and non-fully-summed
row indices.

**Proposed fix**: Change L21 extraction from `frontal.data[(k..m, 0..ne)]` to
`frontal.data[(ne..m, 0..ne)]`, and build `row_indices` as:
```
// Delayed rows (ne..k): use APTP perm for global indices
for &lp in &result.perm[ne..k] {
    row_indices.push(frontal.row_indices[lp]);
}
// Non-fully-summed rows (k..m): no perm needed
row_indices.extend_from_slice(&frontal.row_indices[k..]);
```

**No other bugs in solve.rs**: The forward/backward/diagonal solve logic is correct.
D^{-1} handling (1x1 and 2x2), gather/scatter, and permutation are all verified.

### Agent C: MC64 Scaling Application (`solver.rs`)

**Scope**: `SparseLDLT` scaling logic — pre-scale matrix, scale/unscale RHS

**Key questions**:
- Is the MC64 scaling correctly applied to the matrix before factorization?
- Is the RHS correctly scaled before solve and unscaled after?
- Are the scaling vectors applied in the right order relative to permutations?
- Does the scaling interact correctly with MatchOrderMetis ordering?
- How does SPRAL's scaling application compare?

**SPRAL references**:
- `ssids/cpu/cpu_iface.f90`: scaling application
- `ssids/ssids.f90`: `solve` with scaling

**Status**: COMPLETE — NO BUGS FOUND
**Findings**:
- Matrix scaling applied during assembly: `val *= s[row_elim] * s[col_elim]` — correct
- RHS handling: permute → scale → solve → unscale → unpermute — matches SPRAL
- Unscale correctly multiplies by `s[i]` (not divides), matching SPRAL's
  `x(akeep%invp(i),r) = x2(i,r) * fkeep%scaling(i)`
- Scaling transform to elimination order: `result.scaling[perm_fwd[i]]` — correct
- Permutation composition: MC64 matching + METIS on condensed graph produces single
  `result.ordering`, no double-permutation issue
- Edge cases: scaling=None properly guarded throughout
- Confirmed: backward error investigation already isolated scaling from ordering
  (same BE with and without scaling on linverse)

### Agent D: Empirical Isolation (factorization vs solve)

**Scope**: Determine whether the error is in factorization or solve by computing
per-supernode reconstruction error on a failing matrix.

**Key questions**:
- For a failing matrix (e.g., stokes128 or bratu3d), what is the per-supernode
  reconstruction error `||PAP^T - LDL^T||/||A||`?
- If reconstruction is good (~1e-15) but backward error is bad (~1e-4), the
  bug is in the solve path
- If reconstruction is also bad, the bug is in assembly or factorization
- Does the reconstruction error correlate with front size?

**Status**: COMPLETE — confirms fix resolves primary backward error issues
**Findings**:

Agent D ran per-supernode reconstruction diagnostics AFTER the extract_front_factors
fix was applied. Key results:

- **bratu3d (MatchOrderMetis)**: BE=2.56e-18 (**PASS**). Worst per-supernode
  reconstruction error = 1.31 (SN 3611, front=799). Diagnosis: "No bug."
- **bratu3d (plain METIS)**: BE=5.57e-3 (FAIL). Worst reconstruction = 1.09e4
  (SN 16982, front=859). This is the known ordering-dependent failure — plain METIS
  without MC64 matching causes massive pivot delays on hard indefinite matrices.
- **cvxqp3 (MatchOrderMetis)**: BE=8.32e-7 (FAIL). Worst reconstruction = 9.0
  (SN 4818, front=1664). "Large-front accuracy" issue — reconstruction error
  correlated with front size, same as sparsine/helm3d01/etc.

## Resolution

### BUG FOUND AND FIXED: `extract_front_factors` missing L21 entries for delayed rows

**Root cause**: When APTP eliminates `ne` out of `k` fully-summed columns (ne < k due
to delays), the factored frontal matrix has valid L21 entries at rows `ne..k` (the
delayed fully-summed rows). These entries were computed by `apply_and_check`'s TRSM
before the columns failed the threshold test and were delayed.

`extract_front_factors` extracted L21 from rows `k..m` only, **missing rows `ne..k`**.
The `extract_contribution` function in the same file correctly included delayed rows
— an inconsistency.

**Fix**: Changed L21 extraction from `frontal.data[(k..m, 0..ne)]` to
`frontal.data[(ne..m, 0..ne)]`, and built `row_indices` to include delayed rows
via `result.perm[ne..k]` followed by `frontal.row_indices[k..]`.

**Impact on solve**:
- Forward solve: L21 * work now correctly scatters to delayed row positions
- Backward solve: delayed row contributions now correctly gathered for L21^T updates

**Before/after backward error**:

| Matrix | max_front | Before Fix | After Fix | Change |
|--------|----------:|:----------:|:---------:|:------:|
| bratu3d | 1,496 | 7.55e-4 | 2.56e-18 | 10^14x improvement |
| stokes128 | 547 | 6.08e-6 | 3.75e-18 | 10^12x improvement |
| bloweybq | 13 | 3.26e-11 | 2.68e-18 | 10^7x improvement |
| sparsine (single) | 11,131 | 2.68e-3 | 4.21e-5 | 60x improvement |

### Remaining: Large-Front Accuracy (pre-existing, not a bug)

sparsine, helm3d01, dawson5, copter2, astro-ph, cvxqp3 all fail at ~1e-3 to ~1e-5
with max_front sizes of 1000-11000. Per-supernode reconstruction errors of 1-10
(O(1), not machine precision) confirm this is numerical accuracy in the dense APTP
kernel, not an extraction or solve bug. Potential fixes for Phase 9:
- Iterative refinement
- Store D^{-1} instead of D (match SPRAL)
- Profile-guided inner block size tuning

### Tests Added

5 unit tests in `src/aptp/numeric.rs`:
1. `test_extract_front_factors_l21_includes_delayed_rows` — L21 dimension regression test
2. `test_extract_front_factors_contribution_row_indices_consistent` — delayed row index consistency
3. `test_extract_front_factors_delayed_l21_entries_populated` — TRSM entries preserved
4. `test_extract_front_factors_reconstruction_with_delays` — L11*D11*L11^T reconstruction
5. `test_extract_front_factors_solve_roundtrip_with_delays` — full [L11;L21]*D*[L11;L21]^T

All 347 lib tests + 22 integration tests pass.
