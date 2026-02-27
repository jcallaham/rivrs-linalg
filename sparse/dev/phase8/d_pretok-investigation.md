# d_pretok Investigation Log

## Problem Statement

After TPP fallback implementation (Phase 8.1d), 64/65 SuiteSparse matrices pass
the backward error threshold (<5e-11). The sole remaining failure is **d_pretok**
(n=182,269) with backward error **2.50e-6**.

SPRAL achieves machine precision (1.30e-18) on the same matrix with identical
MC64 matching+scaling and condensed METIS ordering.

## Key Comparison Data (SPRAL vs Us)

| Metric | SPRAL | Ours |
|--------|-------|------|
| Backward error (test RHS) | 1.30e-18 | 2.50e-6 |
| 2x2 pivots | 77,468 | 29,336 |
| Delays (after TPP) | 27 | 403 |
| `not_first_pass` | 2 | — |

The 2x2 pivot count discrepancy (77K vs 29K) is the smoking gun — SPRAL accepts
far more 2x2 pivots than we do. This leads to a worse Schur complement cascade
and ultimately poor solution accuracy.

## Confirmed Non-Issues

- MC64 matching+scaling: identical between solvers (verified via spral_solve_comparison)
- Condensed METIS ordering: identical (same comparison)
- Backward error formula: both use compatible norms
- block_ldlt 2x2 determinant test: mathematically identical to SPRAL's test_2x2
- 1x1 fallback threshold check: mathematically impossible to fail with complete pivoting (SPRAL has exit(1) there)
- D storage convention: we store D, SPRAL stores D^{-1}, but both compute L21 = A21 * L11^{-T} * D^{-1} correctly
- apply_and_check: threshold scan matches SPRAL's check_threshold (|L_entry| >= 1/u)
- Default threshold: both use u = 0.01

## Bisection Strategy

### Level 1: Dense kernel vs multifrontal assembly

Export a frontal matrix before factoring → factor with both our kernel and SPRAL's
ldlt_app_factor → compare results. If they match, issue is in assembly.

### Level 2: Per-supernode stats

Dump per-supernode (front_size, k, nelim, delays, 2x2_count) to identify which
supernodes diverge first.

## Findings

### Per-Supernode Stats

- 50,902 total supernodes, overwhelmingly small fronts (k=1-4)
- 403 delays all occur in small fronts (max front size ~few hundred)
- The vast majority of supernodes have k < inner_block_size (32)

### Root Cause: SPRAL Uses TPP for Small Blocks

**SPRAL's `Block::factor()` (ldlt_app.cxx:994-1020)**:
```cpp
if(ncol() < INNER_BLOCK_SIZE || !is_aligned(aval_)) {
    // Use TPP for partial/unaligned blocks
    cdata_[i_].nelim = ldlt_tpp_factor(nrow(), ncol(), ...);
} else {
    // Use complete pivoting for full aligned blocks
    block_ldlt<T, INNER_BLOCK_SIZE>(...);
}
```

For d_pretok with ~50K small supernodes (k < 32), SPRAL uses **TPP as the primary
factorization method** for the vast majority. We were using **complete pivoting
(`factor_block_diagonal`)** for ALL blocks regardless of size.

**Critical behavioral difference:**
- **Complete pivoting** (`block_ldlt`/`factor_block_diagonal`): Find global max entry
  first. If on diagonal → take 1x1 (no 2x2 attempted). Only tries 2x2 when the
  off-diagonal is the global maximum.
- **TPP** (`ldlt_tpp_factor`): Iterate columns, try 2x2 FIRST for each column →
  accepts 2x2 pivots whenever the determinant test passes, even when a good 1x1 exists.

For indefinite matrices, 2x2 pivots pair positive/negative eigenvalues and produce
better-conditioned factors. TPP accepts more 2x2 pivots because it tries 2x2 first
on every column search, while complete pivoting only tries 2x2 when forced by the
maximum location.

### Fix

Modified `aptp_factor_in_place` dispatch:
```rust
let mut result = if num_fully_summed < options.inner_block_size {
    // Small front: use TPP as primary method (matching SPRAL)
    tpp_factor_as_primary(a.rb_mut(), num_fully_summed, options)?
} else if num_fully_summed > options.outer_block_size {
    two_level_factor(a.rb_mut(), num_fully_summed, options)?
} else {
    factor_inner(a.rb_mut(), num_fully_summed, options)?
};
```

### Results

| Metric | Before | After | SPRAL |
|--------|--------|-------|-------|
| Backward error | 2.50e-6 | 7.21e-19 | 1.30e-18 |
| 2x2 pivots | 29,336 | 64,993 | 77,468 |
| Delays | 403 | — | 27 |

**65/65 SuiteSparse matrices now pass** the backward error threshold.

The remaining gap in 2x2 count (65K vs 77K) is from large fronts (k >= 32) where
we still use complete pivoting for inner blocks within `factor_inner`. SPRAL also
uses TPP for partial blocks within its two-level loop. This could be addressed in
a future optimization but is not needed for correctness.
