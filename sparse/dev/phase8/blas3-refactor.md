# BLAS-3 Refactoring Plan for Two-Level APTP

## Status

**Deferred from Phase 8.1.** The two-level outer loop and complete pivoting
kernel are implemented and correct (reconstruction < 1e-12 on all dense tests),
but the inner loop uses BLAS-2 (rank-1/rank-2 Schur updates) instead of blocked
BLAS-3 operations. This document captures what needs to change and why.

## Background

The two-level APTP algorithm (Duff, Hogg & Lopez 2020, Algorithm 3.1) is
designed to exploit BLAS-3 by decomposing each outer block into three phases:

1. **Factor** (BLAS-2): Factor the ib x ib diagonal block with complete pivoting
2. **Apply** (BLAS-3 TRSM): Solve L21 from the panel below, then threshold-check
3. **Update** (BLAS-3 GEMM): Schur complement A22 -= L21 D11 L21^T

The BLAS-3 payoff comes from steps 2 and 3, which operate on rectangular blocks
(panel height x ib) with high arithmetic intensity and good cache utilization.
BLAS-2 rank-1/rank-2 updates are memory-bandwidth-limited for large matrices.

### Current Architecture (BLAS-2)

`factor_inner` processes ib-sized sub-blocks but applies the Schur complement to
ALL remaining rows immediately after each pivot (via `try_1x1_pivot` /
`try_2x2_pivot`). Each pivot does a rank-1 or rank-2 update to the entire
trailing submatrix. Threshold checking happens at every individual pivot.

This means `apply_and_check`, `update_trailing`, and `update_delayed` are
implemented but unused. The three-function decomposition exists in code but the
control flow bypasses it.

### Why It Was Deferred

The block-level backup/restore approach was causing correctness issues during
initial development: when `factor_inner` only processes the diagonal block and
then `apply_and_check` solves the panel, a threshold failure requires restoring
the entire block and re-factoring with failed columns delayed. The backup was
corrupting already-Schur-complemented values in the trailing matrix, leading to
reconstruction errors of O(1).

The per-pivot approach sidesteps this by checking each pivot immediately and only
undoing one column's worth of work on failure. This is correct and tested but
leaves BLAS-3 performance on the table.

## What Needs to Change

### 1. Limit `factor_inner` to the Diagonal Block

Currently `factor_inner` receives the full remaining submatrix and applies Schur
updates to all rows below the diagonal block. It should instead only process the
ib x ib diagonal block, producing:
- L11 (unit lower triangular, ib x ib)
- D11 (1x1 and 2x2 pivots via MixedDiagonal)
- Permutation within the block
- nelim (number of successful pivots, may be < ib)

The rows below the diagonal block should be untouched by this step.

### 2. Apply Step (TRSM + Threshold Check)

After factoring the diagonal block, solve for the panel L21:

```
L21 = A21 P (L11 D11)^{-T}
```

This is a TRSM call on a (remaining_rows x nelim) panel. In SPRAL this is
`apply_pivot` in `ldlt_app.cxx:332`, which calls `host_trsm` (SIDE_RIGHT,
FILL_MODE_LWR, OP_T, DIAG_UNIT) followed by per-column D scaling.

The faer equivalent:
```rust
faer::linalg::triangular_solve::solve_unit_lower_triangular_in_place(
    l11.as_ref().transpose(),  // L11^T on the right
    panel.as_mut(),            // A21 * P → L21
    Par::Seq,
);
// Then scale each column by D^{-1} (handle 1x1 and 2x2 blocks)
```

After TRSM, scan L21 entries against the APTP threshold bound:
- 1x1 pivots: |L_ij| <= u^{-1} (default u = 0.01, so bound = 100)
- 2x2 pivots: |L_ij| <= sqrt(2)

The first column where ANY row fails determines the effective nelim for this
block. This is the `apply_and_check` function already implemented in
`src/aptp/factor.rs:941`.

### 3. Block-Level Backup/Restore on Threshold Failure

If threshold check reduces nelim:

1. Restore the diagonal block from backup (the `BlockBackup` struct already
   exists at `src/aptp/factor.rs:878`)
2. The trailing matrix has NOT been modified yet (GEMM hasn't run), so no
   trailing restore is needed — this is the key insight that makes block-level
   backup/restore safe
3. Re-factor the diagonal block with the failed columns excluded (swapped to
   end of the block)
4. Re-apply TRSM with the reduced nelim

SPRAL implements this via `restore_if_required` in `ldlt_app.cxx:1488` and
the `CopyBackup`/`PoolBackup` classes at `ldlt_app.cxx:416-900`.

### 4. Update Step (GEMM)

Once all threshold checks pass (nelim is final):

```
W = L21 * D11           (scale panel columns by D entries)
A22 -= W * L21^T        (symmetric rank-nelim update)
```

The faer equivalent:
```rust
// Build W = L21 * D11 in workspace
// Then:
faer::linalg::matmul::matmul(
    trailing.as_mut(),     // A22
    Accum::Add,            // beta = 1
    w.as_ref(),            // L21 * D11
    l21.as_ref().transpose(), // L21^T
    -1.0,                  // alpha = -1
    Par::Seq,
);
```

This is the `update_trailing` function already implemented at
`src/aptp/factor.rs:1073`.

### 5. Delayed Column Handling Across Blocks

When columns fail the threshold check in one outer block, they are delayed to
the next block. SPRAL tracks this per-column via `ColumnData` and handles the
update of delayed columns from earlier blocks in the `UpdateNT` and `UpdateTN`
tasks (Duff et al. 2020, Algorithm 3.1 lines 10-13).

The `update_delayed` function at `src/aptp/factor.rs:1142` is the placeholder
for this. It needs to apply the Schur complement from the current block's
factorization to any previously-delayed columns that now appear in the current
block's range.

## Correct Ordering

The ordering that avoids backup corruption:

```
for each outer block j:
    1. backup(diagonal_block[j])
    2. factor(diagonal_block[j])           → nelim_j, L_jj, D_jj, P_j
    3. for each panel block:
         apply_perm_and_backup(panel[i,j])
         apply_pivot(panel[i,j])           → TRSM, threshold check
         update_passed(nelim_j)            → atomic min reduction
    4. adjust(nelim_j)                     → finalize, handle 2x2 boundary
    5. for each trailing block:
         restore_if_required(block)        → undo failed columns
         update(block)                     → GEMM: A -= L*D*L^T
```

Steps 1-4 may iterate if nelim_j is reduced by threshold failures. Step 5 only
runs after nelim_j is final, so the trailing matrix is never corrupted.

## SPRAL Reference

The complete BLAS-3 implementation lives in:
- `references/spral/src/ssids/cpu/kernels/ldlt_app.cxx` — `run_elim_pivoted()`
  (line 1273): task-based outer loop with Factor/Apply/Adjust/Update phases
- `references/spral/src/ssids/cpu/kernels/ldlt_app.cxx` — `apply_pivot()`
  (line 332): TRSM + D scaling for panel solve
- `references/spral/src/ssids/cpu/kernels/block_ldlt.hxx` — inner block
  factorization (complete pivoting at the ib x ib leaf level)
- `references/spral/src/ssids/cpu/kernels/ldlt_app.cxx` — `Block::update()`
  uses `host_gemm` for the trailing matrix GEMM

Key constants: `INNER_BLOCK_SIZE = 32` (line 41), outer `BLOCK_SIZE` is a
template parameter (typically 256).

SPRAL uses OpenMP tasks for parallelism across blocks within a single frontal
matrix. Each task type (Factor, ApplyN, ApplyT, UpdateNN, UpdateNT, UpdateTN)
has explicit `depend` clauses that enforce the correct ordering without
serializing independent blocks.

## Academic References

- **Duff, Hogg & Lopez (2020)**: "A new sparse LDLT solver using a posteriori
  threshold pivoting." Algorithm 3.1 defines the two-level blocked APTP with
  Factor/ApplyN/ApplyT/UpdateNN/UpdateNT/UpdateTN kernels. Section 3 describes
  the threshold checking and backup/restore strategy.
  File: `references/ssids/duff2020.md`

- **Hogg & Scott (2013)**: "A new pivoting strategy for solving sparse
  symmetric definite systems." Algorithm 4.1 defines complete pivoting within
  inner blocks (the leaf-level kernel used by Factor).
  File: `references/ssids/hogg2013.md`

- **Duff & Reid (1983)**: "The multifrontal solution of indefinite sparse
  symmetric linear equations." Foundation for the multifrontal framework.
  File: `references/ssids/duff1983.md`

## faer API Surface

The BLAS-3 operations map to faer as follows:

| Operation | faer function | Module |
|-----------|--------------|--------|
| TRSM (panel solve) | `solve_unit_lower_triangular_in_place` | `faer::linalg::triangular_solve` |
| GEMM (trailing update) | `matmul` | `faer::linalg::matmul` |
| Parallelism control | `Par::Seq` / `Par::rayon(n)` | `faer::Par` |

All accept `Par` for threading control (transparent composition with faer).

## Existing Code to Reuse

These functions in `src/aptp/factor.rs` were written for BLAS-3 and are
currently unused (`#[allow(dead_code)]`):

| Function | Line | Purpose |
|----------|------|---------|
| `apply_and_check` | 941 | TRSM + D scaling + threshold scan |
| `update_trailing` | 1073 | W = L21*D11, GEMM A22 -= W*L21^T |
| `update_delayed` | 1142 | Update previously-delayed columns |
| `BlockBackup` | 878 | Save/restore diagonal block region |

## Expected Performance Impact

With BLAS-2, two-level and single-level perform similarly (confirmed by
benchmarks). BLAS-3 should show significant improvement for large fronts
(n > 256) where the GEMM dominates. The crossover point depends on cache
size and BLAS implementation quality.

Benchmark sizes to test: 256, 512, 1024, 2048. The `kernel/two_level` and
`kernel/single_level` benchmark groups in `benches/solver_benchmarks.rs`
already measure this.
