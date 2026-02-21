# Implementation Plan: Two-Level APTP Factorization

**Branch**: `017-two-level-aptp` | **Date**: 2026-02-16 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/017-two-level-aptp/spec.md`

## Summary

Replace the single-level (column-by-column) dense APTP kernel with a two-level blocked implementation that uses BLAS-3 operations (TRSM/GEMM) for large frontal matrices. The outer loop processes blocks of nb=256 columns; within each block, an inner loop processes sub-blocks of ib=32 columns using APTP, with the innermost ib×ib diagonal blocks factored using complete pivoting (Algorithm 4.1, Duff et al. 2020). Per-column backup is replaced with per-block backup. The multifrontal caller (`AptpNumeric::factor`) is unchanged — dispatching between single-level (front ≤ nb) and two-level (front > nb) happens inside the kernel.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (dense LA: matmul, triangular_solve, Mat/MatMut/MatRef), metis-sys 0.3.x (ordering, existing)
**Storage**: In-memory dense matrices (Mat<f64>) per supernode; MixedDiagonal for D
**Testing**: cargo test (unit + integration), cargo test -- --ignored (full SuiteSparse), criterion benchmarks
**Target Platform**: Linux x86_64 (primary), any platform faer supports
**Project Type**: Single Rust library crate
**Performance Goals**: Two-level faster than single-level on fronts > 256; identify crossover point via benchmarks
**Constraints**: Reconstruction error < 1e-12; backward error no regression vs Phase 7 baseline; sequential only (parallelism deferred to Phase 8.2)
**Scale/Scope**: 91% of SuiteSparse supernodal matrices have max_front > 256; top 10 fronts range 3,249–11,125

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Reconstruction < 1e-12 and backward error no-regression are primary acceptance criteria. Performance optimization must not compromise accuracy. |
| II. Clean Room | PASS | Algorithm from Duff, Hogg & Lopez (2020) paper + SPRAL (BSD-3) reference. No HSL or GPL sources. All references documented in spec. |
| III. TDD | PASS | Tests exist for single-level kernel (935 lines in factor.rs + multifrontal.rs + solve.rs). Two-level must pass all existing tests plus new block-boundary tests. Complete pivoting gets its own unit tests. |
| IV. Documentation | PASS | Academic references cited in spec. New functions will include rustdoc with paper section references and SPRAL equivalents. |
| V. Numerical Stability | PASS | Complete pivoting at innermost level is provably stable (u=0.25 equivalence, Duff2020 Section 4). APTP a posteriori checking maintained at outer level. |
| VI. Structured Development | PASS | Phase 8.1 follows Phase 7 completion. All Phase 7 exit criteria met. Plan document will be updated. |
| VII. Code Quality | PASS | Result types unchanged (AptpFactorResult); faer types at boundary; in-place operations; Result<T,E> error handling preserved. |

No violations. Gate passed.

## Project Structure

### Documentation (this feature)

```text
specs/017-two-level-aptp/
├── spec.md              # Feature specification (complete)
├── plan.md              # This file
├── research.md          # Phase 0: algorithm research & decisions
├── data-model.md        # Phase 1: type/struct changes
├── quickstart.md        # Phase 1: development guide
├── contracts/
│   └── api-surface.md   # Phase 1: public API changes
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code (repository root)

```text
src/aptp/
├── factor.rs            # PRIMARY: Two-level APTP kernel + complete pivoting
│                        #   - Refactor aptp_factor_in_place with dispatch logic
│                        #   - New: two_level_factor (outer block loop)
│                        #   - New: factor_inner (inner APTP on nb×nb diagonal)
│                        #   - New: complete_pivoting_factor (Algorithm 4.1 for ib×ib)
│                        #   - New: apply_and_check (TRSM-based L21 computation + threshold check)
│                        #   - New: update_trailing (GEMM-based Schur complement)
│                        #   - New: BlockBackup (per-block backup/restore)
│                        #   - Existing: swap_symmetric, extract_l (unchanged)
├── numeric.rs           # MINIMAL CHANGE: pass extended options if AptpOptions changes
├── solver.rs            # MINIMAL CHANGE: expose block size options in FactorOptions
├── diagonal.rs          # UNCHANGED
├── pivot.rs             # UNCHANGED
├── solve.rs             # UNCHANGED
├── symbolic.rs          # UNCHANGED
├── ordering.rs          # UNCHANGED
├── matching.rs          # UNCHANGED
├── inertia.rs           # UNCHANGED
├── perm.rs              # UNCHANGED
└── mod.rs               # MINIMAL: re-export any new public types

tests/
├── multifrontal.rs      # EXTEND: add block-boundary test cases
├── solve.rs             # EXTEND: verify no regression with two-level kernel
└── (factor.rs embedded) # EXTEND: complete pivoting unit tests, two-level tests
```

**Structure Decision**: All changes are concentrated in `src/aptp/factor.rs` with the two-level kernel implementation. The multifrontal caller (`numeric.rs`) and solver API (`solver.rs`) require only minimal option-plumbing changes. No new files needed — all new functions are added to `factor.rs`.

## Design Decisions

### D1: Dispatch location — inside kernel, not at caller

The two-level dispatch (front > nb → two-level, else single-level) happens inside `aptp_factor_in_place`, not in `AptpNumeric::factor`. This keeps the multifrontal code unchanged and the kernel self-contained.

**Rationale**: The multifrontal factorization already passes the full frontal matrix to the kernel. The kernel knows its own dimension and can dispatch accordingly. Pushing dispatch up would require the caller to understand internal blocking details.

### D2: Inner APTP uses complete pivoting at ib×ib leaves with blocked updates

The inner Factor phase (processing nb×nb diagonal blocks) loops over ib-sized sub-blocks. For each sub-block: (1) `complete_pivoting_factor` processes the entire ib×ib diagonal block at once (Algorithm 4.1 — no column-by-column processing within the sub-block), (2) an inner Apply (TRSM) updates the remaining columns of the nb block below the diagonal, (3) an inner Update (GEMM) applies the Schur complement within the nb block. The existing `try_1x1_pivot`/`try_2x2_pivot` functions are NOT used within `factor_inner` — complete pivoting replaces them entirely at this level.

**Rationale**: Algorithm 4.1 processes all ib columns simultaneously (searching the entire remaining submatrix for the maximum entry), which is fundamentally different from the column-by-column threshold approach. The ib×ib block fits in L1 cache (32×32 doubles = 8KB), making full-matrix search negligible in cost. This matches SPRAL's `block_ldlt()` at the innermost level of its Factor kernel.

### D3: Complete pivoting as new function, not refactored try_1x1/try_2x2

Complete pivoting (Algorithm 4.1) is fundamentally different from the current threshold-check approach: it searches the entire remaining submatrix for the largest entry, then chooses 1x1 or 2x2 based on determinant conditions. This is a new function `complete_pivoting_factor` that takes an ib×ib block and returns the same `AptpFactorResult` structure. The existing `try_1x1_pivot`/`try_2x2_pivot` functions become unused after the two-level kernel replaces the single-level main loop (they may be retained temporarily for comparison testing or removed entirely).

**Rationale**: Algorithm 4.1 has different control flow (search entire matrix, then decide pivot type) vs current code (try diagonal first, then partner search on failure). A single new function is cleaner than conditional logic inside existing functions. Within `factor_inner`, all ib columns are processed at once by `complete_pivoting_factor`; the remaining columns within the nb block are updated by inner Apply (TRSM) and inner Update (GEMM), not by column-by-column pivoting.

### D4: Per-block backup scope

Before each outer block j is factored, the following regions are backed up:
- The diagonal block A[j*nb..(j+1)*nb, j*nb..(j+1)*nb]
- Each sub-diagonal block A[i*nb..(i+1)*nb, j*nb..(j+1)*nb] for i > j
- Each super-diagonal block (transpose, for previously-failed columns) if applicable

On failure (nelim < nb), the failed columns (nelim..nb of the block) are restored from backup. Successfully factored columns (0..nelim) are retained.

**Rationale**: This matches SPRAL's CopyBackup approach. Backing up only the affected column of blocks (not the entire matrix) keeps memory overhead proportional to the block size, not the matrix size. PoolBackup (dynamic allocation) is deferred — CopyBackup is simpler and sufficient for sequential execution.

### D5: Single-level code replacement strategy

Once two-level is validated, the single-level main loop in `aptp_factor_in_place` will be replaced by the two-level implementation. For fronts ≤ nb, the two-level loop executes exactly one iteration with a single block equal to the full matrix — this IS the single-level case. No separate code path is maintained.

**Rationale**: User explicitly requested no separate single-level code path. The two-level code with nblk=1 degenerates to single-level behavior. The inner APTP (within the single block) reuses the existing column-by-column logic.

### D6: AptpOptions extension

Add `outer_block_size: usize` (default 256) and `inner_block_size: usize` (default 32) to `AptpOptions`. These are passed through from `FactorOptions` → `AptpOptions` → kernel. The dispatch threshold equals `outer_block_size`.

**Rationale**: Making block sizes configurable allows benchmarking different values and tuning for different architectures. Using AptpOptions (existing type) avoids introducing a new TwoLevelAptpOptions struct — the two-level algorithm IS the APTP algorithm, just with blocking parameters.

## Algorithm Architecture

### Two-Level APTP Flow (per frontal matrix)

```
aptp_factor_in_place(a: m×m, p: num_fully_summed, options)
│
├─ if p ≤ options.outer_block_size:
│    // Single block = single-level behavior
│    factor_inner(a[0..p, 0..m], p, options)
│    // (reuses existing loop with complete pivoting at leaf)
│
└─ else:
     nblk = ceil(p / nb)
     global_nelim = 0
     global_perm = [0..m]

     for j in 0..nblk:
       col_start = global_nelim + j * nb  // adjusted for prior delays
       block_cols = min(nb, remaining_fully_summed)

       // 1. BACKUP: Save A[col_start:, col_start:col_start+block_cols]
       backup = BlockBackup::create(a, col_start, block_cols, m)

       // 2. FACTOR: Inner APTP on diagonal block
       //    Uses complete pivoting at ib×ib leaves
       block_result = factor_inner(
           a[col_start..col_start+block_cols, col_start..col_start+block_cols],
           block_cols, options
       )
       block_nelim = block_result.num_eliminated

       // 3. APPLY: L21 = A21 * (L11 * D11)^{-T} via TRSM
       //    Check a posteriori: reduce block_nelim if any |L21[i,j]| > 1/u
       block_nelim = apply_and_check(
           a, col_start, block_nelim, block_cols, m,
           &block_result.d, options.threshold
       )

       // 4. ADJUST: If block_nelim < block_cols, restore failed portion
       if block_nelim < block_cols:
           backup.restore_failed(a, col_start, block_nelim, block_cols, m)

       // 5. UPDATE (trailing submatrix): A22 -= L21 * D * L21^T via GEMM
       update_trailing(a, col_start, block_nelim, m, &block_result.d)

       // 6. UPDATE (previously-delayed columns): UpdateNT/UpdateTN
       update_delayed(a, col_start, block_nelim, delayed_range, &block_result.d)

       // 7. Accumulate into global result
       global_nelim += block_nelim
       // Track delayed columns from this block

     // Final: permute all delayed columns to end
     // Build AptpFactorResult from accumulated blocks
```

### Complete Pivoting (Algorithm 4.1) — Innermost Level

```
complete_pivoting_factor(a: ib×ib block, small: f64) -> AptpFactorResult
│
├─ for each uneliminated column:
│    1. Find (t, m) = argmax |a[i,j]| over all uneliminated entries
│    2. if |a[t,m]| < small: mark all remaining as zero pivots; done
│    3. if t == m: use as 1×1 pivot (diagonal maximum)
│    4. else:
│       Δ = a[m,m] * a[t,t] - a[m,t]²
│       if |Δ| ≥ 0.5 * |a[m,t]|²:
│           use (t,m) as 2×2 pivot
│       else:
│           use max(|a[m,m]|, |a[t,t]|) as 1×1 pivot
│    5. Swap pivot row/col to current position
│    6. Compute L column(s) and D entry/block
│    7. Apply Schur complement update
│
└─ Return: D, permutation, num_eliminated (always == ib for complete pivoting,
          unless matrix is near-singular)
```

### BLAS-3 Operations

**Apply phase (TRSM)**:
```
L21 = A[below, block_cols] * inv(L11^T)    // triangular solve
L21 = L21 * inv(D11)                        // diagonal scaling
```
Using faer: `solve_unit_lower_triangular_in_place` for L11^{-T}, then scale by D11^{-1}.

**Update phase (GEMM)**:
```
A[trailing, trailing] -= L21 * D11 * L21^T
```
Using faer: compute `W = L21 * D11` in workspace, then `matmul(A_trailing, W, L21^T, -1.0, 1.0)`.

## Complexity Tracking

> No constitution violations — table not needed.

## Post-Design Constitution Re-Check

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | All existing tests must pass. New complete pivoting tests added. Reconstruction < 1e-12 verified on all test matrices. |
| II. Clean Room | PASS | Algorithm 4.1 from Duff et al. 2020; backup strategy inspired by SPRAL CopyBackup (BSD-3). |
| III. TDD | PASS | Complete pivoting unit tests written before implementation. Block-boundary tests before two-level integration. |
| IV. Documentation | PASS | All new functions cite Duff et al. 2020 section numbers. SPRAL equivalents noted. |
| V. Numerical Stability | PASS | Complete pivoting provably stable (u=0.25). A posteriori check maintained at outer level. |
| VI. Structured Development | PASS | Phase 8.1 builds on validated Phase 7 solver. |
| VII. Code Quality | PASS | Same output types, same error handling. New code follows existing patterns. |

Gate passed. Design is constitution-compliant.
