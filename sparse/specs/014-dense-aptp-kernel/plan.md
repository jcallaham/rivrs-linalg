# Implementation Plan: Dense APTP Factorization Kernel

**Branch**: `014-dense-aptp-kernel` | **Date**: 2026-02-14 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/014-dense-aptp-kernel/spec.md`

## Summary

Implement the A Posteriori Threshold Pivoting (APTP) algorithm for dense symmetric indefinite matrices — the core numerical kernel for the SSIDS multifrontal solver. The kernel factors a dense frontal matrix in place using an optimistic 1x1 pivot strategy with a posteriori stability checking, falling back to 2x2 Bunch-Kaufman pivots or column delay when stability bounds are violated. Uses faer's dense BLAS for Schur complement updates and triangular solves. Stores the D factor in Phase 2's MixedDiagonal type and returns delayed columns for the caller (Phase 6 multifrontal assembly) to propagate up the elimination tree.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (dense matrix types, matmul, triangular solve), Phase 2 types (MixedDiagonal, PivotType, Block2x2, Inertia)
**Storage**: N/A (in-memory dense matrices only)
**Testing**: cargo test (unit + integration), NumericalValidator (reconstruction error < 10^-12), random matrix generators (test-util feature), hand-constructed + SuiteSparse CI-subset matrices
**Target Platform**: Linux (development), cross-platform Rust
**Project Type**: Single Rust library crate (existing sparse/ project)
**Performance Goals**: Correctness-first; no explicit performance targets for Phase 5. Uses faer BLAS from day one for cache-optimal operations. Two-level blocking optimization deferred to Phase 9.1.
**Constraints**: Reconstruction error < 10^-12, stability bound |l_ij| < 1/threshold enforced, no heap allocation in hot path (where practical)
**Scale/Scope**: Dense matrices from 1x1 to ~500x500 (frontal matrix sizes in typical multifrontal factorization). Single new file (~500-800 lines) + tests.

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Evidence |
|-----------|--------|----------|
| I. Correctness First | PASS | Reconstruction error < 10^-12 as primary validation. SC-001 through SC-008 define measurable correctness targets. TDD workflow: tests before implementation. |
| II. Clean Room | PASS | Algorithm from Duff, Hogg & Lopez (2020) — academic paper. SPRAL (BSD-3) consulted for patterns. No HSL or GPL code. All references documented in spec. |
| III. Test-Driven Development | PASS | Constitution-mandated test categories all addressed: reconstruction tests, unit tests (hand-constructed), property-based (inertia, statistics invariants), edge cases (1x1, zero, all-delayed). |
| IV. Algorithm Documentation | PASS | Academic references cited in spec and research.md. Rustdoc will include algorithm description, SPRAL equivalents, complexity analysis. |
| V. Numerical Stability | PASS | APTP with threshold pivoting, 2x2 BK fallback, delayed columns. Stability bound enforced. Inertia reported. Singularity detection via `small` parameter. |
| VI. Structured Development | PASS | Phase 5 follows Phase 4 (complete). Exit criteria defined. No skipped phases. |
| VII. Code Quality | PASS | Uses faer types at boundary. Result types for errors. No .unwrap() in library code. Rustdoc on all public items. |

**Post-design re-check**: All principles remain satisfied. The in-place API design (research R3) supports both correctness testing (convenience wrapper) and performance (Phase 6 integration). No complexity violations.

## Project Structure

### Documentation (this feature)

```text
specs/014-dense-aptp-kernel/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0: faer APIs, SPRAL patterns, design decisions
├── data-model.md        # Phase 1: Entity definitions and relationships
├── quickstart.md        # Phase 1: Build/test/usage guide
├── contracts/
│   └── api.md           # Phase 1: Public API contract
├── checklists/
│   └── requirements.md  # Specification quality checklist
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code

```text
src/
├── aptp/
│   ├── mod.rs            # MODIFIED: add `pub mod factor;` + re-exports
│   ├── factor.rs         # NEW: APTP factorization kernel
│   │   ├── AptpOptions, AptpFallback (configuration)
│   │   ├── AptpFactorResult (in-place result)
│   │   ├── AptpFactorization (convenience result)
│   │   ├── AptpStatistics, AptpPivotRecord (diagnostics)
│   │   ├── aptp_factor_in_place() (core algorithm)
│   │   ├── aptp_factor() (convenience wrapper)
│   │   └── internal: try_1x1_pivot, try_2x2_pivot, select_2x2_partner,
│   │       update_schur_1x1, update_schur_2x2, extract_l
│   ├── diagonal.rs       # UNCHANGED: MixedDiagonal (Phase 2)
│   ├── pivot.rs          # UNCHANGED: PivotType, Block2x2 (Phase 2)
│   ├── inertia.rs        # UNCHANGED: Inertia (Phase 2)
│   ├── perm.rs           # UNCHANGED: perm_from_forward (Phase 2)
│   ├── symbolic.rs       # UNCHANGED: AptpSymbolic (Phase 3)
│   ├── ordering.rs       # UNCHANGED: METIS ordering (Phase 4.1)
│   └── matching.rs       # UNCHANGED: MC64 matching (Phase 4.2)
├── error.rs              # POSSIBLY MODIFIED: add factorization-specific error variants if needed
└── ...                   # All other modules unchanged
```

**Structure Decision**: Single new file `src/aptp/factor.rs` added to the existing `aptp` module. This follows the established pattern where each Phase adds a module file to `src/aptp/`. No new directories or restructuring needed.

## Algorithm Design

### Column-by-Column APTP (Single-Level)

Phase 5 implements a **single-level** (non-blocked) version of the APTP algorithm. The full two-level blocked version (Algorithm 3.1 with outer block size nb and inner block size ib) is deferred to Phase 9.1.

**Processing order**: Columns are processed left-to-right (k = 0, 1, ..., p-1 where p = num_fully_summed).

**For each column k**:

1. **Compute L column**: For all rows i > k in the active set:
   ```
   l_ik = (a_ik - sum_{s < k, s eliminated} l_is * d_s * l_ks) / d_kk
   ```
   In practice, the Schur complement updates from prior pivots have already been applied to column k, so the formula simplifies to: `l_ik = a_ik / d_kk` (where a_ik is the already-updated value).

2. **Check stability**: `max_i |l_ik| < 1/threshold`

3. **On pass (1x1 pivot)**:
   - Record `PivotType::OneByOne` in MixedDiagonal
   - Update trailing submatrix: `a_ij -= l_ik * d_kk * l_jk` for all i, j > k
   - Advance to column k+1

4. **On fail (fallback)**:
   - **BunchKaufman**: Select partner column (largest off-diagonal entry in column k among uneliminated columns). Attempt 2x2 pivot with columns k and partner.
     - If 2x2 passes: record `PivotType::TwoByTwo`, swap partner column adjacent, update trailing submatrix with rank-2 update, advance by 2
     - If 2x2 fails: mark both columns as `PivotType::Delayed`
   - **Delay**: Mark column k as `PivotType::Delayed` immediately

5. **After all columns processed**: Permute delayed columns to end of ordering. Return result.

### Schur Complement Update

For a **1x1 pivot** at column k with diagonal value d_kk:
```
a_ij -= l_ik * d_kk * l_jk   for all active i, j > k
```
This is a symmetric rank-1 update. Use `faer::linalg::matmul::matmul_with_conj` with the L column as both lhs and rhs (transposed), alpha = -d_kk.

For a **2x2 pivot** at columns k, k+1 with block D_22:
```
a_ij -= [l_ik, l_i(k+1)] * D_22 * [l_jk, l_j(k+1)]^T   for all active i, j > k+1
```
This is a symmetric rank-2 update. Precompute `W = L_cols * D_22` (m x 2), then `A -= W * L_cols^T` using matmul.

### 2x2 Partner Selection

When a 1x1 pivot at column k fails, select the partner for a 2x2 attempt:
- Search the remaining uneliminated fully-summed columns for the one with the largest absolute off-diagonal entry `|a_ik|` (i != k)
- This follows the Bunch-Kaufman strategy: maximize the off-diagonal element of the 2x2 block
- The stability test for the 2x2 block (from Algorithm 4.1): `|det(D_22)| >= 0.5 * |a_21|^2`

### Delayed Column Handling

Delayed columns remain in the matrix but are not eliminated. After all passes:
1. Collect indices of all delayed columns
2. Permute the column ordering so eliminated columns come first, delayed columns last
3. The contribution block (rows/cols beyond num_fully_summed) has been updated by all successful pivots' Schur complement updates — delayed columns' contributions are NOT subtracted (they weren't eliminated)

### In-Place Storage Layout

After factorization, the input matrix contains:
```
For eliminated column k (in permuted order):
  a[k, k]    = d_kk (1x1 pivot value) — stored in MixedDiagonal, not in matrix
  a[i, k]    = l_ik for i > k (L factor entries)
  a[k, j]    = l_jk for j < k (L factor entries, stored as column below diagonal)

For 2x2 pivot columns k, k+1:
  a[k:k+2, k:k+2] = D_22 entries — stored in MixedDiagonal
  a[i, k], a[i, k+1] = L entries for i > k+1

For delayed column k:
  a[*, k] = partially updated column (Schur complement from prior pivots applied)

Contribution block (rows/cols >= num_fully_summed):
  Updated with Schur complement from all successful pivots
```

## Testing Strategy

### Test Categories (Constitution III compliance)

1. **Reconstruction tests** (primary oracle):
   - `||A - P^T L D L^T P|| / ||A|| < 10^-12` for all test matrices
   - For partial factorizations: verify eliminated portion only
   - For full factorizations: verify entire matrix

2. **Unit tests** (hand-constructed matrices with known factorizations):
   - 2x2 diagonal matrix (trivial 1x1 pivots)
   - 3x3 positive definite (all 1x1, no delays)
   - 4x4 indefinite requiring 2x2 pivot
   - 5x5 arrow matrix with known LDL^T
   - Matrix forcing all-delayed (zero matrix, near-singular)
   - 1x1 matrix (trivial case)

3. **Property-based tests** (structural invariants):
   - `num_1x1 + 2*num_2x2 + num_delayed == num_fully_summed` always
   - Permutation is valid (bijection on 0..n)
   - L has unit diagonal for eliminated columns
   - D has correct PivotType for each column
   - Inertia from D matches expected (for known matrices)

4. **Edge case tests**:
   - 1x1 matrix (single element)
   - Zero matrix (all delayed)
   - Matrix with one near-zero diagonal (tests `small` threshold)
   - Threshold = 0.01 (default), 0.5 (strict), edge values
   - Both fallback strategies (BunchKaufman, Delay) on same matrix

5. **Random matrix stress tests** (via generators):
   - 100+ random symmetric positive definite (all 1x1, no delays)
   - 100+ random symmetric indefinite (mixed pivots)
   - Various sizes: 5x5, 20x20, 50x50, 100x100

6. **Integration tests with existing test infrastructure**:
   - Hand-constructed matrices from test-data/ (15 matrices, load as dense)
   - SuiteSparse CI-subset (10 matrices, extract dense submatrices)

### Reconstruction Test Formula for Partial Factorizations

For a full factorization (num_fully_summed = n, num_eliminated = q):
```
Error = ||A[perm, perm] - L[:, :q] * D[:q, :q] * L[:, :q]^T|| / ||A||
```
where L is unit lower triangular with q eliminated columns and the remaining n-q columns are delayed (uneliminated).

For partial factorization (frontal matrix m x m, num_fully_summed = p, num_eliminated = q):
- The first q columns are factored: L, D entries valid
- Columns q..p are delayed: will be handled by parent node
- Rows/cols p..m are the contribution block: updated in place

The reconstruction check verifies the eliminated portion. The contribution block update is verified separately by comparing the updated trailing submatrix against the expected Schur complement.

## faer API Usage

### Key imports

```rust
use faer::{Mat, MatRef, MatMut, Col, ColMut, Par};
use faer::linalg::matmul::matmul_with_conj;
use faer::linalg::matmul::rank_update;
use faer::Accum;
use faer::Conj;
use faer::perm::Perm;
```

### Schur complement (rank-1 update for 1x1 pivot)

```rust
// After eliminating column k with d_kk:
// A[k+1:, k+1:] -= l_k * d_kk * l_k^T
// where l_k = a[k+1:, k] (already computed L column)

let l_col = a.rb().col(k).subrows(k + 1, m - k - 1);
let mut trailing = a.rb_mut().submatrix_mut(k + 1, k + 1, m - k - 1, m - k - 1);

// rank-1 symmetric update
rank_update(trailing, Accum::Add, l_col, Conj::No, l_col.transpose(), Conj::No, &(-d_kk), Par::Seq);
```

### Schur complement (rank-2 update for 2x2 pivot)

```rust
// After eliminating columns k, k+1 with D_22:
// A[k+2:, k+2:] -= L_2 * D_22 * L_2^T
// where L_2 = a[k+2:, k:k+2] (two L columns)

let l_cols = a.rb().submatrix(k + 2, k, m - k - 2, 2);  // (m-k-2) x 2

// Precompute W = L_2 * D_22  (a (m-k-2) x 2 matrix)
// Then: A -= W * L_2^T
matmul_with_conj(trailing, Accum::Add, w, Conj::No, l_cols.transpose(), Conj::No, -1.0, Par::Seq);
```

## Complexity Tracking

> No constitution violations to justify. All design decisions align with established principles.

| Decision | Rationale |
|----------|-----------|
| In-place + convenience dual API | In-place required for Phase 6; convenience for testing. Not over-engineering — both are thin. |
| `small` parameter (separate from threshold) | Matches SPRAL; prevents catastrophic division by near-zero. Distinct purpose from stability bound. |
| Column-by-column (no blocking) | Phase 5 scope. Blocking adds complexity without benefit at small matrix sizes. Phase 9.1 adds two-level. |
