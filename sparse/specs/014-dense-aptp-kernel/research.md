# Research: Dense APTP Factorization Kernel

**Feature Branch**: `014-dense-aptp-kernel`
**Date**: 2026-02-14

## R1: faer 0.22 Dense BLAS APIs

### Decision: Use faer's `matmul_with_conj` for Schur complement, `solve_unit_lower_triangular_in_place_with_conj` for triangular solve

### Rationale

faer provides high-performance dense operations with the exact semantics needed:

**Matrix Multiplication** (`faer::linalg::matmul`):
- `matmul_with_conj(dst, beta, lhs, conj_lhs, rhs, conj_rhs, alpha, par)` — general GEMM
- `Accum::Add` with `alpha = -1.0` gives `dst -= lhs * rhs` (Schur complement update)
- `Accum::Replace` gives `dst = alpha * lhs * rhs` (fresh computation)
- `Par::Seq` for Phase 5 (parallelism deferred to Phase 9)

**Triangular Solve** (`faer::linalg::triangular_solve`):
- `solve_unit_lower_triangular_in_place_with_conj(L, conj, rhs, par)` — solves L*X = B in place
- Unit lower triangular is exactly what APTP produces (diagonal of L is 1)

**Submatrix Views** (`MatMut::submatrix_mut`):
- `mat.submatrix_mut(row_start, col_start, nrows, ncols)` — mutable subview for in-place updates
- Essential for updating trailing submatrix without copies

**Rank-1 Update** (`faer::linalg::matmul::rank_update`):
- Available for column-by-column Schur complement: `S -= l_k * d_kk * l_k^T`
- More efficient than matmul for single-column updates

### Alternatives Considered
- **Manual loops for Schur complement**: Rejected — faer's matmul uses SIMD and cache-optimal blocking
- **nalgebra**: Rejected — faer is the established dependency, already used for sparse infrastructure
- **BLAS via FFI**: Rejected — faer provides pure Rust equivalent, avoids FFI overhead and dependency

## R2: SPRAL Dense Kernel Patterns (BSD-3 Reference)

### Decision: Follow SPRAL's in-place factorization with backup/restore pattern, adapted for single-level (no blocking) Phase 5

### Rationale

SPRAL's CPU kernel (`ldlt_app.cxx`) provides a mature reference for the APTP algorithm:

**In-place storage pattern**: L factors stored in the lower triangle of the input matrix, D stored separately. This avoids allocating a separate L matrix and is cache-friendly. Phase 5 adopts this pattern.

**Backup/restore for fail-in-place**: SPRAL backs up each block before optimistic factorization, then restores on failure. For Phase 5's single-level approach (column-by-column), the backup is a single column vector, making it lightweight.

**2x2 pivot selection** (from SPRAL's `test_2x2`):
```
detscale = 1/|a21|
detpiv = (a11 * detscale) * a22 - |a21|
accept if |detpiv| >= |a21| / 2
```
This matches Algorithm 4.1's condition: `|Δ| >= ½|a_mt|²`.

**Threshold checking**: SPRAL checks `u * max(|L entries|) < 1.0`, equivalent to `max(|L entries|) < 1/u`. Default `u = 0.01`.

**Singularity detection**: SPRAL uses a separate `small` parameter (default 1e-20) to detect effectively zero pivots, distinct from the stability threshold.

### Alternatives Considered
- **Out-of-place factorization (separate L matrix)**: Rejected — doubles memory for large frontal matrices, unnecessary copy
- **faer's Bunch-Kaufman as primary algorithm**: Rejected — faer's BK is not APTP (no fail-in-place, no delayed columns). Used as reference for 2x2 pivot mechanics only.
- **Complete pivoting for all blocks**: Rejected for Phase 5 — complete pivoting searches entire remaining matrix (O(n²) per pivot). Column-by-column with BK fallback is sufficient. Complete pivoting is an option for the inner-block kernel in Phase 9.1's two-level blocking.

## R3: In-place vs Out-of-place API

### Decision: In-place factorization as core API, with convenience wrapper for testing

### Rationale

The core function mutates the input matrix in place (L stored in lower triangle, contribution block updated). This is essential for Phase 6 integration where frontal matrices are large and copying is expensive.

For standalone testing (Phase 5 validation), a convenience function copies the input, calls the in-place version, and extracts L into a separate `Mat<f64>`.

**Core API**: `aptp_factor_in_place(a: MatMut<f64>, num_fully_summed: usize, options: &AptpOptions) -> Result<AptpFactorResult>`
**Convenience**: `aptp_factor(a: MatRef<f64>, options: &AptpOptions) -> Result<AptpFactorization>`

### Alternatives Considered
- **Only out-of-place**: Rejected — Phase 6 needs in-place for performance
- **Only in-place**: Rejected — testing convenience matters; extracting L from a mutated matrix for reconstruction checks is error-prone
- **Builder pattern**: Rejected as premature — simple options struct suffices for Phase 5

## R4: Column Permutation Tracking

### Decision: Track column permutation as a `Vec<usize>` forward mapping within the factorization result

### Rationale

When 2x2 Bunch-Kaufman pivots are used, columns are swapped to bring the partner column adjacent. The factorization must track these swaps so that the reconstruction P^T A P = L D L^T can be verified.

SPRAL stores a `perm[m]` array mapping original column indices to factored column positions. Phase 5 does the same, using a `Vec<usize>` that can be converted to a `faer::perm::Perm<usize>` via `perm_from_forward()` (from Phase 2).

Delayed columns are permuted to the end of the column ordering (after all eliminated columns), following SPRAL's convention.

### Alternatives Considered
- **No permutation tracking (1x1 only)**: Rejected — 2x2 pivots require column swaps
- **faer Perm directly**: The factorization builds the permutation incrementally; storing as Vec<usize> during construction, then converting to Perm at the end is cleaner

## R5: Singularity Detection Threshold

### Decision: Include a separate `small` parameter (default 1e-20) for detecting effectively zero pivots

### Rationale

SPRAL distinguishes between:
- **Stability threshold** (`u = 0.01`): controls growth factor bound |l_ij| < 1/u
- **Singularity threshold** (`small = 1e-20`): detects pivots that are effectively zero

A pivot with |d_kk| < small is classified as zero (Delayed with PivotType::Delayed) rather than attempting division. This prevents catastrophic cancellation and infinite L entries.

The paper's Algorithm 4.1 uses δ for this purpose: "if |a_tm| < δ then all remaining pivots are zero."

### Alternatives Considered
- **Use machine epsilon**: Too strict — would flag many valid small pivots as zero
- **Relative threshold (small * ||A||)**: More robust but adds complexity. Deferred to hardening (Phase 10). Use absolute threshold for Phase 5, matching SPRAL.
