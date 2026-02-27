# Implementation Plan: Triangular Solve & Solver API

**Branch**: `016-triangular-solve-api` | **Date**: 2026-02-16 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/016-triangular-solve-api/spec.md`

## Summary

Implement Phase 7 of the SSIDS solver: the per-supernode triangular solve (forward, diagonal, backward substitution) and the user-facing `SparseLDLT` solver API wrapping the analyze → factor → solve pipeline. This is the "working solver" milestone that transforms the multifrontal factorization from Phases 2-6 into an end-to-end solver validated by backward error on the full test matrix suite.

The core solve is a free function `aptp_solve()` that traverses the assembly tree performing gather/solve/scatter operations per supernode using faer's dense triangular solve and matmul kernels. `SparseLDLT` wraps this with permutation, optional MC64 scaling, and a three-phase (analyze/factor/solve) user API.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (dense LA, sparse types, MemStack, triangular solve, matmul), metis-sys 0.3.x (METIS ordering), serde/serde_json (existing)
**Storage**: In-memory dense matrices per supernode (`Mat<f64>`), `MixedDiagonal` for D storage
**Testing**: `cargo test` (unit + integration), `cargo test -- --ignored --test-threads=1` (full SuiteSparse)
**Target Platform**: Linux (x86_64), portable Rust
**Project Type**: Single Rust library crate
**Performance Goals**: Correctness-focused (Phase 7). No latency targets. MemStack-based solve (no heap allocation in solve hot path).
**Constraints**: Sequential only (parallelism deferred to Phase 8). Single-column RHS (batch deferred to Phase 8).
**Scale/Scope**: 15 hand-constructed matrices + 10 SuiteSparse CI + 67 full SuiteSparse. Matrices up to ~350K dimension.

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | **PASS** | Backward error < 10^-10 as primary validation. Full SuiteSparse test suite. Per-supernode index mapping unit tested before integration. |
| II. Clean Room | **PASS** | Solve algorithm from academic papers (Duff+2020, Hogg+2016, Liu 1992). SPRAL (BSD-3) consulted for scaling/rank-deficiency patterns. No HSL source. |
| III. TDD | **PASS** | Per-supernode gather/scatter/solve unit tests BEFORE integration. Backward error tests on all matrix classes. Edge case tests (rank-deficient, empty, simplicial). |
| IV. Documentation | **PASS** | All public types/functions get rustdoc with algorithm references, SPRAL equivalents, complexity. |
| V. Numerical Stability | **PASS** | APTP pivoting from Phase 5/6. Zero pivot handling (SPRAL convention). MixedDiagonal D-solve with 2x2 Cramer's rule. |
| VI. Structured Development | **PASS** | Phase 6 complete and merged. Phase 7 builds on existing validated types. |
| VII. Code Quality | **PASS** | Result types for all errors. faer type composition. MemStack for solve workspace. |

No violations. No complexity tracking needed.

## Project Structure

### Documentation (this feature)

```text
specs/016-triangular-solve-api/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0: faer API research, algorithm design
├── data-model.md        # Phase 1: type/entity model
├── quickstart.md        # Phase 1: usage examples
├── contracts/
│   ├── solve-api.md     # SparseLDLT public API contract
│   └── internal-api.md  # aptp_solve + helper function contracts
└── tasks.md             # Phase 2 output (from /speckit.tasks)
```

### Source Code (repository root)

```text
src/
├── lib.rs                    # Add SolverPhase::Solve (existing enum)
├── error.rs                  # Add SolveBeforeFactor variant
├── validate.rs               # Add sparse_backward_error()
└── aptp/
    ├── mod.rs                # Re-export new types
    ├── solve.rs              # NEW: aptp_solve() + per-supernode helpers
    ├── solver.rs             # NEW: SparseLDLT, AnalyzeOptions, FactorOptions, SolverOptions
    ├── numeric.rs            # MODIFY: add Option<&[f64]> scaling to factor()
    └── diagonal.rs           # MODIFY: handle zero pivots in solve_in_place()

tests/
├── solve.rs                  # NEW: per-supernode unit tests + end-to-end backward error
├── solve_suitesparse.rs      # NEW: full SuiteSparse backward error (#[ignore])
├── hand_constructed.rs       # MODIFY: add backward error validation
└── multifrontal.rs           # EXISTING: Phase 6 reconstruction tests (unchanged)
```

**Structure Decision**: Single Rust library crate. New functionality in two new files (`solve.rs`, `solver.rs`) under the existing `aptp/` module, plus modifications to four existing files. Test files follow the established pattern of separate integration test files per phase.

## Key faer APIs Used

| Operation | faer API | Module |
|-----------|----------|--------|
| L11 * y = b (forward) | `solve_unit_lower_triangular_in_place_with_conj` | `faer::linalg::triangular_solve` |
| L11^T * z = y (backward) | `solve_unit_upper_triangular_in_place_with_conj` (with `l11.transpose()`) | `faer::linalg::triangular_solve` |
| tmp = L21 * y (scatter) | `matmul_with_conj` (Accum::Replace, alpha=1.0) | `faer::linalg::matmul` |
| x -= L21^T * tmp (gather update) | `matmul_with_conj` (Accum::Add, alpha=-1.0) | `faer::linalg::matmul` |
| Sparse A*x (backward error) | `sparse_dense_matmul` (two passes for symmetric: A_lower + A_lower^T - diag) | `faer::sparse::linalg::matmul` |
| Workspace sizing | `StackReq`, `MemBuffer`, `MemStack` | `faer::dyn_stack` |
| Temporary matrices | `temp_mat_zeroed`, `temp_mat_scratch` | `faer::linalg` |
| Permutation | `PermRef`, `Perm` | `faer::perm` |

## Algorithm Design

### Core Solve: `aptp_solve(symbolic, numeric, rhs, stack)`

Given P^T A P = L D L^T from multifrontal factorization, solve in the **permuted** coordinate system (P already applied by caller):

**Forward solve** (postorder traversal, s = 0..n_supernodes):
```
For each supernode s with ne = num_eliminated, r = |row_indices|:
  1. Gather: local[i] = rhs[col_indices[i]], i = 0..ne
  2. Solve: L11 * y = local (in-place, unit lower triangular)
  3. Write back: rhs[col_indices[i]] = local[i]
  4. Scatter (if r > 0):
     a. tmp = L21 * local     (matmul: r×ne · ne×1 → r×1)
     b. rhs[row_indices[j]] -= tmp[j], j = 0..r
```

**Diagonal solve** (any order, s = 0..n_supernodes):
```
For each supernode s:
  1. Gather: local[i] = rhs[col_indices[i]], i = 0..ne
  2. d11.solve_in_place(&mut local)  (handles 1x1, 2x2, and zero pivots)
  3. Write back: rhs[col_indices[i]] = local[i]
```

**Backward solve** (reverse postorder, s = n_supernodes-1..0):
```
For each supernode s:
  1. Gather: local[i] = rhs[col_indices[i]], i = 0..ne
  2. Gather (if r > 0): tmp[j] = rhs[row_indices[j]], j = 0..r
  3. Update (if r > 0): local -= L21^T * tmp  (matmul: ne×r · r×1 → ne×1, alpha=-1)
  4. Solve: L11^T * z = local (in-place, via l11.transpose() + unit upper triangular)
  5. Write back: rhs[col_indices[i]] = local[i]
```

**Why no explicit `local_perm` during solve**: `col_indices` is constructed during `extract_front_factors` as `frontal.row_indices[local_perm[i]]`, so it already maps factored positions to global permuted indices. The APTP reordering is baked into `col_indices` — no runtime permutation application needed.

**Workspace**: Per supernode, need `ne + r` f64 entries (local + tmp buffers). Since supernodes are processed sequentially, total workspace = `max(ne + r)` over all supernodes ≤ `max_front_size`. Computed from `AptpNumeric::stats().max_front_size`.

### SparseLDLT Wrapper

`SparseLDLT::solve_in_place()` wraps `aptp_solve` with:
1. **Permute**: `rhs_perm[perm_fwd[i]] = rhs[i]` (or via `faer::perm`)
2. **Scale** (if scaling present): `rhs_perm[i] *= scaling[i]` (scaling in elimination order)
3. **Core solve**: `aptp_solve(symbolic, numeric, &mut rhs_perm, stack)`
4. **Unscale**: `rhs_perm[i] *= scaling[i]` (symmetric: same S applied)
5. **Unpermute**: `rhs[i] = rhs_perm[perm_fwd[i]]`

**Scaling coordinate system**: MC64 scaling from `match_order_metis()` is in original matrix order. During `SparseLDLT::analyze()`, scaling is permuted to elimination order: `elim_scaling[i] = orig_scaling[perm_fwd[i]]`. This matches SPRAL's `fkeep%scaling(i) = scaling(akeep%invp(i))`.

### Modifications to Existing Code

1. **`AptpNumeric::factor()` — add scaling parameter**: Add `scaling: Option<&[f64]>` parameter. In `scatter_original_entries`, when scaling is Some, apply `val *= scaling[perm_inv[orig_row]] * scaling[perm_inv[orig_col]]` (scaling in elimination order applied to permuted indices). The APTP kernel and all other logic remain unchanged.

2. **`MixedDiagonal::solve_in_place()` — handle zero pivots**: Replace `debug_assert!(d != 0.0)` with graceful handling: if `d == 0.0` for 1x1 pivot, set `x[col] = 0.0`. For 2x2 blocks with `det == 0.0`, set both `x[col] = 0.0` and `x[partner] = 0.0`. This follows SPRAL's `action=true` convention.

3. **`validate.rs` — add sparse backward error**: New function `sparse_backward_error(a, x, b)` using `sparse_dense_matmul` instead of `a.to_dense()`. The existing `backward_error` is preserved for small matrices.

4. **`build_supernode_info` visibility**: Change from `pub(crate)` to confirm it's accessible from `solve.rs` (already `pub(crate)`, so accessible within `aptp` module).

## Complexity Tracking

No constitution violations. No complexity tracking needed.
