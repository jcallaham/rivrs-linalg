# Implementation Plan: Direct GEMM into Contribution Buffer

**Branch**: `024-direct-gemm-contrib` | **Date**: 2026-02-24 (revised) | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/024-direct-gemm-contrib/spec.md`

## Summary

Eliminate the dominant factorization bottleneck (40.1% of factor time on c-71) by restructuring the BLAS-3 blocking loop to defer the NFS×NFS Schur complement computation. Per-block trailing updates are restricted to the fully-summed region and cross-terms; a single GEMM after the loop computes the entire NFS×NFS Schur complement directly into a pre-allocated contribution buffer. This eliminates both the O(n²) copy and the per-supernode allocation, aligning with SPRAL's zero-copy architecture (`factor.hxx:92-103`).

Key research finding: Phase 9.1d proved that eliminating allocation overhead alone nets zero improvement (pool eliminated sys time 3.1s→0.5s but factor time was unchanged). The copy itself — 48 MB per large supernode — is the real cost. The only way to eliminate it is to avoid writing the Schur complement to the wrong place in the first place.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (dense LA, CSC, `tri_matmul`, `matmul`), rayon 1.x (parallel tree traversal), serde/serde_json (diagnostic export)
**Storage**: N/A (in-memory computation)
**Testing**: cargo test (358 unit tests + 65 SuiteSparse integration tests)
**Target Platform**: Linux x86_64 (primary), cross-platform (Rust)
**Project Type**: Single Rust crate (library)
**Performance Goals**: ExtractContr near zero; contribution GEMM visible as separate sub-phase. No fixed ratio target.
**Constraints**: No `unsafe` code; all 65 SuiteSparse matrices must pass backward error <5e-11
**Scale/Scope**: ~2 files modified (factor.rs + numeric.rs), estimated ~300-500 lines changed

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Correctness validated via existing 65-matrix test suite + backward error thresholds. The Schur complement is mathematically identical whether computed incrementally or in a single GEMM — only FP operation ordering changes (FR-008: relaxed to backward error tolerance). |
| II. Clean Room | PASS | Architecture follows SPRAL's `factor.hxx:92-103` pattern (BSD-3, freely consultable). No restricted code. |
| III. TDD | PASS | Existing tests serve as regression suite. New tests for: deferred GEMM correctness (small matrices with known Schur complements), delayed column interaction, contribution buffer lifecycle. |
| IV. Documentation | PASS | Academic references cited in spec.md. SPRAL source files documented. |
| V. Numerical Stability | PASS | APTP kernel, TRSM, threshold checking unchanged. The deferred GEMM computes the same Schur complement via a different BLAS call ordering. Backward error tolerance (not bit-identical) is the correctness standard. |
| VI. Structured Development | PASS | Phase 9.1e follows 9.1c/9.1d per ssids-plan.md. |
| VII. Code Quality | PASS | No `unsafe` code. Internal API changes only. |

**Post-design re-check**: All gates still pass. The deferred GEMM changes the FP operation ordering in the Schur complement but not the mathematical result. Backward error tolerance is the appropriate correctness standard (Constitution III numerical accuracy standards: backward error < 10^-10).

## Project Structure

### Documentation (this feature)

```text
specs/024-direct-gemm-contrib/
├── plan.md              # This file
├── spec.md              # Feature specification (revised)
├── research.md          # Phase 0: architectural decisions (revised)
├── data-model.md        # Phase 1: modified data structures
├── quickstart.md        # Phase 1: verification guide
├── contracts/
│   └── internal-api.md  # Phase 1: modified internal API contracts
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # Phase 2 output (/speckit.tasks)
```

### Source Code (modified files)

```text
src/aptp/
├── factor.rs            # MODIFIED: update_trailing restricted to FS region + cross-terms;
│                        #   new deferred contribution GEMM after blocking loop;
│                        #   aptp_factor_in_place returns contribution data in buffer
├── numeric.rs           # MODIFIED: FactorizationWorkspace (add contrib_buffer),
│                        #   extract_contribution becomes index-only (+ small delayed copy),
│                        #   extend_add/extend_add_mapped (return recycled buffer),
│                        #   factor_single_supernode (thread buffer, call deferred GEMM),
│                        #   factor_tree_levelset (swap orchestration)
└── solver.rs            # UNCHANGED: Public SparseLDLT API

examples/
└── profile_matrix.rs    # MINOR: Add ContribGEMM sub-phase display
```

**Structure Decision**: Single Rust crate, existing `src/aptp/` module. Changes span two files: `factor.rs` (GEMM restructuring) and `numeric.rs` (buffer management + orchestration). No new files or modules.

## Design Artifacts

### Research (Phase 0)

See [research.md](research.md) for six resolved decisions:
- **R1**: Deferred contribution GEMM — skip NFS×NFS per-block, single final GEMM to contribution buffer. 9.1d data proves pre-allocation alone is insufficient.
- **R2**: Buffer sized to `max_front × max_front` from symbolic analysis
- **R3**: Delayed columns — NFS×NFS from deferred GEMM, small delayed portion copied from workspace
- **R4**: ContributionBlock stays owned (`Mat<f64>`) for parallel path; swap protocol for zero-allocation steady state
- **R5**: Deferred GEMM executes while L21/D are in workspace, before workspace reuse
- **R6**: Custom allocators (AppendAlloc, BuddyAllocator) not needed

### Data Model (Phase 1)

See [data-model.md](data-model.md) for:
- `FactorizationWorkspace` extended with `contrib_buffer: Mat<f64>`
- Buffer lifecycle: final GEMM writes to buffer → buffer moved into ContributionBlock → parent's extend-add recycles buffer back to workspace
- Per-block trailing update decomposition (regions 1, 2, 3)

### API Contracts (Phase 1)

See [contracts/internal-api.md](contracts/internal-api.md) for modified internal contracts:
- C1-C2: Workspace construction/capacity (add contrib_buffer)
- C3: `update_trailing` (restricted to FS region + cross-terms)
- C4: New deferred contribution GEMM function
- C5: `extract_contribution` (index-only + small delayed copy)
- C6-C7: `extend_add`/`extend_add_mapped` (return recycled buffer)
- C8: `factor_single_supernode` (thread buffer, invoke deferred GEMM)

### Implementation Approach

The implementation follows a layered approach:

**Layer 1 — Workspace extension**: Add `contrib_buffer` to `FactorizationWorkspace`, allocate in `new()`, grow in `ensure_capacity()`. Sized to `max_front × max_front`.

**Layer 2 — Trailing update restriction**: Modify `update_trailing` to perform two GEMMs instead of one: (1) a lower-triangular GEMM on `A[ts..p, ts..p]` (FS×FS, region 1) and (2) a rectangular GEMM on `A[p..m, ts..p]` (NFS×FS cross-term, region 2). The NFS×NFS region `A[p..m, p..m]` is no longer touched during the blocking loop.

**Layer 3 — Deferred contribution GEMM**: Add a new function `compute_contribution_gemm` called from `factor_single_supernode` in `numeric.rs`, after `aptp_factor_in_place` returns (not inside `factor_inner`). This function copies the assembled NFS×NFS values from the workspace into `contrib_buffer`, then applies the rank-`ne` symmetric update in-place (`contrib_buffer -= L21_NFS * D * L21_NFS^T`). The function is self-contained: it reads all inputs from the workspace and produces the final Schur complement in `contrib_buffer`.

**Layer 4 — Index-only extraction**: Rewrite `extract_contribution` to build only `row_indices` and `num_delayed`. For the NFS×NFS portion, the data is already in `contrib_buffer` (from Layer 3). For the small delayed portion (rows `ne..k`), copy from the workspace. Move the buffer into `ContributionBlock`.

**Layer 5 — Extend-add buffer return + orchestration**: Modify `extend_add`/`extend_add_mapped` to take ownership and return the consumed buffer. In the factorization loop, swap the returned buffer back into the workspace. Thread-local buffers for the parallel path.

**Layer 6 — Validation**: Run full test suite + baseline comparison. Verify ExtractContr near zero in diagnostic output. Verify new ContribGEMM sub-phase appears.

### Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| Splitting `update_trailing` into FS+cross-term GEMMs introduces subtle correctness bugs | Medium | High | Regression tests on all 65 matrices + hand-constructed matrices with known Schur complements. Small-matrix exact-value tests where Schur complement can be computed analytically. |
| The deferred GEMM's FP ordering produces backward errors above 5e-11 on some matrices | Low | Medium | Mathematically identical computation; only rounding order changes. Monitor backward error deltas in baseline comparison. |
| `two_level_factor` and `tpp_factor_as_primary` paths need different deferred GEMM integration | Medium | Medium | All three paths produce the same workspace layout (L21 + D). The deferred GEMM is applied uniformly after whichever path runs. Verify with matrices that exercise each path. |
| Parallel path buffer management across thread boundaries | Low | Medium | Rust ownership prevents data races at compile time. Thread-local workspace pattern already proven in Phase 8.2. |
| Per-block GEMM splitting hurts performance on matrices with large `p` (many fully-summed columns) | Low | Low | For large `p`, the FS×FS GEMM is substantial but the NFS region is small — the deferred GEMM is cheap. The split is only costly when both FS and NFS are large, which is uncommon. |
