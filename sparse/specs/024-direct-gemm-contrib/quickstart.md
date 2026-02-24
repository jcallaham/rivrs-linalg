# Quickstart: Direct GEMM into Contribution Buffer

**Feature**: 024-direct-gemm-contrib
**Date**: 2026-02-24 (revised)

## What This Feature Does

Restructures the BLAS-3 blocking loop in the multifrontal factorization so that the Schur complement (NFS×NFS region) is computed in a single GEMM after the loop, writing directly into a pre-allocated contribution buffer. This eliminates the O(n²) post-factorization copy that accounts for 40.1% of factor time on c-71, aligning with SPRAL's zero-copy architecture.

## Files That Change

| File | Change Description |
|------|-------------------|
| `src/aptp/factor.rs` | `update_trailing` restricted to FS region + cross-terms (skips NFS×NFS). New deferred contribution GEMM function after blocking loop. |
| `src/aptp/numeric.rs` | `FactorizationWorkspace` gains `contrib_buffer`. `extract_contribution` becomes index-only (+ small delayed-column copy). `extend_add`/`extend_add_mapped` return recycled buffer. `factor_single_supernode` orchestrates deferred GEMM + buffer swap. |
| `examples/profile_matrix.rs` | New `ContribGEMM` sub-phase in diagnostic output. |

## Files That Don't Change

| File | Why Not |
|------|---------|
| `src/aptp/solver.rs` | Public `SparseLDLT` API is unchanged |
| `src/aptp/symbolic.rs` | Symbolic analysis is unchanged |
| `src/aptp/solve.rs` | Solve path is unchanged |

## How to Verify

```bash
# Unit tests (should all pass — backward error within tolerance)
cargo test

# Full SuiteSparse suite (65 matrices, backward error < 5e-11)
cargo test -- --ignored --test-threads=1

# Performance measurement (requires diagnostic feature)
# Look for: ExtractContr near zero, new ContribGEMM sub-phase
cargo run --example profile_matrix --features diagnostic --release -- c-71

# Baseline comparison
cargo run --example baseline_collection --features diagnostic --release -- --ci-only
```

## Key Design Decision

The BLAS-3 blocking loop's per-block trailing update is **restructured**: each block updates only the FS×FS region and NFS×FS cross-terms, skipping the NFS×NFS region entirely. After the loop, a single GEMM computes the NFS×NFS Schur complement directly into the contribution buffer. This is the same approach SPRAL uses (`factor.hxx:92-103`), where `host_gemm` writes to `node.contrib` rather than the frontal workspace.

Phase 9.1d proved that pre-allocating the buffer alone (without changing the GEMM target) nets zero improvement — the copy itself is the bottleneck, not the allocation. See `research.md` R1 for the full analysis.
