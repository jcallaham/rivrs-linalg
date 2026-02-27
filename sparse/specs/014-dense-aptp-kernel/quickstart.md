# Quickstart: Dense APTP Factorization Kernel

**Feature Branch**: `014-dense-aptp-kernel`
**Date**: 2026-02-14

## Prerequisites

- Rust 1.87+ (edition 2024)
- faer 0.22 (already in Cargo.toml)
- All Phase 2 types available (MixedDiagonal, PivotType, Block2x2, Inertia)
- Test infrastructure from Phases 0.4/1.1 (NumericalValidator, test generators)

## New Files

| File | Purpose |
|------|---------|
| `src/aptp/factor.rs` | Core APTP factorization kernel: types, in-place algorithm, convenience wrapper |

## Modified Files

| File | Change |
|------|--------|
| `src/aptp/mod.rs` | Add `pub mod factor;` and re-exports for new public types |

## Build & Test

```bash
cd /workspace/rivrs-linalg/sparse

# Build
cargo build

# Run unit tests
cargo test aptp::factor

# Run full test suite (includes hand-constructed and CI matrices)
cargo test

# Run ignored tests (full SuiteSparse collection)
cargo test -- --ignored --test-threads=1
```

## Key Design Decisions

1. **In-place core API**: `aptp_factor_in_place(MatMut, num_fully_summed, options)` — L stored in lower triangle of mutated input. Essential for Phase 6 integration.

2. **Convenience wrapper**: `aptp_factor(MatRef, options)` — copies input, calls in-place, extracts L. For testing and small matrices.

3. **Single-level (column-by-column)**: No blocking in Phase 5. Two-level blocking deferred to Phase 9.1.

4. **Separate singularity threshold**: `AptpOptions.small` (default 1e-20) detects zero pivots, distinct from stability threshold `threshold` (default 0.01).

5. **Phase 2 types for D storage**: `MixedDiagonal` stores D factor, `PivotType` classifies each column.

## Algorithm Reference

Duff, Hogg & Lopez (2020), "A New Sparse LDL^T Solver Using A Posteriori Threshold Pivoting", SIAM J. Sci. Comput. 42(4). Available at `/workspace/rivrs-linalg/references/ssids/duff2020.md`.

SPRAL CPU kernel (BSD-3): `/opt/references/spral/src/ssids/cpu/kernels/ldlt_app.cxx`

## Validation

Reconstruction error `||A - P^T L D L^T P|| / ||A|| < 10^-12` verified using `NumericalValidator` from test infrastructure.
