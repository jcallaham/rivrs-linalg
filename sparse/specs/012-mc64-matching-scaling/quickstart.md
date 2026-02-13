# Quickstart: MC64 Matching & Scaling

**Feature**: 012-mc64-matching-scaling
**Branch**: `012-mc64-matching-scaling`

## Prerequisites

- Rust 1.87+ (edition 2024)
- Working `cargo build` in `sparse/` directory
- Test data extracted (see `test-data/README.md`)

## Build & Test

```bash
cd sparse/

# Build
cargo build

# Run unit tests (fast, hand-constructed matrices)
cargo test matching

# Run CI integration tests (SuiteSparse CI subset, ~10 matrices)
cargo test mc64 -- --test-threads=1

# Run full SuiteSparse validation (67 matrices, slower)
cargo test mc64 -- --ignored --test-threads=1

# Lint and format
cargo fmt --check
cargo clippy
```

## Key Files

### Implementation
- `src/aptp/matching.rs` — MC64 matching algorithm (core implementation)
- `src/aptp/mod.rs` — Module exports (`mc64_matching`, `Mc64Result`, `Mc64Job`)

### Tests
- `src/aptp/matching.rs` (inline `#[cfg(test)]` module) — Unit tests for internal functions
- `tests/mc64_matching.rs` — Integration tests (hand-constructed + SuiteSparse)

### References
- `/workspace/rivrs-linalg/references/ssids/duff2001.md` — Core MC64 algorithm (Algorithm MPD)
- `/workspace/rivrs-linalg/references/ssids/duff2005.md` — MC64SYM symmetric adaptation
- `/opt/references/spral/src/scaling.f90` — SPRAL Hungarian matching implementation
- `/opt/references/spral/src/match_order.f90` — SPRAL combined matching+ordering (`mo_split()` for future condensation)
- `/opt/references/spral/tests/scaling.f90` — SPRAL matching test patterns

### Design Artifacts
- `specs/012-mc64-matching-scaling/spec.md` — Feature specification
- `specs/012-mc64-matching-scaling/plan.md` — Implementation plan
- `specs/012-mc64-matching-scaling/research.md` — Research decisions
- `specs/012-mc64-matching-scaling/contracts/api.md` — API contract

## Quick Validation

After implementation, verify scaling properties match SPRAL's criteria:

1. All entries of scaled matrix satisfy `|s_i * a_ij * s_j| <= 1.0`
2. Row maxima of scaled matrix satisfy `max_j |s_i * a_ij * s_j| >= 0.75`
3. Matched diagonal entries satisfy `|s_i * a_{i,σ(i)} * s_{σ(i)}| ≈ 1.0`
4. All scaling factors are positive and finite
5. Matching decomposes into singletons + 2-cycles only (no longer cycles)
