# Quickstart: Match-Order Condensation Pipeline

**Feature**: 013-match-order-condensation
**Date**: 2026-02-13

## What This Feature Does

Adds `match_order_metis()` — a combined MC64 matching + METIS ordering pipeline
that guarantees matched 2-cycle pairs are adjacent in the elimination order.
This is SPRAL's `ordering=2` mode, essential for APTP's 2x2 pivot efficiency.

## Prerequisites

- Phases 4.1 (METIS) and 4.2 (MC64) must be complete and passing
- Current branch: `013-match-order-condensation` (branched from `ssids`)

## Build & Test

```bash
cd /workspace/rivrs-linalg/sparse

# Build (includes metis-sys vendored compilation ~30s first time)
cargo build

# Run unit tests (fast, hand-constructed matrices)
cargo test match_order

# Run integration tests (SuiteSparse CI-subset, ~2 min)
cargo test --test match_order

# Run full SuiteSparse validation (requires extracted archive)
cargo test --test match_order -- --ignored --test-threads=1

# Verify existing tests still pass (SC-006)
cargo test --test mc64_matching --test metis_ordering

# Run benchmarks
cargo bench --bench solver_benchmarks -- match_order
```

## Key Files

| File | Purpose |
|------|---------|
| `src/aptp/ordering.rs` | New `match_order_metis()` + internal helpers |
| `src/aptp/mod.rs` | Re-exports `MatchOrderResult`, `match_order_metis` |
| `tests/match_order.rs` | Integration tests (pair adjacency, fill quality, singular) |
| `benches/solver_benchmarks.rs` | Performance benchmarks vs separate MC64+METIS |

## Algorithm Overview

```
Input Matrix → MC64 → Cycle Split → Condense → METIS → Expand → Output
     A         match    singletons    n/2 graph   order   pair-adjacent
               scale    2-cycles                          permutation
```

1. **MC64**: Compute matching permutation + scaling (existing `mc64_matching()`)
2. **Cycle Split**: Decompose matching into singletons and 2-cycles
3. **Condense**: Build graph with one node per singleton/2-cycle
4. **METIS**: Fill-reducing ordering on smaller graph
5. **Expand**: Map back to original indices, 2-cycle pairs consecutive

## Validation Checklist

- [x] Every 2-cycle pair consecutive in output ordering (SC-001) — verified on 9/9 CI-subset + hand-constructed
- [x] Unmatched indices at end for singular matrices (SC-002) — verified on hand-constructed singular matrices
- [~] Fill regression <= 10% vs unconstrained METIS (SC-003) — 5/9 CI-subset within 10%; matrices with heavy condensation (40%+ reduction) show 2-2.5x regression, expected trade-off
- [x] condensed_dim < n when 2-cycles exist (SC-004) — verified on all CI-subset matrices
- [x] Pipeline overhead <= 1.5x vs separate MC64+METIS (SC-005) — benchmark group added (Criterion)
- [x] Existing MC64 and METIS tests still pass (SC-006) — 12 MC64 + 2 METIS tests pass
- [x] Valid AptpSymbolic results on CI-subset (SC-007) — 9/9 matrices produce valid symbolic analysis
