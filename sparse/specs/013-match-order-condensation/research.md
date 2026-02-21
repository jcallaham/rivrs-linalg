# Research: Match-Order Condensation Pipeline

**Feature**: 013-match-order-condensation
**Date**: 2026-02-13

## Overview

No NEEDS CLARIFICATION items existed in the Technical Context — all technology
choices are established from Phases 4.1 and 4.2. Research focused on algorithm
details from the SPRAL reference implementation and API design decisions.

## Research Findings

### R1: SPRAL's mo_split Algorithm (Cycle Splitting)

**Decision**: Follow SPRAL's `mo_split` cycle-walking algorithm exactly.

**Rationale**: The algorithm is simple, well-tested, and handles all edge cases
(singletons, 2-cycles, longer cycles, unmatched nodes). It uses a single `iwork`
array with a clear state convention:

| iwork value | Meaning |
|-------------|---------|
| -2 | Unmatched by MC64 |
| -1 | Singleton (self-matched) |
| 0 | Not yet processed |
| > 0 | Matched with node iwork[i] |

**Cycle-splitting logic for longer cycles**: Walk the cycle via `cperm`. At each
step, pair node `j` with `cperm(j)`, then advance to `cperm(cperm(j))`. If the
cycle has odd length, the last node becomes a singleton. This produces only
singletons and 2-cycles from any matching.

**Interface with Mc64Result**: SPRAL's `mo_split` receives the raw matching
arrays directly from `mo_match`. In our design, `Mc64Result.matching` is a
`Perm<usize>` where `build_singular_permutation()` assigns unmatched rows to
arbitrary free columns, making `fwd[i]==i` unreliable for distinguishing
singletons from unmatched indices. To resolve this, `Mc64Result` exposes an
`is_matched: Vec<bool>` field, and `split_matching_cycles()` takes
`(matching_fwd, is_matched, n)` — checking `is_matched[i]` before processing
each index.

**Alternatives considered**:
- Custom cycle decomposition optimizing for specific 2-cycle patterns → Rejected:
  SPRAL's approach is simpler and equivalent for our use case.
- Preserving longer cycles as higher-order blocks → Rejected: APTP only uses
  1x1 and 2x2 pivots, so longer cycles have no downstream utility.
- Reconstructing matched status from permutation structure alone → Rejected:
  `build_singular_permutation()` creates a valid permutation for unmatched indices
  by assigning them to free columns, so `fwd[i]==i` does not reliably indicate
  a singleton. The `is_matched` field is already computed internally in
  `mc64_matching()` and costs nothing to expose.

### R2: Condensed Graph Construction (Deduplication Strategy)

**Decision**: Use SPRAL's marker-array deduplication technique.

**Rationale**: When merging edges from both members of a 2-cycle pair into one
condensed column, duplicates arise naturally. SPRAL reuses `iwork` as a per-column
marker: `iwork[krow] = current_column_index` means "row krow already added to this
column." This avoids explicit sorting or hash-set operations.

**Key detail**: The condensed graph is built as a full symmetric graph first, then
filtered to lower-triangular before passing to METIS. SPRAL's METIS wrapper
(`metis5_wrapper.F90:half_to_full_drop_diag`) expands lower-triangular to full
symmetric internally. Our `extract_adjacency()` already handles this expansion
(builds full symmetric from upper-triangular input), so we can either:

(a) Build lower-triangular condensed graph → pass to modified METIS call, or
(b) Build full symmetric condensed graph → pass directly to existing `extract_adjacency()` pattern.

**Decision**: Option (b) — build full symmetric condensed graph and pass to an
internal METIS wrapper that mirrors `extract_adjacency()`. This reuses established
patterns and avoids introducing a second adjacency construction code path.

**Alternatives considered**:
- Hash-set deduplication → Rejected: More allocation, slower for sparse graphs.
- Sort-and-unique → Rejected: O(nnz log nnz) vs O(nnz) for marker approach.

### R3: METIS Input for Condensed Graph

**Decision**: Build CSR adjacency from the condensed graph and call
`METIS_NodeND` directly, reusing the same FFI pattern as `metis_ordering()`.

**Rationale**: The condensed graph is a valid symmetric sparse graph. METIS
doesn't care that nodes represent super-nodes — it only sees adjacency structure.
We need to:
1. Build condensed CSC (full symmetric, diagonal excluded)
2. Convert to CSR (METIS input format) — same as `extract_adjacency()`
3. Call `METIS_NodeND` with the condensed dimension
4. Map METIS output back through new_to_old

**Key constraint**: METIS uses `i32` indices. The condensed dimension is at most n
(typically ~n/2), so if the original dimension fits in i32, the condensed dimension
always fits.

### R4: Ordering Expansion Algorithm

**Decision**: Follow SPRAL's two-step expansion (lines 376-395).

**Rationale**: METIS returns `order[i] = position of condensed node i`. The
expansion:
1. Build inverse: `inv_order[position] = condensed_node`
2. Walk positions 1..ncomp: for each condensed node, emit its original index(es)
   into the final ordering. Matched pairs get consecutive positions.

This guarantees pair adjacency by construction — it's not a post-hoc property
but inherent to the expansion algorithm.

### R5: Where to Place New Code

**Decision**: Add `match_order_metis()` and supporting types to
`src/aptp/ordering.rs`.

**Rationale**: The condensed match-order pipeline is fundamentally an ordering
algorithm that composes MC64 (from `matching.rs`) with METIS (already in
`ordering.rs`). Placing it in `ordering.rs`:
- Follows the existing pattern (ordering algorithms live in `ordering.rs`)
- Allows reuse of `extract_adjacency()` and METIS FFI infrastructure
- Keeps the module focused on its responsibility (computing orderings)

The cycle decomposition and condensed graph construction are internal helper
functions, not public API. Only `match_order_metis()` and `MatchOrderResult`
are public.

**Alternatives considered**:
- New `match_order.rs` module → Rejected: Would duplicate METIS FFI setup and
  create an artificial split. The function is ~100-150 lines of new code, not
  enough to warrant a separate module.
- Add to `matching.rs` → Rejected: The output is an ordering, not a matching.
  `matching.rs` is about bipartite matching; ordering is the downstream consumer.

### R6: MatchOrderResult Design

**Decision**: Return a struct containing ordering (`Perm<usize>`), scaling
(`Vec<f64>`), match count, and diagnostics.

**Rationale**: Callers need both the ordering (for symbolic analysis) and the
scaling (for numeric factorization). Bundling them avoids the caller having to
run MC64 separately. The diagnostics (condensed dimension, cycle counts) enable
validation and performance analysis without separate profiling.

**Design**: `MatchOrderResult` wraps `Mc64Result` fields plus the ordering:
```
MatchOrderResult {
    ordering: Perm<usize>,      // Fill-reducing ordering with pair adjacency
    scaling: Vec<f64>,          // MC64 scaling factors (pass-through)
    matched: usize,             // Number of matched entries
    condensed_dim: usize,       // Dimension of condensed graph
    singletons: usize,          // Number of singleton nodes
    two_cycles: usize,          // Number of 2-cycle pairs
}
```

### R7: Test Strategy

**Decision**: Four test categories targeting distinct validation goals.

1. **Structural invariant tests** (hand-constructed matrices):
   - Pair adjacency: every 2-cycle pair occupies consecutive positions (SC-001)
   - Unmatched-at-end: unmatched indices have highest positions (SC-002)
   - Valid permutation: output is a permutation of {0..n-1} (FR-007)
   - Cycle decomposition correctness: known matchings decompose correctly (FR-001, FR-002)
   - Round-trip: condensed_dim < n when 2-cycles exist (SC-004)

2. **Fill quality comparison** (SuiteSparse CI-subset):
   - One-sided: `condensed_nnz <= unconstrained_nnz * 1.10` (SC-003)
   - Run both `match_order_metis()` and `metis_ordering()` on same matrices
   - Compare `AptpSymbolic::analyze()` predicted nnz(L)

3. **Integration tests** (SuiteSparse CI-subset + full):
   - Pipeline produces valid `AptpSymbolic` results (SC-007)
   - Existing MC64 and METIS tests still pass (SC-006)

4. **Performance benchmarks** (Criterion):
   - Time `match_order_metis()` vs `mc64_matching() + metis_ordering()` (SC-005)
   - Report condensed dimension / original dimension ratio

**Alternatives considered**:
- Property-based testing (proptest) for random matchings → Deferred: SPRAL's
  algorithm is deterministic and well-understood. Hand-constructed cases with
  known cycle structure are more informative.
