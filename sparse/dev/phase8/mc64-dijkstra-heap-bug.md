# MC64 Dijkstra Heap Bug: Investigation and Fix

## Problem Statement

After implementing the MC64 matching/scaling pipeline (Phase 4.2), integration tests
on production-scale SuiteSparse matrices showed **dual feasibility violations** in
the Hungarian algorithm's scaling output. Specifically, the property
`|s_i * a_ij * s_j| <= 1` for all stored entries was violated on several matrices:

| Matrix | n | Matched | Max Violation |
|--------|---|---------|---------------|
| TSOPF_FS_b39_c7 | 28,216 | 28,216/28,216 (full) | 5.0 |
| d_pretok | 182,730 | 182,730/182,730 (full) | 0.98 |
| turon_m | 189,924 | 189,924/189,924 (full) | 0.15 |
| astro-ph | 16,706 | 15,678/16,706 (singular) | 73.0 |

An `enforce_scaling_bound` post-processing function was introduced as a workaround,
iteratively capping scaling factors. While effective, SPRAL's reference implementation
has no such workaround, indicating a bug in our Dijkstra implementation.

## Investigation Methodology

### 1. SPRAL as Ground Truth

Compiled SPRAL's `hungarian_scale_sym` and `hungarian_match` directly from the BSD-3
source at `/opt/references/spral/src/scaling.f90` using gfortran. Created a suite of
Fortran test drivers:

- `spral_scaling_test.f90` — full scaling pipeline, checks `|s*a*s| <= 1 + 5e-14`
- `spral_scaling_dump.f90` — dumps scaling factors and matching
- `spral_dump_duals.f90` — dumps raw dual variables (u, v, cmax, matching)
- `spral_match_costgraph.f90` — runs SPRAL's `hungarian_match` on an externally-provided cost graph
- `spral_greedy_compare.f90` — runs only SPRAL's greedy initialization phase

### 2. Modular Isolation

Created Rust example binaries to export intermediate state for cross-comparison:

- `export_cost_graph` — exports our internal cost graph in SPRAL-compatible format
- `export_for_spral` — exports matrix as lower-triangle CSC for SPRAL consumption
- `compare_scaling` — side-by-side comparison of our vs SPRAL scaling factors
- `check_raw_violations` — surveys all matrices for violations without `enforce_scaling_bound`
- `dump_scaling` — dumps our MC64 output (scaling, matching, violation stats)
- `check_duplicates` — verifies cost graph has no duplicate entries after dedup

### 3. Per-Augmentation Tracing

Added `MC64_TRACE_DUALS` environment variable instrumentation to check dual feasibility
after every single Dijkstra augmentation, pinpointing the exact augmentation where the
first violation appeared.

## Areas Confirmed Correct

### Cost Graph Construction (CONFIRMED CORRECT)

**Test**: Exported our cost graph for TSOPF and fed it to SPRAL's `hungarian_match`.

**Result**: SPRAL achieved dual feasibility of 8.88e-16 on our cost graph.

**Conclusion**: `build_cost_graph` (c[i,j] = col_max_log[j] - ln|a[i,j]|, with dedup
of mirrored entries from full symmetric CSC) produces the correct cost matrix.

### Greedy Initial Matching (CONFIRMED CORRECT)

**Test**: Compared our `greedy_initial_matching` output with SPRAL's
`hungarian_init_heurisitic` on the same cost graph.

**Results**:
- Both match exactly 21,172 / 28,216 columns on TSOPF
- Both leave exactly 7,044 columns for Dijkstra augmentation
- Greedy dual infeasibility: 0.0 in both implementations
- Row dual statistics match: positive=7,093, zero=21,123, min=0.0, max=6.149

**Note**: The two implementations match *different* columns (our first unmatched column
is 64, SPRAL's is 20). This is expected — the Phase 2 length-2 augmenting path scan
is path-dependent. Both are valid greedy matchings. Column order does not affect
correctness of subsequent Dijkstra augmentations.

### Post-Hoc Column Dual Computation (CONFIRMED CORRECT)

`v[j] = c[matched_row, j] - u[matched_row]` matches SPRAL line 1160-1167.
This is standard complementary slackness — correct by construction given correct
row duals u[i].

### Scaling Symmetrization (CONFIRMED CORRECT)

`scaling[i] = exp((u[i] + v[i] - col_max_log[i]) / 2)` matches SPRAL's formula
(scaling.f90 line 751). The Duff-Pralet log-space correction for singular matrices
also matches SPRAL lines 777-797.

### Augmentation Path Traceback (CONFIRMED CORRECT)

The path traceback via `pr[j]`/`out[j]` and the matching update loop match SPRAL
lines 1124-1137 exactly. `isp` stores the edge index (not row index), and `jperm`
is updated along the path.

### Dual Update Formula (CONFIRMED CORRECT)

`u[i] += d[i] - csp` for finalized rows matches SPRAL lines 1140-1143. The update
is applied to exactly the set of finalized rows (those extracted from Q1 and processed).

### Reduced Cost Computation (CONFIRMED CORRECT)

`vj = dq0 - cost[jperm[j]] + u[q0]` and `dnew = vj + cost[idx] - u[i]` match
SPRAL lines 1081 and 1086. The `jperm[j]` O(1) matched-edge cost lookup was
implemented correctly.

### Q1 Skip Guard (CONFIRMED CORRECT)

Our `label[i] == LABEL_Q1` check (skip rows already at optimal distance in current
batch) correctly implements SPRAL's `l(i) >= low` guard (line 1098). Combined with
the finalized check (`label[i] == LABEL_VISITED` / SPRAL `l(i) >= up` at line 1084),
the skipping logic is equivalent.

### Condensation Pipeline (CONFIRMED CORRECT)

The `extract_matched_subgraph` + re-matching pipeline for structurally singular
matrices matches SPRAL's approach. The re-matching now uses the same fixed Dijkstra.

## Root Cause: Lazy Heap Deletion

### The Bug

Our Dijkstra implementation used Rust's `std::collections::BinaryHeap` (a max-heap
used with `Reverse<>` for min-heap behavior). This heap does **not** support
deletion by key or decrease-key. When a row's distance decreased and it moved from
the heap to Q1, the old heap entry with the higher distance remained as a
**stale entry** (lazy deletion pattern).

SPRAL uses an **indexed binary min-heap** with O(log n) position-tracked operations:
- `l(i)` tracks each row's position in the heap (1-based, 0 = not in heap)
- `heap_delete(lpos, ...)` removes a row from the heap when moving to Q1
- `heap_update(i, ...)` performs in-place decrease-key

The critical code difference is at SPRAL line 1101-1102:
```fortran
lpos = l(i)
if (lpos .ne. 0) call heap_delete(lpos, qlen, m, Q, D, L)
low = low - 1
q(low) = i
l(i) = low
```

Our old code simply pushed to Q1 without removing from the heap:
```rust
q1.push(i);
label[i] = LABEL_Q1;
// Old heap entry remains — "lazy deletion"
```

### Why Lazy Deletion Fails

Although our code had filtering logic for stale entries (checking `label[i]` and
`d[i]` when popping), the interaction between stale entries and the Q1 extraction
loop caused rows to be processed in subtly wrong order, leading to **1,400 fewer
rows being finalized** than SPRAL (14,114 vs 15,514 on TSOPF). Rows that were
never finalized retained stale (too-high) dual values, causing post-hoc column
duals to violate `u[i] + v[j] <= c[i,j]` for non-matched edges.

The first violation appeared at augmentation 4 (root_col=68), only 4 augmentations
into the 7,044 needed. The specific violating edge was (69, 14127) where:
- u[69] = 0.0 (row 69 was never finalized — stale initial dual)
- v[14127] = 1.885 (derived from row 68's updated dual: u[68] = -1.885)
- SPRAL had u[69] = -4.175 (row 69 WAS finalized in SPRAL's run)

### Why Standard Lazy Deletion Arguments Don't Apply

In textbook Dijkstra with lazy deletion, correctness is maintained because each node
is processed at most once (via a "visited" check). The Hungarian algorithm's Dijkstra
variant has additional complexity:

1. **Batch extraction**: All rows with d[i] == dmin are extracted together into Q1.
   Stale entries at the same distance as legitimate entries can cause incorrect
   batch composition.

2. **Heap-to-Q1 migration**: When a row's distance drops to <= dmin mid-batch,
   it should be atomically moved from heap to Q1. Without explicit deletion,
   the heap's size (`qlen`) is wrong, affecting subsequent heap operations.

3. **Cross-augmentation state**: The `l[i]`/position array serves dual purposes —
   tracking heap membership AND Q1/finalized membership. Lazy deletion breaks
   this unified state tracking.

## The Fix

### Indexed Binary Min-Heap

Implemented SPRAL's exact heap as three inline functions operating on shared arrays:

```rust
fn heap_update_inline(idx, q, d, pos, qlen)  // Insert or decrease-key
fn heap_pop_inline(q, d, pos, qlen) -> usize  // Extract minimum
fn heap_delete_inline(pos0, q, d, pos, qlen)  // Delete at position
```

These operate on:
- `q[0..qlen]` — heap array (0-indexed, maps to SPRAL's `Q(1:qlen)`)
- `pos[i]` — 1-based position of row i in heap (0 = not in heap)
- `d[i]` — external distance array used for comparisons

The `l[i]` array now correctly encodes:
- `0` = not touched
- `1..qlen` = in heap (position in heap array)
- `low..up-1` = in Q1 (current minimum-distance batch)
- `up..n` = finalized (processed and dual-updated)

### Persistent DijkstraState

Added `DijkstraState` struct that persists across augmentations, matching SPRAL's
allocation pattern. Arrays `d`, `l`, `jperm`, `pr`, `out`, `q` are allocated once
and selectively reset after each augmentation (SPRAL lines 1144-1153), avoiding
O(n) re-initialization per augmentation.

### Removed enforce_scaling_bound

The `enforce_scaling_bound` workaround function was deleted along with its two
unit tests. With correct dual feasibility, the scaling bound is guaranteed by
the Hungarian algorithm's LP duality properties — no post-processing needed.

## Verification

| Test | Result |
|------|--------|
| TSOPF dual infeasibility | 8.88e-16 (was 3.66) |
| Rows with u < 0 (TSOPF) | 15,514 (matches SPRAL exactly) |
| Unit tests (`cargo test --lib -- matching`) | 32/32 pass |
| MC64 CI subset (`test_mc64_suitesparse_ci_subset`) | Pass |
| Full SuiteSparse MC64 (`test_mc64_suitesparse_full`) | 67/67 pass |
| Solver unit tests (`cargo test --test solve`) | 22/22 pass |
| All library tests (`cargo test --lib`) | 353/353 pass |

## Key Lesson

When porting algorithms that use indexed priority queues (common in graph algorithms,
network flow, shortest-path problems), **do not substitute lazy-deletion heaps**.
The indexed heap's position tracking serves as unified state that the algorithm
depends on for correctness, not just performance. Rust's `BinaryHeap` is suitable
for simple Dijkstra but not for variants where elements must be explicitly moved
between data structures mid-algorithm.

## Files Modified

| File | Change |
|------|--------|
| `src/aptp/matching.rs` | Replaced `BinaryHeap` Dijkstra with indexed heap port of SPRAL. Added `DijkstraState`, `heap_update_inline`, `heap_delete_inline`, `heap_pop_inline`. Removed `enforce_scaling_bound`, `check_dual_feasibility_raw`, `OrderedFloat`, trace instrumentation. |
