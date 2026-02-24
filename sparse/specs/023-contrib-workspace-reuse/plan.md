# Implementation Plan: Contribution Workspace Reuse

**Branch**: `023-contrib-workspace-reuse` | **Date**: 2026-02-23 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/023-contrib-workspace-reuse/spec.md`

## Summary

Eliminate per-supernode heap allocation for contribution blocks in the multifrontal factorization loop. Phase 9.1c profiling showed 73% of c-71 factor time is spent in `extract_contribution` (40.1%) and `extend_add` (33.3%). R8 diagnostic testing confirmed the dominant cost is first-touch page faults during faer's `Mat::zeros` zero-write pass on freshly mapped memory (~15-25% of factor time), compounded by the column-by-column data copy (~10-15%). Buffer reuse provides 1.7-2.3x empirical speedup on extraction; direct extend-add eliminates the copy entirely.

Two complementary approaches:
1. **Contribution buffer pool** (Approach 1): Free-list of reusable `Mat<f64>` on `FactorizationWorkspace`. Eliminates per-supernode allocation for all supernodes.
2. **Direct extend-add from frontal workspace** (Approach 2): For the last child of each parent in sequential DFS postorder, skip the contribution copy entirely and read the Schur complement in-place from the dual frontal buffer.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22, rayon 1.x, serde/serde_json (diagnostic export)
**Storage**: N/A (in-memory buffers; no persistent storage changes)
**Testing**: `cargo test` (358 unit tests), `cargo test -- --ignored --test-threads=1` (65 SuiteSparse matrices)
**Target Platform**: Linux (x86_64), development in Docker devcontainer
**Project Type**: Single Rust library crate
**Performance Goals**: c-71 factor time ≤ 1.5× SPRAL (from 2.48×); sys time < 10% (from 32%); dTLB misses < 64B (from 644B)
**Constraints**: 2× frontal workspace memory overhead (two max_front × max_front buffers); no correctness regressions (backward error < 5e-11 for all 65 matrices); no unsafe code

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | FR-009 requires all 65 SuiteSparse matrices pass with backward error < 5e-11. Performance optimization must not compromise numerical accuracy. Reconstruction tests and backward error tests remain the primary validation. |
| II. Clean Room | PASS | SPRAL source (BSD-3) freely consulted for workspace management patterns. Academic references (Liu 1992, Duff 2020) cited. No HSL or GPL code consulted. |
| III. TDD | PASS | Regression tests before refactoring: encode current correct behavior before modifying factorization loop. New unit tests for contribution pool, direct extend-add, dual-buffer swap. Property tests for buffer lifecycle invariants. |
| IV. Documentation | PASS | Academic references cited in spec and research.md. Rustdoc for new public API surface. SPRAL routine cross-references. |
| V. Numerical Stability | PASS | Optimization changes memory management, not numerical computation. Same APTP kernel, same pivot decisions, same factored values. Bit-exact results expected. |
| VI. Structured Development | PASS | Phase 9.1d follows Phase 9.1c in the ssids-plan.md roadmap. Exit criteria are measurable (SC-001 through SC-006). |
| VII. Code Quality | PASS | Extends existing FactorizationWorkspace pattern. No unsafe code. Result types for errors. Feature-gated diagnostics. |

**Post-Phase 1 re-check**: No violations introduced. The dual-buffer approach adds memory (2× frontal workspace) but this is documented as an accepted tradeoff. No new dependencies. No new unsafe code.

## Project Structure

### Documentation (this feature)

```text
specs/023-contrib-workspace-reuse/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Phase 0: SPRAL patterns, Liu/Duff references, design decisions
├── data-model.md        # Phase 1: Modified entities, state transitions
├── quickstart.md        # Phase 1: Build/test/validate instructions
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code

```text
src/
└── aptp/
    └── numeric.rs       # Primary file: FactorizationWorkspace, factor_tree_levelset,
                         #   factor_single_supernode, extract_contribution, extend_add,
                         #   ContributionBlock, new DFS sequential path

examples/
├── profile_matrix.rs    # Updated sub-phase timing for new code paths
└── baseline_collection.rs  # Updated baseline format
```

**Structure Decision**: All changes are internal to `src/aptp/numeric.rs`. No new modules or files needed. The optimization modifies the factorization loop internals and workspace types within the existing module.

## Design

### Architecture Overview

The optimization has three layers, each independently valuable:

**Layer 1: Contribution Buffer Pool (Approach 1)**
- Add `contribution_pool: Vec<Mat<f64>>` to `FactorizationWorkspace`
- `extract_contribution` takes a buffer from pool instead of `Mat::zeros`
- After `extend_add` consumes a `ContributionBlock`, its `Mat<f64>` is returned to pool
- Eliminates first-touch page fault overhead for ALL supernodes (sequential and parallel paths) by keeping buffer pages physically resident
- Works with the existing wave-based traversal — no control flow changes

**Layer 2: DFS Sequential Path + Dual Buffers (Approach 2)**
- Add `frontal_data_alt: Mat<f64>` to `FactorizationWorkspace`
- New `factor_tree_dfs` function: iterative DFS postorder traversal for sequential path
- Parent's frontal matrix assembled in the alternate buffer while children factor in the primary buffer
- Last child: direct extend-add from primary buffer into alternate buffer (zero-copy)
- Earlier children: use contribution pool (Layer 1)
- Buffer roles swap at each parent/child transition

**Layer 3: Parallel Path Pool Extension**
- Per-thread contribution pool via existing `Cell<FactorizationWorkspace>` pattern
- No direct extend-add (contributions cross thread boundaries)
- Wave-based traversal unchanged

### Detailed Design: Layer 1 — Contribution Buffer Pool

**Pool operations** on `FactorizationWorkspace`:

- `take_contribution_buffer(size: usize) -> Mat<f64>`: Takes a buffer from pool that is ≥ size×size. If none available, allocates fresh. Zeroes the used triangle.
- `return_contribution_buffer(buf: Mat<f64>)`: Returns buffer to pool for reuse.

**Integration with existing code**:

`extract_contribution` currently:
```
let mut data = Mat::zeros(size, size);  // ← BOTTLENECK: allocation + page faults
for j in 0..size { /* column-by-column copy */ }
ContributionBlock { data, row_indices, num_delayed }
```

With pool:
```
let mut data = workspace.take_contribution_buffer(size);  // ← reuse existing buffer
// zero only the triangle we'll write (fast: pages already resident)
for j in 0..size { data.col_as_slice_mut(j)[j..size].fill(0.0); }
for j in 0..size { /* column-by-column copy */ }
ContributionBlock { data, row_indices, num_delayed }
```

**Consumption and return**: After `extend_add` consumes a `ContributionBlock`:
```
// In factor_single_supernode, after extend_add:
if let Some(cb) = opt_cb {
    extend_add(&mut frontal, &cb, ...);
    workspace.return_contribution_buffer(cb.data);  // recycle
}
```

This requires changing `ContributionBlock` consumption to return the `Mat<f64>` instead of dropping it:
- `ContributionBlock::into_parts(self) -> (Mat<f64>, Vec<usize>, usize)` method — returns all three fields, allowing the caller to recycle the `Mat<f64>` via `return_contribution_buffer` while retaining access to `row_indices` and `num_delayed` if needed

**Pool sizing**: Starts empty; grows to accommodate the working set. For sequential processing, the maximum simultaneous live contributions equals the maximum number of children of any supernode (typically < 10). Pool buffers are sized to max_front but only the needed subregion is zeroed/written.

### Detailed Design: Layer 2 — DFS Sequential Path + Dual Buffers

**New traversal function**: `factor_tree_dfs` replaces `factor_tree_levelset` for the sequential path (`Par::Seq`). The wave-based `factor_tree_levelset` is retained for the parallel path.

**DFS postorder algorithm** (iterative, using explicit stack):

```
stack = [(root, ENTER)]
while stack not empty:
  (s, state) = stack.pop()
  match state:
    ENTER:
      stack.push((s, PROCESS))
      for each child c of s (in reverse order):
        stack.push((c, ENTER))
    PROCESS:
      // All children of s have been processed
      // Collect child contributions and factor s
```

**Dual-buffer integration**:

Two `Mat<f64>` buffers: `buf[0]` and `buf[1]`. A `current_buf: usize` index tracks which buffer is active (being factored into).

When processing supernode s:
1. s factors in `buf[current_buf]`
2. If s has a parent p:
   - If s is the LAST child of p: direct extend-add from `buf[current_buf]` trailing submatrix into `buf[1 - current_buf]` (parent's accumulating buffer)
   - If s is not the last child: extract contribution into pool buffer, extend-add pool buffer into `buf[1 - current_buf]`, return pool buffer
3. When all children of p are done: scatter p's original entries into `buf[1 - current_buf]`, then `current_buf = 1 - current_buf` (swap), factor p in new `buf[current_buf]`

**Parent pre-allocation timing**:
- Pre-allocate parent's frontal region in `buf[1 - current_buf]` when the FIRST child of that parent begins processing
- Use estimated front size (without delays) for initial allocation
- After all children are factored and delayed columns are known, verify size. If delays caused growth beyond estimate, grow the buffer (rare with MatchOrderMetis)
- Scatter of original entries happens AFTER all children, when delays are known and the front size is exact

**Last-child detection**:
- During DFS, maintain `remaining_children[s]` counter (same as current wave-based approach)
- When `remaining_children[parent] == 0` after processing child, the just-processed child was the last
- Alternatively, compare against `children[parent].last()` to identify the designated "last child"

**Direct extend-add function**:

New function `extend_add_from_frontal`:
```
fn extend_add_from_frontal(
    parent: &mut FrontalMatrix,     // parent's frontal in alt buffer
    child_frontal_data: &Mat<f64>,  // child's factored frontal in primary buffer
    child_m: usize,                 // child's front size
    child_ne: usize,                // child's eliminated count
    child_k: usize,                 // child's fully-summed count
    child_row_indices: &[usize],    // child's frontal row indices
    child_perm: &[usize],          // child's APTP permutation
    parent_g2l: &[usize],          // parent's global-to-local
)
```

This function reads the trailing `(m - ne) × (m - ne)` submatrix of `child_frontal_data` and maps it into `parent` using the same index arithmetic as `extract_contribution` + `extend_add` combined — but without the intermediate copy.

**Row index computation for direct extend-add**:
1. Delayed columns: positions `ne..k` in the child's frontal, mapped through `child_perm` and `child_row_indices`
2. Pattern rows: positions `k..m` in the child's frontal, directly from `child_row_indices`
3. Map each global index to parent local index via `parent_g2l`
4. Add entries from child's lower triangle at `child_frontal_data[ne+i, ne+j]` to `parent[parent_li, parent_lj]`

### Detailed Design: Layer 3 — Parallel Path Pool Extension

**Changes to parallel path**:
- `FactorizationWorkspace` already includes contribution pool (Layer 1)
- Thread-local `Cell<FactorizationWorkspace>` pattern unchanged
- Each rayon worker's pool starts empty, warms up independently
- `extract_contribution` uses `workspace.take_contribution_buffer()` instead of `Mat::zeros`
- Consumed contributions return buffers to the thread-local pool

**No direct extend-add in parallel path**: Contributions must cross thread boundaries as owned `ContributionBlock` values. The `contributions[s]` vector stores owned values accessed by different threads across waves.

**No traversal changes**: Wave-based `factor_tree_levelset` retained for parallel path. Only the allocation mechanism changes (pool vs. fresh allocation).

### Testing Strategy

**Regression tests (before any changes)**:
- Collect baseline results for all 65 SuiteSparse matrices: backward error, factor time, sub-phase timing
- Snapshot hand-constructed matrix factorizations with exact expected values

**Unit tests for new components**:
- `test_contribution_pool_take_return`: Pool take/return lifecycle, buffer reuse, correct sizing
- `test_contribution_pool_empty_allocates`: Empty pool allocates fresh buffer
- `test_contribution_pool_reuses_buffer`: Returned buffer is reused on next take
- `test_extend_add_from_frontal_no_delays`: Direct extend-add matches extract+extend_add for zero-delay case
- `test_extend_add_from_frontal_with_delays`: Direct extend-add matches extract+extend_add with delayed columns
- `test_dual_buffer_single_child`: Single-child chain uses direct extend-add throughout
- `test_dual_buffer_multi_child`: Multi-child parent: pool for earlier children, direct for last
- `test_dfs_vs_levelset_results`: DFS and wave-based traversals produce identical FrontFactors

**Integration tests**:
- All 65 SuiteSparse matrices pass with backward error < 5e-11
- Parallel and sequential paths produce identical backward errors
- No matrix regresses in factor time by > 5%

**Property-based tests**:
- Buffer pool invariant: pool.len() + live_contributions == total_allocated (no leaks)
- Dual-buffer invariant: at any point, exactly one buffer is "active" and one is "accumulating"
- DFS postorder: all children processed before parent; all supernodes processed exactly once

### Performance Validation

**Primary benchmark**: c-71 (6,350 supernodes, max front ~2,475)
- Factor time: target ≤ 1.5× SPRAL (from 2.48×)
- `perf stat`: sys time < 10%, dTLB-load-misses < 64B
- Sub-phase timing: extract_contrib near zero for direct extend-add path

**Secondary benchmarks**: Full SuiteSparse suite via `baseline_collection`
- Compare against Phase 9.1c baseline
- Verify no regressions > 5%
- Verify simplicial matrices (dixmaanl, bloweybq, etc.) benefit from pool

## Algorithm References

| Reference | Sections Used | File |
|-----------|---------------|------|
| Liu (1992) | §6 Algorithm 6.1 (supernodal multifrontal, stack-based contribution management), §7.4 (stack storage optimization, child ordering for minimum workspace) | `references/ssids/liu1992.md` |
| Duff, Hogg & Lopez (2020) | §5 (two-tier allocators: stack for factors, buddy for contributions; OS page cache exhaustion) | `references/ssids/duff2020.md` |
| Hogg, Duff & Scott (2016) | §2.3 (factorize phase with contribution on temporary stack), §3.1 (assembly mapping arrays) | `references/ssids/hogg2016.md` |
| Duff & Reid (1983/1984) | Foundational multifrontal method, assembly tree traversal | `references/ssids/duff1984.md` |

### SPRAL Source (BSD-3)

| File | Lines | What We Consult |
|------|-------|----------------|
| `BuddyAllocator.hxx` | 18-148 | Variable-size buddy allocator for contribution blocks |
| `NumericNode.hxx` | 37-51 | `alloc_contrib()` / `free_contrib()` lifecycle |
| `factor_cpu.cxx` | 73-168 | `factor_subtree` with workspace reuse |
| `assemble.hxx` | 143-443 | `assemble_pre()` / `assemble_post()` two-phase assembly |
| `NumericSubtree.hxx` | 75-81, 267-340 | Per-thread workspace init, assembly path |
