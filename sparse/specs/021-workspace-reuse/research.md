# Research: Workspace Reuse & Per-Supernode Allocation Optimization

**Feature**: 021-workspace-reuse
**Date**: 2026-02-22

## R1: Which allocations can be hoisted out of the per-supernode loop?

### Decision

Allocations fall into three tiers: **hoistable workspace** (reusable across supernodes), **per-supernode outputs** (must be freshly owned), and **per-block temporaries** (reusable within a single APTP kernel call).

### Hoistable Workspace (Tier 1 — per-supernode loop in numeric.rs)

These are allocated in `factor_single_supernode` and freed when the function returns. Dimensions are bounded by `max_front_size` from symbolic analysis.

| Buffer | Current Location | Size | Can Hoist? |
|--------|-----------------|------|-----------|
| `frontal.data` (Mat) | numeric.rs:757 | m x m | YES — use pre-allocated max_front x max_front buffer |
| `frontal.row_indices` (Vec) | numeric.rs:739 | m entries | YES — pre-allocate max_front capacity |
| `delayed_cols` (Vec) | numeric.rs:727 | <= k entries | YES — pre-allocate max_k capacity |

### Per-Supernode Outputs (CANNOT hoist — stored in AptpNumeric/ContributionBlock)

| Buffer | Location | Why Not Hoistable |
|--------|----------|-------------------|
| `l11` (Mat) | numeric.rs:1197 | Stored in FrontFactors, used in solve |
| `l21` (Mat) | numeric.rs:1269 | Stored in FrontFactors, used in solve |
| `d11` (MixedDiagonal) | numeric.rs:1235 | Stored in FrontFactors, used in solve |
| `local_perm` (Vec) | numeric.rs:1281 | Stored in FrontFactors |
| `col_indices` (Vec) | numeric.rs:1289 | Stored in FrontFactors |
| `row_indices` (Vec) | numeric.rs:1296 | Stored in FrontFactors |
| `contribution.data` (Mat) | numeric.rs:1330 | Propagated to parent via extend_add |
| `contribution.row_indices` (Vec) | numeric.rs:1340 | Propagated to parent |

### Per-Block Temporaries Inside APTP Kernel (Tier 2 — factor.rs)

These are allocated per block iteration inside `factor_inner`. Pre-allocating them at the start of each `aptp_factor_in_place` call (or hoisting into the workspace struct) eliminates per-block allocation overhead.

| Buffer | Location | Size | Notes |
|--------|----------|------|-------|
| `backup` (Mat) | factor.rs:1654 | (m-k) x ib | Per block. Could pre-allocate m x ib once |
| `l11_copy` (Mat) | factor.rs:1287 | bn x bn | Per apply_and_check. Pre-allocate ib x ib |
| `w` in update_trailing (Mat) | factor.rs:1397 | ts x ne | Per block. Pre-allocate m x ib |
| `l21_copy` (Mat) | factor.rs:1426 | ts x ne | Per block. Pre-allocate m x ib |
| `l_blk`, `w_blk` (Mat) | factor.rs:1520-1526 | nf x ne | Per cross-term. Pre-allocate ib x ib |
| `l_panel`, `w_panel` (Mat) | factor.rs:1543-1548 | ts x ne | Per cross-term. Pre-allocate m x ib |
| `panel_perm_buf` | factor.rs:1626 | ib | ALREADY reused across blocks |
| `row_perm_buf` | factor.rs:1627 | ib | ALREADY reused across blocks |
| `col_order_buf` | factor.rs:1628 | ib | ALREADY reused across blocks |

### APTP Kernel Per-Call Outputs (returned to caller)

| Buffer | Location | Why Not Hoistable |
|--------|----------|-------------------|
| `col_order` (Vec) | factor.rs:1617 | Returned as AptpFactorResult.perm |
| `d` (MixedDiagonal) | factor.rs:1618 | Returned as AptpFactorResult.d |
| `pivot_log` (Vec) | factor.rs:1620 | Returned as AptpFactorResult.pivot_log |

### Rationale

The key insight is the **output vs workspace** distinction. Roughly half the per-supernode allocations (frontal matrix, row indices, delayed columns, per-block backup/copies) are temporary workspace that can be pre-allocated and reused. The other half (L11, L21, D11, contribution blocks, APTP result vectors) are final outputs that must be freshly owned because they persist past the supernode's processing.

### Alternatives Considered

- **Arena allocation for outputs**: Allocate all FrontFactors from a single arena. This would require unsafe code and complex lifetime management. Deferred — the workspace reuse for temporaries should close most of the gap first.
- **Compact factor storage**: Store L factors in CSC format instead of dense Mat. This changes the solve path significantly and is out of scope.

---

## R2: Contribution block lifecycle and in-place feasibility

### Decision

In-place contribution block handling (US2) is **partially feasible**. The contribution block copy in `extract_contribution` (numeric.rs:1330-1335) can be eliminated by transferring ownership of a subregion, but cannot use a view into the frontal workspace because the workspace is reused for the next supernode.

### Findings

1. **Current flow**: `extract_contribution` copies the trailing (m-ne) x (m-ne) submatrix of the frontal matrix into a new owned `Mat<f64>`. The child's frontal workspace is then dropped.

2. **Why views are infeasible**: The child's frontal workspace goes out of scope at the end of `factor_single_supernode` (or would be overwritten with workspace reuse). The parent processes the contribution block later, possibly much later in the tree traversal. The contribution must survive independently.

3. **Alternative: owned transfer instead of copy**: Instead of zeroing a new Mat and copying element-by-element, the contribution block can potentially be extracted by taking ownership of the frontal workspace's trailing submatrix. This requires restructuring how the frontal data is stored — either as separate fully-summed and trailing regions, or by extracting a contiguous subblock.

4. **Extend-add access pattern**: `extend_add` reads contribution blocks sequentially via `&ContributionBlock` references. It processes children one at a time. No simultaneous access to multiple contribution blocks is needed.

5. **Level-set ordering**: Children complete before parents (postorder traversal). A child's contribution is extracted and stored in `contributions[c]`. The parent later takes it via `contributions[c].take()`.

### Alternatives Considered

- **Separate frontal layout** (fully-summed region + trailing region as independent Mats): Eliminates the copy but increases assembly complexity. Could be explored as a refinement.
- **Memory pool for contribution blocks**: Pre-allocate a pool of Mats sized to common contribution block dimensions. Returns to pool after extend-add instead of deallocating. Medium complexity.

---

## R3: SPRAL workspace patterns

### Decision

Adopt a simplified version of SPRAL's workspace strategy: pre-allocate workspace buffers per factorization (sequential) or per thread (parallel), sized to max front dimensions. Do not implement buddy allocator or block pool — the Rust allocator (jemalloc-style) with pre-sized buffers is sufficient for Phase 9.1b.

### SPRAL's Approach

1. **Stack allocator for factors**: Never deallocated during factorization. Uses `calloc()` for zero-on-allocation.
2. **Buddy allocator for contribution blocks and workspaces**: Variable-size allocation with merge-on-free. Thread-safe via OpenMP locks.
3. **Per-thread workspace**: Each thread gets fixed-size workspace (SSIDS_PAGE_SIZE ~1-2MB), reused across all nodes.
4. **Small-leaf subtree optimization**: L2-cache-resident subtrees with zero allocation overhead. Pre-calculated factor sizes.
5. **Lazy zeroing**: GEMM treats destination as zero via beta=0 rather than explicit memset.

### What to Adopt

- **Per-thread workspace** (high priority): Match SPRAL's per-thread workspace pattern. Pre-allocate one `FactorizationWorkspace` per thread, reuse across supernodes via Cell-based move (matching existing `g2l` pattern).
- **Lazy zeroing** (medium priority): For the frontal matrix workspace, zero only the lower triangle of the m x m subregion actually used, not the full max_front x max_front buffer.

### What to Defer

- **Buddy allocator**: Contribution blocks and factor outputs require fresh allocations (they're final outputs, not workspace). A buddy allocator would only help if we pooled contribution block storage, which adds complexity for modest gains.
- **Small-leaf subtrees**: This is essentially the simplicial fast path, which was deferred from scope.
- **Stack allocator for factors**: Rust's ownership model handles factor storage naturally. No need for a custom allocator.

### References

- Duff, Hogg, Lopez (2020), `/workspace/rivrs-linalg/references/ssids/duff2020.md` lines 246-289
- Liu (1992), `/workspace/rivrs-linalg/references/ssids/liu1992.md` lines 718-794
- SPRAL source: `BlockPool.hxx:16-82`, `BuddyAllocator.hxx:18-148`, `NumericSubtree.hxx:75-81`

---

## R4: Thread-local workspace for parallel factorization

### Decision

Thread-local `FactorizationWorkspace` via `Cell<FactorizationWorkspace>` in a `thread_local!` block, matching the existing `G2L_BUF` pattern in `factor_tree_levelset`.

### Rationale

- The `g2l` thread-local pattern (numeric.rs:937-963) is proven correct with nested parallelism (faer Par::rayon within rayon par_iter)
- Cell-based move semantics avoid RefCell borrow panics during work-stealing
- Each worker thread allocates its workspace lazily on first use, then reuses across all level-set waves
- Workspace sized to max front across all supernodes (conservative but simple; subtree-specific sizing would require symbolic analysis per subtree)

### Alternative Considered

- **Per-subtree workspace sizing**: Size each worker's workspace to the max front in its assigned subtree rather than the global max. This saves memory when subtrees have very different front sizes. Adds complexity to the scheduling logic. Could be a follow-up optimization if memory pressure is observed.

---

## R5: Zeroing strategy for reused workspace

### Decision

Use partial zeroing: zero only the lower triangle of the m x m subregion used for each supernode, where m is the supernode's front size. The rest of the pre-allocated buffer is untouched.

### Rationale

- The current code uses `Mat::zeros(m, m)` which zeros the full m x m region
- With workspace reuse, we have a max_front x max_front buffer. Zeroing the full buffer on every supernode would be O(max_front^2) even for small fronts
- Only the m x m subregion is used. Zeroing the lower triangle (m*(m+1)/2 entries) is sufficient because:
  - `scatter_original_entries_multi` writes specific entries
  - `extend_add` adds to specific entries
  - `aptp_factor_in_place` operates on the full m x m subregion
- After factorization, the used region must be left in a clean state for the next supernode. Two options:
  1. Zero before use (current approach, adapted to subregion)
  2. Clean up after use (restore to zero). More complex but avoids redundant zeroing if previous supernode had same or larger m
- **Choose option 1** (zero before use): simpler, predictable cost, matches SPRAL's `calloc()` pattern

### SPRAL Reference

Duff et al. (2020), `/workspace/rivrs-linalg/references/ssids/duff2020.md` line 256: "the memory associated with the fully summed columns should be zeroed before being assembled... the contribution block should be zeroed as well... it is done directly since the _gemm() operations can be told to treat memory as zero before the calculation."
