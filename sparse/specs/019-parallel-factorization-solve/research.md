# Research: Parallel Factorization & Solve (Phase 8.2)

**Date**: 2026-02-21
**Branch**: `019-parallel-factorization-solve`

## 1. faer's `Par` Enum — API and Composition

### Decision
Use faer's `Par` enum (`Par::Seq` / `Par::Rayon(NonZeroUsize)`) as the parallelism control surface. Add a `par: Par` field to `FactorOptions` and `SolverOptions`.

### Findings

**Par definition** (faer 0.22, `src/lib.rs`):
```rust
pub enum Par {
    Seq,
    #[cfg(feature = "rayon")]
    Rayon(NonZeroUsize),
}
```
- `Par::Rayon` is feature-gated in faer behind `rayon` (enabled by default)
- Helper: `Par::rayon(0)` = use all available threads; `Par::rayon(n)` = use n threads
- `par.degree()` returns thread count as `usize`

**How faer passes Par to BLAS operations**:
- `matmul`, `matmul_with_conj`, triangular solve, etc. all accept `par: Par` as final parameter
- Par::Seq runs a direct loop; Par::Rayon uses atomic job counter + `spindle::for_each` with rayon par_iter

**Nested parallelism via `join_raw`** (faer `src/utils/mod.rs`):
```rust
pub fn join_raw(op_a: impl FnOnce(Par), op_b: impl FnOnce(Par), par: Par)
```
- `Par::Seq` → sequential execution
- `Par::Rayon(n)` with n=1 → both ops get `Par::Seq`
- `Par::Rayon(n)` with n>1 → thread budget split between two operations
- Used internally by faer for recursive triangular operations

**Key implication**: faer's parallel BLAS shares Rayon's global thread pool. When rivrs-sparse uses Rayon for tree-level parallelism and faer uses Rayon for intra-node BLAS, they share the same pool — no oversubscription at the pool level, but tasks compete for worker slots.

### Alternatives Considered
- Custom `Parallelism` enum: rejected (violates transparent composition principle)
- Global thread pool config: rejected (per-call control is more flexible)
- Separate `ParallelOptions` struct: rejected (over-engineered; `Par` is sufficient)

## 2. Rayon Tree-Parallel Patterns

### Decision
Use `rayon::scope` for tree-level scheduling. Each supernode spawns its children as Rayon tasks within a scope; the scope barrier ensures all children complete before the parent processes.

### Findings

**rayon::join**:
- Takes exactly 2 closures; blocking but work-stealing
- Stack-allocated, ~100-500ns overhead
- Suitable for binary splits only

**rayon::scope** (recommended for N-ary trees):
```rust
rayon::scope(|s| {
    for &child in &children[node] {
        s.spawn(move |_| factor_node(child, ...));
    }
}); // Implicit barrier — all children complete before parent
```
- Handles arbitrary child counts naturally
- Heap-allocated per spawn (~1-5µs per task) — acceptable for tree-level granularity
- Maps directly to SPRAL's `#pragma omp taskgroup` semantics

**rayon::scope_fifo**: Same as scope but FIFO scheduling (breadth-first). Useful if cache locality matters for sequential-like traversal. Consider for solve phase.

### Composability with faer
- faer's BLAS operations use the same global Rayon pool
- No risk of deadlock from nested parallelism (work-stealing prevents it)
- Thread budget is implicitly shared — when tree tasks occupy pool threads, fewer are available for intra-node BLAS
- Mitigation: front size threshold (256) ensures only large fronts use parallel BLAS

### SPRAL Reference (OpenMP)
From `NumericSubtree.hxx`: SPRAL uses `#pragma omp taskgroup` with `depend` clauses for task ordering. The `depend(inout: this_lcol)` / `depend(in: parent_lcol)` pattern enforces child→parent ordering. Rayon's scope barrier achieves the same effect structurally.

### Alternatives Considered
- Recursive binary `rayon::join`: rejected (awkward for N children, requires binary splitting)
- Level-set scheduling with `par_iter`: rejected (explicit barriers between levels are less efficient than scope)
- Task graph with atomic dependency counters: rejected (unnecessary complexity; scope barrier is sufficient)

## 3. Current Sequential Factorization Structure

### Key Shared Mutable State (parallelization impact)

| State | Current Pattern | Parallel Strategy |
|-------|----------------|-------------------|
| `contributions: Vec<Option<ContributionBlock>>` | Index-based: `contributions[c].take()` consumed by parent | Return-value pattern: `factor_supernode()` returns `(FrontFactors, Option<ContributionBlock>, PerSupernodeStats)` as owned values. `factor_subtree()` collects children's returned contributions after scope barrier, then passes to parent assembly. No shared mutable Vec needed; no `unsafe` code. |
| `global_to_local: Vec<usize>` | Built/reset per supernode | Each parallel supernode needs its own; allocate per-task |
| `front_factors_vec: Vec<FrontFactors>` | Sequential push | Pre-allocate `Vec<Option<FrontFactors>>`; each supernode writes at its index |
| `per_sn_stats: Vec<PerSupernodeStats>` | Sequential push | Pre-allocate `Vec<Option<PerSupernodeStats>>`; each supernode writes at its index |
| `stats: FactorizationStats` | Incremental accumulation | Thread-local accumulators + final merge, or compute from per_sn_stats |
| Frontal matrix allocation | Per-supernode, dropped after extraction | No change needed (each task allocates its own) |

### Assembly Tree Properties
- **Postorder**: `parent_index > child_index` always (faer's supernodal layout)
- **Children map**: `build_children_map()` returns `Vec<Vec<usize>>` — each supernode knows its children
- **Contribution flow**: child creates `ContributionBlock`, parent consumes via `.take()`
- **No sibling dependencies**: children of the same parent are fully independent

### APTP Kernel BLAS-3 Sites (factor.rs)
Two `Par::Seq` call sites in `factor_inner()`:

1. **TRSM** (line 1278): `solve_unit_lower_triangular_in_place(l11, panel, Par::Seq)`
   - Dimensions: L11 is `b×b` (b = inner_block_size = 32), panel is `r×b`
   - Parallel benefit: moderate (panel rows can be large for big fronts)

2. **GEMM** (line 1415): `tri_matmul::matmul_with_conj(a22, ..., w, ..., l21, ..., Par::Seq)`
   - Dimensions: A22 is `s×s`, W is `s×b`, L21 is `s×b` where s = trailing_size
   - Parallel benefit: HIGH (this is the dominant cost — cubic in trailing_size)

### Solve Data Dependencies
- **Forward solve** (postorder): Each supernode's scatter (`rhs[row_indices[j]] -= ...`) writes to indices that overlap with descendants' col_indices → must respect postorder
- **Diagonal solve**: No inter-supernode dependencies (each supernode reads/writes only its own col_indices)
- **Backward solve** (reverse postorder): Each supernode's gather-update reads `rhs[row_indices[j]]` → must respect reverse postorder

## 4. Threading Par Through the Dense Kernel

### Decision
Add `par: Par` to `AptpOptions`. Thread it through `aptp_factor_in_place()` → `two_level_factor()` / `factor_inner()` → `apply_and_check()` + `update_trailing()`.

### Call Chain
```
AptpNumeric::factor(options with par)
  → aptp_factor_in_place(a, k, options)  // dispatcher
    → tpp_factor_as_primary(...)          // k < 32: no BLAS-3, Par ignored
    → two_level_factor(...)               // k > 256: outer loop
      → factor_inner(block, ...)          // inner BLAS-3 loop
        → apply_and_check(...)            // TRSM with par
        → update_trailing(...)            // GEMM with par
    → factor_inner(...)                   // 32 <= k <= 256
      → (same as above)
```

### Threshold Logic
- Front dimension < 256: pass `Par::Seq` to kernel regardless of options.par
- Front dimension >= 256: pass `options.par` to kernel
- The threshold is applied at the `AptpNumeric::factor()` level per-supernode, not inside the kernel

## 5. Solve Parallelism Architecture

### Decision
Three strategies, from simple to complex:

1. **Diagonal solve**: `par_iter` over all supernodes (embarrassingly parallel)
2. **Forward/backward solve**: Recursive `rayon::scope` matching factorization tree traversal
3. **Intra-node dense ops**: Pass `par` to per-supernode triangular solve and matmul

### Current Par::Seq Sites in solve.rs
- Line 147: `solve_unit_lower_triangular_in_place_with_conj(..., Par::Seq)` (forward)
- Line 170: `matmul_with_conj(..., Par::Seq)` (forward scatter)
- Line 250: `matmul_with_conj(..., Par::Seq)` (backward gather-update)
- Line 261: `solve_unit_upper_triangular_in_place_with_conj(..., Par::Seq)` (backward)

### RHS Vector Safety
Forward/backward solve read and write shared `rhs: &mut [f64]`:
- Forward: reads `rhs[col_indices]`, writes `rhs[col_indices]` and `rhs[row_indices]`
- Backward: reads `rhs[col_indices]` and `rhs[row_indices]`, writes `rhs[col_indices]`

For tree-level parallel solve: sibling supernodes at the same tree level may have overlapping `row_indices` (both share a parent's column space). This means:
- **Siblings in forward solve**: both scatter to potentially overlapping row_indices → data race
- **Solution**: use scope barrier per tree level (same as factorization — process siblings first, then parent)

### Workspace Allocation
Current: two `Vec<f64>` buffers reused across all supernodes (sequential).
Parallel: each parallel task needs its own workspace. Options:
- Per-task allocation (simple, small overhead for tree-level)
- Thread-local workspace pools

## 6. Rayon Dependency Integration

### Decision
Rayon is an unconditional dependency. faer already enables Rayon by default.

### Cargo.toml Changes
```toml
[dependencies]
rayon = "1"  # New: tree-level parallelism
```

No feature gating needed. faer's `rayon` feature is default-enabled, so `Par::Rayon` is always available.

### faer's Rayon Feature
faer's `Cargo.toml` has `rayon` in its default features. The `Par::Rayon` variant is compiled when faer has `rayon` enabled. Since rivrs-sparse uses `faer = "0.22"` (default features), `Par::Rayon` is always available.

## 7. Determinism Analysis

### Decision
Require bitwise determinism at fixed thread count (FR-005). Relax to backward-error equivalence across thread counts (FR-006).

### Why Bitwise Determinism at Fixed Thread Count is Achievable
1. Assembly tree schedule is fixed (postorder, not work-stealing-order-dependent)
2. `rayon::scope` barrier ensures all children complete before parent — same as sequential
3. faer's BLAS operations produce deterministic results for a given Par configuration
4. Pivot decisions depend only on matrix values, not scheduling order

### Why Bitwise Identity Across Thread Counts is NOT Guaranteed
1. faer's parallel GEMM uses atomic job counters for work distribution — same computation, but accumulation order may differ between Par::Seq and Par::Rayon(4)
2. Floating-point addition is not associative — different summation order = different ULP-level results
3. Different ULP-level results can cause different pivot decisions at threshold boundaries
