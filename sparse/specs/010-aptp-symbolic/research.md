# Research: APTP Symbolic Analysis

**Feature**: 010-aptp-symbolic
**Date**: 2026-02-10

## Research Questions & Findings

### RQ-1: faer Symbolic Cholesky API Surface

**Question**: What is the exact public API of `factorize_symbolic_cholesky` and `SymbolicCholesky<I>` that AptpSymbolic will compose with?

**Findings**:

`factorize_symbolic_cholesky` signature:
```rust
pub fn factorize_symbolic_cholesky<I: Index>(
    A: SymbolicSparseColMatRef<'_, I>,
    side: Side,
    ord: SymmetricOrdering<'_, I>,
    params: CholeskySymbolicParams<'_>,
) -> Result<SymbolicCholesky<I>, FaerError>
```

`SymmetricOrdering` variants:
- `Amd` (default) — approximate minimum degree
- `Identity` — no reordering
- `Custom(PermRef<'a, I>)` — user-supplied permutation

`SymbolicCholesky<I>` public methods:
- `nrows()`, `ncols()` — matrix dimensions
- `perm() -> Option<PermRef<'_, I>>` — fill-reducing permutation
- `len_val() -> usize` — predicted nonzero count (L + D values storage)
- `raw() -> &SymbolicCholeskyRaw<I>` — access inner simplicial or supernodal variant
- `factorize_numeric_ldlt_scratch<T>()` — workspace requirement for numeric factorization
- `factorize_numeric_ldlt<'out, T>()` — perform numeric factorization
- `factorize_numeric_intranode_lblt<'out, T>()` — BK within supernodes
- `solve_in_place_scratch<T>()` — workspace requirement for solve

`SymbolicCholeskyRaw<I>` is an enum:
- `Simplicial(SymbolicSimplicialCholesky<I>)` — for sparse matrices
- `Supernodal(SymbolicSupernodalCholesky<I>)` — for denser/larger matrices

**faer source locations** (faer 0.22.6, `/opt/references/faer-rs/faer/src/sparse/linalg/`):
- `factorize_symbolic_cholesky`: `cholesky.rs:4192`
- `prefactorize_symbolic_cholesky`: `cholesky.rs:584`
- `SymbolicCholesky<I>` impl block: `cholesky.rs:3556`
- `SymbolicSupernodalCholesky<I>` impl block: `cholesky.rs:1554`
- `SymmetricOrdering` enum: `cholesky.rs:510`
- `CholeskySymbolicParams`: `cholesky.rs:3529`
- `EliminationTreeRef`: `cholesky.rs:535`
- `postorder` (in qr module): `qr.rs:302`

**Decision**: Use `factorize_symbolic_cholesky` (high-level), not low-level primitives. Store result directly.
**Rationale**: Transparent composition principle. faer handles the simplicial-vs-supernodal decision internally. We add only APTP-specific metadata.

### RQ-2: Supernodal Structure Accessibility

**Question**: Can we access the assembly tree (parent pointers, postorder, descendant counts) from faer's supernodal structure?

**Findings**:

**Critical limitation**: The assembly tree fields in `SymbolicSupernodalCholesky` are `pub(crate)`:
- `supernode_postorder: Vec<I>` — NOT publicly accessible
- `supernode_postorder_inv: Vec<I>` — NOT publicly accessible
- `descendant_count: Vec<I>` — NOT publicly accessible

Similarly, `SymbolicSimplicialCholesky.etree: Vec<I>` is private (no `pub`).

**Public methods available on SymbolicSupernodalCholesky**:
- `n_supernodes()` — count of supernodes
- `supernode_begin() -> &[I]` — start column of each supernode
- `supernode_end() -> &[I]` — end column of each supernode
- `col_ptr_for_row_idx() -> &[I]`, `col_ptr_for_val() -> &[I]` — column pointers
- `row_idx() -> &[I]` — row indices
- `supernode(s: usize) -> SymbolicSupernodeRef` — per-supernode accessor
  - `start() -> usize` — starting column
  - `pattern() -> &[I]` — row indices for this supernode

**However**, faer DOES expose a standalone `prefactorize_symbolic_cholesky` function:
```rust
pub fn prefactorize_symbolic_cholesky<'out, I: Index>(
    etree: &'out mut [I::Signed],
    col_counts: &mut [I],
    A: SymbolicSparseColMatRef<'_, I>,
    stack: &mut MemStack,
) -> EliminationTreeRef<'out, I>
```

And `EliminationTreeRef` has `into_inner() -> &[I::Signed]` which gives the parent pointer array.

There is also a public `postorder` function in `faer::sparse::linalg::qr`:
```rust
pub fn postorder<I: Index>(post: &mut [I], etree: EliminationTreeRef<'_, I>, stack: &mut MemStack)
```

**Decision**: Two-path approach:
1. Primary: Use `factorize_symbolic_cholesky` (high-level) for the core symbolic result
2. Supplementary: Also call `prefactorize_symbolic_cholesky` to obtain the elimination tree parent pointer array, which is cheap (almost-linear time) and gives us direct access to the tree structure

**Rationale**: The etree is needed by Phase 6 (multifrontal) for assembly tree traversal and by the pivot buffer heuristic. Computing it separately adds negligible cost (~O(nnz) vs the O(nnz * alpha) of the full symbolic factorization) and avoids forking faer or relying on private internals.

**Alternatives considered**:
- Reconstruct assembly tree from supernode patterns alone → fragile, non-trivial to get right
- Fork faer to expose private fields → violates transparent composition principle
- Use only low-level primitives → rebuilds what `factorize_symbolic_cholesky` already does

### RQ-3: Simplicial vs Supernodal Selection

**Question**: Does `factorize_symbolic_cholesky` always produce supernodal output, or can it produce simplicial? How do we handle both?

**Findings**:

`SymbolicCholeskyRaw<I>` is an enum — the result is **one or the other**, never both. The decision is automatic based on `supernodal_flop_ratio_threshold` in `CholeskySymbolicParams`:
- If the estimated BLAS speedup ratio exceeds the threshold → supernodal
- Otherwise → simplicial

For small hand-constructed matrices (n ≤ ~30), faer may choose simplicial. For SuiteSparse matrices, supernodal is more likely.

`CholeskySymbolicParams` fields:
```rust
pub struct CholeskySymbolicParams<'a> {
    pub amd_params: amd::Control,                         // dense threshold, aggressive absorption
    pub supernodal_flop_ratio_threshold: SupernodalThreshold, // BLAS speedup threshold
    pub supernodal_params: SymbolicSupernodalParams<'a>,  // relaxation thresholds
}
```

**Decision**: AptpSymbolic must handle both variants via pattern matching on `raw()`. The public API should abstract over this — callers should not need to know whether the result is simplicial or supernodal.

**Rationale**: Small test matrices may be simplicial; production matrices will typically be supernodal. Forcing supernodal on tiny matrices wastes memory and may degrade performance. The dual-variant handling is the standard faer pattern.

### RQ-4: Pivot Buffer Estimation Heuristic

**Question**: What is a reasonable heuristic for estimating APTP delayed-pivot buffer requirements?

**Findings**:

The ssids-plan.md proposes "10% buffer per column for delayed pivots." This means for each column in the factor, allocate 10% extra space for entries that may be delayed from descendant supernodes.

From the academic references:
- Hogg et al. (2016), Section 2.4: Delayed pivots propagate entries up the assembly tree. A delayed column from child supernode `c` adds its entire row pattern to parent supernode `p`.
- The number of delayed pivots depends on the matrix's numerical properties (not just structure), so any pre-factorization estimate is inherently heuristic.
- SPRAL's approach: allocate based on structure, then dynamically resize if needed during factorization.

**Decision**: Implement the 10% heuristic as the initial approach, applied per-supernode (not per-column) for supernodal mode. Each supernode's buffer estimate = ceil(0.10 * column_count_in_supernode). For simplicial mode, apply per-column: buffer[j] = ceil(0.10 * nnz_in_column_j_of_L). Store as `Vec<usize>` with one entry per supernode (supernodal) or per column (simplicial).

**Rationale**: 10% is a conservative starting point. The estimate will be validated empirically during Phase 5-6 numeric factorization. The heuristic is cheap to compute and easy to tune later (it's just a multiplier).

**Alternatives considered**:
- Fixed absolute buffer per column → doesn't scale with fill-in
- No buffer estimation (always dynamic allocation) → Phase 5-6 would need to handle reallocation in the hot loop
- Analytical bound from inertia prediction → inertia isn't known at symbolic phase

### RQ-5: Error Handling Strategy

**Question**: How should errors from faer's `factorize_symbolic_cholesky` be mapped to our `SparseError`?

**Findings**:

`factorize_symbolic_cholesky` returns `Result<SymbolicCholesky<I>, FaerError>`. The `FaerError` enum has limited variants — primarily allocation failures and internal errors.

Pre-validation (dimension mismatch, non-square, invalid custom permutation) should happen before calling faer.

**Decision**:
- Validate inputs (matrix dimensions, permutation validity) → `SparseError::DimensionMismatch` or `SparseError::InvalidInput`
- Map `FaerError` → `SparseError::AnalysisFailure` with the faer error message as context
- Structural singularity is not directly detected during symbolic analysis (it manifests at the numeric phase via zero pivots), so `SparseError::StructurallySingular` is not used here

### RQ-6: Module Placement

**Question**: Where should AptpSymbolic live in the codebase?

**Findings**:

Current aptp module (`src/aptp/`) contains Phase 2 data structures: PivotType, Block2x2, MixedDiagonal, Inertia, perm_from_forward. These are all "numeric-phase" data structures.

AptpSymbolic is an "analysis-phase" structure that precedes all numeric operations.

**Decision**: Add `src/aptp/symbolic.rs` to the existing aptp module, exported as `aptp::AptpSymbolic`. Keep `SymbolicStatistics` in the same file.

**Rationale**: AptpSymbolic is part of the APTP algorithm (it's the "analyze" step of the three-phase analyze→factorize→solve pipeline). Keeping it in `src/aptp/` maintains cohesion with the other APTP types it will interact with. A separate `src/symbolic/` module would scatter related types across the codebase.

**Alternatives considered**:
- `src/symbolic/mod.rs` → separates analysis from the APTP types it feeds into
- `src/aptp/analysis.rs` → fine but `symbolic.rs` matches faer's naming convention
