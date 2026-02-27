# Implementation Plan: MC64 Matching & Scaling

**Branch**: `012-mc64-matching-scaling` | **Date**: 2026-02-13 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `specs/012-mc64-matching-scaling/spec.md`

## Summary

Implement the MC64 weighted bipartite matching algorithm (Duff & Koster 2001) with symmetric scaling (Duff & Pralet 2005, MC64SYM) for preprocessing sparse symmetric indefinite matrices. MC64 is genuinely new code — faer has no matching or scaling functionality. The algorithm uses Dijkstra-based shortest augmenting paths on a logarithmic cost graph to find a maximum-product matching, then derives symmetric scaling factors from the dual variables. This improves diagonal dominance, reducing delayed pivots in subsequent APTP factorization.

MC64 matching and METIS ordering operate independently: MC64 produces scaling (for Phase 5 numeric factorization) and matching metadata (pivot classification), while METIS produces fill-reducing ordering (for Phase 3 symbolic analysis). SPRAL-style graph condensation (`mo_split()`) is deferred as a lightweight follow-up.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (sparse matrix types, permutations), std::collections::BinaryHeap (Dijkstra priority queue)
**Storage**: N/A (in-memory data structures only)
**Testing**: cargo test + hand-constructed matrices + SuiteSparse CI subset (10 matrices) + full SuiteSparse (67 matrices, `#[ignore]`)
**Target Platform**: Linux (development), cross-platform (library)
**Project Type**: Single Rust library crate (existing)
**Performance Goals**: Complete on all CI subset matrices (largest ~100K dimension) without timeout. O(n(tau + n) log n) worst case.
**Constraints**: No new external crate dependencies. Pure Rust implementation. Clean room from academic references + SPRAL (BSD-3).
**Scale/Scope**: ~600-900 lines implementation + ~400-600 lines tests. Single module file `src/aptp/matching.rs` following `ordering.rs` pattern.

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Validation via SPRAL's scaling property criteria (entries <= 1, row max >= 0.75, unit diagonal for matched entries). Hand-constructed + SuiteSparse tests. |
| II. Clean Room | PASS | Implemented from Duff & Koster (2001) and Duff & Pralet (2005) academic papers. SPRAL (BSD-3) consulted for patterns and test approach. No HSL source code. |
| III. TDD | PASS | Tests encode scaling property constraints from SPRAL's `tests/scaling.f90`. Red-green-refactor cycle. |
| IV. Documentation | PASS | Algorithm references cited. Rustdoc with academic attribution. |
| V. Numerical Stability | PASS | Log-domain cost formulation avoids overflow. Duff-Pralet correction handles structural singularity. Scaling factors returned in linear domain with log-domain option. |
| VI. Structured Development | PASS | Phase 4.2 follows completed Phase 4.1 (METIS). Exit criteria defined. |
| VII. Code Quality | PASS | Result-based errors, faer types at boundary, no unwrap in production code. |

No violations. All gates pass.

## Project Structure

### Documentation (this feature)

```text
specs/012-mc64-matching-scaling/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0: research decisions
├── data-model.md        # Phase 1: entity definitions
├── quickstart.md        # Phase 1: build & test guide
├── contracts/
│   └── api.md           # Phase 1: API contract
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code

```text
src/aptp/
├── mod.rs               # Add: pub mod matching; pub use matching::{mc64_matching, Mc64Result, Mc64Job};
├── matching.rs           # NEW: MC64 matching algorithm (~600-900 lines)
├── ordering.rs          # Existing: METIS ordering (Phase 4.1)
├── symbolic.rs          # Existing: AptpSymbolic (Phase 3)
├── diagonal.rs          # Existing: MixedDiagonal (Phase 2)
├── perm.rs              # Existing: perm_from_forward (Phase 2)
├── pivot.rs             # Existing: PivotType, Block2x2 (Phase 2)
└── inertia.rs           # Existing: Inertia (Phase 2)

tests/
├── mc64_matching.rs      # NEW: Integration tests
├── metis_ordering.rs    # Existing: reference for test patterns
└── ...                  # Other existing test files
```

**Structure Decision**: Single new file `src/aptp/matching.rs` following the established pattern from `ordering.rs`. No new directories, no new dependencies. Integration tests in `tests/mc64_matching.rs`.

## Algorithm Design

### Core Algorithm: Dijkstra-Based Augmenting Paths (Duff & Koster 2001, Algorithm MPD)

**Cost matrix formulation**:
```
c[i,j] = log(max_k |a[k,j]|) - log|a[i,j]|   for a[i,j] != 0
c[i,j] = infinity                               for a[i,j] = 0
```
All costs are non-negative. Minimizing the assignment cost maximizes the diagonal product.

**Dual variables and reduced costs**:
```
Extremality: u[i] + v[j] <= c[i,j]     for all edges
             u[i] + v[j] = c[i,j]      for matched edges
Reduced cost: c_bar[i,j] = c[i,j] - u[i] - v[j] >= 0
```

**Algorithm outline**:
1. Build cost graph from input matrix (expand upper-triangle to full bipartite graph)
2. Compute initial dual variables via column-minimum and row-minimum passes
3. Build greedy initial matching from zero-reduced-cost edges (~80% cardinality)
4. For each unmatched column j0:
   a. Run Dijkstra from j0 on the reduced-cost graph
   b. If augmenting path found: augment matching and update duals
   c. If no path exists: column j0 is part of the structurally singular block
5. Extract matching permutation and compute symmetric scaling from duals

**Dijkstra details** (per augmentation, from root column j0):
- Maintain distance array `d[i]` from j0 to each row node i
- Priority queue (BinaryHeap) for frontier nodes
- For each extracted node i: follow matched edge to column m[i], scan column m[i]'s neighbors
- Stop when `lsap <= lsp` (shortest augmenting path found) or queue empty

**Path reconstruction**: Parent pointers `pr[j]` trace back from the augmenting path endpoint to root j0. Alternate matched/unmatched edges along the path to augment.

**Dual variable update** after augmentation:
```
u'[i] = u[i] + d[i] - lsap    for visited row nodes i
v'[j] = c[i,j] - u'[i]        for each matched edge (i,j)
```

### Symmetric Scaling: MC64SYM (Duff & Pralet 2005)

**Symmetrization**: From row duals u and column duals v:
```
scaling[i] = exp(-(u[i] + v[i]) / 2)
```

**Property 4.1 guarantee**: The scaled matrix `D A D` (where `D = diag(scaling)`) satisfies:
- All entries: `|d_i * a_ij * d_j| <= 1`
- Matched diagonal: `|d_i * a_{i,σ(i)} * d_{σ(i)}| = 1`
- Row/column maxima = 1

### Structural Singularity: Duff-Pralet Correction

When matching cardinality < n:
1. Identify matched indices I (rows with `row_match[i].is_some()`)
2. Extract submatrix A[I, I] (structurally nonsingular by Property 4.2)
3. Re-run MC64 on submatrix to get proper weighted matching + duals
4. Compute scaling for matched indices from submatrix duals
5. For unmatched index i: `scaling[i] = 1.0 / max_k |a[i,k] * scaling[k]|` over matched k
6. Convention: if max is 0 (isolated row), `scaling[i] = 1.0`
7. Place unmatched indices at end of matching permutation

### Input Handling: Upper-Triangle Expansion

The input matrix is in faer's upper-triangular CSC format. MC64 needs the full bipartite graph:
1. Scan upper triangle: for each stored entry (i, j) with i < j, add edges (i,j) and (j,i)
2. Diagonal entries: add edge (i,i) for each stored diagonal
3. Build full CSC cost graph with `c[i,j] = log(col_max_j) - log|a[i,j]|`

This mirrors the adjacency extraction in `ordering.rs::extract_adjacency` but retains numeric values for cost computation.

## Testing Strategy

### Layer 1: Unit Tests (inline in `matching.rs`)

Internal function tests with small hand-constructed examples:
- `build_cost_graph`: Verify cost matrix for a known 3x3 matrix
- `greedy_initial_matching`: Verify cardinality and dual feasibility on a 4x4 matrix
- `dijkstra_augment`: Verify augmenting path on a small bipartite graph
- `symmetrize_scaling`: Verify `exp(-(u+v)/2)` computation against known duals
- `duff_pralet_correction`: Verify unmatched scaling on a known singular matrix
- Matching cycle structure: verify only singletons + 2-cycles for symmetric input

### Layer 2: Integration Tests (`tests/mc64_matching.rs`)

SPRAL-style scaling property validation on hand-constructed matrices:
- **Scaling constraint**: `|s_i * a_ij * s_j| <= 1.0` for all entries
- **Row maximum constraint**: `max_j |s_i * a_ij * s_j| >= 0.75` for each row
- **Matching validity**: matched pairs correspond to nonzero entries, no duplicate matches
- **Matching cardinality**: `matched == n` for nonsingular matrices
- **Cycle structure**: only singletons + 2-cycles

Test matrices:
- Diagonal matrix (trivial: identity matching, unit scaling)
- Tridiagonal indefinite (singletons + 2-cycles)
- Arrow matrix (indefinite variant)
- Badly-scaled matrix (entries spanning many orders of magnitude)
- Matrix with zero diagonal entries
- Well-conditioned PD matrix (should produce near-identity matching)
- 1x1 and 2x2 trivial cases

### Layer 3: SuiteSparse CI Integration (`tests/mc64_matching.rs`, regular tests)

Run MC64 on 10 CI subset matrices:
- Verify full matching (matched == n) for nonsingular matrices
- Verify all scaling properties (constraints 1-4 from SPRAL)
- Verify scaling factors positive and finite
- Verify cycle structure (singletons + 2-cycles only)

### Layer 4: Full SuiteSparse Validation (`tests/mc64_matching.rs`, `#[ignore]`)

Run MC64 on all 67 SuiteSparse matrices:
- Same validation as Layer 3
- Also verify diagonal dominance improvement metric
- Log 2-cycle counts per matrix (useful metadata for future condensation)

### Cross-validation with METIS (Layer 3)

Verify independent composition:
- Compute MC64 scaling on matrix
- Compute METIS ordering on same matrix (independently)
- Feed METIS ordering into AptpSymbolic::analyze
- Verify symbolic analysis succeeds (scaling doesn't affect symbolic phase)

## Implementation Sequence

### Step 1: Module Skeleton + Types

Create `src/aptp/matching.rs` with:
- Public types: `Mc64Result`, `Mc64Job`
- Public function signature: `mc64_matching` (returns `Err(AnalysisFailure)` stub)
- Internal types: `CostGraph`, `MatchingState`
- Module registration in `mod.rs`
- Basic compilation test

### Step 2: Cost Graph Construction

Implement `build_cost_graph`:
- Expand upper-triangle to full CSC
- Compute column maxima
- Compute logarithmic costs
- Unit tests with known 3x3 matrix

### Step 3: Greedy Initialization

Implement `greedy_initial_matching`:
- Initial dual variables from min-cost passes
- Greedy matching from zero-reduced-cost edges
- Secondary matching (2-pass rearrangement from Duff & Koster 2001)
- Unit tests: verify feasibility of dual variables, ~80% cardinality

### Step 4: Dijkstra Augmenting Path

Implement `dijkstra_augment`:
- Priority queue (BinaryHeap) with distance tracking
- Path construction via parent pointers
- Augmentation along alternating path
- Dual variable update
- Unit tests: verify augmentation on small graph

### Step 5: Main Matching Loop

Wire steps 2-4 into `mc64_matching`:
- Input validation
- Trivial cases (n=0, n=1, diagonal matrix)
- Main loop: greedy init → augment unmatched columns
- Integration tests with hand-constructed matrices

### Step 6: Symmetric Scaling

Implement `symmetrize_scaling`:
- `D[i] = exp(-(u[i] + v[i]) / 2)`
- Unit tests against known dual values
- Integration tests: verify SPRAL scaling properties

### Step 7: Structural Singularity

Implement Duff-Pralet correction:
- Detect partial matching (matched < n)
- Extract nonsingular submatrix
- Re-run MC64 on submatrix
- Apply scaling correction for unmatched indices
- Tests with structurally singular matrix

### Step 8: SuiteSparse Validation

- Integration tests on CI subset (10 matrices)
- Full validation on 67 matrices (`#[ignore]`)
- Cross-validation with METIS ordering
- Log cycle structure statistics

### Step 9: Documentation & Cleanup

- Rustdoc with algorithm description and academic references
- Module-level documentation
- Update `ssids-log.md` with Phase 4.2 completion
- `cargo fmt`, `cargo clippy`, final test run

## Complexity Tracking

> No constitution violations to justify.

## Post-Phase 1 Constitution Re-Check

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | SPRAL-derived scaling property tests as primary validation. Hand-constructed + SuiteSparse coverage. |
| II. Clean Room | PASS | Duff & Koster (2001) + Duff & Pralet (2005) as primary sources. SPRAL (BSD-3) for test patterns. |
| III. TDD | PASS | Tests written before each implementation step. Scaling property checks from SPRAL `tests/scaling.f90`. |
| IV. Documentation | PASS | Algorithm references, rustdoc, data-model.md, contracts/api.md all prepared. |
| V. Numerical Stability | PASS | Log-domain costs, Duff-Pralet singularity correction, scaling factor validation. |
| VI. Structured Development | PASS | Phase 4.2 follows Phase 4.1. Clear step sequence with exit criteria per step. |
| VII. Code Quality | PASS | Single module, Result-based errors, faer types at boundary, no new dependencies. |
