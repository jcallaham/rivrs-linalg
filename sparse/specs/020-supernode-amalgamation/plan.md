# Implementation Plan: Supernode Amalgamation

**Branch**: `020-supernode-amalgamation` | **Date**: 2026-02-21 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/020-supernode-amalgamation/spec.md`

## Summary

Implement a SPRAL-style supernode amalgamation pass that merges small parent-child supernode pairs after faer's symbolic analysis. This addresses the primary performance bottleneck on Schenk optimization matrices (c-71: 24.5x slower than SPRAL, c-big: 11.1x) caused by ~35K narrow supernodes where SPRAL has ~8K. The amalgamation is a post-processing transformation of `Vec<SupernodeInfo>` in `numeric.rs`, using SPRAL's two-condition merge predicate (structural match OR both-nodes-small with nemin=32). No changes to faer, the dense APTP kernel, or the solve path are required.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22, rayon 1.x (existing)
**Storage**: N/A (in-memory transformation of supernode metadata)
**Testing**: cargo test (unit + SuiteSparse integration), criterion benchmarks, baseline collection
**Target Platform**: Linux (primary), macOS (secondary)
**Project Type**: Single Rust library crate
**Performance Goals**: c-71/c-big within 2x of SPRAL factor time; no regression >10% on any other matrix; amalgamation pass <5% of symbolic analysis time
**Constraints**: All safe Rust (no unsafe). Must maintain backward error <5e-11 on all 65 SuiteSparse matrices. Must not modify faer.
**Scale/Scope**: ~600-800 LOC new module, ~50-100 LOC modifications to existing modules. ~15-20 new tests.

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Backward error <5e-11 maintained as hard constraint (FR-010). Full SuiteSparse regression suite required before merge. |
| II. Clean Room | PASS | Algorithm from SPRAL (BSD-3) `core_analyse.f90:806-853` and Liu (1992). No HSL or GPL sources consulted. |
| III. TDD | PASS | Tests written before implementation: unit tests for merge predicate, integration tests for supernode count reduction, regression tests for backward error. |
| IV. Documentation | PASS | Academic references cited in spec. Module will document SPRAL source and Liu (1992) as algorithm origins. |
| V. Numerical Stability | PASS | Amalgamation is a symbolic transformation — does not affect numerical computation paths. APTP kernel and solve are unchanged. Correctness validated by existing backward error tests. |
| VI. Structured Development | PASS | Phase 9.1a per ssids-plan.md. Phase 8 exit criteria met. |
| VII. Code Quality | PASS | New module follows existing patterns (SupernodeInfo, children map). Result type for errors. No unsafe. |

**Post-design re-check**: All gates still pass. No new dependencies, no changes to numerical code paths.

## Project Structure

### Documentation (this feature)

```text
specs/020-supernode-amalgamation/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0 research decisions
├── data-model.md        # Data model documentation
├── quickstart.md        # Developer quickstart
├── checklists/
│   └── requirements.md  # Spec quality checklist
└── tasks.md             # Phase 2 output (/speckit.tasks command)
```

### Source Code

```text
src/aptp/
├── amalgamation.rs      # NEW — amalgamation pass
├── numeric.rs           # MODIFIED — calls amalgamation after build_supernode_info()
├── solver.rs            # MODIFIED — adds nemin to AnalyzeOptions/SolverOptions
├── mod.rs               # MODIFIED — declares amalgamation module
└── (all other files unchanged)
```

**Structure Decision**: Single new module `amalgamation.rs` in the existing `src/aptp/` directory. This follows the established pattern where each algorithmic concern has its own module (factor.rs, solve.rs, ordering.rs, matching.rs, etc.).

## Design

### Architecture Overview

The amalgamation pass is a pure transformation: `Vec<SupernodeInfo> → Vec<SupernodeInfo>` that runs between `build_supernode_info()` and `factor_tree_levelset()` in `AptpNumeric::factor()`.

```
AptpNumeric::factor()
    │
    ├── build_supernode_info(symbolic)          // existing: ~35K supernodes for c-71
    │
    ├── amalgamate(supernodes, nemin)           // NEW: ~8K supernodes for c-71
    │
    ├── build_children_map(&supernodes)         // existing: rebuilt from merged parents
    │
    └── factor_tree_levelset(supernodes, ...)   // existing: unchanged
```

### Core Algorithm: `amalgamate()`

**Input**: `Vec<SupernodeInfo>` (fundamental supernodes from faer), `nemin: usize`
**Output**: `Vec<SupernodeInfo>` (amalgamated supernodes, fewer entries)

**Algorithm** (adapted from SPRAL `find_supernodes` + `do_merge`, `core_analyse.f90:618-641`):

1. Build children lists from parent pointers (same as `build_children_map`)
2. Process parents in ascending order (postorder ensures children processed first)
3. For each parent, iterate over ALL children and check `do_merge`:
   - Condition (a): `nelim(parent) == 1 AND cc(parent) == cc(child) - 1` (structural match)
   - Condition (b): `nelim(parent) < nemin AND nelim(child) < nemin` (both small)
4. If merge: absorb child's columns and pattern into parent, mark child as deleted
5. After all merges: compact the `Vec<SupernodeInfo>`, renumber parent pointers

**Merge operation** (adapted from SPRAL `merge_nodes`, `core_analyse.f90:827-853`):

When merging child C into parent P:
- `P.col_begin = min(P.col_begin, C.col_begin)` (extend range downward)
- `P.col_end = max(P.col_end, C.col_end)` (extend range upward)
- `P.pattern = sorted_union(P.pattern, C.pattern)` minus columns in `[P.col_begin, P.col_end)`
- `P.nelim += C.nelim` (tracked separately for merge decisions)
- C's children become P's children (reparent)
- C is marked deleted

**Column count (cc)**: In SPRAL, `cc(node)` is the total column count including pattern rows. For our `SupernodeInfo`: `cc = (col_end - col_begin) + pattern.len()`.

### Key Design Decisions (from research.md)

1. **No column renumbering**: faer's supernodes are already contiguous. Non-contiguous merges (nemin condition) produce ranges that span other supernodes' columns, but this is safe — `scatter_original_entries` only scatters the supernode's own columns, and extra columns in the frontal matrix are zero-initialized.

2. **Pattern union**: Merged pattern = `sorted_union(parent.pattern, child.pattern) - merged_columns`. The fill entries are implicit zeros in the frontal matrix.

3. **Postorder preservation**: Merging child into parent preserves postorder because we only delete children (lower-indexed) and absorb into parents (higher-indexed). The compacted vec maintains ascending order.

4. **nemin on AnalyzeOptions**: `nemin: usize` field with default 32. `nemin = 1` disables amalgamation. Propagated through `SolverOptions`.

### Integration Points

**Propagation path**: `AnalyzeOptions.nemin → SparseLDLT.nemin → AptpNumeric::factor(nemin) → amalgamate(supernodes, nemin)`

This mirrors the existing pattern for `scaling`: set during analyze, stored on
`SparseLDLT`, passed to `AptpNumeric::factor()` as a parameter.

**solver.rs changes** (~30 LOC):
- Add `pub nemin: usize` to `AnalyzeOptions` (default 32) and `SolverOptions` (default 32)
- Add `nemin: usize` field to `SparseLDLT` struct (set from `AnalyzeOptions` during `analyze()` / `analyze_with_matrix()`)
- In `SparseLDLT::factor()`: pass `self.nemin` to `AptpNumeric::factor()`
- In `SparseLDLT::solve_full()`: propagate `options.nemin` to `analyze_opts`

**numeric.rs changes** (~10 LOC):
- Add `nemin: usize` parameter to `AptpNumeric::factor()` signature
- After `let supernodes = build_supernode_info(symbolic);` (line ~405), call `amalgamate(supernodes, nemin)`
- The rest of `factor()` is unchanged — it already works with `Vec<SupernodeInfo>`

**mod.rs changes** (~2 LOC):
- Add `pub(crate) mod amalgamation;`

### Edge Cases

1. **All supernodes already large**: Merge loop runs but no merges trigger. Cost: one pass over the vec (O(n_supernodes)).
2. **Simplicial matrices (n supernodes)**: Amalgamation merges aggressively, but these matrices are Phase 9.1b's target. Amalgamation helps somewhat but doesn't address the per-supernode allocation overhead.
3. **Deep chains**: A chain of small supernodes (s1 → s2 → s3 → ... → sk) where all have <nemin cols. Processing in postorder: s1 merges into s2, then s2 (now larger) may or may not merge into s3 depending on the accumulated nelim.
4. **Multiple children merging into one parent**: Parent accumulates nelim from multiple children. After absorbing first child, parent's nelim increases, which may prevent subsequent merges if `nelim >= nemin`.

### Testing Strategy

**Unit tests** (amalgamation.rs):
- `test_no_merges_large_supernodes`: All supernodes >nemin → output identical to input
- `test_structural_match_merge`: Parent with 1 col, cc matches child → merge
- `test_nemin_merge`: Both parent and child <nemin → merge
- `test_nemin_no_merge_one_large`: One node ≥nemin → no merge
- `test_pattern_union`: Merged pattern is correct sorted union minus fully-summed cols
- `test_parent_reparenting`: Child's children become parent's children after merge
- `test_chain_merge`: Chain of small supernodes merges correctly
- `test_bushy_tree_merge`: Parent with many small children merges some/all
- `test_nemin_1_no_merge`: nemin=1 disables all merges
- `test_postorder_preserved`: Output supernodes maintain postorder invariant

**Integration tests** (existing test files, extended):
- c-71 supernode count <12K with default nemin
- c-big supernode count reduction proportional to c-71
- All 65 SuiteSparse matrices backward error <5e-11
- All 65 matrices factor time regression <10%
- Backward error identical with nemin=1 vs pre-amalgamation baseline

**Benchmark tests** (baseline_collection example):
- Before/after baseline comparison for c-71, c-big, and full suite

## Complexity Tracking

No constitution violations to justify.
