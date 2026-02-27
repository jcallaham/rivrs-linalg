# Implementation Plan: Multifrontal Numeric Factorization

**Branch**: `015-multifrontal-factorization` | **Date**: 2026-02-15 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/015-multifrontal-factorization/spec.md`

## Summary

Implement the multifrontal numeric factorization loop for sparse symmetric indefinite matrices using APTP (A Posteriori Threshold Pivoting). The factorization traverses the assembly tree in postorder, assembling dense frontal matrices from original sparse entries and child contribution blocks, then factoring each front using Phase 5's `aptp_factor_in_place()` kernel. Per-supernode factors (L11, D11, L21) are stored in an `AptpNumeric` result struct. Delayed pivots propagate to parent supernodes as additional fully-summed columns.

**Key architectural insight**: The Phase 5 kernel's Schur complement updates propagate to ALL trailing rows when the entire frontal matrix is passed (verified in research). This means L21 and F22 are computed implicitly — no separate TRSM/GEMM steps are needed in Phase 6.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22 (dense LA, sparse types, permutations), Phase 2/3/5 modules (existing)
**Storage**: N/A (in-memory dense matrices only)
**Testing**: cargo test + criterion benchmarks. Reconstruction error validation via existing `validate.rs`
**Target Platform**: Linux (x86_64, aarch64)
**Project Type**: Single Rust library crate (existing `sparse/`)
**Performance Goals**: Correctness-first. Sequential factorization. BLAS-3 optimization deferred to Phase 8.1
**Constraints**: Reconstruction error < 10^-12 for all test matrices. No new external dependencies
**Scale/Scope**: Factor sparse matrices up to ~1.6M dimension (SuiteSparse collection). Single new source file (~1000-1500 lines) + test file

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Evidence |
|-----------|--------|----------|
| I. Correctness First | PASS | Reconstruction error < 10^-12 as primary validation. SC-001 through SC-007 define measurable correctness criteria |
| II. Clean Room | PASS | Algorithm from Duff & Reid 1983, Liu 1992, Hogg et al. 2016/2020 (all BSD-3 or academic). SPRAL consulted as reference. No HSL source |
| III. TDD | PASS | Test strategy defined: hand-constructed matrices, CI subset, full SuiteSparse. Tests written before implementation |
| IV. Documentation | PASS | Algorithm references cited with file paths. Rustdoc required for all public APIs |
| V. Numerical Stability | PASS | APTP with threshold pivoting, 2x2 Bunch-Kaufman fallback, delayed pivot propagation. Inertia validation |
| VI. Structured Development | PASS | Phase 6 builds on completed Phases 2, 3, 5. Exit criteria defined |
| VII. Code Quality | PASS | Uses faer types at boundary, Result<T,E> for errors, type-enforced correctness |

**Post-design re-check**: All principles remain satisfied. The "pass entire frontal matrix" strategy is simpler than separate TRSM/GEMM, reducing code complexity without sacrificing correctness.

## Project Structure

### Documentation (this feature)

```text
specs/015-multifrontal-factorization/
├── spec.md              # Feature specification
├── plan.md              # This file
├── research.md          # Phase 0: research decisions
├── data-model.md        # Phase 1: entity definitions
├── quickstart.md        # Phase 1: build/test/usage guide
├── contracts/
│   ├── numeric_api.md   # Public API contract (AptpNumeric)
│   └── internal_api.md  # Internal API contracts (SupernodeInfo, FrontalMatrix, etc.)
└── tasks.md             # Phase 2 output (created by /speckit.tasks)
```

### Source Code

```text
src/aptp/
├── mod.rs              # MODIFIED: add pub mod numeric; re-export new public types
├── numeric.rs          # NEW: multifrontal factorization (~1000-1500 lines)
│   ├── SupernodeInfo           (pub(crate))
│   ├── FrontalMatrix           (pub(crate))
│   ├── ContributionBlock       (pub(crate))
│   ├── FrontFactors            (pub)
│   ├── FactorizationStats      (pub)
│   ├── AptpNumeric             (pub)
│   ├── build_supernode_info()  (pub(crate))
│   ├── assemble_frontal()      (pub(crate))
│   ├── extend_add()            (pub(crate))
│   └── extract_front_factors() (pub(crate))
├── symbolic.rs         # EXISTING: AptpSymbolic (may need minor additions)
├── factor.rs           # EXISTING: aptp_factor_in_place (no changes expected)
├── diagonal.rs         # EXISTING: MixedDiagonal (no changes expected)
├── pivot.rs            # EXISTING: PivotType, Block2x2 (no changes expected)
├── inertia.rs          # EXISTING: Inertia (no changes expected)
├── perm.rs             # EXISTING: perm_from_forward (no changes expected)
├── ordering.rs         # EXISTING: METIS ordering (no changes expected)
└── matching.rs         # EXISTING: MC64 matching (no changes expected)

tests/
├── hand_constructed.rs # MODIFIED: add multifrontal reconstruction tests
├── suitesparse_ci.rs   # MODIFIED: add multifrontal factorization tests
└── multifrontal.rs     # NEW: dedicated multifrontal tests
    ├── assembly correctness tests
    ├── extend-add tests
    ├── delayed pivot propagation tests
    ├── single-supernode-matches-dense tests
    └── full reconstruction tests
```

**Structure Decision**: Single new module `src/aptp/numeric.rs` within the existing `aptp` module hierarchy. This follows the established pattern (symbolic.rs, factor.rs) and keeps the APTP implementation cohesive. One new test file `tests/multifrontal.rs` for dedicated multifrontal tests, plus additions to existing integration tests.

## Algorithm Overview

### Factorization Loop

```
function multifrontal_factor(symbolic, matrix, options):
    // Precompute
    supernodes = build_supernode_info(symbolic)
    children = build_children_map(supernodes)
    perm = symbolic.perm()
    global_to_local = array of size n, initialized to SENTINEL
    contributions = array of Option<ContributionBlock>, size n_supernodes

    // Postorder traversal (index order = postorder)
    for s in 0..n_supernodes:
        // 1. Determine frontal matrix structure
        child_delayed = collect delayed columns from children[s]
        fully_summed = supernode_cols(s) + child_delayed
        k = fully_summed.len()
        frontal_rows = fully_summed ∪ supernodes[s].pattern
        m = frontal_rows.len()

        // 2. Build global-to-local mapping
        for (i, &global) in frontal_rows.iter().enumerate():
            global_to_local[global] = i

        // 3. Assemble frontal matrix
        frontal = Mat::zeros(m, m)
        scatter_original_entries(frontal, matrix, perm, frontal_rows, global_to_local)
        for c in children[s]:
            extend_add(frontal, contributions[c], global_to_local)
            contributions[c] = None  // deallocate

        // 4. Factor
        result = aptp_factor_in_place(frontal, k, options)?
        ne = result.num_eliminated

        // 5. Extract and store
        front_factors[s] = extract_front_factors(frontal, result, supernodes[s])
        accumulate_stats(stats, result, m)

        // 6. Prepare contribution for parent
        if supernodes[s].parent.is_some() and ne < m:
            contributions[s] = extract_contribution(frontal, result, frontal_rows, k)

        // 7. Cleanup global-to-local
        for &global in frontal_rows:
            global_to_local[global] = SENTINEL

    return AptpNumeric { front_factors, stats, n }
```

### Assembly: Scatter Original Entries

For each column j in the supernode's column range `[col_begin..col_end)`:
1. Map j through the fill-reducing permutation: `perm_j = perm_fwd[j]`
2. For each nonzero in column `perm_j` of the original matrix:
   - Get row index: `perm_r = row_index`
   - Map through inverse permutation: `r = perm_inv[perm_r]`
   - Map to local position: `local_r = global_to_local[r]`, `local_c = global_to_local[j]`
   - If both in frontal: `frontal[local_r, local_c] += value`
   - For symmetric: also set `frontal[local_c, local_r] += value` if needed (lower triangle)

### Assembly: Extend-Add

For each child contribution block:
1. For each entry at local position (i, j) in the contribution:
   - Map to global: `gi = contribution.row_indices[i]`, `gj = contribution.row_indices[j]`
   - Map to parent local: `li = global_to_local[gi]`, `lj = global_to_local[gj]`
   - Add: `parent_frontal[li, lj] += contribution.data[i, j]`
2. Only process lower triangle (i >= j)

### Extraction: Front Factors

After `aptp_factor_in_place()` on the (m x m) frontal matrix with num_fully_summed = k:
1. ne = result.num_eliminated
2. L11 = copy of frontal[0..ne, 0..ne] (lower triangle)
3. D11 = result.d (first ne entries of MixedDiagonal)
4. L21 = copy of frontal[k..m, 0..ne]
5. local_perm = result.perm[0..k].to_vec()
6. col_indices = local_perm[0..ne] mapped through frontal_row_indices
7. row_indices = frontal_row_indices[k..m].to_vec()

### Extraction: Contribution Block

1. size = m - ne
2. data = copy of frontal[ne..m, ne..m]
3. row_indices = [delayed_global_indices..., frontal_row_indices[k..m]...]
4. num_delayed = k - ne

## Key Implementation Notes

### Permutation Tracking

There are three levels of permutation:
1. **Fill-reducing permutation** (from `AptpSymbolic::perm()`): maps original matrix indices to permuted indices. Applied globally during scatter.
2. **Supernode-to-frontal mapping** (`global_to_local`): maps permuted global indices to local frontal matrix positions. Computed per-supernode.
3. **APTP local permutation** (`AptpFactorResult.perm`): maps post-pivot positions to original frontal-local positions. Stored per-supernode in `FrontFactors.local_perm`.

### MixedDiagonal Truncation

`aptp_factor_in_place()` returns a `MixedDiagonal` of size `m` (the full frontal matrix dimension). For FrontFactors, we need only the first `ne` entries. The MixedDiagonal needs to be truncated or a new one built from the first `ne` pivot entries.

### Swap Effects on F21

When the APTP kernel swaps columns within the fully-summed region (delayed columns swapped to end), the corresponding F21 column entries are also permuted (via `swap_symmetric`). The L21 block at `frontal[k..m, 0..ne]` has columns in APTP-permuted order. The `local_perm` maps these back to original column positions.

### Memory Management

- **Frontal matrices**: Allocated per-supernode, deallocated after extraction
- **Contribution blocks**: Stored in `Vec<Option<ContributionBlock>>`, deallocated after extend-add into parent
- **FrontFactors**: Persists for the lifetime of `AptpNumeric` (needed by Phase 7)
- **global_to_local**: Single allocation of size n, reused across supernodes (cleaned up per-iteration)

## Testing Strategy

### Unit Tests (in `src/aptp/numeric.rs` or `tests/multifrontal.rs`)

1. **build_supernode_info**: Verify correct extraction for both supernodal and simplicial symbolic analysis
2. **assemble_frontal (leaf)**: Hand-constructed matrix, verify scattered entries match expected dense block
3. **assemble_frontal (with children)**: Verify extend-add correctly merges contribution blocks
4. **extract_front_factors**: Verify L11, D11, L21 extraction from known factored matrix
5. **extract_contribution**: Verify contribution block and row indices
6. **extend_add**: Two supernodes with known structure, verify parent assembly matches dense factorization
7. **Single-supernode regression**: Matrix that is one supernode — must match Phase 5's dense APTP exactly
8. **Delayed pivot propagation**: Construct matrix where child delays, parent succeeds
9. **All-delayed front**: Front where every column is delayed — no factors stored, all propagated
10. **Simplicial path**: Small matrix where faer chooses simplicial, verify correctness

### Integration Tests

11. **Hand-constructed matrices**: All 15 matrices, reconstruction error < 10^-12
12. **CI SuiteSparse subset**: All 10 matrices, reconstruction error < 10^-12
13. **Dense equivalence**: For small matrices (n < 100), verify multifrontal matches dense APTP
14. **Inertia validation**: All hand-constructed matrices with known inertia
15. **Statistics accuracy**: Verify pivot counts match expected for small matrices

### Ignored Tests (full SuiteSparse)

16. **Full SuiteSparse**: All 67 matrices, reconstruction error < 10^-12

## Complexity Tracking

No constitution violations. No complexity justifications needed.
