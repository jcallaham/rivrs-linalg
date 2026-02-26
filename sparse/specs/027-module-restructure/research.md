# Research: Module Restructure

**Date**: 2026-02-26
**Feature**: 027-module-restructure

## Decision 1: Which files move to `ordering/`

**Decision**: Move `aptp/ordering.rs`, `aptp/matching.rs`, `aptp/perm.rs`

**Rationale**: These three files form a self-contained ordering/preprocessing pipeline with no dependencies on APTP-specific types (PivotType, MixedDiagonal, etc.):
- `ordering.rs`: depends on `matching::{Mc64Job, mc64_matching}` (sibling) and `crate::error`
- `matching.rs`: depends only on `crate::error`
- `perm.rs`: depends on `crate::error` and `crate::validate`

**Alternatives considered**:
- Moving only `ordering.rs` + `matching.rs` (not `perm.rs`): Rejected because `perm_from_forward()` is used by `factor.rs` for permutation conversion — it's a general utility, not APTP-specific.
- Moving `symbolic.rs` too: Rejected because `AptpSymbolic` is tightly coupled to the symmetric solver's symbolic analysis (wraps faer's `SymbolicCholesky`, stores APTP-specific buffer estimates).

## Decision 2: File renaming within `ordering/`

**Decision**: Rename `ordering.rs` → `metis.rs` within the new `ordering/` module.

**Rationale**: Inside an `ordering/` module, having a file called `ordering.rs` is redundant and confusing. The file's content is METIS-specific: `metis_ordering()`, `match_order_metis()`, `MatchOrderResult`. The name `metis.rs` is self-documenting.

**Alternatives considered**:
- Keep as `ordering.rs`: Confusing — `ordering/ordering.rs`
- Name it `nested_dissection.rs`: More general but less discoverable; METIS is the well-known brand

## Decision 3: Internal module path within `symmetric/`

**Decision**: Keep flat structure (no `aptp/` subdirectory within `symmetric/`).

**Rationale**: The roadmap proposed `symmetric/aptp/` containing `kernel.rs`, `pivot.rs`, `diagonal.rs`. This adds a nesting level with no practical benefit — there's only one factorization strategy (APTP). If a second strategy is added (e.g., Bunch-Kaufman without APTP), restructuring at that point is cheap. For now, flat is simpler.

**Alternatives considered**:
- `symmetric/aptp/` subdirectory: Premature — no second consumer exists
- `symmetric/kernel/` subdirectory: Considered, but factor.rs + pivot.rs + diagonal.rs don't form a coherent sub-API

## Decision 4: `multifrontal/` extraction

**Decision**: Do NOT extract. Keep all multifrontal machinery (workspace, assembly maps, tree traversal, amalgamation) inside `symmetric/`.

**Rationale**: Per roadmap Section 4's "extract on second consumer" principle. The interfaces are tightly coupled to symmetric factorization specifics (FrontalMatrix shapes, MixedDiagonal storage, PivotType tracking). Designing generic interfaces without a concrete second use case would likely produce wrong abstractions.

## Dependency Graph Analysis

### Files that need import updates (internal `use crate::aptp::*`)

| File | Current import | New import |
|------|---------------|------------|
| `symmetric/factor.rs` | `super::perm::perm_from_forward` | `crate::ordering::perm_from_forward` |
| `symmetric/solver.rs` | `super::ordering::{match_order_metis, metis_ordering}` | `crate::ordering::{match_order_metis, metis_ordering}` |
| `io/reference.rs` | `crate::aptp::Inertia` | `crate::symmetric::Inertia` |
| `testing/mc64_validation.rs` | `crate::aptp::Mc64Result` + `crate::aptp::matching::count_cycles` | `crate::ordering::Mc64Result` + `crate::ordering::matching::count_cycles` |
| `ordering/metis.rs` | `super::matching::{Mc64Job, mc64_matching}` | `super::matching::{Mc64Job, mc64_matching}` (no change — still siblings) |

### Files that need import updates (external `rivrs_sparse::aptp::*`)

| File | Current import | New import |
|------|---------------|------------|
| `tests/mc64_matching.rs` | `rivrs_sparse::aptp::matching::count_cycles` + `rivrs_sparse::aptp::{Mc64Job, Mc64Result, mc64_matching}` | `rivrs_sparse::ordering::matching::count_cycles` + `rivrs_sparse::ordering::{Mc64Job, Mc64Result, mc64_matching}` |
| `tests/match_order.rs` | `rivrs_sparse::aptp::{AptpSymbolic, MatchOrderResult, Mc64Job, match_order_metis, mc64_matching, metis_ordering}` | Split: ordering types from `rivrs_sparse::ordering::*`, AptpSymbolic from `rivrs_sparse::symmetric::*` |
| `tests/metis_ordering.rs` | `rivrs_sparse::aptp::{AptpSymbolic, metis_ordering}` | `rivrs_sparse::ordering::metis_ordering` + `rivrs_sparse::symmetric::AptpSymbolic` |
| `tests/symbolic_analysis.rs` | `rivrs_sparse::aptp::AptpSymbolic` | `rivrs_sparse::symmetric::AptpSymbolic` |
| `tests/symbolic_analysis_full.rs` | `rivrs_sparse::aptp::{AptpSymbolic, metis_ordering}` | Split between ordering and symmetric |
| `tests/multifrontal.rs` | `rivrs_sparse::aptp::{...}` | Split: AptpNumeric/AptpOptions/AptpSymbolic → symmetric, metis_ordering → ordering |
| `tests/solve.rs` | `rivrs_sparse::aptp::{SparseLDLT, ...}` | `rivrs_sparse::symmetric::{SparseLDLT, ...}` |
| `tests/hand_constructed.rs` | `rivrs_sparse::aptp::{SparseLDLT, ...}` | `rivrs_sparse::symmetric::{SparseLDLT, ...}` |
| `tests/suitesparse_ci.rs` | `rivrs_sparse::aptp::{SparseLDLT, ...}` | `rivrs_sparse::symmetric::{SparseLDLT, ...}` |
| `tests/aptp_data_structures.rs` | `rivrs_sparse::aptp::{...}` | `rivrs_sparse::symmetric::{...}` |
| `tests/property.rs` | `rivrs_sparse::aptp::{SparseLDLT, ...}` | `rivrs_sparse::symmetric::{SparseLDLT, ...}` |
| `tests/adversarial.rs` | `rivrs_sparse::aptp::{SparseLDLT, ...}` | `rivrs_sparse::symmetric::{SparseLDLT, ...}` |
| All examples | `rivrs_sparse::aptp::{SparseLDLT, ...}` | `rivrs_sparse::symmetric::{SparseLDLT, ...}` |
| `benches/solver_benchmarks.rs` | `rivrs_sparse::aptp::{...}` | Split between ordering and symmetric |
| `comparisons/src/spral_benchmark.rs` | `rivrs_sparse::aptp::{...}` | `rivrs_sparse::symmetric::{...}` |

### Documentation files that need updates

- `sparse/CLAUDE.md` — Source Code Layout section, module references throughout
- `docs/ssids-log.md` — Add restructure entry
- `docs/ssids-plan.md` — Update if it references `aptp/` paths
- `docs/commands.md` — Update if it references `aptp/` paths
- Root `CLAUDE.md` — May reference `aptp/` in sparse section
