# Implementation Plan: Module Restructure

**Branch**: `027-module-restructure` | **Date**: 2026-02-26 | **Spec**: [spec.md](spec.md)
**Input**: Feature specification from `/specs/027-module-restructure/spec.md`

## Summary

Restructure the `rivrs-sparse` module hierarchy by extracting a top-level `ordering/` module from `aptp/` (containing METIS ordering, MC64 matching, and permutation utilities), renaming `aptp/` to `symmetric/`, and updating all imports across the crate. This is a mechanical refactor with no behavioral changes — all existing tests must pass with bit-identical results.

## Technical Context

**Language/Version**: Rust 1.87+ (edition 2024)
**Primary Dependencies**: faer 0.22, metis-sys 0.3, rayon, serde (unchanged)
**Storage**: N/A
**Testing**: cargo test (546+ tests), cargo test -- --ignored (65 SuiteSparse)
**Target Platform**: Linux (primary), cross-platform Rust
**Project Type**: Single Rust library crate
**Performance Goals**: Zero performance regression (this is a pure refactor)
**Constraints**: All 546+ tests pass, zero clippy warnings, cargo doc succeeds
**Scale/Scope**: ~30 source files, ~15 test files, 5 examples, 1 benchmark, 1 comparison driver

## Constitution Check

*GATE: Must pass before Phase 0 research. Re-check after Phase 1 design.*

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Correctness First | PASS | Pure refactor — no algorithm changes, no numerical behavior changes |
| II. Clean Room | PASS | No new code from external sources; just file moves and import rewrites |
| III. TDD | PASS | All existing tests preserved and updated; no new algorithmic code to TDD |
| IV. Documentation | PASS | CLAUDE.md and docs updated to reflect new structure |
| V. Numerical Stability | PASS | No numerical code changes |
| VI. Structured Development | PASS | Follows phased plan; documented in ssids-log |
| VII. Code Quality | PASS | Improved module organization; Rust import hygiene maintained |

## Project Structure

### Documentation (this feature)

```text
specs/027-module-restructure/
├── plan.md              # This file
├── spec.md              # Feature specification
├── research.md          # Dependency graph analysis
└── checklists/
    └── requirements.md  # Quality checklist
```

### Source Code Changes

```text
src/
├── lib.rs                          # UPDATE: mod aptp → mod ordering + mod symmetric
├── error.rs                        # UNCHANGED
├── validate.rs                     # UNCHANGED
│
├── io/                             # MINOR UPDATES
│   ├── mod.rs                      # UNCHANGED
│   ├── mtx.rs                      # UNCHANGED
│   ├── registry.rs                 # UNCHANGED
│   └── reference.rs                # UPDATE: crate::aptp::Inertia → crate::symmetric::Inertia
│
├── ordering/                       # NEW MODULE (extracted from aptp/)
│   ├── mod.rs                      # NEW: module declarations + pub use re-exports
│   ├── metis.rs                    # MOVED from aptp/ordering.rs, update: super:: → super::
│   ├── matching.rs                 # MOVED from aptp/matching.rs (minimal changes)
│   └── perm.rs                     # MOVED from aptp/perm.rs (minimal changes)
│
├── symmetric/                      # RENAMED from aptp/
│   ├── mod.rs                      # UPDATE: remove ordering/matching/perm, add crate::ordering imports
│   ├── pivot.rs                    # UNCHANGED
│   ├── diagonal.rs                 # UNCHANGED
│   ├── inertia.rs                  # UNCHANGED
│   ├── symbolic.rs                 # UNCHANGED
│   ├── factor.rs                   # UPDATE: super::perm → crate::ordering::perm
│   ├── numeric.rs                  # UNCHANGED
│   ├── solve.rs                    # UNCHANGED
│   ├── solver.rs                   # UPDATE: super::ordering → crate::ordering
│   └── amalgamation.rs             # UNCHANGED
│
├── profiling/                      # UNCHANGED
├── debug/                          # UNCHANGED
├── testing/
│   └── mc64_validation.rs          # UPDATE: crate::aptp → crate::ordering
└── benchmarking/                   # UNCHANGED
```

**Structure Decision**: Single Rust library crate with two top-level domain modules (`ordering`, `symmetric`) plus general infrastructure modules. No workspace changes.

## Implementation Steps

### Step 1: Create `ordering/` module (file moves)

1. Create `src/ordering/` directory
2. Copy `src/aptp/ordering.rs` → `src/ordering/metis.rs`
3. Copy `src/aptp/matching.rs` → `src/ordering/matching.rs`
4. Copy `src/aptp/perm.rs` → `src/ordering/perm.rs`
5. Create `src/ordering/mod.rs` with module declarations and re-exports

**ordering/mod.rs content**:
```rust
//! Fill-reducing orderings and matrix preprocessing.
//!
//! This module provides reusable ordering and scaling algorithms for sparse matrices:
//! - METIS nested dissection ordering ([`metis_ordering`])
//! - MC64 weighted bipartite matching and scaling ([`mc64_matching`])
//! - Combined matching + ordering pipeline ([`match_order_metis`])
//! - Permutation utilities ([`perm_from_forward`])

mod matching;
mod metis;
mod perm;

pub use matching::{Mc64Job, Mc64Result, mc64_matching};
pub use metis::{MatchOrderResult, match_order_metis, metis_ordering};
pub use perm::perm_from_forward;

// Re-export matching module for semi-public items (count_cycles used in tests)
pub use self::matching;
```

### Step 2: Update internal imports in moved files

**ordering/metis.rs** (was aptp/ordering.rs):
- Change `super::matching::{Mc64Job, mc64_matching}` → `super::matching::{Mc64Job, mc64_matching}` (stays the same — still sibling modules within ordering/)
- Verify no other `super::` references to aptp-specific modules

**ordering/matching.rs** (was aptp/matching.rs):
- Verify `crate::error` imports still work (they will — crate path unchanged)
- No `super::` references to other aptp modules expected

**ordering/perm.rs** (was aptp/perm.rs):
- Verify `crate::error` and `crate::validate` imports still work

### Step 3: Rename `aptp/` → `symmetric/`

1. Rename directory: `src/aptp/` → `src/symmetric/`
2. Update `src/lib.rs`: `mod aptp` → `pub mod ordering; pub mod symmetric;`
3. Remove the three moved files from `symmetric/` (ordering.rs, matching.rs, perm.rs)

### Step 4: Update `symmetric/mod.rs`

- Remove `mod ordering;`, `mod matching;`, `mod perm;` declarations
- Remove ordering-related `pub use` lines
- Add re-exports that pull ordering types through symmetric for backward compat (or not — clean break)
- Keep all other `pub use` re-exports unchanged

### Step 5: Update internal cross-module imports

**symmetric/factor.rs**:
- `super::perm::perm_from_forward` → `crate::ordering::perm_from_forward`

**symmetric/solver.rs**:
- `super::ordering::{match_order_metis, metis_ordering}` → `crate::ordering::{match_order_metis, metis_ordering}`

**io/reference.rs**:
- `crate::aptp::Inertia` → `crate::symmetric::Inertia`

**testing/mc64_validation.rs**:
- `crate::aptp::Mc64Result` → `crate::ordering::Mc64Result`
- `crate::aptp::matching::count_cycles` → `crate::ordering::matching::count_cycles`

### Step 6: Update integration tests

All `tests/*.rs` files: replace `rivrs_sparse::aptp::` with appropriate module:
- Ordering types → `rivrs_sparse::ordering::`
- Solver types → `rivrs_sparse::symmetric::`

Files to update: `aptp_data_structures.rs`, `mc64_matching.rs`, `mc64_profiling.rs`, `match_order.rs`, `metis_ordering.rs`, `multifrontal.rs`, `symbolic_analysis.rs`, `symbolic_analysis_full.rs`, `solve.rs`, `hand_constructed.rs`, `suitesparse_ci.rs`, `property.rs`, `adversarial.rs`

### Step 7: Update examples, benchmarks, comparisons

All `examples/*.rs`: replace `rivrs_sparse::aptp::` → `rivrs_sparse::symmetric::`
`benches/solver_benchmarks.rs`: split imports between ordering and symmetric
`comparisons/src/spral_benchmark.rs`: replace `rivrs_sparse::aptp::` → `rivrs_sparse::symmetric::`

### Step 8: Verify and fix

1. `cargo build` — fix any remaining import errors
2. `cargo test` — all tests must pass
3. `cargo clippy --all-targets --features diagnostic` — zero warnings
4. `cargo fmt --check` — formatting clean
5. `cargo doc --no-deps` — documentation builds

### Step 9: Update documentation

1. Update `sparse/CLAUDE.md` Source Code Layout section
2. Add entry to `dev/ssids-log.md`
3. Update any `aptp/` path references in `dev/` files
4. Update root `CLAUDE.md` if it references `aptp/`

## Complexity Tracking

No constitution violations. This is a pure refactor with no new algorithmic code, no new dependencies, and no behavioral changes.
