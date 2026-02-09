# Research: APTP Data Structures

**Feature**: 009-aptp-data-structures
**Date**: 2026-02-08

## R1: faer Perm<usize> Construction API

**Decision**: Use `Perm::new_checked(forward, inverse, dim)` where forward/inverse are `Box<[usize]>`.

**Rationale**: faer 0.22's `Perm<usize>` (with default shape `N = usize`) accepts `Box<[usize]>` for both arrays. The `new_checked` method validates that both arrays are valid permutations and are inverses of each other (panics on violation). With `N = usize`, `Idx<usize, usize>` collapses to plain `usize`.

**Alternatives considered**:
- `PermRef::new_checked(&[usize], &[usize], usize)` — borrowed view, doesn't own data. Useful for passing around but not for construction from ordering output which produces owned arrays.
- Custom permutation wrapper — rejected per transparent composition principle.

**Key API surface**:
- `Perm::new_checked(fwd: Box<[usize]>, inv: Box<[usize]>, dim: usize) -> Perm<usize>`
- `perm.inverse() -> PermRef<'_, usize>` — zero-copy swap of forward/inverse references
- `perm.as_ref() -> PermRef<'_, usize>` — borrow
- `perm.arrays() -> (&[usize], &[usize])` — access forward and inverse arrays
- No built-in composition operator — must implement manually via `composed[i] = p2[p1[i]]`

## R2: 2x2 Block Eigenvalue Sign Determination

**Decision**: Use trace/determinant to classify eigenvalue signs without computing eigenvalues.

**Rationale**: For inertia computation, we only need the *signs* of eigenvalues, not their values. The trace (a + c) and determinant (ac - b²) of a 2x2 symmetric block [[a, b], [b, c]] fully determine the sign classification:

| Condition | Eigenvalue signs | Inertia contribution |
|-----------|-----------------|---------------------|
| det > 0, trace > 0 | both positive | (+2, 0, 0) |
| det > 0, trace < 0 | both negative | (0, +2, 0) |
| det < 0 | one positive, one negative | (+1, +1, 0) |
| det = 0, trace > 0 | one positive, one zero | (+1, 0, +1) |
| det = 0, trace < 0 | one negative, one zero | (0, +1, +1) |
| det = 0, trace = 0 | both zero | (0, 0, +2) |

This avoids square root computation and is branch-free for the common cases.

**Alternatives considered**:
- Compute actual eigenvalues via quadratic formula — unnecessary overhead, and requires sqrt.
- Use faer's eigendecomposition — overkill for 2x2 blocks.
- Sylvester's criterion — only works for positive definiteness, not full inertia.

## R3: 2x2 Symmetric System Solve

**Decision**: Use Cramer's rule (analytical inverse) for 2x2 block solve.

**Rationale**: For a 2x2 symmetric system [[a, b], [b, c]] * [x1, x2]^T = [r1, r2]^T:
```
det = ac - b²
x1 = (c * r1 - b * r2) / det
x2 = (a * r2 - b * r1) / det
```

This is exact in floating point (no iterative refinement needed for 2x2) and requires only 5 multiplications, 2 subtractions, and 2 divisions.

**Alternatives considered**:
- LDL^T factorization of the 2x2 block — unnecessary complexity for a 2x2 system.
- faer's dense solve — would require constructing Mat objects for a 2x2 system.

## R4: Permutation Composition

**Decision**: Implement composition manually since faer doesn't provide it.

**Rationale**: faer's Perm type has no `Mul` operator or `compose` method. Composition is straightforward: `composed_fwd[i] = p2_fwd[p1_fwd[i]]` for applying p1 then p2. The inverse is computed from the composed forward array.

This is only needed in tests (acceptance scenario US3-2) and Phase 4 ordering, not in the hot path.

## R5: Inertia Struct Relocation

**Decision**: Move `Inertia` from `io/reference.rs` to the new APTP module, re-export from original location.

**Rationale**: `Inertia` is a domain concept (eigenvalue sign classification) used by both:
1. Reference factorization loading (`io/reference.rs`) — existing usage
2. MixedDiagonal inertia computation — new usage

The struct has no I/O-specific logic (no Deserialize dependency needed for the struct itself — serde derives can remain). Moving it to the APTP module and re-exporting preserves all existing imports via `use crate::io::reference::Inertia`.

**Impact**: `io/reference.rs` adds one `pub use` line. `validate.rs` import path unchanged (still imports from `crate::io::reference`). Zero behavioral changes.
