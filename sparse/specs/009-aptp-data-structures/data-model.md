# Data Model: APTP Data Structures

**Feature**: 009-aptp-data-structures
**Date**: 2026-02-08

## Entity Definitions

### PivotType

Classification of a single column's pivot decision during APTP factorization.

**Variants**:

| Variant | Fields | Description |
|---------|--------|-------------|
| OneByOne | (none) | Standard 1x1 scalar pivot |
| TwoByTwo | `partner: usize` | 2x2 Bunch-Kaufman pivot; `partner` is the index of the paired column |
| Delayed | (none) | Column failed APTP stability check; deferred to ancestor node |

**Invariants**:
- For TwoByTwo at column `i` with `partner = j`: column `j` must also be TwoByTwo with `partner = i`
- Partners are always adjacent: `|i - j| == 1`
- The lower-indexed column "owns" the Block2x2 data

**Derives**: Debug, Clone, Copy, PartialEq, Eq

### Block2x2

Storage for a single 2x2 symmetric diagonal block.

**Fields**:

| Field | Type | Description |
|-------|------|-------------|
| `first_col` | `usize` | Index of the first (lower-indexed) column of the 2x2 block |
| `a` | `f64` | Top-left element: D[i, i] |
| `b` | `f64` | Off-diagonal element: D[i, i+1] = D[i+1, i] |
| `c` | `f64` | Bottom-right element: D[i+1, i+1] |

**Computed properties**:
- Determinant: `a * c - b * b`
- Trace: `a + c`

**Invariants**:
- Represents the symmetric matrix `[[a, b], [b, c]]`
- `first_col + 1` must be within matrix dimension

**Derives**: Debug, Clone, Copy, PartialEq

### MixedDiagonal

The complete D factor in P^T A P = L D L^T, supporting mixed 1x1 and 2x2 blocks.

**Fields**:

| Field | Type | Visibility | Description |
|-------|------|------------|-------------|
| `pivot_map` | `Vec<PivotType>` | private | Per-column pivot classification, length = n |
| `diag_1x1` | `Vec<f64>` | private | Diagonal values for 1x1 pivots (length = n; entries at 2x2/delayed columns are unused) |
| `blocks_2x2` | `Vec<Block2x2>` | private | 2x2 blocks (one per pair; owned by the lower-indexed column) |
| `n` | `usize` | private | Matrix dimension |

**Construction pattern**: Incremental, column-by-column
1. `MixedDiagonal::new(n)` — all columns initially Delayed (unset)
2. `set_1x1(col, value)` — marks column as 1x1, stores diagonal value
3. `set_2x2(first_col, block)` — marks two columns as 2x2, stores block

**Query methods**:
- `dimension() -> usize` — returns n
- `get_pivot_type(col) -> PivotType` — pivot classification for column
- `get_1x1(col) -> f64` — diagonal value (debug-asserts column is 1x1)
- `get_2x2(first_col) -> &Block2x2` — block data (debug-asserts column is 2x2 owner)
- `num_delayed() -> usize` — count of Delayed columns
- `num_1x1() -> usize` — count of 1x1 pivots
- `num_2x2_pairs() -> usize` — count of 2x2 pivot pairs

**Operations**:
- `solve_in_place(x: &mut [f64])` — solves D x = b in place
  - Debug-asserts: no delayed columns, no singular pivots
  - 1x1: `x[i] /= d[i]`
  - 2x2: Cramer's rule on the 2x2 block
- `compute_inertia() -> Inertia` — eigenvalue sign counts
  - 1x1: sign of d[i]
  - 2x2: trace/determinant classification (see research.md R2)
  - Delayed: not counted (debug-asserts none exist, or counts as zero)

### Inertia (relocated)

Eigenvalue sign classification of a symmetric matrix. Currently defined in `io/reference.rs`, relocated to the APTP module.

**Fields** (unchanged):

| Field | Type | Description |
|-------|------|-------------|
| `positive` | `usize` | Count of positive eigenvalues |
| `negative` | `usize` | Count of negative eigenvalues |
| `zero` | `usize` | Count of zero eigenvalues |

**Methods** (unchanged):
- `dimension() -> usize` — positive + negative + zero

**Derives** (unchanged): Debug, Clone, PartialEq, Eq, Serialize, Deserialize

### perm_from_forward (function)

Standalone utility function, not a struct.

**Signature**: `perm_from_forward(fwd: Vec<usize>) -> Result<Perm<usize>, SparseError>`

**Algorithm**:
1. Validate: check for duplicates and out-of-bounds (reuse `validate_permutation`)
2. Compute inverse: `inv[fwd[i]] = i` for all i
3. Convert to boxed slices: `fwd.into_boxed_slice()`, `inv.into_boxed_slice()`
4. Construct: `Perm::new_checked(fwd_box, inv_box, n)`

**Error cases**:
- Duplicate indices → `SparseError::InvalidInput`
- Out-of-bounds indices → `SparseError::InvalidInput`

## Relationships

```
MixedDiagonal
  ├── contains Vec<PivotType>     (one per column)
  ├── contains Vec<f64>           (1x1 diagonal values)
  ├── contains Vec<Block2x2>     (2x2 block data)
  └── produces Inertia            (via compute_inertia)

perm_from_forward
  └── produces Perm<usize>        (faer type, not custom)

Inertia
  ├── used by MixedDiagonal::compute_inertia (new)
  ├── used by io::reference::ReferenceFactorization (existing)
  └── used by validate::check_inertia (existing)
```

## State Transitions

### MixedDiagonal Column States

```
Delayed (initial) ──set_1x1()──→ OneByOne (terminal)
Delayed (initial) ──set_2x2()──→ TwoByTwo (terminal, both columns)
```

All transitions are one-way. Once a column is set as 1x1 or 2x2, it cannot be changed. Setting a column that is already set is a debug-assert violation (programmer error).
