# Quickstart: APTP Data Structures

**Feature**: 009-aptp-data-structures
**Date**: 2026-02-08

## What This Feature Adds

New `aptp` module in `rivrs-sparse` providing data structures for the APTP (A Posteriori Threshold Pivoting) algorithm:

- **PivotType** — classify columns as 1x1, 2x2 (Bunch-Kaufman), or delayed
- **Block2x2** — store 2x2 symmetric diagonal blocks
- **MixedDiagonal** — the D factor in LDL^T with mixed block sizes, plus solve and inertia
- **perm_from_forward** — bridge ordering output to faer's `Perm<usize>`
- **Inertia** — relocated from I/O module to its proper domain home

## Usage Examples

### Constructing a MixedDiagonal

```rust
use rivrs_sparse::aptp::{MixedDiagonal, Block2x2, PivotType, Inertia};

// For a 6x6 matrix: columns 0,1 are a 2x2 block; columns 2,3,4,5 are 1x1
let mut diag = MixedDiagonal::new(6);

// Set a 2x2 block at columns 0-1
diag.set_2x2(Block2x2 { first_col: 0, a: 2.0, b: 0.5, c: -3.0 });

// Set 1x1 pivots at columns 2-5
diag.set_1x1(2, 4.0);
diag.set_1x1(3, -1.0);
diag.set_1x1(4, 7.0);
diag.set_1x1(5, 2.0);

// Query
assert_eq!(diag.num_delayed(), 0);
assert_eq!(diag.num_2x2_pairs(), 1);
assert_eq!(diag.num_1x1(), 4);
assert_eq!(diag.get_pivot_type(0), PivotType::TwoByTwo { partner: 1 });
assert_eq!(diag.get_pivot_type(2), PivotType::OneByOne);
```

### Solving D x = b

```rust
// b = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
let mut x = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
diag.solve_in_place(&mut x);
// x now contains the solution to D x = b
```

### Computing Inertia

```rust
let inertia = diag.compute_inertia();
// For the example above: 2x2 block [[2, 0.5], [0.5, -3]] has det < 0 → one +, one -
// 1x1 pivots: 4.0 (+), -1.0 (-), 7.0 (+), 2.0 (+)
// Total: positive=4, negative=2, zero=0
assert_eq!(inertia, Inertia { positive: 4, negative: 2, zero: 0 });
```

### Constructing Permutations from Ordering Output

```rust
use rivrs_sparse::aptp::perm_from_forward;

// AMD ordering produces a forward permutation
let fwd = vec![3, 1, 4, 0, 2];
let perm = perm_from_forward(fwd).expect("valid permutation");

// Use with faer operations
let (fwd_arr, inv_arr) = perm.as_ref().arrays();
assert_eq!(inv_arr[3], 0); // inverse of forward[0]=3 is 0
```

## Build & Test

```bash
cd sparse/
cargo build          # Builds with new aptp module
cargo test           # Runs all tests including new APTP tests
cargo clippy         # Lint check
cargo doc --open     # View documentation
```

## Module Structure

```
src/
├── lib.rs              # adds: pub mod aptp
├── aptp/
│   ├── mod.rs          # re-exports: PivotType, Block2x2, MixedDiagonal, Inertia, perm_from_forward
│   ├── pivot.rs        # PivotType enum, Block2x2 struct
│   ├── diagonal.rs     # MixedDiagonal struct with solve and inertia
│   ├── inertia.rs      # Inertia struct (relocated from io/reference.rs)
│   └── perm.rs         # perm_from_forward function
├── io/
│   └── reference.rs    # adds: pub use crate::aptp::Inertia (re-export)
└── ... (existing modules unchanged)
```
