# Quickstart: Triangular Solve & Solver API

**Feature**: 016-triangular-solve-api
**Date**: 2026-02-16

## One-Shot Solve

The simplest way to solve Ax = b:

```rust
use faer::sparse::{SparseColMat, Triplet};
use faer::Col;
use rivrs_sparse::aptp::{SparseLDLT, SolverOptions};

// Build a sparse symmetric matrix (lower triangle)
let triplets = vec![
    Triplet::new(0, 0, 4.0),
    Triplet::new(1, 0, 1.0), Triplet::new(1, 1, 3.0),
    Triplet::new(2, 0, 0.5), Triplet::new(2, 2, -2.0),
];
let a = SparseColMat::try_new_from_triplets(3, 3, &triplets).unwrap();

// Right-hand side
let b = Col::from_fn(3, |i| [5.5, 4.0, -1.5][i]);

// Solve in one call
let x = SparseLDLT::solve_full(&a, &b, &SolverOptions::default()).unwrap();

// Verify: backward error should be < 1e-10
```

## Three-Phase API (Reusable Analysis)

For solving multiple systems with the same sparsity pattern:

```rust
use faer::dyn_stack::{MemBuffer, MemStack};
use rivrs_sparse::aptp::{
    SparseLDLT, AnalyzeOptions, FactorOptions, OrderingStrategy,
};

// 1. Analyze (once per sparsity pattern)
let options = AnalyzeOptions {
    ordering: OrderingStrategy::Metis,
};
let mut solver = SparseLDLT::analyze(a.symbolic(), &options).unwrap();

// 2. Factor
let factor_opts = FactorOptions::default();
solver.factor(&a, &factor_opts).unwrap();

// Check factorization diagnostics
if let Some(stats) = solver.stats() {
    println!("Pivots: {} 1x1, {} 2x2, {} zero",
        stats.total_1x1_pivots, stats.total_2x2_pivots, stats.zero_pivots);
}
if let Some(inertia) = solver.inertia() {
    println!("Inertia: +{} -{} 0{}", inertia.positive, inertia.negative, inertia.zero);
}

// 3. Solve (reuse factorization for multiple RHS)
let scratch = solver.solve_scratch(1);
let mut mem = MemBuffer::new(scratch);
let mut stack = MemStack::new(&mut mem);

let x1 = solver.solve(&b1, &mut stack).unwrap();
let x2 = solver.solve(&b2, &mut stack).unwrap();
```

## Refactoring (Same Pattern, Different Values)

When the sparsity pattern is unchanged but numeric values change (e.g., interior point iteration):

```rust
// First system
solver.factor(&a1, &factor_opts).unwrap();
let x1 = solver.solve(&b1, &mut stack).unwrap();

// Second system (same sparsity, different values)
solver.refactor(&a2, &factor_opts).unwrap();
let x2 = solver.solve(&b2, &mut stack).unwrap();
```

## MC64 Scaling (Hard Indefinite Problems)

For hard indefinite matrices, MC64 matching + scaling improves pivot quality:

```rust
let options = AnalyzeOptions {
    ordering: OrderingStrategy::MatchOrderMetis,
};
let mut solver = SparseLDLT::analyze(a.symbolic(), &options).unwrap();

// Scaling is applied automatically during factor() and solve()
solver.factor(&a, &FactorOptions::default()).unwrap();
let x = solver.solve(&b, &mut stack).unwrap();
// x is in the original (unscaled) coordinate system
```

## Detecting Rank Deficiency

Zero pivots indicate a rank-deficient matrix. The solver produces a solution with zeroed components at zero-pivot positions (SPRAL's default behavior):

```rust
solver.factor(&singular_matrix, &factor_opts).unwrap();

if let Some(stats) = solver.stats() {
    if stats.zero_pivots > 0 {
        println!("Warning: matrix is rank-deficient ({} zero pivots)", stats.zero_pivots);
        // Solution will have zeros at rank-deficient positions
    }
}

// Solve still succeeds — check solution quality via backward error
let x = solver.solve(&b, &mut stack).unwrap();
```

## Error Handling

```rust
// Solve before factor → clear error
let analyzer_only = SparseLDLT::analyze(a.symbolic(), &options).unwrap();
let err = analyzer_only.solve(&b, &mut stack);
assert!(matches!(err, Err(SparseError::SolveBeforeFactor { .. })));

// Dimension mismatch → clear error
let wrong_b = Col::from_fn(999, |i| i as f64);
let err = solver.solve(&wrong_b, &mut stack);
assert!(matches!(err, Err(SparseError::DimensionMismatch { .. })));
```
