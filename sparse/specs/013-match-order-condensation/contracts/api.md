# API Contract: Match-Order Condensation Pipeline

**Feature**: 013-match-order-condensation
**Date**: 2026-02-13

## Public API

### `match_order_metis()`

Combined matching-ordering pipeline that guarantees matched pair adjacency.

```rust
/// Compute a combined MC64 matching + METIS ordering with cycle condensation.
///
/// This implements SPRAL's `match_order_metis` pipeline:
/// 1. MC64 matching → scaling + matching permutation
/// 2. Cycle splitting → decompose matching into singletons and 2-cycles
/// 3. Condensed graph → fuse 2-cycle pairs into single super-nodes
/// 4. METIS ordering → fill-reducing ordering on condensed graph
/// 5. Expansion → map back to original indices with pair adjacency
///
/// Matched 2-cycle pairs are guaranteed to occupy consecutive positions
/// in the output ordering, making them natural candidates for 2x2 pivots
/// in APTP factorization.
///
/// # Arguments
///
/// * `matrix` — Sparse symmetric matrix in upper-triangular CSC format.
///   Must be square. Numeric values required (used by MC64).
///
/// # Returns
///
/// * `Ok(MatchOrderResult)` — Ordering, scaling, and diagnostics.
/// * `Err(SparseError::NotSquare)` — Matrix is not square.
/// * `Err(SparseError::InvalidInput)` — Zero dimension or non-finite entries.
/// * `Err(SparseError::AnalysisFailure)` — METIS or MC64 internal failure.
///
/// # Algorithm References
///
/// - Hogg & Scott (2013), HSL_MC80 — cycle condensation approach
/// - Duff & Koster (2001), Algorithm MPD — MC64 matching
/// - Duff & Pralet (2005) — MC64SYM symmetric scaling
/// - SPRAL `match_order.f90` (BSD-3) — reference implementation
///
/// # SPRAL Equivalent
///
/// Corresponds to `spral_match_order::match_order_metis` with
/// `options%ordering = 2` in SSIDS.
pub fn match_order_metis(
    matrix: &SparseColMat<usize, f64>,
) -> Result<MatchOrderResult, SparseError>;
```

**Preconditions**:
- Matrix must be square
- Matrix must be in upper-triangular CSC format (faer convention)
- All entries must be finite (no NaN/Inf)
- Dimension must fit in i32 (METIS constraint)

**Postconditions**:
- `result.ordering` is a valid permutation of {0..n-1}
- Every 2-cycle pair (i,j) in the MC64 matching occupies consecutive positions in `result.ordering`
- Unmatched indices (if any) occupy the last `n - result.matched` positions
- `result.scaling` has length n with all positive entries
- `result.condensed_dim <= n`; strictly `< n` when 2-cycles exist

**Error behavior**: Returns `Err` on invalid input; never panics in release builds.

### `MatchOrderResult`

```rust
/// Result of the combined matching-ordering pipeline.
///
/// Contains both the fill-reducing ordering (for symbolic analysis) and
/// the MC64 scaling factors (for numeric factorization), plus diagnostics
/// about the condensation process.
pub struct MatchOrderResult {
    /// Fill-reducing permutation with matched pair adjacency guarantee.
    /// Use with `SymmetricOrdering::Custom(result.ordering.as_ref())`.
    pub ordering: Perm<usize>,

    /// MC64 symmetric scaling factors (linear domain).
    /// Apply as `A_scaled[i,j] = scaling[i] * A[i,j] * scaling[j]`.
    pub scaling: Vec<f64>,

    /// Number of matched entries from MC64.
    /// Equals n for structurally nonsingular matrices.
    pub matched: usize,

    /// Dimension of the condensed graph passed to METIS.
    /// Strictly less than n when 2-cycles exist.
    pub condensed_dim: usize,

    /// Number of singleton nodes (self-matched).
    pub singletons: usize,

    /// Number of 2-cycle pairs in the decomposed matching.
    pub two_cycles: usize,
}
```

## Usage Patterns

### Basic usage (recommended for indefinite matrices)

```rust
use rivrs_sparse::aptp::{match_order_metis, AptpSymbolic};
use faer::sparse::linalg::cholesky::SymmetricOrdering;

let matrix = load_matrix("indefinite_problem.mtx")?;
let result = match_order_metis(&matrix)?;

// Use ordering for symbolic analysis
let symbolic = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Custom(result.ordering.as_ref()),
)?;

// Store scaling for Phase 5 numeric factorization
let scaling = result.scaling;
```

### Comparison with separate MC64 + METIS

```rust
// Option A: Condensed pipeline (pair adjacency guaranteed)
let condensed = match_order_metis(&matrix)?;

// Option B: Separate components (no pair adjacency guarantee)
let mc64 = mc64_matching(&matrix, Mc64Job::MaximumProduct)?;
let metis_perm = metis_ordering(matrix.symbolic())?;
// Naive composition: mc64.matching * metis_perm
// WARNING: matched pairs may be scattered in the ordering
```

### Inspecting diagnostics

```rust
let result = match_order_metis(&matrix)?;
println!("Original dim: {}", matrix.nrows());
println!("Condensed dim: {} ({:.0}% reduction)",
    result.condensed_dim,
    (1.0 - result.condensed_dim as f64 / matrix.nrows() as f64) * 100.0);
println!("Cycle breakdown: {} singletons, {} 2-cycles",
    result.singletons, result.two_cycles);
```

## Internal API (not public)

These functions are implementation details within `ordering.rs`:

```rust
// Decompose matching into singletons, 2-cycles, unmatched.
// The is_matched slice (from Mc64Result.is_matched) distinguishes
// true singletons (matched, fwd[i]==i) from unmatched indices
// (not matched, assigned to arbitrary free columns by build_singular_permutation).
fn split_matching_cycles(
    matching_fwd: &[usize],
    is_matched: &[bool],
    n: usize,
) -> CycleDecomposition;

// Build condensed CSR adjacency (xadj, adjncy) for METIS
fn build_condensed_adjacency(
    matrix: SymbolicSparseColMatRef<'_, usize>,
    decomp: &CycleDecomposition,
) -> Result<(Vec<i32>, Vec<i32>), SparseError>;

// Expand condensed METIS ordering to full-size permutation.
// Computes inverse of condensed_order internally.
fn expand_ordering(
    condensed_order: &[i32],
    decomp: &CycleDecomposition,
    n: usize,
) -> Perm<usize>;
```
