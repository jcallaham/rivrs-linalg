# Quickstart: Multifrontal Numeric Factorization

**Feature**: 015-multifrontal-factorization
**Module**: `src/aptp/numeric.rs`

## Build

```bash
cd /workspace/rivrs-linalg/sparse
cargo build
```

No new dependencies required. Uses existing faer 0.22 + Phase 2/3/5 modules.

## Test

```bash
# Unit tests (fast — hand-constructed matrices)
cargo test aptp::numeric

# Integration tests (CI subset — 10 SuiteSparse matrices)
cargo test --test suitesparse_ci

# Full SuiteSparse test set (slow — 67 matrices, requires test-data archive)
cargo test -- --ignored --test-threads=1

# Run with reconstruction error details
cargo test aptp::numeric -- --nocapture
```

## Usage

```rust
use faer::sparse::SparseColMat;
use rivrs_sparse::aptp::{
    AptpSymbolic, AptpNumeric, AptpOptions, AptpFallback,
};
use faer::sparse::linalg::cholesky::SymmetricOrdering;

// 1. Load or construct a sparse symmetric matrix
let matrix: SparseColMat<usize, f64> = /* ... */;

// 2. Symbolic analysis (Phase 3)
let symbolic = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Amd,
)?;

// 3. Numeric factorization (Phase 6 — this feature)
let options = AptpOptions {
    threshold: 0.01,
    small: 1e-20,
    fallback: AptpFallback::BunchKaufman,
};
let numeric = AptpNumeric::factor(&symbolic, &matrix, &options)?;

// 4. Inspect factorization statistics
let stats = numeric.stats();
println!("Total pivots: {} 1x1, {} 2x2", stats.total_1x1_pivots, stats.total_2x2_pivots);
println!("Delayed: {}", stats.total_delayed);
println!("Max front size: {}", stats.max_front_size);
```

## File Layout

```
src/aptp/
├── mod.rs          # Add: pub mod numeric; re-export AptpNumeric, FrontFactors, FactorizationStats
├── numeric.rs      # NEW: Multifrontal factorization (primary implementation file)
├── symbolic.rs     # Existing: AptpSymbolic (Phase 3)
├── factor.rs       # Existing: aptp_factor_in_place (Phase 5)
├── diagonal.rs     # Existing: MixedDiagonal (Phase 2)
├── pivot.rs        # Existing: PivotType, Block2x2 (Phase 2)
├── inertia.rs      # Existing: Inertia (Phase 2)
├── perm.rs         # Existing: perm_from_forward (Phase 2)
├── ordering.rs     # Existing: METIS ordering (Phase 4.1)
└── matching.rs     # Existing: MC64 matching (Phase 4.2)

tests/
├── hand_constructed.rs   # Existing: update with multifrontal reconstruction tests
├── suitesparse_ci.rs     # Existing: update with multifrontal factorization tests
└── multifrontal.rs       # NEW: dedicated multifrontal unit/integration tests
```

## Key Design Decisions

1. **Pass entire frontal matrix to Phase 5 kernel**: `aptp_factor_in_place(frontal, k, options)` where k = num_fully_summed. The kernel's Schur complement updates propagate to all trailing rows (L21 and F22 are computed implicitly).

2. **Unified supernode abstraction**: Both supernodal and simplicial decompositions are mapped to `SupernodeInfo` structs, so the factorization loop has a single code path.

3. **Postorder is trivial**: faer guarantees parent supernode index > child index. Iterate `0..n_supernodes` for bottom-up traversal.

4. **Contribution block includes delayed columns**: The trailing `(m - ne) x (m - ne)` submatrix is passed to the parent as-is. Delayed columns are the first `(k - ne)` rows/columns.
