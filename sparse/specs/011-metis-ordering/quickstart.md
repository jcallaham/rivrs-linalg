# Quickstart: METIS Nested Dissection Ordering

**Feature**: 011-metis-ordering

## Usage

```rust
use rivrs_sparse::aptp::{AptpSymbolic, metis_ordering};
use rivrs_sparse::io::mtx::load_matrix_market;
use faer::sparse::linalg::cholesky::SymmetricOrdering;

// Load a symmetric sparse matrix
let matrix = load_matrix_market("test-data/suitesparse-ci/bmwcra_1.mtx")?;

// Compute METIS fill-reducing ordering
let perm = metis_ordering(matrix.symbolic())?;

// Use the ordering in symbolic analysis
let symbolic = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Custom(perm.as_ref()),
)?;

// Inspect results
println!("Predicted nnz(L): {}", symbolic.predicted_nnz());
println!("Statistics: {}", symbolic.statistics());
```

## Comparing METIS vs AMD

```rust
use faer::sparse::linalg::cholesky::SymmetricOrdering;

// AMD ordering (built into faer)
let sym_amd = AptpSymbolic::analyze(matrix.symbolic(), SymmetricOrdering::Amd)?;

// METIS ordering
let perm = metis_ordering(matrix.symbolic())?;
let sym_metis = AptpSymbolic::analyze(
    matrix.symbolic(),
    SymmetricOrdering::Custom(perm.as_ref()),
)?;

println!("AMD  predicted nnz(L): {}", sym_amd.predicted_nnz());
println!("METIS predicted nnz(L): {}", sym_metis.predicted_nnz());
// METIS typically produces 2-10x less fill for matrices with geometric structure
```

## Build Requirements

- Rust 1.87+ (edition 2024)
- C compiler (gcc, clang, or MSVC) — needed to compile vendored METIS source
- No system library installation required
