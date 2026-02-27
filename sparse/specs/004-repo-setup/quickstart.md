# Quickstart: 004-repo-setup

**Date**: 2026-02-06

## Prerequisites

- Rust 1.87+ (edition 2024)
- The repository cloned with `test-data/hand-constructed/` and `test-data/suitesparse-ci/` present (standard clone includes these)

## Build & Test

```bash
cd sparse/

# Build
cargo build

# Run all tests (hand-constructed + CI-subset)
cargo test

# Run just the hand-constructed matrix validation tests
cargo test --test hand_constructed

# Run just the CI-subset SuiteSparse tests
cargo test --test suitesparse_ci

# Run benchmarks
cargo bench

# Lint checks
cargo clippy --all-targets -- -D warnings
cargo fmt --check

# Documentation build
cargo doc --no-deps
```

## Using the Test Infrastructure in New Tests

```rust
use rivrs_sparse::io::registry::load_test_matrix;
use rivrs_sparse::validate;

#[test]
fn test_my_solver_on_arrow_matrix() {
    let test = load_test_matrix("arrow-5-pd")
        .expect("registry error")
        .expect("matrix not on disk");

    // Use the sparse matrix
    let a = test.matrix.as_ref();
    assert_eq!(a.nrows(), 5);
    assert_eq!(a.ncols(), 5);

    // Check reference factorization if available
    if let Some(ref reference) = test.reference {
        let error = validate::reconstruction_error(a, reference);
        assert!(error < 1e-12, "reconstruction error: {error:.2e}");
    }
}
```

## Module Layout After Implementation

```
src/
├── lib.rs          # Crate root, re-exports
├── error.rs        # SparseError (extended with IO variants)
├── io.rs           # pub mod mtx, reference, registry
├── io/
│   ├── mtx.rs      # Matrix Market parser
│   ├── reference.rs # JSON factorization loader
│   └── registry.rs  # Test matrix registry
└── validate.rs     # Reconstruction error, backward error, inertia check

tests/
├── common/
│   └── mod.rs      # Shared test utilities
├── hand_constructed.rs  # Integration tests for 15 hand-constructed matrices
└── suitesparse_ci.rs    # Integration tests for 10 CI-subset matrices

benches/
└── matrix_loading.rs    # Criterion benchmarks
```

## New Dependencies (dev-only)

```toml
[dependencies]
serde = { version = "1", features = ["derive"] }
serde_json = "1"

[dev-dependencies]
# (existing: approx, criterion, rand, rand_distr)
```

Note: `serde` and `serde_json` are regular dependencies (not dev-only) because the IO and validation modules are part of the library crate (used by both tests and future solver code).
