# Quickstart: Core Test Infrastructure

**Branch**: `005-test-infrastructure` | **Date**: 2026-02-06

## Prerequisites

- Rust 1.87+ (edition 2024)
- Existing `rivrs-sparse` crate with Phase 0.4 complete
- Test data directory with metadata.json and matrix files

## Setup

The test infrastructure is gated behind the `test-util` Cargo feature flag:

```toml
# In Cargo.toml
[features]
test-util = []

[dev-dependencies]
rivrs-sparse = { path = ".", features = ["test-util"] }
```

## Usage Examples

### 1. Validate a hand-constructed matrix using NumericalValidator

```rust
use rivrs_sparse::testing::{NumericalValidator, load_test_cases, TestCaseFilter};

#[test]
fn validate_arrow_5_pd() {
    let cases = load_test_cases(&TestCaseFilter::hand_constructed())
        .expect("failed to load test cases");

    let validator = NumericalValidator::new(); // default: recon < 1e-12, berr < 1e-10

    for case in &cases {
        if let Some(ref reference) = case.reference {
            let result = validator.check_reconstruction(&case.matrix, reference);
            assert!(
                result.passed,
                "{}: reconstruction error {:.2e} exceeds {:.2e}",
                case.name, result.value, result.threshold
            );
        }
    }
}
```

### 2. Test a solver implementation with SolverTest trait

```rust
use rivrs_sparse::testing::{SolverTest, SolverTestCase, TestResult, load_test_cases, TestCaseFilter};

struct MySolver;

impl SolverTest for MySolver {
    fn test_analyze(&self, case: &SolverTestCase) -> TestResult {
        // Call your analyze implementation, validate results
        todo!()
    }

    fn test_factor(&self, case: &SolverTestCase) -> TestResult {
        // Call your factorize implementation, validate results
        todo!()
    }

    fn test_solve(&self, case: &SolverTestCase) -> TestResult {
        // Call your solve implementation, validate results
        todo!()
    }

    fn test_roundtrip(&self, case: &SolverTestCase) -> TestResult {
        // Run full pipeline, validate final result
        todo!()
    }
}

#[test]
fn test_my_solver_on_hand_constructed() {
    let solver = MySolver;
    let cases = load_test_cases(&TestCaseFilter::hand_constructed().require_reference())
        .expect("failed to load");

    for case in &cases {
        let result = solver.test_roundtrip(&case);
        assert!(result.passed, "{}: {}", case.name, result.diagnostics.join("; "));
    }
}
```

### 3. Generate random test matrices

```rust
use rivrs_sparse::testing::{RandomMatrixConfig, generate_random_symmetric};
use rand::thread_rng;

#[test]
fn solver_handles_random_pd_matrices() {
    let mut rng = thread_rng();
    let config = RandomMatrixConfig {
        size: 100,
        target_nnz: 500,
        positive_definite: true,
    };

    let matrix = generate_random_symmetric(&config, &mut rng)
        .expect("generation failed");

    assert_eq!(matrix.nrows(), 100);
    assert_eq!(matrix.ncols(), 100);
    // matrix is guaranteed PD by construction (diagonal dominance)
}
```

### 4. Generate structured test matrices

```rust
use rivrs_sparse::testing::{generate_arrow, generate_tridiagonal, generate_banded};
use rand::thread_rng;

#[test]
fn test_arrow_matrices() {
    let mut rng = thread_rng();
    let arrow = generate_arrow(50, true, &mut rng);  // 50x50 PD arrow
    assert_eq!(arrow.nrows(), 50);
    // Arrow has nnz = 2*(n-1) + 1 = 99 (for symmetric full storage: 3*(n-1) + 1)
}
```

### 5. Filter test cases by category

```rust
use rivrs_sparse::testing::{load_test_cases, TestCaseFilter};

#[test]
fn test_hard_indefinite_matrices() {
    let filter = TestCaseFilter::all()
        .with_category("hard-indefinite")
        .ci_only();

    let cases = load_test_cases(&filter).expect("failed to load");
    // Only hard-indefinite matrices in CI subset
    for case in &cases {
        assert!(case.properties.indefinite);
    }
}
```

### 6. Custom tolerances for ill-conditioned matrices

```rust
use rivrs_sparse::testing::NumericalValidator;

let relaxed = NumericalValidator::new()
    .with_reconstruction_tol(1e-6)
    .with_backward_error_tol(1e-4);

// Use relaxed tolerances for ill-conditioned test cases
```

## Running Tests

```bash
# Run all tests (includes test infrastructure via dev-dependencies)
cargo test

# Run only tests that use the harness
cargo test --test hand_constructed
cargo test --test suitesparse_ci

# Run benchmarks
cargo bench
```

## Key Conventions

1. **Default tolerances** match the constitution: reconstruction < 1e-12, backward error < 1e-10
2. **Missing matrices** are silently skipped (gitignored SuiteSparse files)
3. **TestResult.diagnostics** always includes matrix name and dimensions for debugging
4. **Generators** use `&mut impl Rng` for reproducibility with seeded RNGs
