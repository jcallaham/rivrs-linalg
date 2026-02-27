# Data Model: Robustness — Testing & Hardening

**Feature**: 026-robustness-hardening
**Date**: 2026-02-25

## Entities

### Perturbation Helpers

Functions that transform a well-formed dense symmetric matrix into an adversarial variant. All operate on `Mat<f64>` in-place.

| Helper | Inputs | Effect | Used By |
|--------|--------|--------|---------|
| `cause_delays` | matrix, rng | Multiplies n/8 random rows and n/8 random entries by 1000 | Torture tests (70% probability) |
| `make_singular` | matrix, col1, col2 | Makes col2 a scaled copy of col1 (creates rank deficiency) | Torture tests (20% probability) |
| `make_dblk_singular` | matrix, block_row, block_size | Makes a specific diagonal block singular | Torture tests (10% probability, APP only) |

**Location**: `src/testing/perturbations.rs` (new module, behind `test-util` feature)

### Proptest Strategies

Custom proptest strategies for generating test inputs with automatic shrinking.

| Strategy | Output Type | Parameters | Used By |
|----------|------------|------------|---------|
| `arb_symmetric_pd` | `Mat<f64>` | size range (5–500) | Property tests |
| `arb_symmetric_indefinite` | `Mat<f64>` | size range (5–500) | Property tests |
| `arb_sparse_symmetric` | `SparseColMat<usize, f64>` | size range (5–500), density range | End-to-end property tests |
| `arb_aptp_options` | `AptpOptions` | threshold range, block sizes | Property tests |

**Location**: `src/testing/strategies.rs` (new module, behind `test-util` feature)

### Torture Test Configuration

Parameters controlling probabilistic torture test generation.

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `num_instances` | usize | 500 | Random instances per configuration |
| `delay_probability` | f64 | 0.70 | Chance of applying `cause_delays` |
| `singular_probability` | f64 | 0.20 | Chance of applying `make_singular` |
| `dblk_singular_probability` | f64 | 0.10 | Chance of applying `make_dblk_singular` |
| `matrix_sizes` | Vec<(usize, usize)> | [(32,32), (64,64), (128,128), (128,48)] | (m, n) pairs to test |
| `backward_error_threshold` | f64 | 5e-11 | Max backward error for non-singular |
| `seed` | Option<u64> | None | Fixed seed for reproducibility |

### Audit Mapping Document

A markdown document recording the two-directional audit.

**Structure**:
```
dev/spral-test-audit.md
├── SPRAL → rivrs mapping table
│   ├── ldlt_app.cxx scenarios (48 deterministic + 2 torture)
│   └── ldlt_tpp.cxx scenarios (23 deterministic + 2 torture)
├── rivrs test evaluation
│   ├── Retained tests (with rationale)
│   └── Removed tests (with rationale)
└── Gap analysis summary
```

## State Transitions

No runtime state transitions — this feature adds test infrastructure only.

## Relationships

```
perturbations.rs ──uses──> faer::Mat<f64>
strategies.rs ──uses──> generators.rs (existing)
strategies.rs ──uses──> perturbations.rs
torture_tests ──uses──> perturbations.rs
torture_tests ──uses──> aptp_factor_in_place / tpp_factor_as_primary
property_tests ──uses──> strategies.rs
property_tests ──uses──> SparseLDLT (end-to-end)
adversarial_tests ──uses──> SparseLDLT
```
