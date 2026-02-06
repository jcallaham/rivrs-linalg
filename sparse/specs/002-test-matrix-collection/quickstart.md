# Quickstart: Test Matrix Collection

**Feature**: `002-test-matrix-collection`

## Overview

The test matrix collection provides a comprehensive set of sparse symmetric matrices for validating the SSIDS solver. Matrices are organized by source and difficulty, with complete metadata for filtering and selection.

## Directory Layout

```
sparse/test-data/
├── metadata.json                    # Complete index (all matrices)
├── hand-constructed/                # Small matrices with exact factorizations
│   ├── arrow-5-pd.mtx
│   ├── arrow-5-pd.json             # Exact L, D, P factorization
│   ├── ...
├── suitesparse/
│   ├── easy-indefinite/            # ~30 well-behaved indefinite matrices
│   ├── hard-indefinite/            # ~18 challenging indefinite matrices
│   └── positive-definite/          # ~19 SPD matrices for fast-path validation
└── scripts/
    ├── download_suitesparse.py     # Download curated SuiteSparse matrices
    └── generate_hand_constructed.py # Generate hand-constructed matrices + factorizations
```

## Setup

### Prerequisites

- Python 3.8+ with `ssgetpy` installed: `pip install ssgetpy`
- Git LFS installed and initialized: `git lfs install`

### First-Time Setup

```bash
cd sparse/

# 1. Generate hand-constructed matrices (always available, no download needed)
python test-data/scripts/generate_hand_constructed.py

# 2. Download SuiteSparse matrices (requires internet)
python test-data/scripts/download_suitesparse.py

# 3. Verify collection
python test-data/scripts/download_suitesparse.py --verify-only
```

### Selective Download

```bash
# Download only small matrices (< 50K rows) for quick testing
python test-data/scripts/download_suitesparse.py --max-rows 50000

# Download only easy indefinite set
python test-data/scripts/download_suitesparse.py --category easy-indefinite

# Download a specific matrix by name
python test-data/scripts/download_suitesparse.py --name "GHS_indef/ncvxqp3"
```

## Using the Collection

### Reading metadata.json

```python
import json

with open("test-data/metadata.json") as f:
    index = json.load(f)

# Find all hard indefinite matrices under 50K rows
hard = [m for m in index["matrices"]
        if m["category"] == "hard-indefinite"
        and m["size"] < 50000]

for m in hard:
    print(f"{m['name']}: {m['size']}x{m['size']}, nnz={m['nnz']}")
```

### Finding Matrices by Property

```python
# Matrices with known high delayed pivots (APTP stress tests)
killers = [m for m in index["matrices"]
           if m["properties"].get("killer_case", False)]

# Matrices from a specific paper
paper_matrices = [m for m in index["matrices"]
                  if any("Hogg" in ref for ref in m.get("paper_references", []))]

# All matrices with exact factorizations
exact = [m for m in index["matrices"]
         if m.get("factorization_path")]
```

### Loading in Rust Tests

```rust
// Future test infrastructure (Phase 0.4) will provide:
use sparse_ldlt_test::TestMatrix;

let matrix = TestMatrix::load("hand-constructed/arrow-10-indef")?;
assert_eq!(matrix.size(), 10);
assert!(matrix.is_indefinite());

let factorization = matrix.exact_factorization().unwrap();
// factorization.l, factorization.d, factorization.permutation
```

## Matrix Categories

| Category | Count | Size Range | Purpose |
|----------|-------|------------|---------|
| Hand-constructed | 15+ | 1–20 | Unit tests, exact factorization verification |
| Easy indefinite | ~30 | 4K–1M | Basic solver validation, regression testing |
| Hard indefinite | ~18 | 10K–390K | APTP stress testing, delayed pivot validation |
| Positive definite | ~19 | 18K–1.6M | Fast-path validation, Cholesky comparison |

## Key "Killer" Matrices

These matrices are known to cause failures in solvers with only static pivoting:

| Matrix | Size | Why It's Hard |
|--------|------|---------------|
| GHS_indef/ncvxqp3 | 75K | Nonconvex QP, high delayed pivots |
| GHS_indef/c-71 | 76.6K | Nonlinear optimization |
| ND/nd6k | 18K | 3D mesh, APTP-sensitive |
| ND/nd12k | 36K | 3D mesh, APTP-sensitive |
| GHS_indef/stokes128 | 49.7K | Stokes flow saddle-point |
| Schenk_IBMNA/c-big | 345K | Nonlinear optimization |

## Git LFS

Matrix Market files (`.mtx`) are tracked by Git LFS. After cloning:

```bash
git lfs pull    # Downloads actual matrix files
```

Large matrices (>50MB) are not stored in the repo. Use the download script to fetch them on demand.
