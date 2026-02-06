# Test Matrix Collection for SSIDS Solver

Comprehensive set of sparse symmetric matrices for validating the SSIDS
(Sparse Symmetric Indefinite Direct Solver) implementation.

## Prerequisites

- Python 3.8+ with pip
- Git LFS (for matrix file tracking)
- uv (recommended for virtual environment management)

## First-Time Setup

```bash
cd sparse/test-data/scripts

# Create and activate virtual environment
uv venv .venv
source .venv/bin/activate  # or: .venv/bin/python for direct invocation

# Install dependencies
uv pip install -r requirements.txt

# Generate hand-constructed matrices (no internet needed)
python generate_hand_constructed.py

# Download SuiteSparse matrices (requires internet)
python download_suitesparse.py

# Validate the collection
python validate_collection.py
```

## Directory Layout

```
sparse/test-data/
├── metadata.json                    # Complete index of all matrices
├── hand-constructed/                # 15+ small matrices with exact LDL^T factorizations
│   ├── arrow-5-pd.mtx              # Matrix Market format
│   ├── arrow-5-pd.json             # Exact L, D, P factorization
│   └── ...
├── suitesparse/
│   ├── easy-indefinite/            # ~30 well-behaved indefinite matrices
│   ├── hard-indefinite/            # ~18 challenging indefinite matrices
│   └── positive-definite/          # ~19 SPD matrices for fast-path validation
└── scripts/
    ├── requirements.txt
    ├── generate_hand_constructed.py
    ├── download_suitesparse.py
    └── validate_collection.py
```

## Selective Download

```bash
# Download only small matrices (< 50K rows) for quick testing
python download_suitesparse.py --max-rows 50000

# Download only one category
python download_suitesparse.py --category easy-indefinite

# Download a specific matrix
python download_suitesparse.py --name "GHS_indef/ncvxqp3"

# Dry-run (list what would be downloaded)
python download_suitesparse.py --dry-run

# Verify existing downloads against curated lists
python download_suitesparse.py --verify-only
```

## Matrix Categories

| Category | Count | Size Range | Purpose |
|----------|-------|------------|---------|
| Hand-constructed | 15+ | 1-20 | Unit tests, exact factorization verification |
| Easy indefinite | ~30 | 4K-416K | Basic solver validation |
| Hard indefinite | ~18 | 10K-345K | APTP stress testing, delayed pivot validation |
| Positive definite | ~19 | 18K-1.6M | Fast-path validation, Cholesky comparison |

## Adding New Matrices

### Hand-Constructed

Add a generator function to `generate_hand_constructed.py`, then re-run:

```bash
python generate_hand_constructed.py
```

The script will generate the `.mtx` and `.json` files, verify the factorization,
and update `metadata.json`.

### SuiteSparse

Add the `"Group/Name"` entry to the appropriate curated list in
`download_suitesparse.py`, then re-run the download.

## Git LFS

Matrix Market files (`.mtx`) are tracked by Git LFS via the root `.gitattributes`:

```
sparse/test-data/**/*.mtx filter=lfs diff=lfs merge=lfs -text
```

After cloning, run `git lfs pull` to download actual matrix files. Large matrices
(>50MB) are not stored in the repo and must be downloaded on demand.
