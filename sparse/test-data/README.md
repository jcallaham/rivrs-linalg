# Test Matrix Collection for SSIDS Solver

Comprehensive set of sparse symmetric matrices for validating the SSIDS
(Sparse Symmetric Indefinite Direct Solver) implementation.

## Directory Layout

```
sparse/test-data/
├── metadata.json                    # Complete index of all 82 matrices
├── hand-constructed/                # 15 small matrices with exact LDL^T factorizations
│   ├── arrow-5-pd.mtx              # Matrix Market format
│   ├── arrow-5-pd.json             # Exact L, D, P factorization
│   └── ...
├── suitesparse-ci/                  # 10 representative matrices (committed to git)
│   ├── easy-indefinite/
│   ├── hard-indefinite/
│   └── positive-definite/
├── suitesparse/                     # Full 67-matrix collection (gitignored, local only)
│   ├── easy-indefinite/             # ~30 well-behaved indefinite matrices
│   ├── hard-indefinite/             # ~18 challenging indefinite matrices
│   └── positive-definite/           # ~19 SPD matrices for fast-path validation
└── scripts/
    ├── requirements.txt
    ├── generate_hand_constructed.py
    ├── download_suitesparse.py
    └── validate_collection.py
```

## Storage Strategy

Matrices are stored in three tiers:

| Tier | Directory | In Git | Size | Used By |
|------|-----------|--------|------|---------|
| Hand-constructed | `hand-constructed/` | Yes (plain git) | 144KB | Unit tests, CI |
| CI subset | `suitesparse-ci/` | Yes (plain git) | ~73MB | CI integration tests |
| Full collection | `suitesparse/` | No (gitignored) | ~4GB | Developer deep testing |

The full SuiteSparse collection is gitignored because it's ~4GB. It's available via:
1. **Docker entrypoint** — automatically extracted from `references/ssids/suitesparse.tar.gz` on container build
2. **Download script** — `python download_suitesparse.py` fetches from SuiteSparse API
3. **Manual archive** — stored at `references/ssids/suitesparse.tar.gz` (outside git)

## First-Time Setup

### In Docker (recommended)

The container entrypoint automatically extracts the full SuiteSparse collection
from `references/ssids/suitesparse.tar.gz`. No manual steps needed.

### Manual setup

```bash
cd sparse/test-data/scripts

# Create and activate virtual environment
python -m venv .venv
source .venv/bin/activate

# Install dependencies
pip install -r requirements.txt

# Generate hand-constructed matrices (no internet needed)
python generate_hand_constructed.py

# Download full SuiteSparse collection (requires internet, slow)
python download_suitesparse.py

# Or download only a subset
python download_suitesparse.py --max-rows 50000
python download_suitesparse.py --category easy-indefinite
python download_suitesparse.py --name "GHS_indef/ncvxqp3"

# Validate the collection
python validate_collection.py
```

## Matrix Categories

| Category | Count | Size Range | Purpose |
|----------|-------|------------|---------|
| Hand-constructed | 15 | 1-20 | Unit tests, exact factorization verification |
| Easy indefinite | 30 | 4K-416K | Basic solver validation |
| Hard indefinite | 18 | 10K-345K | APTP stress testing, delayed pivot validation |
| Positive definite | 19 | 18K-1.6M | Fast-path validation, Cholesky comparison |

## CI Subset (10 matrices)

A representative subset committed to git for CI without SuiteSparse API dependency:

| Matrix | Category | Size | Disk |
|--------|----------|------|------|
| `t2dal` | easy-indefinite | 4.3K | 0.6MB |
| `bloweybq` | easy-indefinite | 10K | 0.4MB |
| `sparsine` | easy-indefinite | 50K | 22.9MB |
| `ncvxqp1` | hard-indefinite | 12.1K | 0.6MB |
| `bratu3d` | hard-indefinite | 27.8K | 1.3MB |
| `cvxqp3` | hard-indefinite | 17.5K | 1.0MB |
| `stokes128` | hard-indefinite (killer) | 49.7K | 8.7MB |
| `ncvxqp3` | hard-indefinite (killer) | 75K | 4.4MB |
| `nd6k` | positive-definite | 18K | 16.6MB |
| `cfd2` | positive-definite | 123K | 16.6MB |

These are identified in `metadata.json` by the `ci_subset: true` flag.

## Adding New Matrices

### Hand-Constructed

Add a generator function to `generate_hand_constructed.py`, then re-run:

```bash
python generate_hand_constructed.py
```

### SuiteSparse

Add the `"Group/Name"` entry to the appropriate curated list in
`download_suitesparse.py`, then re-run the download.
