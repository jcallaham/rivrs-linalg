# Contract: SPRAL Driver CLI

**Type**: Command-line interface
**Component**: `spral_driver` (C program)

## Interface

```
spral_driver <matrix.mtx> [options]
```

### Arguments

| Argument | Type | Required | Description |
|----------|------|----------|-------------|
| `<matrix.mtx>` | file path | Yes | Path to Matrix Market file (.mtx) |

### Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `--output <file>` | file path | stdout | Write JSON result to file instead of stdout |
| `--posdef` | boolean flag | false | Assume matrix is positive definite |
| `--ordering <N>` | integer | 1 | Ordering strategy (1=METIS) |
| `--scaling <N>` | integer | 0 | Scaling strategy (0=none) |
| `--u <val>` | float | 0.01 | Pivot tolerance |
| `--nemin <N>` | integer | 32 | Supernode amalgamation parameter |

### Output

Writes a JSON object to stdout (or `--output` file) conforming to the SPRALResult schema defined in `data-model.md`.

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success (result file written, even if factorization had warnings) |
| 1 | Matrix file could not be read |
| 2 | Matrix is not symmetric |
| 3 | SPRAL analysis failed with unrecoverable error |
| 4 | SPRAL factorization failed with unrecoverable error |
| 5 | Memory allocation failure |

### Example

```bash
# Run on a single matrix, output to stdout
./spral_driver test-data/hand-constructed/arrow-5-pd.mtx

# Run on a matrix, write result to file
./spral_driver test-data/suitesparse-ci/bcsstk14/bcsstk14.mtx \
    --output test-data/spral-reference/suitesparse-ci/bcsstk14.json
```

---

# Contract: Run Script

**Type**: Shell script
**Component**: `run-spral.sh`

## Interface

```
run-spral.sh <tier> [--matrix <name>]
```

### Arguments

| Argument | Values | Description |
|----------|--------|-------------|
| `<tier>` | `hand-constructed`, `suitesparse-ci`, `suitesparse` | Which tier of test matrices to process |
| `--matrix <name>` | matrix name | Run only on a single named matrix (optional) |

### Behavior

1. Reads `metadata.json` to find matrices in the specified tier
2. Creates output directory `test-data/spral-reference/<tier>/`
3. Runs `spral_driver` on each matrix
4. Writes `config.json` with SPRAL build configuration
5. Reports progress and any failures to stderr

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | All matrices processed (some may have factorization failures recorded in results) |
| 1 | Driver binary not found or not executable |
| 2 | metadata.json not found |

---

# Contract: Summary Generator

**Type**: Shell script
**Component**: `generate-summary.sh`

## Interface

```
generate-summary.sh [--output <file>]
```

### Behavior

1. Reads all JSON result files from `test-data/spral-reference/*/`
2. Combines with matrix metadata from `metadata.json`
3. Classifies each matrix by difficulty
4. Writes summary JSON to `test-data/spral-reference/summary.json` (or `--output` path)

### Difficulty Classification

| Category | Criteria |
|----------|----------|
| `trivial` | n < 100 AND backward_error < 1e-14 |
| `easy` | backward_error < 1e-12 AND num_delay == 0 |
| `moderate` | backward_error < 1e-10 OR num_delay > 0 |
| `hard` | backward_error >= 1e-10 OR num_two > 0 |
| `failed` | status != success |

---

# Contract: Reproducibility Verifier

**Type**: Shell script
**Component**: `verify-reproducibility.sh`

## Interface

```
verify-reproducibility.sh <tier>
```

### Behavior

1. Runs `spral_driver` on all matrices in the specified tier
2. Saves results to a temporary directory
3. Runs again and saves to a second temporary directory
4. Compares structural fields for exact match
5. Compares numerical fields within 1e-14 relative tolerance
6. Reports pass/fail for each matrix

### Exit Codes

| Code | Meaning |
|------|---------|
| 0 | All matrices reproducible |
| 1 | One or more matrices failed reproducibility check |
