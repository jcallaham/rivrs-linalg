# External Solver Comparisons

This directory contains tools for benchmarking rivrs-sparse against external
reference solvers. Each solver has a standalone Fortran or C driver that runs
as a subprocess — no FFI linking is required. This keeps the rivrs-sparse crate
free of build-time dependencies on external solver libraries.

## Supported Solvers

| Solver | License | Driver | Status |
|--------|---------|--------|--------|
| SPRAL SSIDS | BSD-3 | `drivers/spral_benchmark.f90` | Working |
| MUMPS | Public domain | `drivers/mumps_benchmark.f90` | Working |
| HSL MA27 | Academic/Commercial | `drivers/ma27_benchmark.f90` | Working (requires user-provided source) |

## Directory Structure

```
comparisons/
├── README.md                          # This file
├── drivers/                           # Fortran/C driver source and build scripts
│   ├── build_spral.sh                 # Builds libspral.a from SPRAL source
│   ├── build_mumps.sh                 # Builds MUMPS driver (requires libmumps-seq-dev)
│   ├── build_ma27.sh                  # Builds MA27 driver (requires user-provided source)
│   ├── build_spral_dense_factor.sh    # Builds dense factor comparison driver
│   ├── spral_benchmark.f90            # SPRAL SSIDS benchmark driver (subprocess)
│   ├── mumps_benchmark.f90            # MUMPS benchmark driver (subprocess)
│   ├── ma27_benchmark.f90             # MA27 benchmark driver (subprocess)
│   ├── spral_full_solve.f90           # SPRAL full solve driver
│   ├── spral_match_order.f90          # SPRAL ordering comparison driver
│   └── spral_dense_factor.cpp         # SPRAL dense APTP comparison
└── src/
    ├── common.rs                      # Shared types, subprocess protocol, formatters
    ├── spral_benchmark.rs             # SPRAL orchestration binary
    ├── mumps_benchmark.rs             # MUMPS orchestration binary
    └── ma27_benchmark.rs              # MA27 orchestration binary
```

## Prerequisites

### SPRAL

1. Ensure SPRAL source is available at `/opt/references/spral/` (or adjust paths
   in the build script).

2. Build the SPRAL library and driver binaries:

```sh
comparisons/drivers/build_spral.sh
```

This produces `/tmp/spral_benchmark` (and other driver binaries).

### MUMPS

MUMPS is available as a system package (sequential version with dummy MPI).

1. Install the library:

```sh
# Debian/Ubuntu
apt-get install libmumps-seq-dev
```

2. Build the MUMPS driver:

```sh
comparisons/drivers/build_mumps.sh
```

This produces `/tmp/mumps_benchmark`.

**Ordering**: MUMPS auto-selects its ordering by default (ICNTL(7)=7). Override
via the `MUMPS_ORDERING` environment variable: `auto`, `metis`, `amd`, `scotch`, `pord`.

### MA27

MA27 source is **not redistributable**. It must be obtained from HSL:

1. Visit https://www.hsl.rl.ac.uk/catalogue/ma27.html
2. Register for a free academic license (or purchase commercial)
3. Download the source package
4. Extract to `/opt/hsl/ma27/` (or set `MA27_SRC` env var)

5. Build the MA27 driver:

```sh
MA27_SRC=/path/to/ma27 comparisons/drivers/build_ma27.sh
```

This produces `/tmp/ma27_benchmark`.

**Ordering note**: MA27 uses its own built-in minimum degree ordering — it does
not accept external orderings (unlike SPRAL/rivrs which use METIS). Timing
comparisons should note this difference.

### METIS

The SPRAL drivers require METIS. The `metis-sys` crate vendored by rivrs-sparse
provides a static `libmetis.a`:

```sh
METIS_LIB=$(find target -name "libmetis.a" | head -1)
```

## Running Comparisons

### SPRAL Benchmark Suite

```sh
# SPRAL-only on CI subset
cargo run --bin spral-comparison --release -- --ci-only

# Side-by-side with rivrs
cargo run --bin spral-comparison --release -- --ci-only --rivrs

# Control thread count (sets OMP_NUM_THREADS for SPRAL, Par::rayon for rivrs)
cargo run --bin spral-comparison --release -- --ci-only --threads 4

# Compare against a previously collected rivrs baseline
cargo run --bin spral-comparison --release -- --ci-only \
  --compare target/benchmarks/baselines/baseline-latest.json

# Filter by category
cargo run --bin spral-comparison --release -- --category hard-indefinite
```

JSON output: `target/benchmarks/spral/`

### MUMPS Benchmark Suite

```sh
# MUMPS-only on CI subset
cargo run --bin mumps-comparison --release -- --ci-only

# Side-by-side with rivrs
cargo run --bin mumps-comparison --release -- --ci-only --rivrs

# Control rivrs thread count (MUMPS sequential is single-threaded)
cargo run --bin mumps-comparison --release -- --ci-only --rivrs --threads 4

# Compare against a previously collected rivrs baseline
cargo run --bin mumps-comparison --release -- --ci-only \
  --compare target/benchmarks/baselines/baseline-latest.json
```

JSON output: `target/benchmarks/mumps/`

### MA27 Benchmark Suite

```sh
# MA27-only on CI subset
cargo run --bin ma27-comparison --release -- --ci-only

# Side-by-side with rivrs
cargo run --bin ma27-comparison --release -- --ci-only --rivrs

# Compare against a previously collected rivrs baseline
cargo run --bin ma27-comparison --release -- --ci-only \
  --compare target/benchmarks/baselines/baseline-latest.json
```

JSON output: `target/benchmarks/ma27/`

### Common CLI Options

All comparison binaries support:

| Flag | Description |
|------|-------------|
| `--ci-only` | Run on CI subset only (10 small matrices) |
| `--rivrs` | Also run rivrs solver for side-by-side comparison |
| `--threads N` | Set thread count for rivrs (and SPRAL via OMP_NUM_THREADS) |
| `--category CAT` | Filter by category: `positive-definite`, `easy-indefinite`, `hard-indefinite` |
| `--compare FILE` | Compare solver timings against a rivrs baseline JSON file |

## Adding a New Solver

To add a comparison against a new solver:

1. Write a Fortran/C driver in `drivers/` that:
   - Reads a matrix from stdin or a file (lower-triangle COO or CSC, 1-indexed)
   - Runs analyze → factor → solve
   - Prints structured results between sentinel markers (see existing drivers)
2. Add a build script in `drivers/`
3. Add a Rust binary in `src/` using `common.rs` shared infrastructure
4. Register the binary in `Cargo.toml`

The subprocess protocol uses sentinel-delimited key-value output, making it easy
to add new solvers without modifying the Rust crate itself.

### Matrix Input Formats

- **CSC** (SPRAL): `n nnz` header, then column pointers (1-indexed), then `row val` pairs
- **COO** (MUMPS, MA27): `n nnz` header, then `row col val` lines (1-indexed, lower triangle)

Use `common::format_spral_input()` or `common::format_lower_coo_text()` from the Rust side.
