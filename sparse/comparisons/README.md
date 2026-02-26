# External Solver Comparisons

This directory contains tools for benchmarking rivrs-sparse against external
reference solvers. Each solver has a standalone Fortran or C driver that runs
as a subprocess — no FFI linking is required. This keeps the rivrs-sparse crate
free of build-time dependencies on external solver libraries.

## Supported Solvers

| Solver | License | Driver | Status |
|--------|---------|--------|--------|
| SPRAL SSIDS | BSD-3 | `drivers/spral_benchmark.f90` | Working |
| MUMPS | Public domain | — | Planned |
| HSL MA27 | Academic/Commercial | — | Planned |

## Directory Structure

```
comparisons/
├── README.md                          # This file
├── drivers/                           # Fortran/C driver source and build scripts
│   ├── build_spral.sh                 # Builds libspral.a from SPRAL source
│   ├── build_spral_dense_factor.sh    # Builds dense factor comparison driver
│   ├── spral_benchmark.f90            # SPRAL SSIDS benchmark driver (subprocess)
│   ├── spral_full_solve.f90           # SPRAL full solve driver
│   ├── spral_match_order.f90          # SPRAL ordering comparison driver
│   └── spral_dense_factor.cpp         # SPRAL dense APTP comparison
└── src/
    └── spral_benchmark.rs             # Rust orchestration binary
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

### METIS

The SPRAL drivers require METIS. The `metis-sys` crate vendored by rivrs-sparse
provides a static `libmetis.a`:

```sh
METIS_LIB=$(find target -name "libmetis.a" | head -1)
```

## Running Comparisons

### SPRAL Benchmark Suite

The Rust orchestration binary (`src/spral_benchmark.rs`) loads matrices from the
SuiteSparse test registry, runs the SPRAL subprocess, and produces side-by-side
comparison tables.

```sh
# Build the comparison binary
cargo build --bin spral-comparison --release

# SPRAL-only on CI subset
cargo run --bin spral-comparison --release -- --ci-only

# Side-by-side with rivrs
cargo run --bin spral-comparison --release -- --ci-only --rivrs

# Control SPRAL thread count (sets OMP_NUM_THREADS for the subprocess)
cargo run --bin spral-comparison --release -- --ci-only --threads 4

# Compare against a previously collected rivrs baseline
cargo run --bin spral-comparison --release -- --ci-only \
  --compare target/benchmarks/baselines/baseline-latest.json

# Filter by category
cargo run --bin spral-comparison --release -- --category hard-indefinite
```

JSON output is written to `target/benchmarks/spral/`.

## Adding a New Solver

To add a comparison against a new solver (e.g., MUMPS):

1. Write a Fortran/C driver in `drivers/` that:
   - Reads a matrix from stdin or a file (lower-triangle CSC, 1-indexed)
   - Runs analyze → factor → solve
   - Prints structured results between sentinel markers (see `spral_benchmark.f90`)
2. Add a build script in `drivers/`
3. Add a Rust module in `src/` with subprocess invocation and result parsing
4. Register the solver in the orchestration binary

The subprocess protocol uses sentinel-delimited key-value output, making it easy
to add new solvers without modifying the Rust crate itself.
