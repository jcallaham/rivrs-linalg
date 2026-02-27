# Quickstart: SPRAL Golden Results Generation

**Feature**: 003-spral-golden-results
**Date**: 2026-02-06

## Prerequisites

1. Development container running (Debian 12 / arm64)
2. SPRAL source at `/opt/references/spral/` (via references symlink)
3. Test matrices from Phase 0.2 at `test-data/` (82 matrices)
4. `metadata.json` listing all matrices with paths and properties

## Step 1: Install Build Dependencies

```bash
# In the development container
apt-get update && apt-get install -y \
    python3 python3-pip ninja-build \
    libopenblas-dev liblapack-dev \
    libmetis-dev libhwloc-dev

pip3 install meson
```

## Step 2: Build SPRAL (CPU-only)

```bash
cd /opt/references/spral

meson setup builddir --prefix=/usr/local \
    --buildtype=release \
    -Dgpu=false \
    -Dopenmp=true \
    -Dexamples=true \
    -Dtests=true

cd builddir
meson compile
meson install  # installs libspral and headers to /usr/local
ldconfig       # update shared library cache
```

Verify: `meson test` should pass SSIDS-related tests.

## Step 3: Build the Driver Program

```bash
cd /workspace/rivrs-linalg/sparse/test-data/scripts/

# Compile the C driver (links against installed SPRAL)
gcc -O2 -o spral-driver spral_driver.c mmio.c \
    -I/usr/local/include \
    -L/usr/local/lib -lspral \
    -lopenblas -lmetis -lhwloc \
    -lgfortran -lm -lgomp -lpthread

# Verify it runs on a small matrix
export OMP_CANCELLATION=TRUE
export OMP_PROC_BIND=TRUE
./spral-driver ../../test-data/hand-constructed/arrow-5-pd.mtx
```

## Step 4: Generate Golden Results

```bash
cd /workspace/rivrs-linalg/sparse

# Run on hand-constructed matrices (Tier 1)
./test-data/scripts/run-spral.sh hand-constructed

# Run on CI subset (Tier 2)
./test-data/scripts/run-spral.sh suitesparse-ci

# Run on full SuiteSparse set (Tier 3 — requires archive extraction)
./test-data/scripts/run-spral.sh suitesparse
```

## Step 5: Verify Reproducibility

```bash
# Run twice on CI subset and diff results
./test-data/scripts/verify-reproducibility.sh suitesparse-ci
```

## Step 6: Generate Summary

```bash
./test-data/scripts/generate-summary.sh
# Outputs test-data/spral-reference/summary.json
```

## Key File Locations

| Artifact | Path |
|----------|------|
| SPRAL build | `/opt/references/spral/builddir/` |
| Driver source | `test-data/scripts/spral_driver.c` |
| Run script | `test-data/scripts/run-spral.sh` |
| Tier 1 results | `test-data/spral-reference/hand-constructed/` |
| Tier 2 results | `test-data/spral-reference/suitesparse-ci/` |
| Tier 3 results | `test-data/spral-reference/suitesparse/` |
| Config record | `test-data/spral-reference/config.json` |
| Summary | `test-data/spral-reference/summary.json` |

## Validation Checklist

- [ ] SPRAL builds and passes its own tests
- [ ] Driver processes a hand-constructed matrix and outputs valid JSON
- [ ] All 15 hand-constructed matrices produce results
- [ ] All 10 CI subset matrices produce results
- [ ] All 67 SuiteSparse matrices produce results (or document failures)
- [ ] Running driver twice on same matrix gives identical output
- [ ] Summary report generated with difficulty categorizations
- [ ] Tier 1 and Tier 2 results committed to git
- [ ] Tier 3 results in gitignored directory
