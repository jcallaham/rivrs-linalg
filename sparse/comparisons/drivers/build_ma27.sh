#!/bin/bash
# Build MA27 benchmark driver for comparison testing.
#
# MA27 source is NOT redistributable — it must be obtained from HSL:
#   https://www.hsl.rl.ac.uk/catalogue/ma27.html
#   (Free for academic use; commercial license available)
#
# Expected source layout (set MA27_SRC to override):
#   $MA27_SRC/
#   ├── ma27ad.f     (or ma27d.f, depending on HSL package)
#   └── (dependency sources if separate)
#
# Produces: /tmp/ma27_benchmark

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DRIVER="$SCRIPT_DIR/ma27_benchmark.f90"

# MA27 source location
# Check environment, then common locations
MA27_SRC="${MA27_SRC:-}"
if [ -z "$MA27_SRC" ]; then
    for candidate in /opt/hsl/ma27 /workspace/rivrs-linalg/references/ma27; do
        if [ -d "$candidate" ]; then
            MA27_SRC="$candidate"
            break
        fi
    done
fi
MA27_SRC="${MA27_SRC:-/opt/hsl/ma27}"

if [ ! -d "$MA27_SRC" ]; then
    echo "ERROR: MA27 source directory not found at: $MA27_SRC"
    echo ""
    echo "MA27 is not redistributable and must be obtained from HSL."
    echo ""
    echo "To get MA27:"
    echo "  1. Visit https://www.hsl.rl.ac.uk/catalogue/ma27.html"
    echo "  2. Register for a free academic license (or purchase commercial)"
    echo "  3. Download the source package"
    echo "  4. Extract to /opt/hsl/ma27/ (or set MA27_SRC env var)"
    echo ""
    echo "Expected files in \$MA27_SRC:"
    echo "  ma27ad.f (or ma27d.f) — the main MA27 source"
    echo ""
    echo "To build after placing the source:"
    echo "  MA27_SRC=/path/to/ma27 comparisons/drivers/build_ma27.sh"
    exit 1
fi

echo "=== Building MA27 benchmark driver ==="
echo "MA27 source: $MA27_SRC"

# Find the main MA27 source file (naming varies by package version)
MA27_MAIN=""
for name in ma27ad.f ma27d.f ma27.f; do
    if [ -f "$MA27_SRC/$name" ]; then
        MA27_MAIN="$MA27_SRC/$name"
        break
    fi
done

if [ -z "$MA27_MAIN" ]; then
    echo "ERROR: Cannot find MA27 Fortran source in $MA27_SRC"
    echo "Looked for: ma27ad.f, ma27d.f, ma27.f"
    ls -la "$MA27_SRC"/ 2>/dev/null || true
    exit 1
fi
echo "MA27 source: $MA27_MAIN"

# Check for dependency files (some HSL packages include fd15 for machine constants)
MA27_DEPS=""
for dep in fd15ad.f fd15d.f; do
    if [ -f "$MA27_SRC/$dep" ]; then
        MA27_DEPS="$MA27_DEPS $MA27_SRC/$dep"
    fi
done

# Build temporary directory
BUILD_DIR="/tmp/ma27_build"
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Compile MA27 source
echo "Compiling MA27..."
gfortran -O2 -fPIC -c "$MA27_MAIN"
if [ -n "$MA27_DEPS" ]; then
    for dep in $MA27_DEPS; do
        echo "Compiling dependency: $(basename $dep)"
        gfortran -O2 -fPIC -c "$dep"
    done
fi

# Compile and link the benchmark driver
echo "Compiling benchmark driver..."
gfortran -O2 -o /tmp/ma27_benchmark \
    "$DRIVER" \
    "$BUILD_DIR"/*.o \
    -lopenblas -lm

echo "=== Done ==="
echo "Driver binary: /tmp/ma27_benchmark"
ls -la /tmp/ma27_benchmark
