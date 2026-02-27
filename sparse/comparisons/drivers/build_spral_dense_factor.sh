#!/bin/bash
# Build the SPRAL dense APTP kernel driver for comparison testing.
#
# Prerequisites:
#   1. Run comparisons/drivers/build_spral.sh first (builds libspral.a)
#   2. Have OpenBLAS installed (libopenblas-dev)
#
# Produces: /tmp/spral_dense_factor
#
# Usage:
#   ./comparisons/drivers/build_spral_dense_factor.sh
#   cargo run --example export_frontal --release
#   /tmp/spral_dense_factor [-v] [-u 0.01] [-b 256] < /tmp/frontal_matrix.txt

set -euo pipefail

SPRAL=/opt/references/spral
SRC=$SPRAL/src
SPRAL_BUILD=/tmp/spral_ssids
TOOLS_DIR="$(cd "$(dirname "$0")" && pwd)"
DRIVER=$TOOLS_DIR/spral_dense_factor.cpp
OUTPUT=/tmp/spral_dense_factor

if [ ! -f "$SPRAL_BUILD/libspral.a" ]; then
    echo "ERROR: $SPRAL_BUILD/libspral.a not found."
    echo "Run comparisons/drivers/build_spral.sh first."
    exit 1
fi

echo "=== Building SPRAL dense factor driver ==="

CXXFLAGS="-O2 -std=c++11 -I$SPRAL/include -I$SRC -I$SPRAL_BUILD"

# Compile the driver
echo "Compiling $DRIVER..."
g++ $CXXFLAGS -c "$DRIVER" -o /tmp/spral_dense_factor.o

# Link against libspral.a and BLAS
echo "Linking..."
g++ -o "$OUTPUT" /tmp/spral_dense_factor.o \
    -L"$SPRAL_BUILD" -lspral \
    -lopenblas -lstdc++ -lm -lpthread

echo "=== Done ==="
echo "Binary: $OUTPUT"
ls -la "$OUTPUT"
echo ""
echo "Usage:"
echo "  cargo run --example export_frontal --release"
echo "  $OUTPUT [-v] [-u 0.01] [-b 256] < /tmp/frontal_matrix.txt"
