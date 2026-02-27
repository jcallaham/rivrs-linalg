#!/bin/bash
# Build MUMPS benchmark driver for comparison testing.
#
# Requires: gfortran, libmumps-seq-dev (sequential MUMPS with dummy MPI)
#
# On Debian/Ubuntu:
#   apt-get install libmumps-seq-dev
#
# Produces: /tmp/mumps_benchmark

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DRIVER="$SCRIPT_DIR/mumps_benchmark.f90"

# Check for MUMPS sequential library
if ! dpkg -l libmumps-seq-dev > /dev/null 2>&1; then
    echo "ERROR: libmumps-seq-dev not found."
    echo ""
    echo "Install with:  apt-get install libmumps-seq-dev"
    echo ""
    echo "This provides the sequential (dummy MPI) version of MUMPS,"
    echo "which is sufficient for single-process benchmarking."
    exit 1
fi

echo "=== Building MUMPS benchmark driver ==="

# Find include directory for MUMPS Fortran module files
MUMPS_INC=""
for dir in /usr/include /usr/include/mumps-seq-* /usr/include/x86_64-linux-gnu/mumps_seq; do
    if [ -f "$dir/dmumps_struc.h" ]; then
        MUMPS_INC="$dir"
        break
    fi
done

if [ -z "$MUMPS_INC" ]; then
    # Try a broader search
    MUMPS_INC=$(find /usr -name "dmumps_struc.h" -printf '%h\n' 2>/dev/null | head -1)
fi

if [ -z "$MUMPS_INC" ]; then
    echo "ERROR: Cannot find dmumps_struc.h. Is libmumps-seq-dev installed?"
    exit 1
fi
echo "MUMPS include: $MUMPS_INC"

# Compile the driver
# Link order matters: MUMPS, then its dependencies
gfortran -O2 -I "$MUMPS_INC" -o /tmp/mumps_benchmark \
    "$DRIVER" \
    -ldmumps_seq -lmumps_common_seq -lpord_seq -lmpiseq_seq \
    -lmetis -lopenblas -lm -lpthread

echo "=== Done ==="
echo "Driver binary: /tmp/mumps_benchmark"
ls -la /tmp/mumps_benchmark
