#!/bin/bash
# Build SPRAL SSIDS library (CPU-only) for comparison testing.
#
# Produces /tmp/spral_ssids/libspral.a and .mod files.
# Requires: gfortran, g++, libopenblas-dev, METIS (vendored via metis-sys)

set -euo pipefail

SPRAL_TOOLS="$(cd "$(dirname "$0")" && pwd)"
SPRAL=/opt/references/spral
SRC=$SPRAL/src
OUT=/tmp/spral_ssids
METIS_INCLUDE=$(find /home/node/.cargo/registry -path '*/metis-sys-*/vendor/metis/include' -type d 2>/dev/null | head -1)
METIS_LIB=$(find /workspace/rivrs-linalg/sparse/target -name "libmetis.a" 2>/dev/null | head -1)

if [ -z "$METIS_LIB" ]; then
    echo "ERROR: libmetis.a not found. Run 'cargo build' first to compile METIS."
    exit 1
fi
echo "METIS lib: $METIS_LIB"
echo "METIS include: $METIS_INCLUDE"

rm -rf "$OUT"
mkdir -p "$OUT"
cd "$OUT"

FFLAGS="-O2 -fPIC -J$OUT -I$OUT"
CXXFLAGS="-O2 -fPIC -std=c++11 -I$SPRAL/include -I$SRC -I$OUT"
FPPFLAGS="-DSPRAL_HAVE_METIS_H=0"

if [ -n "$METIS_INCLUDE" ]; then
    CXXFLAGS="$CXXFLAGS -I$METIS_INCLUDE"
fi

echo "=== Compiling Fortran modules ==="

# Leaf modules
gfortran $FFLAGS -c $SRC/blas_iface.f90
gfortran $FFLAGS -c $SRC/lapack_iface.f90
gfortran $FFLAGS -c $SRC/random.f90
gfortran $FFLAGS -c $SRC/matrix_util.f90
gfortran $FFLAGS -c $SRC/core_analyse.f90
gfortran $FFLAGS -c $SRC/cuda/cuda_nocuda.f90
gfortran $FFLAGS -c $SRC/pgm.f90
gfortran $FFLAGS -c $SRC/hw_topology/hw_topology.f90

# Preprocessing needed
gfortran $FFLAGS $FPPFLAGS -c $SRC/metis5_wrapper.F90

# Depends on matrix_util
gfortran $FFLAGS -c $SRC/scaling.f90

# Depends on metis_wrapper, scaling
gfortran $FFLAGS -c $SRC/match_order.f90

# Depends on matrix_util, random
gfortran $FFLAGS -c $SRC/rutherford_boeing.f90

# SSIDS modules in dependency order
gfortran $FFLAGS -c $SRC/ssids/datatypes.f90
gfortran $FFLAGS -c $SRC/ssids/inform.f90
gfortran $FFLAGS -c $SRC/ssids/contrib.f90
gfortran $FFLAGS -c $SRC/ssids/subtree.f90 -o ssids_subtree.o
gfortran $FFLAGS -c $SRC/ssids/cpu/cpu_iface.f90
gfortran $FFLAGS -c $SRC/ssids/cpu/subtree.f90 -o cpu_subtree.o
gfortran $FFLAGS -c $SRC/ssids/gpu/subtree_no_cuda.f90
gfortran $FFLAGS -c $SRC/ssids/akeep.f90
gfortran $FFLAGS -c $SRC/ssids/profile_iface.f90
gfortran $FFLAGS $FPPFLAGS -c $SRC/ssids/anal.F90
gfortran $FFLAGS $FPPFLAGS -c $SRC/ssids/fkeep.F90
gfortran $FFLAGS -c $SRC/ssids/contrib_free.f90
gfortran $FFLAGS -c $SRC/ssids/ssids.f90

echo "=== Compiling C++ modules ==="

g++ $CXXFLAGS -c $SRC/compat.cxx
g++ $CXXFLAGS -c $SRC/omp.cxx
g++ $CXXFLAGS -c $SRC/hw_topology/guess_topology.cxx
g++ $CXXFLAGS -c $SRC/ssids/profile.cxx
g++ $CXXFLAGS -c $SRC/ssids/cpu/ThreadStats.cxx
g++ $CXXFLAGS -c $SRC/ssids/cpu/SymbolicSubtree.cxx
g++ $CXXFLAGS -c $SRC/ssids/cpu/NumericSubtree.cxx
g++ $CXXFLAGS -c $SRC/ssids/cpu/kernels/wrappers.cxx
g++ $CXXFLAGS -c $SRC/ssids/cpu/kernels/cholesky.cxx
g++ $CXXFLAGS -c $SRC/ssids/cpu/kernels/ldlt_app.cxx
g++ $CXXFLAGS -c $SRC/ssids/cpu/kernels/ldlt_nopiv.cxx
g++ $CXXFLAGS -c $SRC/ssids/cpu/kernels/ldlt_tpp.cxx

echo "=== Building static library ==="

ar rcs libspral.a *.o

echo "=== Compiling SPRAL driver binaries ==="

# Benchmark driver (SPRAL does its own ordering)
gfortran -O2 -fopenmp -I $OUT -o /tmp/spral_benchmark \
  $SPRAL_TOOLS/spral_benchmark.f90 \
  -Wl,--whole-archive $OUT/libspral.a -Wl,--no-whole-archive \
  $METIS_LIB -lopenblas -lstdc++ -lm -lgomp
echo "Benchmark driver: /tmp/spral_benchmark"

# Full solve driver (user-supplied ordering, for spral_solve_comparison)
gfortran -O2 -fopenmp -I $OUT -o /tmp/spral_full_solve \
  $SPRAL_TOOLS/spral_full_solve.f90 \
  -Wl,--whole-archive $OUT/libspral.a -Wl,--no-whole-archive \
  $METIS_LIB -lopenblas -lstdc++ -lm -lgomp
echo "Full solve driver: /tmp/spral_full_solve"

echo "=== Done ==="
echo "Library: $OUT/libspral.a"
echo "Modules: $OUT/*.mod"
ls -la $OUT/libspral.a /tmp/spral_benchmark /tmp/spral_full_solve
