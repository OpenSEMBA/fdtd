#!/usr/bin/env bash
# Build with Intel oneAPI
set -euo pipefail

BUILD_TYPE="Release"
ENABLE_MPI="OFF"
ENABLE_MTLN="ON"
ENABLE_HDF="ON"
ENABLE_DOUBLE="OFF"
BUILD_DIR="build-intel"
JOBS=$(nproc)

# --------------------------------------------------------------------------- #
# Compiler environment setup (Intel oneAPI)
# --------------------------------------------------------------------------- #
INTEL_ENV_SCRIPT="/opt/intel/oneapi/setvars.sh"
# setvars picks up the latest Intel MPI (2021.16) which lacks compiler
# wrappers; explicitly source 2021.14 which provides mpiifx/mpiicx/mpiicpx.
MPI_VARS="/opt/intel/oneapi/mpi/2021.14/env/vars.sh"

echo "Sourcing Intel oneAPI environment: $INTEL_ENV_SCRIPT"
set +u
source "$INTEL_ENV_SCRIPT" --force || true
if [[ -f "$MPI_VARS" ]]; then
    source "$MPI_VARS" || true
fi
set -u

export FC="ifx"
export CC="icx"
export CXX="icpx"
echo "Using Intel MPI compiler wrappers: FC=$FC CC=$CC CXX=$CXX"

# --------------------------------------------------------------------------- #
# CMake configure
# --------------------------------------------------------------------------- #
echo ""
echo "--- CMake configure ---"
cmake -S . -B "$BUILD_DIR" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
    -DSEMBA_FDTD_ENABLE_MPI="$ENABLE_MPI" \
    -DSEMBA_FDTD_ENABLE_HDF="$ENABLE_HDF" \
    -DSEMBA_FDTD_ENABLE_MTLN="$ENABLE_MTLN" \
    -DSEMBA_FDTD_ENABLE_DOUBLE_PRECISION="$ENABLE_DOUBLE"

# --------------------------------------------------------------------------- #
# CMake build
# --------------------------------------------------------------------------- #
echo ""
echo "--- CMake build (jobs=$JOBS) ---"
cmake --build "$BUILD_DIR" -j "$JOBS"

echo ""
echo "Build finished successfully."
