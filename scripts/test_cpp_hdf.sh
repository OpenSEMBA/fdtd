#!/usr/bin/env bash
# C++ build with HDF5 output enabled; runs @pytest.mark.hdf integration tests.
# See doc/cpp_test_migration.md.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/cpp_build_hdf}"

cd "$ROOT"

echo "==> Configuring HDF-enabled C++ build in $BUILD_DIR"
cmake -S . -B "$BUILD_DIR" \
  -DSEMBA_FDTD_BUILD_CXX=ON \
  -DSEMBA_FDTD_ENABLE_MTLN=OFF \
  -DSEMBA_FDTD_ENABLE_MPI=OFF \
  -DSEMBA_FDTD_ENABLE_HDF=ON \
  -DSEMBA_FDTD_COMPONENTS_LIB=OFF \
  -DSEMBA_FDTD_OUTPUTS_LIB=ON \
  -DSEMBA_FDTD_MAIN_LIB=ON \
  -DSEMBA_FDTD_EXECUTABLE=ON \
  -DSEMBA_FDTD_ENABLE_TEST=ON

echo "==> Building semba-fdtd-cpp, semba-outputs, cpp_tests"
cmake --build "$BUILD_DIR" -j --target semba-fdtd-cpp semba-outputs cpp_tests

export SEMBA_FDTD_ENABLE_MTLN=OFF
export SEMBA_FDTD_ENABLE_MPI=OFF
export SEMBA_FDTD_ENABLE_HDF=ON
export SEMBA_EXE="$BUILD_DIR/bin/semba-fdtd-cpp"

echo "==> GoogleTest: XdmfH5"
"$BUILD_DIR/bin/cpp_tests" --gtest_filter=XdmfH5.*

echo "==> pytest: HDF marker"
PYTHONPATH=. pytest test/pyWrapper/ -m hdf -x -v

echo "All HDF C++ gate tests passed."
