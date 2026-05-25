#!/usr/bin/env bash
# Minimal no-MTLN full-system C++ test gate (Tier 0 + Tier 1).
# See doc/cpp_test_migration.md for scope and exclusions.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/cpp_build_nomtln}"

cd "$ROOT"

echo "==> Configuring slim C++ build in $BUILD_DIR"
cmake -S . -B "$BUILD_DIR" \
  -DSEMBA_FDTD_BUILD_CXX=ON \
  -DSEMBA_FDTD_ENABLE_MTLN=OFF \
  -DSEMBA_FDTD_ENABLE_MPI=OFF \
  -DSEMBA_FDTD_ENABLE_HDF=OFF \
  -DSEMBA_FDTD_COMPONENTS_LIB=OFF \
  -DSEMBA_FDTD_OUTPUTS_LIB=OFF \
  -DSEMBA_FDTD_MAIN_LIB=OFF \
  -DSEMBA_FDTD_EXECUTABLE=ON \
  -DSEMBA_FDTD_ENABLE_TEST=ON

echo "==> Building semba-fdtd-cpp and cpp_tests"
cmake --build "$BUILD_DIR" -j --target semba-fdtd-cpp cpp_tests

export SEMBA_FDTD_ENABLE_MTLN=OFF
export SEMBA_FDTD_ENABLE_MPI=OFF
export SEMBA_FDTD_ENABLE_HDF=OFF
export SEMBA_EXE="$BUILD_DIR/bin/semba-fdtd-cpp"

echo "==> Tier 0: pytest (cpp_slim marker)"
PYTHONPATH=. pytest test/pyWrapper/ -m cpp_slim -x -v

echo "==> Tier 1: SystemInitSolver GoogleTest"
"$BUILD_DIR/bin/cpp_tests" --gtest_filter=SystemInitSolver.*

echo "==> Tier 1: full non-MTLN GoogleTest suite"
"$BUILD_DIR/bin/cpp_tests" --gtest_filter=-mtln.*

echo "All minimal no-MTLN gate tests passed."
