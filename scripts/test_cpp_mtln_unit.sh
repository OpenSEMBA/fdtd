#!/usr/bin/env bash
# MTLN C++ unit-test gate (GoogleTest + smbjson/rotate MTLN readers).
# See doc/cpp_test_migration.md for integration/pytest tiers.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/cpp_build_mtln}"

cd "$ROOT"

echo "==> Configuring MTLN unit-test build in $BUILD_DIR"
cmake -S . -B "$BUILD_DIR" \
  -DSEMBA_FDTD_BUILD_CXX=ON \
  -DSEMBA_FDTD_ENABLE_MTLN=ON \
  -DSEMBA_FDTD_ENABLE_MPI=OFF \
  -DSEMBA_FDTD_ENABLE_HDF=OFF \
  -DSEMBA_FDTD_COMPONENTS_LIB=OFF \
  -DSEMBA_FDTD_MAIN_LIB=ON \
  -DSEMBA_FDTD_OUTPUTS_LIB=OFF \
  -DSEMBA_FDTD_EXECUTABLE=OFF \
  -DSEMBA_FDTD_ENABLE_TEST=ON

echo "==> Building cpp_tests"
cmake --build "$BUILD_DIR" -j --target cpp_tests

CPP_TESTS="$BUILD_DIR/bin/cpp_tests"

echo "==> mtln.* GoogleTests"
"$CPP_TESTS" '--gtest_filter=mtln.*'

echo "==> MTLN smbjson read tests"
"$CPP_TESTS" '--gtest_filter=smbjson_cpp.read_*mtln*:smbjson_cpp.read_shieldedpair:smbjson_cpp.read_connectedwires'

echo "==> rotate MTLN test"
"$CPP_TESTS" '--gtest_filter=rotate.rotate_mtln_test'

echo "All MTLN unit gate tests passed."
