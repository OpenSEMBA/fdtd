#!/usr/bin/env bash
# MTLN Tier 1 gate: unit tests + semba-fdtd-cpp + mtlnProblem:true pytest.
# See doc/cpp_test_migration.md. Do not use pytest -m mtln over all pyWrapper (Tier 2 blocked).

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/cpp_build_mtln}"

cd "$ROOT"

echo "==> Configuring MTLN slim build in $BUILD_DIR"
cmake -S . -B "$BUILD_DIR" \
  -DSEMBA_FDTD_BUILD_CXX=ON \
  -DSEMBA_FDTD_ENABLE_MTLN=ON \
  -DSEMBA_FDTD_ENABLE_MPI=OFF \
  -DSEMBA_FDTD_ENABLE_HDF=OFF \
  -DSEMBA_FDTD_COMPONENTS_LIB=OFF \
  -DSEMBA_FDTD_MAIN_LIB=OFF \
  -DSEMBA_FDTD_OUTPUTS_LIB=OFF \
  -DSEMBA_FDTD_EXECUTABLE=ON \
  -DSEMBA_FDTD_ENABLE_TEST=ON

echo "==> Building cpp_tests and semba-fdtd-cpp"
cmake --build "$BUILD_DIR" -j --target cpp_tests semba-fdtd-cpp

CPP_TESTS="$BUILD_DIR/bin/cpp_tests"

echo "==> MTLN unit gate (mtln.* + smbjson/rotate)"
"$CPP_TESTS" '--gtest_filter=mtln.*'
"$CPP_TESTS" '--gtest_filter=smbjson_cpp.read_*mtln*:smbjson_cpp.read_shieldedpair:smbjson_cpp.read_connectedwires:smbjson_cpp.read_paul_*'
"$CPP_TESTS" '--gtest_filter=rotate.rotate_mtln_test'

export SEMBA_FDTD_ENABLE_MTLN=ON
export SEMBA_FDTD_ENABLE_MPI=OFF
export SEMBA_FDTD_ENABLE_HDF=OFF
export SEMBA_EXE="$BUILD_DIR/bin/semba-fdtd-cpp"

echo "==> Tier 1: pytest (mtlnProblem:true standalone)"
PYTHONPATH=. pytest test/pyWrapper/test_mtln_standalone.py -m mtln_standalone -v

echo "All MTLN Tier 1 slim gate tests passed."
