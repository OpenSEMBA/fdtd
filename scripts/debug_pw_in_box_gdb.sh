#!/usr/bin/env bash
# GDB helper for pw-in-box TF/SF + Mur debugging.
# See doc/cpp_planewave_mur_mapping.md

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT/cpp_build_dbg}"
CASE_DIR="$ROOT/testData/cases/planewave"

if [[ ! -x "$BUILD_DIR/bin/semba-fdtd-cpp" ]]; then
  echo "Debug binary not found. Building..."
  cmake -S "$ROOT" -B "$BUILD_DIR" \
    -DSEMBA_FDTD_BUILD_CXX=ON \
    -DSEMBA_FDTD_ENABLE_MTLN=OFF \
    -DSEMBA_FDTD_ENABLE_MPI=OFF \
    -DSEMBA_FDTD_ENABLE_HDF=OFF \
    -DSEMBA_FDTD_COMPONENTS_LIB=OFF \
    -DSEMBA_FDTD_OUTPUTS_LIB=OFF \
    -DSEMBA_FDTD_MAIN_LIB=OFF \
    -DSEMBA_FDTD_EXECUTABLE=ON \
  -DCMAKE_BUILD_TYPE=Debug
  cmake --build "$BUILD_DIR" -j --target semba-fdtd-cpp
fi

cd "$CASE_DIR"
echo "GDB in $CASE_DIR"
echo "Breakpoints: evolucion, advancePlaneWaveE, applyMurH"
exec gdb -ex 'set pagination off' \
  -ex 'break semba_fdtd.cpp:485' \
  -ex 'break semba_fdtd.cpp:838' \
  -ex 'break semba_fdtd.cpp:632' \
  -ex 'break semba_fdtd.cpp:1160' \
  -ex "run --args $BUILD_DIR/bin/semba-fdtd-cpp -i pw-in-box.fdtd.json" \
  "$BUILD_DIR/bin/semba-fdtd-cpp"
