# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**semba-fdtd** is an open-source Finite-Difference Time-Domain (FDTD) electromagnetic solver written primarily in Fortran. It supports MPI cluster processing, OpenMP parallelization, CPML/Mur boundary conditions, dispersive and anisotropic materials, multiconductor transmission line (MTLN) solver with SPICE coupling via ngspice, and near-to-far field transformations.

## Build Commands

**First time setup (required):**
```bash
git submodule init
git submodule update
```

**Configure and build:**
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

**Key CMake options:**
- `-DSEMBA_FDTD_ENABLE_MPI=ON` — distributed cluster support
- `-DSEMBA_FDTD_ENABLE_HDF=ON` — HDF5 output (ON by default)
- `-DSEMBA_FDTD_ENABLE_MTLN=ON` — transmission line solver (ON by default)
- `-DSEMBA_FDTD_ENABLE_SMBJSON=ON` — JSON input parser (ON by default)
- `-DSEMBA_FDTD_ENABLE_DOUBLE_PRECISION=ON` — 8-byte reals (OFF by default)
- `-DSEMBA_FDTD_ENABLE_TEST=ON` — compile unit tests (ON by default)

**Binary output:** `./build/bin/semba-fdtd`

## Running Tests

**C++/Fortran unit tests (GoogleTest):**
```bash
./build/bin/fdtd_tests
```

**Python integration tests:**
```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt

# Fortran full suite (set feature flags to match the binary)
export SEMBA_EXE=$PWD/build/bin/semba-fdtd
export SEMBA_FDTD_ENABLE_MTLN=ON
export SEMBA_FDTD_ENABLE_HDF=ON
export SEMBA_FDTD_ENABLE_MPI=ON
PYTHONPATH=. pytest test/ --durations=20

# C++ migration gates (build semba-fdtd-cpp first, then pytest by marker)
cmake -S . -B cpp_build_nomtln -DSEMBA_FDTD_BUILD_CXX=ON \
  -DSEMBA_FDTD_ENABLE_MTLN=OFF -DSEMBA_FDTD_ENABLE_MPI=OFF -DSEMBA_FDTD_ENABLE_HDF=OFF \
  -DSEMBA_FDTD_MAIN_LIB=ON -DSEMBA_FDTD_EXECUTABLE=ON -DSEMBA_FDTD_ENABLE_TEST=ON
cmake --build cpp_build_nomtln -j --target semba-fdtd-cpp
export SEMBA_EXE=$PWD/cpp_build_nomtln/bin/semba-fdtd-cpp
export SEMBA_FDTD_ENABLE_MTLN=OFF SEMBA_FDTD_ENABLE_MPI=OFF SEMBA_FDTD_ENABLE_HDF=OFF
PYTHONPATH=. pytest test/pyWrapper/ -m cpp_migration

# HDF / MTLN: same pattern with cpp_build_hdf or cpp_build_mtln and -m hdf or
# test/pyWrapper/test_mtln_standalone.py -m mtln_standalone
```

Test markers are defined in `pytest.ini` (`cpp_migration`, `hdf`, `mtln_standalone`, `mtln`, `mpi`, …).

Unit test source: Fortran/C++ GoogleTest under `test/fortran/` (`mtln/`, `smbjson/`, etc.); C++ solver tests under `test/cpp/`; Python integration tests under `test/pyWrapper/`.

## Architecture

### Language & Build
- Primary language: Fortran (free-form, ~49K+ lines)
- C++ solver in `src_cpp/` (mirrors `src/` layout); C++ also used for unit tests (GoogleTest)
- Python used for integration tests and the `pyWrapper/` interface
- Build system: CMake 3.15+

### Library Dependency Chain

The project compiles into layered static libraries linked into the final executable:

```
semba-types          (FDTD/NFDE/MTLN/conformal type definitions)
    └── semba-reports      (error reporting, XDMF snapshot I/O)
        └── smbjson        (JSON input parser — optional)
        └── conformal      (conformal mapping module)
        └── semba-components  (all physics: PML/Mur BCs, dispersive materials,
                               plane waves, nodal sources, far-field, MTLN wires)
            └── mtlnsolver    (MTLN circuit solver + ngspice interface — optional)
            └── semba-outputs (MPI comm, observation probes, VTK/XDMF/HDF5 output)
                └── semba-main   (time-stepping, preprocessing/postprocessing, launcher)
                    └── semba-fdtd  (executable entry point)
```

### Execution Flow

1. `src/main/launcher.F90` — entry point, creates `semba_fdtd_t`
2. `src/main/semba_fdtd.F90` — main module:
   - `init()`: load input (`.fdtd.json` via smbjson, or legacy `.fdtd` NFDE format)
   - `launch()`: run the time-stepping loop
   - `end()`: finalize and write outputs
3. Time-step loop in `src/main/timestepping.F90`:
   - Update E-fields → apply materials, boundary conditions, wire coupling
   - Update H-fields → apply MTLN/SPICE if enabled
   - Sample observation probes, write snapshots

### Key Source Directories

- `src/main/` — core solver, time-stepping, preprocessing, geometry, main types
- `src/conformal/` — conformal mapping (staircase reduction)
- `src/mtln/` — MTLN circuit/transmission-line solver and ngspice coupling
- `src/json_parser/` — `.fdtd.json` input format parser
- `src/wires/` — wire/thin-wire models
- `external/` — submodules: `json-fortran`, `fhash`, `googletest`, `ngspice`, `lapack`

### Input/Output

- **Input**: `.fdtd.json` (primary — see `doc/fdtdjson.md`) or legacy `.fdtd` NFDE format
- **Output**: ASCII probe `.dat` files, XDMF+HDF5 movies/snapshots, VTK (Paraview)
- Test data and example cases live under `testData/`

### Optional Features and Conditional Compilation

Many modules are only compiled when their CMake flag is enabled. The smbjson parser, MTLN solver, and HDF5 output are all conditionally compiled. MPI support wraps communication in `src/main/mpicomm.F90` and is activated via the `SEMBA_FDTD_ENABLE_MPI` flag.

## Platform Notes

- **Linux**: Install `libhdf5-dev libopenmpi-dev`; set `-DHDF5_ROOT=<path>` if using precompiled HDF5
- **Windows**: Requires Intel OneAPI Base Kit + HPC Kit; use Ninja generator (`-G Ninja`)
- **WSL2**: See `doc/development.md` for detailed setup with VSCode

## Contributing

From `CONTRIBUTING.md`: PRs must pass both unit tests and Python integration tests. AI-generated code is allowed but the contributor is responsible for its correctness. New functionality should include corresponding tests.
