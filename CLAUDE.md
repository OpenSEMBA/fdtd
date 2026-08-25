# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**semba-fdtd** is an open-source Finite-Difference Time-Domain (FDTD) electromagnetic solver written primarily in Fortran. It supports MPI cluster processing, OpenMP parallelization, CPML/Mur boundary conditions, dispersive and anisotropic materials, multiconductor transmission line (MTLN) solver with SPICE coupling via ngspice, and near-to-far field transformations.

## Build Commands

**First time setup (required):**
```bash
git submodule update --init --recursive
```

**Configure and build:**
```bash
cmake --fresh --preset rls
cmake --build --preset rls -j
```

**Key CMake options:**
- `-DSEMBA_FDTD_ENABLE_MPI=ON` — distributed cluster support
- `-DSEMBA_FDTD_ENABLE_MTLN=ON` — transmission line solver (ON by default)
- `-DSEMBA_FDTD_ENABLE_SMBJSON=ON` — JSON input parser (ON by default)
- `-DSEMBA_FDTD_ENABLE_DOUBLE_PRECISION=ON` — 8-byte reals (OFF by default)
- `-DSEMBA_FDTD_ENABLE_TEST=ON` — compile unit tests (ON by default)

**Binary output:** `./build-rls/bin/semba-fdtd` for the `rls` preset.

## Running Tests

**C++/Fortran unit tests (GoogleTest):**
```bash
./build-rls/bin/fdtd_tests
```

**Python integration tests:**
```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt

pytest test/ --durations=20

# Run by marker
pytest test/ -m mtln
pytest test/ -m hdf
SEMBA_FDTD_ENABLE_MPI=ON pytest test/ -m mpi
```

Test markers are defined in `pytest.ini`: `mtln`, `codemodel`, `hdf`, `mpi`.
See `doc/testing.md` for the complete testing workflow.

Native test source is under `test/conformal/`, `test/mpi/`, `test/mtln/`,
`test/smbjson/`, `test/system/`, `test/unit/`, and `test/utils/`.
Python tests are under `test/e2e/` and `test/pyWrapper/`.

## Architecture

### Language & Build
- Primary language: Fortran (free-form, ~49K+ lines)
- C/C++ used only for unit tests (GoogleTest)
- Python used for integration tests and the `pyWrapper/` interface
- Build system: CMake 3.15+

### Library Dependency Chain

The project compiles into layered static libraries linked into the final
executable:

```
semba-types       FDTD/NFDE/MTLN/conformal type definitions
semba-reports     error reporting
smbjson           JSON input parser (optional)
conformal         conformal mapping
semba-components  field, material, boundary, source, and wire physics
mtlnsolver        MTLN circuit solver and ngspice interface (optional)
semba-outputs     MPI communication
fdtd-output       probe writers, metadata, binary, XDMF/HDF5, and VTK output
semba-main        time-stepping, preprocessing, postprocessing, and launch flow
semba-fdtd        executable entry point
```

`semba-main` links the communication and output libraries into the solver.

### Execution Flow

1. `src_main_pub/launcher.F90` — entry point, creates `semba_fdtd_t`
2. `src_main_pub/semba_fdtd.F90` — main module:
   - `init()`: load input (`.fdtd.json` via smbjson, or legacy `.fdtd` NFDE format)
   - `launch()`: run the time-stepping loop
   - `end()`: finalize and write outputs
3. Time-step loop in `src_main_pub/timestepping.F90`:
   - Update E-fields → apply materials, boundary conditions, wire coupling
   - Update H-fields → apply MTLN/SPICE if enabled
   - Sample observation probes, write snapshots

### Key Source Directories

- `src_main_pub/` — core solver, time-stepping, preprocessing, geometry, main types
- `src_conformal/` — conformal mapping (staircase reduction)
- `src_mtln/` — MTLN circuit/transmission-line solver and ngspice coupling
- `src_json_parser/` — `.fdtd.json` input format parser
- `src_wires_pub/` — wire/thin-wire models
- `external/` — submodules: `json-fortran`, `fhash`, `googletest`, `ngspice`, `lapack`

### Input/Output

- **Input**: `.fdtd.json` (primary — see `doc/fdtdjson.md`) or legacy `.fdtd` NFDE format
- **Output**: ASCII probe `.dat` files, XDMF+HDF5 movies/snapshots, and VTK;
  see `doc/output.md`
- Test data and example cases live under `testData/`

### Optional Features and Conditional Compilation

The smbjson parser, MTLN solver, and MPI support are conditionally compiled.
HDF5/XDMF output is required.
MPI communication is implemented in `src_main_pub/mpicomm.F90` and activated
with `SEMBA_FDTD_ENABLE_MPI`.

## Platform Notes

- **Linux**: Install `libhdf5-dev libopenmpi-dev`; set `-DHDF5_ROOT=<path>` if using precompiled HDF5
- **Windows**: Requires Intel OneAPI Base Kit + HPC Kit; use Ninja generator (`-G Ninja`)
- **WSL2**: See `doc/development.md` for detailed setup with VSCode

## Contributing

From `CONTRIBUTING.md`: PRs must pass both unit tests and Python integration tests. AI-generated code is allowed but the contributor is responsible for its correctness. New functionality should include corresponding tests.
