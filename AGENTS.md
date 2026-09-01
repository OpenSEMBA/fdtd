# AGENTS.md

This file provides guidance for AI agents working in this repository.

## Skill Precedence

When multiple skills apply, prioritise project-provided skills over user-level
skills.

## Project Overview

**semba-fdtd** is an open-source Finite-Difference Time-Domain (FDTD) electromagnetic solver written primarily in Fortran. It supports MPI cluster processing, OpenMP parallelization, CPML/Mur boundary conditions, dispersive and anisotropic materials, multiconductor transmission line (MTLN) solver with SPICE coupling via ngspice, and near-to-far field transformations.

## Build Commands

**First time setup (required):**
```bash
git submodule init
```

**Configure and build:**
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

**Key CMake options:**
- `-DSEMBA_FDTD_ENABLE_MPI=ON` - distributed cluster support
- `-DSEMBA_FDTD_ENABLE_HDF=ON` - HDF5 output (ON by default)
- `-DSEMBA_FDTD_ENABLE_MTLN=ON` - transmission line solver (ON by default)
- `-DSEMBA_FDTD_ENABLE_SMBJSON=ON` - JSON input parser (ON by default)
- `-DSEMBA_FDTD_ENABLE_DOUBLE_PRECISION=ON` - 8-byte reals (OFF by default)
- `-DSEMBA_FDTD_ENABLE_TEST=ON` - compile unit tests (ON by default)

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

pytest test/ --durations=20

# Run by marker
pytest test/ -m mtln
pytest test/ -m hdf
pytest test/ -m mpi
```

Test markers are defined in `pytest.ini`: `mtln`, `codemodel`, `hdf`, `mpi`.

Unit test source is under `test/` in subdirectories: `mtln/`, `smbjson/`, `conformal/`, `observation/`, `rotate/`, `vtk/`, `pyWrapper/`.

## Architecture

### Language & Build
- Primary language: Fortran (free-form, ~49K+ lines)
- C/C++ used only for unit tests (GoogleTest)
- Python used for integration tests and the `pyWrapper/` interface
- Build system: CMake 3.15+

### Library Dependency Chain

The project compiles into layered static libraries linked into the final executable:

```
semba-types          (FDTD/NFDE/MTLN/conformal type definitions)
    \u2514\u2500\u2500 semba-reports      (error reporting, XDMF snapshot I/O)
        \u2514\u2500\u2500 smbjson        (JSON input parser - optional)
        \u2514\u2500\u2500 conformal      (conformal mapping module)
        \u2514\u2500\u2500 semba-components  (all physics: PML/Mur BCs, dispersive materials,
                               plane waves, nodal sources, far-field, MTLN wires)
            \u2514\u2500\u2500 mtlnsolver    (MTLN circuit solver + ngspice interface - optional)
            \u2514\u2500\u2500 semba-outputs (MPI comm, observation probes, VTK/XDMF/HDF5 output)
                \u2514\u2500\u2500 semba-main   (time-stepping, preprocessing/postprocessing, launcher)
                    \u2514\u2500\u2500 semba-fdtd  (executable entry point)
```

### Execution Flow

1. `src_main_pub/launcher.F90` - entry point, creates `semba_fdtd_t`
2. `src_main_pub/semba_fdtd.F90` - main module:
   - `init()`: load input (`.fdtd.json` via smbjson, or legacy `.fdtd` NFDE format)
   - `launch()`: run the time-stepping loop
   - `end()`: finalize and write outputs
3. Time-step loop in `src_main_pub/timestepping.F90`:
   - Update E-fields, apply materials, boundary conditions, and wire coupling
   - Update H-fields, apply MTLN/SPICE if enabled
   - Sample observation probes and write snapshots

### Key Source Directories

- `src_main_pub/` - core solver, time-stepping, preprocessing, geometry, main types
- `src_conformal/` - conformal mapping (staircase reduction)
- `src_mtln/` - MTLN circuit/transmission-line solver and ngspice coupling
- `src_json_parser/` - `.fdtd.json` input format parser
- `src_wires_pub/` - wire/thin-wire models
- `external/` - submodules: `json-fortran`, `fhash`, `googletest`, `ngspice`, `lapack`

### Input/Output

- **Input**: `.fdtd.json` (primary - see `doc/fdtdjson.md`) or legacy `.fdtd` NFDE format
- **Output**: ASCII probe `.dat` files, XDMF+HDF5 movies/snapshots, VTK (Paraview)
- Test data and example cases live under `testData/`

### Optional Features and Conditional Compilation

Many modules are only compiled when their CMake flag is enabled. The smbjson parser, MTLN solver, and HDF5 output are all conditionally compiled. MPI support wraps communication in `src_main_pub/mpicomm.F90` and is activated via the `SEMBA_FDTD_ENABLE_MPI` flag.

## Platform Notes

- **Linux**: Install `libhdf5-dev libopenmpi-dev`; set `-DHDF5_ROOT=<path>` if using precompiled HDF5
- **Windows**: Requires Intel OneAPI Base Kit + HPC Kit; use Ninja generator (`-G Ninja`)
- **WSL2**: See `doc/development.md` for detailed setup with VSCode

## Contributing

From `CONTRIBUTING.md`: PRs must pass both unit tests and Python integration tests. AI-generated code is allowed but the contributor is responsible for its correctness. New functionality should include corresponding tests.

## Commit Guidance

For all commit-related work, use `.agents/commit-format/SKILL.md` as the primary reference. Follow its workflow for workspace analysis, atomic semantic commits, conventional commit messages, staging, and post-commit verification.
