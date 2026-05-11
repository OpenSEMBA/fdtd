# AGENTS.md

This file provides guidance for AI coding agents (Claude, GitHub Copilot, ChatGPT, etc.) when working with the **semba-fdtd** repository.

## Project Overview

**semba-fdtd** is an open-source Finite-Difference Time-Domain (FDTD) electromagnetic solver written primarily in Fortran. Key capabilities include:

- MPI cluster processing and OpenMP parallelization
- CPML/Mur boundary conditions
- Dispersive and anisotropic materials
- Multiconductor transmission line (MTLN) solver with SPICE coupling via ngspice
- Near-to-far field transformations
- Wire/thin-wire models and plane-wave sources

## Build System

### First-time Setup (REQUIRED)
```bash
git submodule init
git submodule update
```

### Build Commands
```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

### Key CMake Options
- `-DSEMBA_FDTD_ENABLE_MPI=ON` — distributed cluster support
- `-DSEMBA_FDTD_ENABLE_HDF=ON` — HDF5 output (ON by default)
- `-DSEMBA_FDTD_ENABLE_MTLN=ON` — transmission line solver (ON by default)
- `-DSEMBA_FDTD_ENABLE_SMBJSON=ON` — JSON input parser (ON by default)
- `-DSEMBA_FDTD_ENABLE_DOUBLE_PRECISION=ON` — 8-byte reals (OFF by default)
- `-DSEMBA_FDTD_ENABLE_TEST=ON` — compile unit tests (ON by default)
- `-DSEMBA_FDTD_ENABLE_CUDA_FORTRAN=ON` — CUDA Fortran GPU path (NVHPC compiler, requires `SEMBA_FDTD_ENABLE_CUF_RUNTIME=1` at runtime)

### CUDA Fortran Runtime Gate

CUDA Fortran execution is opt-in at runtime to avoid crashes on nodes without accessible CUDA devices.

```bash
export SEMBA_FDTD_ENABLE_CUF_RUNTIME=1
```

If this variable is not set, CUDA Fortran builds fall back to CPU execution at runtime.

**Binary output:** `./build/bin/semba-fdtd`

## Running Tests

### Unit Tests (C++/Fortran — GoogleTest)
```bash
./build/bin/fdtd_tests
```

### Python Integration Tests
```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
pytest test/ --durations=20
```

Test markers: `mtln`, `codemodel`, `hdf`, `mpi`

Unit tests live under `test/` in subdirectories: `mtln/`, `smbjson/`, `conformal/`, `observation/`, `rotate/`, `vtk/`, `pyWrapper/`.

## Architecture

### Language & Build
- **Primary language:** Fortran (free-form, ~49K+ lines)
- **C/C++:** Only for unit tests (GoogleTest)
- **Python:** Integration tests and `pyWrapper/` interface
- **Build system:** CMake 3.15+

### Library Dependency Chain

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

1. **`src_main_pub/launcher.F90`** — entry point, creates `semba_fdtd_t`
2. **`src_main_pub/semba_fdtd.F90`** — main module:
   - `init()`: load input (`.fdtd.json` via smbjson, or legacy `.fdtd` NFDE format)
   - `launch()`: run the time-stepping loop
   - `end()`: finalize and write outputs
3. **Time-step loop** in `src_main_pub/timestepping.F90`:
   - Update E-fields → apply materials, boundary conditions, wire coupling
   - Update H-fields → apply MTLN/SPICE if enabled
   - Sample observation probes, write snapshots

### Key Source Directories

| Directory | Purpose |
|-----------|---------|
| `src_main_pub/` | Core solver, time-stepping, preprocessing, geometry, main types |
| `src_conformal/` | Conformal mapping (staircase reduction) |
| `src_mtln/` | MTLN circuit/transmission-line solver and ngspice coupling |
| `src_json_parser/` | `.fdtd.json` input format parser |
| `src_wires_pub/` | Wire/thin-wire models |
| `src_pyWrapper/` | Python interface |
| `external/` | Submodules: `json-fortran`, `fhash`, `googletest`, `ngspice`, `lapack` |
| `testData/` | Test data and example cases |
| `doc/` | Documentation (fdtdjson.md, development.md, tutorials) |

### Input/Output

- **Input:** `.fdtd.json` (primary — see `doc/fdtdjson.md`) or legacy `.fdtd` NFDE format
- **Output:** ASCII probe `.dat` files, XDMF+HDF5 movies/snapshots, VTK (Paraview)

## Coding Conventions

### Fortran Style
- Free-form source code
- Follow the existing indentation and naming conventions in surrounding code
- Use meaningful module and variable names
- Always check for conditional compilation flags before modifying optional features

### Testing Requirements
- **Every PR must pass both unit tests and Python integration tests**
- When adding new functionality, include corresponding tests
- Unit tests: C++/Fortran under `test/` using GoogleTest
- Integration tests: Python under `test/` using pytest

### Conditional Compilation
Many modules are only compiled when their CMake flag is enabled. Be aware of:
- `SEMBA_FDTD_ENABLE_MPI` — MPI support (wraps communication in `src_main_pub/mpicomm.F90`)
- `SEMBA_FDTD_ENABLE_HDF` — HDF5 output
- `SEMBA_FDTD_ENABLE_MTLN` — MTLN solver + ngspice
- `SEMBA_FDTD_ENABLE_SMBJSON` — JSON input parser
- `SEMBA_FDTD_ENABLE_DOUBLE_PRECISION` — 8-byte reals

Always check `CMakeLists.txt` and use `#ifdef` guards when adding code that depends on optional features.

## Contributing Guidelines

From `CONTRIBUTING.md`:

1. **Keep changes focused** — small, well-scoped commits
2. **Follow existing code style** — match surrounding Fortran/Python conventions
3. **Update documentation** — add/update docs in `doc/` when behavior changes
4. **Add tests** — include tests for new functionality
5. **Write clear commit messages** — in English, describe the "why" not just the "what"
6. **Open PRs against `dev`** (or `main` per current workflow)

### AI-Assisted Contributions

AI coding agents may be used to assist with contributions. However:

- **Human author is responsible for correctness** of all AI-generated code
- PRs with **significant AI-generated content** require review by **at least two human reviewers**
- Clearly state in the PR description that AI assistance was used and describe the extent
- Ensure you have reviewed, understood, and taken responsibility for all AI-generated content

## Platform-Specific Notes

### Linux
```bash
sudo apt install libhdf5-dev libopenmpi-dev
```
Set `-DHDF5_ROOT=<path>` if using precompiled HDF5.

### Windows
Requires Intel OneAPI Base Kit + HPC Kit. Use Ninja generator (`-G Ninja`).

### WSL2
See `doc/development.md` for detailed setup with VSCode.

## Debugging

### VSCode Launch Configuration (`.vscode/launch.json`)
```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Fortran Launch (GDB)",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceRoot}/build/bin/semba-fdtd",
            "miDebuggerPath": "gdb",
            "args": ["-i", "shieldingEffectiveness.fdtd.json"],
            "stopAtEntry": false,
            "cwd": "${workspaceRoot}/tmp_cases/sgbcShieldingEffectiveness/"
        }
    ]
}
```

### MPI Debugging
1. Use `mpirun` to launch the parallel executable
2. Attach GDB to a running process using the "Attach" configuration
3. May require: `echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope`

## Common Tasks for AI Agents

### When modifying Fortran code:
1. Understand the module's purpose and dependencies
2. Check for conditional compilation guards
3. Follow existing code style (indentation, naming, structure)
4. Add tests if the change affects behavior
5. Update documentation if the API or input format changes

### When adding new features:
1. Determine which library layer the feature belongs to
2. Add corresponding unit tests and/or integration tests
3. Update `doc/` if the feature changes user-facing behavior
4. Document any new CMake options in CLAUDE.md and AGENTS.md

### When fixing bugs:
1. Reproduce the issue with existing test cases
2. Add a regression test that captures the bug
3. Fix the issue with minimal changes
4. Verify all tests pass

### When working with input files:
1. Primary format: `.fdtd.json` (schema in `doc/fdtdjson.md`)
2. Legacy format: `.fdtd` (NFDE format)
3. Test cases in `testData/`

## Useful Commands

```bash
# Initialize submodules (first time only)
git submodule init && git submodule update

# Build
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j

# Run unit tests
./build/bin/fdtd_tests

# Run Python tests
pytest test/ -v

# Run tests by marker
pytest test/ -m mtln
pytest test/ -m hdf
pytest test/ -m mpi

# Clean build
rm -rf build && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j
```

## Resources

- **Main README:** `README.md`
- **Contributing:** `CONTRIBUTING.md`
- **Development setup:** `doc/development.md`
- **JSON input format:** `doc/fdtdjson.md`
- **MTLN documentation:** `doc/mtln.md`
- **Docker setup:** `doc/docker.md`
- **Tutorial:** `doc/tutorials/veritasium/veritasium.md`

## Skills

Specialized workflows for common roles. Each skill file contains detailed guidance, key file references, and step-by-step workflows.

| Skill | File | Purpose |
|-------|------|---------|
| **Adding Unit Tests** | `.agents/skills/unit-tests/SKILL.md` | GoogleTest patterns (Fortran/C++ glue), pytest integration, test data management |
| **EMC Engineer** | `.agents/skills/emc-engineer/SKILL.md` | Writing `.fdtd.json` inputs, launching simulations, analyzing probe outputs with pyWrapper |
| **HPC Engineer** | `.agents/skills/hpc-engineer/SKILL.md` | MPI domain decomposition, OpenMP tuning, compiler optimization, profiling, GPU roadmap |
| **Code Reviewer** | `.agents/skills/code-reviewer/SKILL.md` | PR review checklist, conditional compilation guards, test coverage verification |
