# Architecture

The project builds a layered Fortran solver with optional parser, output, conformal, and MTLN components. CMake assembles these layers into static libraries and links them into the final `semba-fdtd` executable.

## Entry point and execution flow

- `src_main_pub/launcher.F90` defines the executable entry point.
- The launcher creates a `semba_fdtd_t` instance and calls `init()`, `launch()`, and `end()`.
- `src_main_pub/semba_fdtd.F90` owns high-level solver setup, input parsing, runtime flags, and orchestration.
- `src_main_pub/timestepping.F90` defines the `solver_t` type and the main stepping operations that advance fields, sources, boundaries, outputs, and optional MTLN behavior.

## Main code areas

- `src_main_pub/`: solver control flow, field stepping, observations, VTK/XDMF output, MPI communication.
- `src_json_parser/`: FDTD JSON parser and related labels/types.
- `src_mtln/`: multiconductor transmission-line solver and ngspice integration points.
- `src_conformal/`: conformal geometry support.
- `src_wires_pub/`: wire models and MTLN wire coupling.

## Build composition

At the CMake level, the main libraries are organized around:

- `semba-types`
- `semba-reports`
- `smbjson` when JSON input is enabled
- `mtlnsolver` and `ngspice_interface` when MTLN is enabled
- `conformal`
- `semba-components`
- `semba-outputs`

Those pieces are linked into the final executable and, when enabled, the test binaries.
