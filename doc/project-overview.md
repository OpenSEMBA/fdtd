# Project overview

`semba-fdtd` solves electromagnetic problems with a finite-difference time-domain core and optional multiconductor transmission-line coupling. The codebase is centered on a Fortran solver, with CMake-based builds, GoogleTest-based unit tests, and Python tooling for integration tests and preprocessing/postprocessing workflows.

## Main capabilities

- MPI-based distributed execution and OpenMP parallelism.
- Structured regular or rectilinear meshes.
- PEC, PMC, Mur, and CPML boundary conditions.
- Lossless, lossy, dispersive, anisotropic, and thin-material models.
- Embedded wires and multiconductor transmission-line networks with optional ngspice coupling.
- Probe-based outputs, near-to-far field transforms, VTK export, and XDMF/HDF5 output.

## Documentation map

- Start with [](getting-started.md) for installation, build, and first-run steps.
- Use [](user-guide.md) for input files, tutorials, outputs, and container workflows.
- Use [](developer-guide.md) for architecture, testing, and contribution workflows.
- Use [](solver-topics.md) for focused solver-domain documentation such as MTLN.

## Repository layout at a glance

- `src_main_pub/`: main solver, launcher, preprocessing, postprocessing, core types.
- `src_json_parser/`: FDTD JSON parser.
- `src_mtln/`: multiconductor transmission-line solver and ngspice interface.
- `src_conformal/`: conformal mapping support.
- `src_wires_pub/`: wire models.
- `test/`: unit and Python integration tests.
- `testData/`: sample cases, inputs, and supporting assets.
- `doc/`: documentation sources, now also used to build this site.
