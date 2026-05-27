# C++ test migration status

Native GoogleTest coverage lives in `cpp_tests`
(`-DSEMBA_FDTD_BUILD_CXX=ON`). Python migration smoke tests target
`cpp_build_nomtln/bin/semba-fdtd-cpp` through `test/pyWrapper/utils.py`
(override with `SEMBA_EXE`).

**Fortran-fidelity rule:** C++ outputs must match the Fortran binary on the
same `.fdtd.json` inputs; pytest compares against Fortran reference behavior,
not ad-hoc tolerances.

## Continuous integration

Workflow: [`.github/workflows/cpp.yml`](../.github/workflows/cpp.yml).

| Job | Script | Scope |
| --- | --- | --- |
| `nomtln-migration` | `./scripts/test_cpp_nomtln.sh` | `pytest -m cpp_migration` plus no-MTLN `cpp_tests` |
| `hdf-output` | `./scripts/test_cpp_hdf.sh` | `XdmfH5` gtest plus `pytest -m hdf` (HDF5 build) |
| `mtln-migration` | `./scripts/test_cpp_mtln.sh` | MTLN GoogleTest plus `pytest -m mtln_standalone` |

Local equivalent:

```bash
./scripts/test_cpp_nomtln.sh
./scripts/test_cpp_hdf.sh      # needs libhdf5-dev; HDF movie output
./scripts/test_cpp_mtln.sh   # needs liblapack-dev; builds ngspice submodule
```

## Active C++ executable gate

The current executable is built through the original-named `semba-main` target:

```bash
cmake -S . -B cpp_build_nomtln \
  -DSEMBA_FDTD_BUILD_CXX=ON \
  -DSEMBA_FDTD_ENABLE_MTLN=OFF \
  -DSEMBA_FDTD_ENABLE_MPI=OFF \
  -DSEMBA_FDTD_ENABLE_HDF=OFF \
  -DSEMBA_FDTD_COMPONENTS_LIB=OFF \
  -DSEMBA_FDTD_OUTPUTS_LIB=OFF \
  -DSEMBA_FDTD_MAIN_LIB=ON \
  -DSEMBA_FDTD_EXECUTABLE=ON \
  -DSEMBA_FDTD_ENABLE_TEST=ON
cmake --build cpp_build_nomtln -j --target semba-fdtd-cpp cpp_tests
```

The full translated `semba-components` and `semba-outputs` libraries are still
not part of this gate because several generated/stubbed modules do not compile
yet. They should be enabled module-by-module as their Fortran translations are
completed.

## Python migration marker

Short executable smoke tests use `@pytest.mark.cpp_migration`.

```bash
export SEMBA_FDTD_ENABLE_MTLN=OFF
export SEMBA_FDTD_ENABLE_MPI=OFF
export SEMBA_FDTD_ENABLE_HDF=OFF
export SEMBA_EXE=$PWD/cpp_build_nomtln/bin/semba-fdtd-cpp
PYTHONPATH=. pytest test/pyWrapper/ -m cpp_migration -x -v
```

Currently marked tests cover wrapper execution, cleanup, and one-step geometry
or map-VTK cases. Broader full-system tests remain in their existing pytest
modules and should be promoted to the gate only when the corresponding C++
physics path is byte-compatible with Fortran.

## GoogleTest migration gate

```bash
./cpp_build_nomtln/bin/cpp_tests --gtest_filter=SystemInitSolver.*
./cpp_build_nomtln/bin/cpp_tests \
  '--gtest_filter=-HollandWire.*:BordersMur.MurCxMatchesFortranCalcMurconstants310:PlanewavePwInBox.ProbeFilesExact_First120Steps:PlanewavePwInBox.ProbeFilesExact_FullRun:PlanewavePwInBox.PeriodicProbeFilesExact_FullRun'
```

Holland, lumped, dielectric, surface-impedance, far-field, and planewave strict
tests live in `cpp_tests` or the full-system pytest modules. They should remain
strict against Fortran probe bytes where reference files exist.

The excluded Mur/planewave entries are known strict parity failures, not removed
tests; run them directly while working on those physics paths.

## MTLN gates

Unit-only MTLN checks:

```bash
./scripts/test_cpp_mtln_unit.sh
```

Standalone MTLN integration checks:

```bash
./scripts/test_cpp_mtln.sh
```

`test_mtln_standalone.py` covers `"mtlnProblem": true` cases. In-grid MTLN
cases still depend on the full translated wire/timestepping path and should stay
separate from non-MTLN migration work.
