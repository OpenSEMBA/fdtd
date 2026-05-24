# C++ test migration status

Native GoogleTest coverage lives in `cpp_tests` (`-DSEMBA_FDTD_BUILD_CXX=ON`). Python Tier-0 smoke tests can target `cpp_build_nomtln/bin/semba-fdtd-cpp` via `test/pyWrapper/utils.py` (override with `SEMBA_EXE`).

**Fortran-fidelity rule:** C++ outputs must match the Fortran binary on the same `.fdtd.json` inputs; pytest compares against Fortran reference behavior, not ad-hoc tolerances.

## Minimal no-MTLN gate (Tier 0 + Tier 1)

Run the full slim gate with one script:

```bash
./scripts/test_cpp_nomtln_slim.sh
```

Or run tiers individually (see below).

## Build profile: `cpp_build_nomtln`

Slim executable (no full `semba-components` / `semba-main` chain):

```bash
cmake -S . -B cpp_build_nomtln \
  -DSEMBA_FDTD_BUILD_CXX=ON \
  -DSEMBA_FDTD_ENABLE_MTLN=OFF \
  -DSEMBA_FDTD_ENABLE_MPI=OFF \
  -DSEMBA_FDTD_ENABLE_HDF=OFF \
  -DSEMBA_FDTD_COMPONENTS_LIB=OFF \
  -DSEMBA_FDTD_OUTPUTS_LIB=OFF \
  -DSEMBA_FDTD_MAIN_LIB=OFF \
  -DSEMBA_FDTD_EXECUTABLE=ON \
  -DSEMBA_FDTD_ENABLE_TEST=ON
cmake --build cpp_build_nomtln -j --target semba-fdtd-cpp cpp_tests
```

Environment for pytest skip logic:

```bash
export SEMBA_FDTD_ENABLE_MTLN=OFF
export SEMBA_FDTD_ENABLE_MPI=OFF
export SEMBA_FDTD_ENABLE_HDF=OFF
export SEMBA_EXE=$PWD/cpp_build_nomtln/bin/semba-fdtd-cpp
```

## Tier 0 Python smoke (`semba-fdtd-cpp`)

Infrastructure + geometry only (1-step runs, `-mapvtk` tagging). Marked `@pytest.mark.cpp_slim`.

| Test | Status |
|------|--------|
| `test_fdtd_set_new_folder_to_run` | pass |
| `test_fdtd_clean_up_after_run` | pass |
| `test_1_volume` | pass (minimal `-mapvtk`) |
| `test_1_line` | pass |
| `test_2_volumes` | pass |
| `test_fill_conformal_vtk_corner` | pass |

Run Tier 0 pytest:

```bash
PYTHONPATH=. pytest test/pyWrapper/ -m cpp_slim -x -v
```

## Tier 1 Core solver parity (GoogleTest)

`SystemInitSolver.MatchesFortranHyHzAfterExPulse` in `cpp_tests` mirrors `test/system/test_init_solver.F90` (Hy/Hz nonzero after Ex pulse + one H step). `test_rank_remapping` remains deferred until `advanceEx` matches Fortran literals (±33.8822708).

```bash
./cpp_build_nomtln/bin/cpp_tests --gtest_filter=SystemInitSolver.*
./cpp_build_nomtln/bin/cpp_tests --gtest_filter=-mtln.*
```

## Tier 2 First time-domain parity (deferred)

| Test | Status | Blocker |
|------|--------|---------|
| `test_planewave_in_box` | **pending** | TF/SF Huygens box + Mur E parity with `planewaves.F90` / `bordersmur.F90`; `applyMurE()` empty in `semba_fdtd.cpp` |

Do **not** include Tier 2 in the minimal gate until `test_planewave_in_box` passes against Fortran reference outputs.

## Explicit exclusions (minimal gate)

| Category | Marker / skip | Why |
|----------|---------------|-----|
| MTLN | `@no_mtln_skip` | Requires MTLN ON |
| MPI | `@no_mpi_skip` | Not in slim build |
| HDF / movies | `@no_hdf_skip`, `@pytest.mark.movie` | Needs `semba-outputs` |
| Wires / lumped / nodal | `@pytest.mark.wires`, `lumped`, `nodal_source` | Needs `semba-components` |
| SGBC / dielectric physics | `@pytest.mark.sgbc`, `dielectric` | Not in slim solver loop |
| Conformal physics | `@pytest.mark.conformal` + probe parity | mapvtk-only tests OK; RCS/delay need full roll-up |

**Note:** `@mtln_skip` tests (classic wires without MTLN) still require `wires.cpp` from `semba-components` — out of scope for the slim build.

## Ported to `cpp_tests`

| Suite | Tests | Notes |
|-------|-------|--------|
| `cells` | 16 | |
| `mesh` | 15 | |
| `idchildtable` | 5 + fhash | |
| `smbjson_cpp` | 17–23 | |
| `conformal` | 34 | |
| `rotate` | 16–17 | |
| `mtln` | 28 | MTLN ON |
| `observation` | 13 | |
| `preprocess` | 6 | |
| **system** | **1** | `init_solver` Hy/Hz check via `semba-fdtd-core` |

## Slim C++ executable architecture

- `semba-fdtd-core`: simplified Yee solver + JSON probes + minimal `-mapvtk` (`mapvtk_writer.cpp`)
- `semba-reports`: `errorreport_core.cpp` (Fortran-aligned signaling)
- `semba-fdtd-cpp`: `launcher.cpp` → `semba_fdtd_t::init/launch/end`

Full `solver_t` / `timestepping.cpp` path remains blocked until `semba-components` / `semba-outputs` compile.

## Known gaps

- **Plane-wave TF/SF:** inject on Huygens box faces (not domain boundaries); port `corrigeondaplanaH` + second-order Mur from Fortran for stable 400-step `pw-in-box`.
- **Spice RC transient tests:** isolated subprocesses on Linux (ngspice shared state).

## Future expansion (after Tier 0 + Tier 1 green)

1. Fix plane-wave TF/SF → add `test_planewave_in_box` to gate (Tier 2)
2. Enable `SEMBA_FDTD_COMPONENTS_LIB=ON` with MTLN still OFF → classic wire tests (`@mtln_skip`)
3. Enable `SEMBA_FDTD_OUTPUTS_LIB=ON` → HDF movie tests
4. Re-enable MTLN as separate build (`cpp_build_mtln`) and test matrix
