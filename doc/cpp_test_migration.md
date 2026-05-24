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

## MTLN unit gate (GoogleTest, `cpp_build_mtln`)

MTLN library + `cpp_tests` only (no `semba-fdtd-cpp` required):

```bash
./scripts/test_cpp_mtln_unit.sh
```

Equivalent manual steps:

```bash
cmake -S . -B cpp_build_mtln \
  -DSEMBA_FDTD_BUILD_CXX=ON \
  -DSEMBA_FDTD_ENABLE_MTLN=ON \
  -DSEMBA_FDTD_ENABLE_MPI=OFF \
  -DSEMBA_FDTD_ENABLE_HDF=OFF \
  -DSEMBA_FDTD_COMPONENTS_LIB=OFF \
  -DSEMBA_FDTD_MAIN_LIB=OFF \
  -DSEMBA_FDTD_OUTPUTS_LIB=OFF \
  -DSEMBA_FDTD_EXECUTABLE=OFF \
  -DSEMBA_FDTD_ENABLE_TEST=ON
cmake --build cpp_build_mtln -j --target cpp_tests
./cpp_build_mtln/bin/cpp_tests '--gtest_filter=mtln.*'
./cpp_build_mtln/bin/cpp_tests '--gtest_filter=smbjson_cpp.read_*mtln*:rotate.rotate_mtln_test'
```

## MTLN Tier 1 slim gate (`mtlnProblem: true`)

One command: unit tests, `semba-fdtd-cpp`, and standalone pytest (8 cases):

```bash
./scripts/test_cpp_mtln_slim.sh
```

This runs GoogleTest (`mtln.*`, smbjson MTLN readers including `read_paul_8_6_square_no_endpoint_wire_generators`), then:

```bash
export SEMBA_FDTD_ENABLE_MTLN=ON
export SEMBA_EXE=$PWD/cpp_build_mtln/bin/semba-fdtd-cpp
PYTHONPATH=. pytest test/pyWrapper/test_mtln_standalone.py -m mtln_standalone -v
```

Manual build (if not using the script): same CMake flags as [`scripts/test_cpp_mtln_slim.sh`](scripts/test_cpp_mtln_slim.sh) — MTLN ON, executable ON, `COMPONENTS_LIB`/`MAIN_LIB` OFF.

**Status (2026-05):** unit gate is green (28 `mtln.*` + smbjson/rotate). **Tier 1 standalone** (`test/pyWrapper/test_mtln_standalone.py`, 8 tests) passes on `cpp_build_mtln` (`corrcoef > 0.999` vs Fortran references). Fixes: `IsGeneratorOnWire` interior nodes only (endpoint sources via Spice only); full `writeNodeDescription` terminations; `connectNodesToNetworkCircuit` + opamp-style `TERMINATION_NETWORK` netlists.

### Tier 1 vs Tier 2 (`-m mtln`)

| Tier | Pytest target | JSON flag | C++ path today | `cpp_build_mtln` |
|------|---------------|-----------|----------------|------------------|
| **1 — standalone** | `test_mtln_standalone.py` | `"mtlnProblem": true` | `semba_fdtd_t` → `solveMTLNProblem()` only (no Yee grid coupling) | **pass** (8/8) |
| **2 — in-grid** | `test_integration.py` / `test_full_system.py` with `@pytest.mark.mtln` or `@no_mtln_skip` | `"mtlnProblem": false` or omitted | Needs `InitWires_mtln` / `AdvanceWiresE_mtln` in full `timestepping` + `semba-components` | **blocked** |

**Tier 2 examples (do not expect pass on slim `semba-fdtd-cpp`):**

| Test | Case | `mtlnProblem` |
|------|------|---------------|
| `test_holland_case_checking_number_of_outputs_*` | `holland/holland1981.fdtd.json` | false |
| `test_shieldedPair`, `test_unshielded_multiwires` | `shieldedPair`, `unshielded_multiwires_berenger` | false |
| `test_bundles_mpi_*` | `mpi/bundles_for_mpi*.fdtd.json` | false |
| `test_holland_mtln_mpi` | `holland/holland1981_unshielded.fdtd.json` | false |
| `test_current_generators_*`, `test_voltage_generators` | `sources/sources_*.fdtd.json` | false |
| `test_multiwire_z_collision_*` (`@no_mtln_skip`, mapvtk) | observation cases | false |

**Unblock Tier 2:** build `cpp_build_mtln2` with `-DSEMBA_FDTD_MAIN_LIB=ON` and `-DSEMBA_FDTD_COMPONENTS_LIB=ON` so `timestepping.cpp` can call the Fortran-equivalent in-grid wire path. Until then, CI should run `./scripts/test_cpp_mtln_slim.sh` (or `-m mtln_standalone`), not `pytest -m mtln` over the whole `test/pyWrapper/` tree.

## Future expansion (after Tier 0 + Tier 1 green)

1. Fix plane-wave TF/SF → add `test_planewave_in_box` to gate (Tier 2)
2. Enable `SEMBA_FDTD_COMPONENTS_LIB=ON` with MTLN still OFF → classic wire tests (`@mtln_skip`)
3. Enable `SEMBA_FDTD_OUTPUTS_LIB=ON` → HDF movie tests
4. ~~MTLN Tier 1~~ done; next: `cpp_build_mtln2` for in-grid `test_full_system.py -m mtln`
