# C++ unit test migration status

Native GoogleTest coverage lives in `cpp_tests` (`-DSEMBA_FDTD_BUILD_CXX=ON`). Fortran `fdtd_tests` remains until CI switches the default gate.

## Ported to `cpp_tests`

| Suite | Tests | Notes |
|-------|-------|--------|
| `cells` | 16 | Expanded vs single Fortran wrapper |
| `mesh` | 15 | |
| `idchildtable` | 5 + `idchildtable_fhash` | |
| `smbjson_cpp` | 23 (MTLN ON) / 17 (MTLN OFF) | Full `readProblemDescription` parity |
| `conformal` | 34 | Geometry, cell_map, filling |
| `rotate` | 16 (OFF) / 17 (ON) | `nfde_rotate_m.h` inline ports |
| `mtln` | 28 | All runnable including dispersive + ngspice/circuit |
| `observation` | 13 | Init/update movie harness |
| `preprocess` | 6 | `searchtag` + tag checks via `preprocess_tags.h` |
| **system** | **deferred** | Needs `SEMBA_FDTD_MAIN_LIB` + working C++ `semba_fdtd_t` solver |

## Build

```bash
# MTLN OFF (default smbjson + conformal + rotate + observation + preprocess)
cmake -S . -B build_cpp -DSEMBA_FDTD_BUILD_CXX=ON -DSEMBA_FDTD_ENABLE_MTLN=OFF
cmake --build build_cpp -j --target cpp_tests
./build_cpp/bin/cpp_tests

# MTLN ON (+ mtln solver tests)
cmake -S . -B build_cpp_mtln -DSEMBA_FDTD_BUILD_CXX=ON -DSEMBA_FDTD_ENABLE_MTLN=ON
cmake --build build_cpp_mtln -j --target cpp_tests
./build_cpp_mtln/bin/cpp_tests
```

## Current counts (MTLN ON)

- **152 pass, 0 skip, 0 fail**

## Known gaps

- **System (2 tests):** deferred — full FDTD time step via C++ executable.
- **Spice RC transient tests:** `spice_tran` and `spice_tran_2` run in isolated subprocesses on Linux to avoid ngspice shared-state interference between back-to-back `setStopTimes` runs with identical node names.

## Still Fortran-only

Python `pytest` integration tests are unchanged and target the Fortran binary.

`fdtd_tests` Fortran wrappers remain for regression until `cpp_tests` is CI-default.
