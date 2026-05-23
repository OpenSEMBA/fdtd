# C++ unit test migration status

Native GoogleTest coverage lives in `cpp_tests` (`-DSEMBA_FDTD_BUILD_CXX=ON`). Fortran `fdtd_tests` remains until CI switches the default gate.

## Ported to `cpp_tests`

| Suite | Tests | Notes |
|-------|-------|--------|
| `cells` | 16 | Expanded vs single Fortran wrapper |
| `mesh` | 15 | |
| `idchildtable` | 5 + `idchildtable_fhash` | |
| `smbjson_cpp` | 23 (MTLN ON) / 17 (MTLN OFF) | Full `readProblemDescription` parity |
| `conformal` | 34 | Geometry, cell_map, filling; 5 contour/filling cases pending C++ geometry parity |
| `rotate` | 16 (OFF) / 17 (ON) | `nfde_rotate_m.h` inline ports |
| `mtln` | 28 registered | 14 runnable + 14 skipped (dispersive/spice until ngspice C++ port) |
| `observation` | 13 | 11 runnable + 2 skipped (init/update movie need full `Observa_m` harness) |
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

## Known gaps

- **Conformal (5 tests):** `geometry_vertex_vertex_contour`, `geometry_vertex_side_contour`, `geometry_areas`, `filling_edge_next_cell`, `filling_closed_corner` — C++ `buildSidesContour` / contour area parity with Fortran.
- **MTLN multipolar (3 tests):** in-test numerical port; may need link to `multipolar_expansion.cpp`.
- **MTLN dispersive + spice (12 tests):** `GTEST_SKIP` until `mtl_bundle` dispersive ctor and `circuit_m` ngspice interface are ported.
- **Observation init/update (2 tests):** `GTEST_SKIP` until `InitObservation` / `UpdateObservation` C++ harness exists.
- **System (2 tests):** deferred — full FDTD time step via C++ executable.

## Still Fortran-only

Python `pytest` integration tests are unchanged and target the Fortran binary.

`fdtd_tests` Fortran wrappers remain for regression until `cpp_tests` is CI-default.
