---
name: unit-tests
description: GoogleTest patterns for semba-fdtd: Fortran/C++ glue with bind(C), testingTools assertion helpers (expect_eq_int, expect_near, checkNear), pytest markers (mtln/hdf/mpi/codemodel), test data management, step-by-step guide for adding new unit and integration tests
---

## When to use

- Adding a new feature or bug fix that needs test coverage
- Creating regression tests for a reported issue
- Modifying existing code and wanting to verify behavior hasn't changed
- Setting up a new test subdirectory under `test/`
- Adding Python integration tests for the solver

## Key Files to Read

### GoogleTest (C++/Fortran)
- `test/CMakeLists.txt` — Top-level test build configuration
- `test/fdtd_tests.cpp` — GoogleTest main entry point
- `test/smbjson/smbjson_testingTools.F90` — Assertion helpers for JSON parser tests
- `test/mtln/mtln_testingTools.F90` — Assertion helpers for MTLN tests
- `test/conformal/conformal_tests.h` — Example C++ glue header mapping Fortran to GoogleTest
- `test/pyWrapper/utils.py` — Python test utilities, pytest markers, FDTD/Probe classes

### Test Data
- `testData/cases/` — Full simulation cases used by integration tests
- `testData/input_examples/` — Minimal `.fdtd.json` files for unit test inputs
- `testData/outputs/` — Reference output data for regression comparison

### Python Integration Tests
- `test/pyWrapper/test_integration.py` — End-to-end solver tests
- `test/pyWrapper/test_full_system.py` — Full system tests with MPI, HDF5, MTLN
- `test/pyWrapper/test_mtln_standalone.py` — MTLN-specific integration tests

## GoogleTest Architecture (Tier 1)

### Directory Structure

Every test subdirectory under `test/` follows this pattern:

    test/<subdir>/
      CMakeLists.txt              # Builds a static test library
      <subdir>_tests.cpp          # 1-line C++ file: #include "<subdir>_tests.h"
      <subdir>_tests.h            # Maps Fortran functions to TEST() macros
      <subdir>_testingTools.F90   # Assertion helpers (if needed)
      test_*.F90                  # Fortran test functions

### Fortran Test Function Pattern

Each test is a `bind(C)` integer function returning `0` on success, non-zero on failure:

    integer function test_my_feature() bind(C) result(err)
       use some_module
       use some_testingTools
       implicit none
       integer :: err
       err = 0
       ! ... test logic ...
       if (.not. expected_condition) then
          err = err + 1
          call testFails(err, "Descriptive failure message")
       end if
    end function

### C++ Glue Header Pattern

The `.h` file declares `extern "C"` Fortran functions and maps them to GoogleTest:

    extern "C" int test_my_feature();

    TEST(subdir, my_feature) {
        EXPECT_EQ(0, test_my_feature());
    }

### Assertion Helpers

**Fortran testingTools modules** provide:

| Function | Purpose |
|----------|---------|
| `expect_eq_int(err, expected, provided, msg)` | Integer comparison |
| `expect_near(a, b, tol)` | Floating point comparison |
| `testFails(err, msg)` | Increment error counter with message |
| `comparePULMatrices(error_cnt, m_line, m_input)` | Matrix comparison |
| `checkNear(target, number, rel_tol)` | Relative tolerance check (multiple overloads) |
| `expect_eq(err, ex, pr, ignoreRegions)` | Full structural comparison of derived types |

### Conditional Compilation

Tests are included conditionally based on CMake flags. Check `test/CMakeLists.txt`:

- `SEMBA_FDTD_ENABLE_MTLN` -> enables `mtln/` and `system/` tests
- `SEMBA_FDTD_ENABLE_SMBJSON` -> enables `smbjson/`, `rotate/`, `vtk/`, `observation/`, `utils/` tests
- `SEMBA_FDTD_ENABLE_MPI` -> excludes `observation/` and `system/` from GoogleTest

**Always guard test code with the same `#ifdef` guards** as the source code it tests.

## Python Integration Tests (Tier 2)

### Directory Structure

Python tests live in `test/pyWrapper/`:

    test/pyWrapper/
      utils.py                      # FDTD/Probe classes, pytest markers, helpers
      test_pyWrapper.py             # pyWrapper class tests
      test_integration.py           # End-to-end solver tests
      test_full_system.py           # Full system tests (MPI, HDF5, MTLN)
      test_mtln_standalone.py       # MTLN-specific tests

### Pytest Markers

Defined in `pytest.ini`:

| Marker | Purpose |
|--------|---------|
| `mtln` | Tests using MTLN features |
| `codemodel` | Tests needing xspice codemodels |
| `hdf` | Tests needing HDF5 |
| `mpi` | Tests to be run with mpirun |

Apply via decorators:

    @pytest.mark.mtln
    def test_my_mtln_feature():
        ...

Skip logic (from `utils.py`):

    mtln_skip = pytest.mark.skipif(
        os.getenv("SEMBA_FDTD_ENABLE_MTLN") == "OFF",
        reason="MTLN not enabled"
    )

### Running Tests

    # All tests
    pytest test/ --durations=20

    # By marker
    pytest test/ -m mtln
    pytest test/ -m hdf
    pytest test/ -m mpi
    pytest test/ -m codemodel

    # Exclude marker
    pytest test/ -m "not codemodel"

### Integration Test Pattern

Integration tests run the actual solver and validate outputs:

    from utils import FDTD, Probe, CASES_FOLDER, OUTPUTS_FOLDER

    @pytest.mark.mtln
    def test_my_case(tmp_path):
        case_dir = CASES_FOLDER / "myCase"
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        solver = FDTD(
            input_filename="myCase.fdtd.json",
            path_to_exe=SEMBA_EXE,
            run_in_folder=str(output_dir)
        )
        solver.cleanUp()
        solver.run()
        assert solver.returncode == 0

        # Validate probe outputs
        probe_files = solver.getSolvedProbeFilenames("myProbe")
        assert len(probe_files) > 0
        probe = Probe(probe_files[0])
        # ... validate data against expected values ...

## How to Add a New Test

### For a new Fortran module:

1. **Create test subdirectory:** `test/<my_module>/`
2. **Write Fortran test functions** following the `bind(C)` integer-returning pattern
3. **Create testingTools module** if custom assertions are needed
4. **Create the `.h` header** mapping Fortran functions to `TEST(group, name)` macros
5. **Create 1-line `.cpp` glue file** with `#include "<my_module>_tests.h"`
6. **Add CMakeLists.txt** building a static library from the test files
7. **Register in `test/CMakeLists.txt`** — add `add_subdirectory(<my_module>)` under the appropriate conditional
8. **Include the header in `fdtd_tests.cpp`** with the matching `#ifdef` guard
9. **Add test data** to `testData/input_examples/` or `testData/cases/`
10. **Run `./build/bin/fdtd_tests`** to verify

### For a new Python integration test:

1. **Add test function** to the appropriate file in `test/pyWrapper/`
2. **Apply correct pytest marker** (`mtln`, `hdf`, `mpi`, `codemodel`)
3. **Use `tmp_path` fixture** for isolated temp directories
4. **Use `FDTD` class** from `utils.py` to launch simulations
5. **Use `Probe` class** to parse and validate `.dat` output files
6. **Compare against reference data** using `np.allclose`, correlation coefficients, or file content checks
7. **Run `pytest test/ -m <marker>`** to verify

## Test Data Management

### When adding test cases:

- Place `.fdtd.json` input files in `testData/input_examples/` for small, focused tests
- Place full simulation cases in `testData/cases/<caseName>/` for integration tests
- Place reference outputs in `testData/outputs/` for regression comparison
- Keep test cases minimal — remove unnecessary cells, steps, or probes while still exercising the code path
- Name files descriptively: `myFeature_<detail>.fdtd.json`

### When modifying test data:

- If you change expected outputs, update `testData/outputs/` accordingly
- If you change a `.fdtd.json` input, verify it against `doc/fdtdjson.md` for correctness
- Run the full test suite after modifying shared test data: `pytest test/ --durations=20`

## Common Gotchas

1. **Fortran test functions must use `bind(C)`** — without this, C++ cannot call them
2. **Return `0` for pass, non-zero for fail** — GoogleTest expects `EXPECT_EQ(0, ...)` in the C++ glue
3. **Conditional compilation must match** — if a feature is `#ifdef CompileWithMTLN` in source, tests must use the same guard
4. **MPI-free builds** — `observation/` and `system/` GoogleTest subdirectories are only compiled when `SEMBA_FDTD_ENABLE_MPI=OFF`
5. **Test data paths** — use `PATH_TO_TEST_DATA` parameter defined in testingTools modules, not hardcoded paths
6. **Python test isolation** — always use `tmp_path` fixture; never write to the repo directory
7. **Tolerance differences** — double precision builds (`CompileWithReal8`) may need looser tolerances than single precision
