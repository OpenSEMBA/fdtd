# Testing

Run all commands in this document from the repository root.

## Prerequisites

Initialise the repository submodules before configuring the project:

```shell
git submodule update --init --recursive
```

The native tests require a working CMake, C/C++ compiler, and Fortran compiler.
The Python tests require Python 3 and the dependencies in `requirements.txt`.
MPI tests additionally require an MPI implementation such as Open MPI.

For platform-specific compiler and library setup, see
[`development.md`](development.md).

## Build

Configure and build a Release tree with native tests enabled:

```shell
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DSEMBA_FDTD_ENABLE_TEST=ON
cmake --build build -j
```

The `build/` directory is intentional. Some Python end-to-end tests look for
the solver at `build/bin/semba-fdtd`.

For a Debug build, change `Release` to `Debug`. The project also provides CMake
presets, for example:

```shell
cmake --preset dbg
cmake --build --preset dbg
```

Preset builds use directories such as `build-dbg/`.

## Native Tests

Run all tests registered with CTest:

```shell
ctest --test-dir build --output-on-failure
```

List the registered tests without running them:

```shell
ctest --test-dir build -N
```

The main native test executable can also be run directly:

```shell
./build/bin/fdtd_tests
```

### Running a Specific Native Test

CTest test names can be listed before running them:

```shell
ctest --test-dir build -N
```

Run one exact CTest test with `-R`:

```shell
ctest --test-dir build -R '^fdtd_unit$' --output-on-failure
```

The expression passed to `-R` is a regular expression. For example, run all
tests whose names contain `output`:

```shell
ctest --test-dir build -R output --output-on-failure
```

Add `-V` to show the complete command and output for a test:

```shell
ctest --test-dir build -V -R '^fdtd_unit$'
```

The native executable uses GoogleTest. List its individual test cases with:

```shell
./build/bin/fdtd_tests --gtest_list_tests
```

Run an individual GoogleTest case with `--gtest_filter`:

```shell
./build/bin/fdtd_tests --gtest_filter='TestSuiteName.TestName'
```

Wildcards can be used in a GoogleTest filter:

```shell
./build/bin/fdtd_tests --gtest_filter='*Conformal*'
```

## Python Tests

Create a virtual environment and install the test dependencies:

```shell
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Run the complete Python test suite:

```shell
pytest test/ --durations=20
```

Equivalent invocation:

```shell
python3 -m pytest test/
```

Run a specific file or test:

```shell
pytest test/pyWrapper/test_integration.py
pytest test/pyWrapper/test_integration.py -k test_name
```

Pytest also supports selecting a test by its full node ID. The node ID has the
form `path/to/test_file.py::test_name`, or includes a class name for class
based tests:

```shell
pytest test/pyWrapper/test_integration.py::test_name
pytest test/pyWrapper/test_integration.py::TestClass::test_name
```

Use `--collect-only` to discover available test names without running them:

```shell
pytest test/pyWrapper/test_integration.py --collect-only -q
```

The `-k` expression can match part of a test name. Use `-vv` for more detail
and `-x` to stop after the first failure:

```shell
pytest test/ -k 'conformal and not mpi' -vv
pytest test/ -k test_name -x
```

## Test Markers

Use markers to select tests for optional features:

```shell
pytest test/ -m mtln
pytest test/ -m hdf
pytest test/ -m mpi
pytest test/ -m "not mpi"
```

Available markers are listed in [`../pytest.ini`](../pytest.ini), including
`mtln`, `hdf`, `mpi`, `spice`, `vtk`, `conformal`, `wires`, and `movie`.

## MPI Tests

Configure and build with MPI enabled:

```shell
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DSEMBA_FDTD_ENABLE_TEST=ON \
  -DSEMBA_FDTD_ENABLE_MPI=ON
cmake --build build -j
```

Run the CTest suite, including MPI tests:

```shell
ctest --test-dir build --output-on-failure
```

Run the Python MPI tests separately:

```shell
pytest test/ -m mpi
```

MPI tests require `mpirun` or `mpiexec` to be available in `PATH`.
