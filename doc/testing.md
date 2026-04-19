# Testing

The repository uses both compiled tests and Python integration tests.

## Unit tests

After building with tests enabled, run:

```bash
./build/bin/fdtd_tests
```

## Python integration tests

Create a virtual environment and install the existing test dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
```

Then run the test suite from the repository root:

```bash
pytest test/ --durations=20
```

Running from the repository root matters because several tests reference `testData/` with relative paths.

## Test markers

The Python suite defines these markers in `pytest.ini`:

- `mtln`
- `codemodel`
- `hdf`
- `mpi`

Examples:

```bash
pytest test/ -m mtln
pytest test/ -m hdf
pytest test/ -m mpi
```

## Containerized test runs

If you prefer Docker-based workflows, [](docker.md) documents how to run both the unit tests and the Python suite inside the provided container setup.
