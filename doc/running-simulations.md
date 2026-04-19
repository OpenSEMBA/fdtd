# Running simulations

## Input format

The solver consumes `.fdtd.json` files. The complete object-by-object reference is documented in [](fdtdjson.md), and sample inputs live in `testData/input_examples/`.

## Basic command

Run the solver from the directory that contains your case files:

```bash
build/bin/semba-fdtd -i CASE_NAME.fdtd.json
```

If you installed a release binary instead of building from source, replace `build/bin/semba-fdtd` with the path to the installed executable.

## Typical workflow

1. Prepare the `.fdtd.json` input and any referenced excitation files.
2. Run the solver.
3. Inspect the generated probe data, VTK/XDMF/HDF5 outputs, or postprocess them with the Python tooling used by the repository tests and examples.

## Examples and preprocessing

The repository contains end-to-end examples in `testData/` and more guided narrative material in [](tutorials/index.md). The Python wrapper examples also illustrate how to generate excitation files and automate postprocessing around the solver.
