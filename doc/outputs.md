# Outputs and postprocessing

`semba-fdtd` supports several output styles, depending on the enabled features and the probes present in the input file.

## Common output forms

- Plain-text probe data, typically written as `.dat` files.
- XDMF + HDF5 outputs for snapshots and movies.
- VTK outputs for visualization workflows such as ParaView.

## How outputs are configured

Outputs are driven by the `probes` section of the `.fdtd.json` input. Probe type, sampling period, field/component selection, and spatial domain determine what gets written.

The full probe and domain reference is part of [](fdtdjson.md).

## Postprocessing workflows

- Use ParaView or other VTK/XDMF-aware tools for field visualization.
- Use Python scripts from `testData/` and `doc/tutorials/` as examples for reading and plotting probe outputs.
- Use the Python wrapper tests and examples when you want to drive solver execution and analysis programmatically.
