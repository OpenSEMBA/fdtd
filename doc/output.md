# Probe Output

## Scalar Text Output

Point, wire, bulk, line, and far-field probes publish formatted `.dat` files
only.
These files are written beside the originating `.fdtd.json` input and are not
placed in probe-specific directories.
They do not have binary artifacts or metadata sidecars.

Time-domain filenames end in `_tm.dat`.
Frequency-domain point-probe filenames end in `_fq.dat`.
A point probe configured for both domains publishes one `.dat` file for each
domain.
Far-field filenames end in `.dat` and are always frequency-domain results.

## Volumetric And Geometry Output

Geometry maps retain their geometry and text artifacts.
Volumetric outputs retain their binary, XDMF, and HDF5 artifacts.

Probe outputs do not create JSON descriptors or a run output manifest.
Binary values retain their established byte order, numeric representation,
record size, component order, and complex-value convention.
