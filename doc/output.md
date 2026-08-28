# Probe Output

## Scalar Text Output

Point, non-MTLN wire, bulk, line, and far-field probes publish formatted
`.dat` files only.
These files are written beside the originating `.fdtd.json` input and are not
placed in probe-specific directories.
They do not have binary artifacts, probe descriptors, or manifest entries.

Time-domain filenames end in `_tm.dat`.
Frequency-domain point-probe filenames end in `_fq.dat`.
A point probe configured for both domains publishes one `.dat` file for each
domain.
Far-field filenames end in `.dat` and are always frequency-domain results.

## Coordinated Output

Geometry maps, time-domain volumetric movies, frequency-domain volumetric
slices, and MTLN outputs retain their applicable descriptors and manifest
entries.
Volumetric outputs retain their binary, XDMF, and HDF5 artifacts.

Binary values use the byte order, numeric representation, record size,
component order, and complex-value convention declared by their descriptors.
