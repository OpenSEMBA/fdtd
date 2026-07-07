# Persistence Contract

## Database Contract

| Item | Contract | Evidence |
| --- | --- | --- |
| Database tables | None for the output subsystem | VERIFIED |
| Database columns | None for the output subsystem | VERIFIED |
| Database indexes | None for the output subsystem | VERIFIED |
| Database foreign keys | None for the output subsystem | VERIFIED |
| Database constraints | None for the output subsystem | VERIFIED |
| Database writes | The output subsystem writes no database records | VERIFIED |
| Database reads | The output subsystem reads no database records | VERIFIED |

## Filesystem Persistence

| Artifact | Creation rule | Update rule | Evidence |
| --- | --- | --- | --- |
| Probe folder | Created at registration for every probe | Owns all files for that probe | VERIFIED |
| Time data file | Created inside the probe folder for time outputs | Appended on flush | VERIFIED |
| Frequency data file | Created inside the probe folder for point frequency outputs | Replaced on flush | VERIFIED |
| Output register | Created only when at least one output exists | Written once at registration | VERIFIED |
| Volumetric binary | Created inside the probe folder at registration | Appended or replaced by output type | VERIFIED |
| Volumetric HDF5 | Created inside the probe folder for volumetric time output | Appended during flush | VERIFIED |
| Volumetric XDMF | Created or replaced inside the probe folder during flush | Rewritten during flush | VERIFIED |
| Geometry map VTU | Created inside the probe folder during map registration | Not updated by time stepping | VERIFIED |
| Geometry metadata | Created inside the probe folder during map registration | Not updated by time stepping | VERIFIED |

## Compatibility-Sensitive Names

| Name | Meaning | Change risk | Evidence |
| --- | --- | --- | --- |
| `_tm.dat` | time-domain data suffix | High | VERIFIED |
| `_fq.dat` | frequency-domain point data suffix | High | VERIFIED |
| `<probe-folder>/<probe-basename>_tm.dat` | time-domain file placement | Critical | VERIFIED |
| `<probe-folder>/<probe-basename>_fq.dat` | frequency-domain file placement | Critical | VERIFIED |
| `.bin` | binary stream output suffix | High | VERIFIED |
| `.h5` | HDF5 output suffix | High | VERIFIED |
| `.xdmf` | XDMF output suffix | High | VERIFIED |
| `.vtu` | unstructured VTK output suffix | High | VERIFIED |
| `coordsX` | HDF5 x-coordinate dataset | High | VERIFIED |
| `coordsY` | HDF5 y-coordinate dataset | High | VERIFIED |
| `coordsZ` | HDF5 z-coordinate dataset | High | VERIFIED |
| `times` | HDF5 time dataset | High | VERIFIED |
| `CurrenDensityX` | HDF5 current-density x dataset | Critical | VERIFIED |
| `CurrenDensityY` | HDF5 current-density y dataset | Critical | VERIFIED |
| `CurrenDensityZ` | HDF5 current-density z dataset | Critical | VERIFIED |
| `ElectricFieldX` | HDF5 electric-field x dataset | High | VERIFIED |
| `ElectricFieldY` | HDF5 electric-field y dataset | High | VERIFIED |
| `ElectricFieldZ` | HDF5 electric-field z dataset | High | VERIFIED |
| `MagneticFieldX` | HDF5 magnetic-field x dataset | High | VERIFIED |
| `MagneticFieldY` | HDF5 magnetic-field y dataset | High | VERIFIED |
| `MagneticFieldZ` | HDF5 magnetic-field z dataset | High | VERIFIED |
| `xVal` | HDF5 frequency-slice x-value dataset | High | VERIFIED |
| `yVal` | HDF5 frequency-slice y-value dataset | High | VERIFIED |
| `zVal` | HDF5 frequency-slice z-value dataset | High | VERIFIED |

## File Row and Record Shapes

| Artifact | Shape | Evidence |
| --- | --- | --- |
| Point time row | time, scalar value | VERIFIED |
| Point frequency row | frequency, real value, imaginary value | VERIFIED |
| Wire current row | time, current, delta voltage, plus voltage, minus voltage, voltage difference | VERIFIED |
| Wire charge row | time, charge | VERIFIED |
| Block row | time, integrated value | VERIFIED |
| Movie binary record | time, x, y, z, x value, y value, z value | VERIFIED |
| Frequency binary record | frequency, x, y, z, x complex, y complex, x complex repeated | VERIFIED |

## Do-Not-Change Warnings

| Concern | Warning | Evidence |
| --- | --- | --- |
| Dataset spelling | `CurrenDensity` is missing the letter `t`; preserve it | VERIFIED |
| Frequency binary record | The third component repeats x instead of writing z | VERIFIED |
| Time flush behaviour | time-domain rows append across flushes | VERIFIED |
| Point frequency flush | frequency-domain point rows replace earlier rows | VERIFIED |
| Probe folder ownership | every probe-generated file is inside its own folder | VERIFIED |
| Coordinate suffix order | non-MPI suffix order is x, y, z | VERIFIED |
| MPI suffix order | MPI suffix order changes by selected axis | VERIFIED |

## Non-Applicable Persistence Areas

| Area | Reason | Evidence |
| --- | --- | --- |
| Soft deletes | no database rows exist | VERIFIED |
| Audit fields | no database rows exist | VERIFIED |
| Foreign keys | no database rows exist | VERIFIED |
| Transaction isolation | no database transaction exists | VERIFIED |
| Database constraints | no database schema exists for this subsystem | VERIFIED |
