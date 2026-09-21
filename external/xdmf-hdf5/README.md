# XDMF/HDF5

`xdmf-hdf5` is a solver-independent Fortran library for producing
XDMF 3 metadata and its associated HDF5 heavy data.
It is developed in-tree initially, but has its own build, tests, package
metadata, and public API so it can later move to a separate repository.

## Scope

Version 0.1 provides:

- Uniform, rectilinear, curvilinear, unstructured, and mixed grids.
- Linear and quadratic XDMF topologies.
- Static, temporal, frequency, and generic parameter collections.
- Scalar, vector, tensor, matrix, and identifier attributes.
- Node, edge, face, cell, and grid-centred data.
- `real32`, `real64`, `int32`, and `int64` heavy data.
- Explicit real/imaginary or magnitude/phase attributes for complex data.

The library does not depend on solver types, OpenMP, or SMBJSON.
HDF5 identifiers and XML serialization are private implementation details.

When CMake discovers an MPI-capable HDF5 Fortran installation, the library also
supports collective scalar-series hyperslab writes.
Set `XDMF_HDF5_ENABLE_MPI=ON` to enable the MPI API.
Collective HDF5 support is detected automatically in either build mode.

## Standalone Build

```sh
cmake -S . -B build -DXDMF_HDF5_BUILD_TESTING=ON
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Example Generator

With `XDMF_HDF5_BUILD_TESTING=ON`, the build produces
`xdmf_hdf5_generate_cases` in the build's binary directory.
It generates the XDMF/HDF5 conformance examples:

```sh
./build/bin/xdmf_hdf5_generate_cases <output-directory> [options]
```

The output directory is created when necessary.
The generator preserves existing files by default and rejects the invocation
before writing when one of its output pairs already exists.
Use `--replace` to explicitly replace the generator's named output pairs:

```sh
./build/bin/xdmf_hdf5_generate_cases ./generated --replace
```

Available options are:

```text
--examples  Generate only the committed examples.
--help      Print command usage and exit.
--replace   Replace generated files that already exist.
```

Consumers link the exported target:

```cmake
find_package(XDMFHdf5 CONFIG REQUIRED)
target_link_libraries(my_writer PRIVATE XDMF::HDF5)
```

## API Example

The writer owns the pair and every HDF5 resource.
All operations return an `xdmf_status_t`; the library never prints or stops the
calling application.

```fortran
use, intrinsic :: iso_fortran_env, only: int64, real64
use xdmf_hdf5_m

type(xdmf_writer_t) :: writer
type(xdmf_options_t) :: options
type(xdmf_status_t) :: status
type(xdmf_grid_id_t) :: grid
type(xdmf_attribute_id_t) :: pressure
real(real64) :: values(24)

options%overwrite = .true.
options%series_kind = XDMF_SERIES_TIME

call writer%create('result', options, status)
call writer%define_uniform_grid('volume', [2_int64, 3_int64, 4_int64], &
  [0.0_real64, 0.0_real64, 0.0_real64], &
  [1.0_real64, 1.0_real64, 1.0_real64], grid, status)
call writer%define_attribute(grid, 'pressure', XDMF_CENTER_NODE, &
  XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., pressure, status)

call writer%begin_step(0.0_real64, status)
call writer%write_attribute(pressure, values, status)
call writer%end_step(status)
call writer%close(status)
```

Applications should check `status%is_error()` after each operation.
`status%message()` contains contextual failure information.
`XDMF_ERROR_CONSISTENCY` means an HDF5 rollback or resource close could not be
confirmed; the writer must then be closed and not reused.

`xdmf_writer_t` is a unique resource owner.
Do not assign or copy an open writer; pass it with `intent(inout)` and close it
explicitly before it leaves scope.

## Collective HDF5

Collective output is available only when the library was built against parallel
HDF5 with a compatible MPI Fortran implementation.
All ranks must call writer creation, definitions, step operations, flush, and
close in the same order.
Only `root_rank` publishes the XDMF document.

```fortran
options%overwrite = .true.
options%series_kind = XDMF_SERIES_TIME
options%collective_io = .true.
options%communicator = MPI_COMM_WORLD
options%root_rank = 0

call writer%create('result', options, status)
! Define the same grid and scalar series attribute on every rank.
call writer%begin_step(time, status)
call writer%write_attribute_hyperslab(field, local_values, &
  local_offset, local_shape, status)
call writer%end_step(status)
```

`local_offset` is zero-based and `local_shape` is in Fortran I/J/K order.
They select a disjoint portion of the globally defined scalar attribute.
Ranks with no cells pass zero for every `local_shape` entry and an empty value
array; they still participate in the collective call.
Compression is intentionally unavailable in collective mode.

## Data Conventions

The Fortran API accepts connectivity with one-based node indices.
The writer validates and converts it to the zero-based indexing required by
XDMF 3.

Attribute values are supplied as a rank-one array in Fortran column-major
order.
The expected shape is fixed when the attribute is defined and is validated on
every write.

Structured dimensions are supplied in I/J/K order.
The generated HDF5 and XDMF metadata expose dimensions in K/J/I order, with
the series axis first and the component axis last where present.

## Reader Compatibility

The conformance suite checks all generated metadata and heavy data directly.
When `pvpython` is available, it also loads representative uniform,
curvilinear, unstructured, mixed, and temporal outputs through ParaView's
XDMF reader.

Ready-to-open examples of those outputs are committed under
[`examples/generated`](examples/generated/README.md).
Keep each `.xdmf` file beside its corresponding `.h5` file when opening it in
ParaView.

Some ParaView releases still use the legacy XDMF2 reader internally.
That reader ignores edge- and face-centred attributes and can fail on valid
higher-rank `Matrix` attributes.
Scalar, vector, tensor, topology, and temporal data remain independently
validated even when a reader lacks support for one XDMF feature.

## Schema Stability

Every HDF5 file carries `schema_name` and `schema_version` root attributes.
Version 1.0 stores immutable grid data under `/grids`, attribute values under
`/attributes`, and series coordinates under `/series/values`.
Display names are kept in XDMF metadata and never become HDF5 paths.

Backward-incompatible storage changes require a new major schema version.
