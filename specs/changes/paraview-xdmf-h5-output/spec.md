# Spec: ParaView-Compatible XDMF/HDF5 Output

## Problem

Simulation users need XDMF/HDF5 volumetric outputs that open reliably in
ParaView 5 without manual conversion, renaming, or guessing dataset meaning.
The current output requirements identify some file names and dataset names,
but they do not fully define the generated file contract from a ParaView user's
point of view.

This feature is needed so time-domain and frequency-domain volumetric probe
outputs have a stable, documented, and testable shape for visual inspection and
post-processing.

## Goals

- A user can open each generated XDMF file in ParaView 5 and see the matching
  HDF5 data without manually editing paths or dataset names.
- Each probe's XDMF/HDF5 files are stored inside that probe's output folder.
- Time-domain selected-point outputs appear as a single continuous time series,
  regardless of how many times the solver flushed during the run.
- Spatial coordinates, time values, and field values have stable names,
  dimensions, and meanings.
- Field and current-density components are visible as node-centred scalar
  attributes in ParaView 5, and full vector quantities are also visible as
  node-centred vector attributes.
- Frequency-domain selected-point outputs expose one inspectable result per
  requested frequency both as a frequency sequence and as separately named
  frequency entries.
- MPI runs produce ParaView-readable merged logical XDMF/HDF5 outputs after
  rank-local output has been collected.
- MPI merged logical outputs order rank-local data by global physical grid
  coordinates, independent of rank number.
- MPI runs let users choose whether rank-local XDMF/HDF5 files are retained
  after the merged logical output is produced.

## Non-Goals

- Changing the numerical field, current-density, or frequency values produced
  by the solver.
- Changing which grid points are selected by a probe.
- Defining non-visualisation binary stream files.
- Supporting ParaView versions older than version 5.
- Requiring users to understand rank-local partial output files for normal
  visualisation.

## User Stories

- As a simulation user, I want to open a probe's XDMF file in ParaView 5 so
  that I can inspect volumetric results without manual conversion. [P1]
- As a simulation user, I want time-domain outputs to appear as one continuous
  time series so that flush cadence does not change what I see. [P1]
- As a simulation user, I want field components to have stable display names so
  that saved ParaView states and post-processing scripts remain reusable. [P1]
- As a simulation user, I want full vector quantities to be available as vector
  attributes so that ParaView 5 vector tools work without manual assembly. [P1]
- As a simulation user running MPI jobs, I want one merged logical output per
  probe so that I do not need to load rank-local pieces manually. [P1]
- As a simulation user, I want frequency-domain outputs to expose requested
  frequencies clearly so that I can inspect frequency response data. [P2]

## Success Criteria

- For every time-domain selected-point probe, the system creates exactly one
  XDMF file and one matching HDF5 file inside the probe folder for the logical
  probe output.
- For every frequency-domain selected-point probe, the system creates exactly
  one XDMF file and one matching HDF5 file inside the probe folder for the
  logical probe output.
- Opening the generated XDMF file in ParaView 5 displays the selected-point
  dataset without missing-file errors, missing-dataset errors, or manual path
  edits.
- The minimum ParaView 5 readability acceptance check is that ParaView 5 opens
  the generated XDMF file without missing-file errors, missing-dataset errors,
  or manual path edits.
- The XDMF file references its HDF5 file by a path that is valid from the XDMF
  file's own folder.
- Time-domain HDF5 data includes three one-dimensional spatial coordinate
  arrays named `coordsX`, `coordsY`, and `coordsZ`.
- The coordinate arrays `coordsX`, `coordsY`, and `coordsZ` have equal length
  `N`, where `N` is the number of selected probe points.
- The coordinate arrays align by index: for every point index `p`, the selected
  point coordinate is `(coordsX[p], coordsY[p], coordsZ[p])`.
- Selected point indexes are ordered by ascending global physical coordinates:
  first x, then y, then z.
- Time-domain HDF5 data includes one one-dimensional time array named `times`.
- The number of time values equals the number of time states visible in
  ParaView 5.
- Time-domain scalar value arrays have dimensions `T x N`, where `T` is the
  number of time states and `N` is the number of selected probe points.
- Time-domain scalar value `values[t, p]` belongs to time state `times[t]` and
  selected point `(coordsX[p], coordsY[p], coordsZ[p])`.
- Time-domain vector value arrays have dimensions `T x N x 3`, with component
  order x, y, z.
- Time-domain vector value `vectors[t, p, 0..2]` belongs to time state
  `times[t]` and selected point `(coordsX[p], coordsY[p], coordsZ[p])`.
- Repeated solver flushes append to the logical time series; they do not create
  duplicate time states, reset earlier time states, or require multiple XDMF
  files for one logical probe.
- Time-domain field values are stored so that ParaView 5 displays requested
  components as node-centred scalar attributes.
- Time-domain full-vector requests also expose one node-centred vector
  attribute assembled from the x, y, and z components.
- Current-density outputs expose x, y, and z components when all components are
  requested, and expose only the requested component when a single component is
  requested.
- Electric-field outputs expose x, y, and z components when all components are
  requested, and expose only the requested component when a single component is
  requested.
- Magnetic-field outputs expose x, y, and z components when all components are
  requested, and expose only the requested component when a single component is
  requested.
- Current-density, electric-field, and magnetic-field full-vector outputs each
  expose one vector attribute in addition to their x, y, and z scalar component
  attributes.
- Component names remain stable across runs for the same probe type and
  requested quantity.
- Only selected probe points are represented in ParaView 5-readable outputs;
  unselected points are omitted rather than zero-filled, NaN-filled, or masked.
- Sparse point outputs remain loadable in ParaView 5 without requiring a full
  rectangular volume.
- For MPI runs, users can open the merged logical XDMF output for each probe
  without manually selecting rank-local files.
- For MPI runs, merged coordinate and value data is ordered by ascending global
  physical coordinates: first x, then y, then z.
- For MPI runs, merged output ordering does not depend on MPI rank number,
  rank-local file discovery order, or output-register entry order.
- For MPI runs, rank-local files may exist, but they are not the primary user
  artifact for ParaView 5 inspection.
- For MPI runs, one global output register identifies the merged logical output
  for each probe.
- For MPI runs, a user-selectable retention setting determines whether
  rank-local XDMF/HDF5 files are kept or removed after a successful merge.
- If a user chooses to retain rank-local XDMF/HDF5 files, the merged logical
  output remains the primary ParaView 5 artifact.
- If a user chooses to remove rank-local XDMF/HDF5 files, removal occurs only
  after the merged logical output has been created successfully.
- Frequency-domain outputs provide a ParaView 5-readable frequency sequence
  containing one state per requested frequency for every selected probe point.
- Frequency-domain scalar value arrays have dimensions `F x N`, where `F` is
  the number of requested frequency states and `N` is the number of selected
  probe points.
- Frequency-domain vector value arrays have dimensions `F x N x 3`, with
  component order x, y, z.
- Frequency-domain value index `f` belongs to the `f`th requested frequency and
  selected point index `p` follows the same coordinate alignment as time-domain
  selected-point output.
- Frequency-domain outputs also provide separately named frequency entries for
  every requested frequency.
- Separately named frequency entries use scientific notation in hertz with the
  prefix `freq_` and suffix `Hz`, for example `freq_1.000e6Hz`.
- Frequency names use a dot as the decimal separator, exactly three digits after
  the decimal point, lowercase `e`, no plus sign for positive exponents, no
  leading zeroes in the exponent, and a minus sign for negative exponents.
- Frequency value `1000000` is named `freq_1.000e6Hz`; frequency value `0.001`
  is named `freq_1.000e-3Hz`.
- Frequency-domain output clearly distinguishes frequency values from physical
  simulation time values in the generated visualisation contract.
- The generated output contract is stable enough that automated tests can
  verify file presence, dataset names, dataset dimensions, XDMF references, and
  ParaView 5 readability assumptions.
- ParaView 5 readability acceptance records the ParaView major version, the
  XDMF file path opened, whether the file opened without missing-file errors,
  whether it opened without missing-dataset errors, and whether manual path
  edits were required.
- A ParaView 5 readability acceptance record passes only when the major version
  is `5`, the opened XDMF path is the generated logical probe XDMF path, no
  missing-file errors occur, no missing-dataset errors occur, and manual path
  edits are not required.
- The minimum readability check does not require visual rendering, saved-state
  reload, or verification that every expected array is visible in the ParaView
  user interface.

## Open Questions

- No open questions remain.
