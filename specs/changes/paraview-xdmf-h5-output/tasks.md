# Tasks: ParaView-Compatible XDMF/HDF5 Output

## Phase 1: Setup & Foundations

- [ ] [T01] [US1] Add failing XDMF API unit tests for sparse selected-point
  coordinates, relative HDF5 references, scalar `T x N` arrays, and vector
  `T x N x 3` arrays in `test/output/test_xdmfAPI.F90`.
- [ ] [T02] [US6] Add failing XDMF API unit tests for frequency sequence
  datasets, `F x N` scalar arrays, `F x N x 3` vector arrays, and
  `freq_1.000e6Hz`/`freq_1.000e-3Hz` names in
  `test/output/test_xdmfAPI.F90`.
- [ ] [T03] [US5] Add failing output utility tests for MPI probe folder names
  that preserve original x/y/z coordinate order and original field prefixes in
  `test/output/test_output_utils.F90`.

## Phase 2: Domain & Application Layer

- [ ] [T04] [US1] Add selected-point output metadata fields for sparse
  coordinate count, coordinate ordering, value-array shape, and vector packing
  in `src_output/outputTypes.F90`.
- [ ] [T05] [US3] Add component-name and vector-attribute naming helpers for
  current density, electric field, and magnetic field in
  `src_output/outputUtils.F90`.
- [ ] [T06] [US6] Add a frequency-entry name formatter that emits
  `freq_<scientific-value>Hz` with three decimal digits, lowercase `e`, no
  positive plus sign, and no exponent leading zeroes in
  `src_output/outputUtils.F90`.
- [ ] [T07] [US5] Change MPI coordinate extension and field-prefix helpers to
  keep original x/y/z coordinate order and original field prefixes in
  `src_output/outputUtils.F90`.

## Phase 3: Infrastructure & Adapters

- [ ] [T08] [US1] Extend `src_output/xdmfAPI.F90` with HDF5 writers for sparse
  `coordsX`, `coordsY`, `coordsZ`, `T x N`, `F x N`, and `* x N x 3` datasets.
- [ ] [T09] [US1] Extend `src_output/xdmfAPI.F90` with XDMF Polyvertex geometry
  writers that reference sparse coordinate arrays by relative HDF5 path from
  the XDMF folder.
- [ ] [T10] [US4] Extend `src_output/xdmfAPI.F90` with scalar and vector
  node-centred attribute writers for selected-point datasets.
- [ ] [T11] [US2] Update time-domain movie output creation, append, and XDMF
  publication to write one logical sparse XDMF/HDF5 pair with continuous
  `times`, scalar component arrays, and full-vector arrays in
  `src_output/movieProbeOutput.F90`.
- [ ] [T12] [US6] Update frequency-domain output creation and publication to
  write one logical sparse XDMF/HDF5 pair with frequency sequence states and
  named frequency entries in `src_output/frequencySliceProbeOutput.F90`.

## Phase 4: MPI Merge & Register

- [ ] [T13] [US5] Add rank-local visualisation output collection and
  deterministic merged selected-point ordering by ascending global x, then y,
  then z in `src_output/output.F90`.
- [ ] [T14] [US5] Publish one global output register entry per merged logical
  probe output instead of rank-local visualisation artifacts in
  `src_output/output.F90`.
- [ ] [T15] [US5] Apply user-selectable rank-local XDMF/HDF5 retention and
  failure-safe cleanup after merged output readability succeeds in
  `src_output/output.F90`.

## Phase 5: Integration & Validation

- [ ] [T16] [US1] Implement and test generated XDMF reference validation for
  missing HDF5 files, missing datasets, and invalid relative paths in
  `src_output/xdmfAPI.F90` and `test/output/test_xdmfAPI.F90`.
- [ ] [T17] [US1] Implement and test ParaView 5 readability acceptance records
  for major version, opened XDMF path, missing-file state, missing-dataset
  state, and manual path edit requirement in `src_output/outputTypes.F90`,
  `src_output/xdmfAPI.F90`, and `test/output/test_xdmfAPI.F90`.
- [ ] [T18] [US2] Add movie-output integration tests covering repeated flushes,
  continuous time-state append, and unchanged prior states in
  `test/output/test_output.F90`.
- [ ] [T19] [US6] Add frequency-output integration tests covering frequency
  sequence states, named entries, and frequency/time distinction in
  `test/output/test_output.F90`.
- [ ] [T20] [US5] Add MPI merge-order and rank-local retention tests for merged
  logical outputs in `test/output/test_output.F90`.

## Phase 6: Documentation & Spec Sync

- [ ] [T21] [US1] Update canonical output feature specs with the accepted
  ParaView-compatible XDMF/HDF5 contract in
  `specs/features/structured-export-contracts/spec.md`.
- [ ] [T22] [US2] Update time-domain volumetric output canonical behaviour in
  `specs/features/volumetric-time-output/spec.md`.
- [ ] [T23] [US6] Update frequency-domain volumetric output canonical behaviour
  in `specs/features/volumetric-frequency-output/spec.md`.

## Dependency Graph

T01 -> T04 -> T08 -> T09 -> T10 -> T11 -> T16 -> T17
T02 -> T06 -> T08 -> T12 -> T19
T03 -> T07 -> T13 -> T14 -> T15 -> T20
T05 -> T10 -> T11
T11 -> T18
T16 -> T21 -> T22 -> T23

## MVP Scope

The minimum deployable slice is T01, T04, T08, T09, T10, T11, T16, and T18.
That slice makes time-domain selected-point movie outputs produce one logical
ParaView-readable sparse XDMF/HDF5 pair with stable coordinates, scalar
components, vector attributes, and repeated-flush append semantics.

Frequency-domain output, MPI merge behaviour, retention policy, and canonical
spec sync follow after the MVP slice.
