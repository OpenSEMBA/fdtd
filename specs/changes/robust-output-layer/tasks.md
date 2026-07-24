# Tasks: Robust Output Layer

## Phase 1: Setup & Foundations

- [x] [T01] [US1] [US3] [US4] Define the versioned probe artifact, metadata, lifecycle, and binary-layout contracts, including explicit byte order, numeric representation, complex-value convention, and required artifacts per probe family in `src_output/outputTypes.F90`.
- [x] [T02] [US1] [US3] [US4] Add contract tests for declared, active, complete, failed, and zero-sample probe outputs in `test/output/test_output.F90` and register them in `test/output/CMakeLists.txt`.
- [x] [T03] [P] [US2] Extend component and rank-boundary ownership tests to prove disjoint complete volumetric coverage and zero-based offsets in `test/output/test_output_decomposition.F90`.
- [x] [T04] [P] [US2] Define serial-versus-distributed equivalence fixtures that compare artifact identities, coordinate coverage, sample counts, and duplicate boundary ownership in `test/output/test_output.F90`.
- [x] [T05] [US1] [US4] Document the metadata schema, relative-artifact-reference rule, lifecycle state transitions, and binary compatibility guarantees in `specs/changes/robust-output-layer/contracts/output-metadata-schema.md`.
- [ ] [T10] [US6] Replace shell-dependent folder and file operations with portable nested-directory creation and atomic replacement adapters in `src_utils/directoryUtils.F90`, `src_utils/CMakeLists.txt`, and `src_utils/filesystemAdapter.c`.

## Phase 2: Domain & Application Layer

- [x] [T06] [US1] [US2] [US4] Add probe-output lifecycle coordination, artifact declaration, terminal-state validation, and a root-owned run manifest model to `src_output/output.F90`.
- [x] [T07] [US2] Add per-probe participant selection and scalar writer ownership using the output contract rather than rank-derived file names in `src_output/output.F90` and `src_output/outputTypes.F90`.
- [x] [T08] [US2] Apply `build_output_partition` to every volumetric output request and retain each rank's global and local bounds in `src_output/output.F90` and `src_output/outputDecomposition.F90`.
- [ ] [T09] [US5] Route point, wire, bulk, far-field, and geometry probe publication through declared artifacts while preserving their existing human-readable payloads in `src_output/pointProbeOutput.F90`, `src_output/wireProbeOutput.F90`, `src_output/bulkProbeOutput.F90`, `src_output/farFieldProbeOutput.F90`, and `src_output/mapVTKOutput.F90`.

## Phase 3: Infrastructure & Adapters

- [ ] [T11] [US1] [US4] Implement initial and final machine-readable probe metadata publication, including relative artifact references and failure state, in `src_output/outputMetadata.F90` and add it to `src_output/CMakeLists.txt`.
- [ ] [T12] [US3] [US4] Implement the portable stream binary writer and metadata-backed layout validation in `src_output/outputBinary.F90`, `src_output/outputTypes.F90`, and `src_output/CMakeLists.txt`.
- [ ] [T13] [US3] Adapt the external visualisation writer behind a volumetric publication interface and link it into the output target in `src_output/outputVisualisation.F90`, `src_output/CMakeLists.txt`, and `external/CMakeLists.txt`.
- [ ] [T14] [US2] Implement the distributed output collective adapter for deterministic probe ownership, collective partition publication, and root aggregation fallback in `src_output/outputCollective.F90` and `src_output/CMakeLists.txt`.

## Phase 4: Integration & Validation

- [ ] [T15] [US3] Integrate the binary writer, visualisation adapter, and metadata lifecycle into time-domain volumetric probes in `src_output/movieProbeOutput.F90` and `src_output/output.F90`.
- [ ] [T16] [US3] Integrate the binary writer, visualisation adapter, and metadata lifecycle into frequency-domain volumetric probes, preserving complex values and correcting component record ordering, in `src_output/frequencySliceProbeOutput.F90` and `src_output/output.F90`.
- [ ] [T17] [US2] [US6] Integrate collective hyperslab publication and root aggregation fallback with the solver output lifecycle in `src_output/output.F90`, `src_output/movieProbeOutput.F90`, and `src_output/frequencySliceProbeOutput.F90`.
- [ ] [T18] [US1] [US2] [US4] Replace the stale per-rank output-request register with an artifact-derived root-owned run manifest in `src_output/output.F90` and remove obsolete registration handling from `src_main_pub/semba_fdtd.F90`.
- [ ] [T19] [US2] Retire the legacy observation dispatch and direct output writers only after the new path is covered, updating `CMakeLists.txt`, `src_main_pub/timestepping.F90`, `src_main_pub/observation.F90`, `src_main_pub/vtk.F90`, `src_main_pub/xdmf.F90`, and `src_main_pub/xdmf_h5.F90`.
- [ ] [T20] [US1] [US3] [US4] [US5] [US6] Add serial integration tests covering every probe family, zero-sample probes, nested paths, spaced paths, metadata descriptors, binary layout, and complete visualisation pairs in `test/output/test_output.F90`, `test/output/test_output_utils.F90`, `test/output/CMakeLists.txt`, and `external/xdmf-hdf5/test/verify_pairs.py`.
- [ ] [T21] [US2] Add distributed integration tests for collective and root-aggregation output modes, asserting serial equivalence and no duplicate interface data, in `test/output/test_output_mpi.F90`, `test/output/CMakeLists.txt`, and `CMakeLists.txt`.
- [ ] [T22] [US2] [US3] [US6] Add the output verification matrix to continuous integration for Linux single-process, Linux distributed collective, Linux distributed aggregation, and Windows single-process builds in `.github/workflows/ubuntu.yml` and `.github/workflows/windows.yml`.

## Dependency Graph

T01 → T02 → T06 → T07 → T09 → T18 → T19 → T20 → T22

T01 → T03 → T08 → T14 → T17 → T21 → T22

T01 → T04 → T21

T01 → T05 → T10 → T09

T10 → T11 → T18

T01 → T06 → T12 → T15 → T20

T01 → T06 → T13 → T15

T01 → T06 → T13 → T16 → T20

T14 → T17

## MVP Scope

The deployable minimum slice is T01, T02, T05 through T07, T09 through T13, T15, T16, T18, and T20.
It establishes a serial-safe output lifecycle, one metadata descriptor per probe, retained human-readable scalar outputs, complete volumetric binary and visualisation artifacts, portable path handling, and serial validation.
T03, T04, T08, T14, T17, T19, T21, and T22 extend that slice to distributed parity, legacy retirement, and the full platform verification matrix.
