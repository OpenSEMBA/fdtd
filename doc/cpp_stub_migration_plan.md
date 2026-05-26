# C++ Stub Migration Plan

This file tracks the C++ migration stubs found under `src_cpp`.  The goal is to
decide, module by module, whether a stub is only a temporary compatibility
artifact or whether it must be translated from the original Fortran source.

## Classification

- `required`: a Fortran implementation exists and the C++ behavior is missing or
  materially simplified.
- `artifact`: the code only exists to make a standalone C++ target compile and
  is currently replaced by another implemented path.
- `defer`: real work, but independent from the non-MTLN full-system failures or
  tied to MPI/MTLN-specific work.

## Execution Order

1. `src_cpp/src_main_pub/semba_fdtd.cpp`: keep stabilizing the standalone C++
   solver used by the Python full-system suite. This is linked through the
   original-named `semba-main` C++ target.
2. `src_cpp/src_main_pub/maloney_nostoch.cpp`: translate the SGBC/Maloney
   surface-impedance module from `src_main_pub/maloney_nostoch.F90`.
3. `src_cpp/src_main_pub/planewaves.cpp`: translate the original TF/SF
   planewave module from `src_main_pub/planewaves.F90`, then remove the
   duplicated simplified planewave logic from the standalone solver.
4. `src_cpp/src_main_pub/observation.cpp` and output writers: make probe/movie
   files byte-for-byte compatible with Fortran where possible.
5. Remaining material, geometry, boundary, wire, conformal, MPI, and MTLN
   modules.

## Stub Inventory

| Area | C++ file and stubs | Fortran source | Class | Decision and next action |
| --- | --- | --- | --- | --- |
| SGBC surface impedance | `src_main_pub/maloney_nostoch.cpp`: `InitSGBCs`, `calc_SGBCconstants`, `AdvanceSGBCE`, `AdvanceSGBCH`, `StoreFieldsSGBCs`, `DestroySGBCs`, `calc_g1g2gm1gm2_compo`, `g1g2`, `gm1gm2`, `g1g2_Dispersive`, `depth`, tridiagonal solvers, `test_stab` | `src_main_pub/maloney_nostoch.F90` | required | Translate from Fortran. Current SGBC full-system failures are consistent with this module being empty or bypassed. Start with scalar helpers and depth logic, then initialization and time advance. |
| Planewaves and TF/SF | `src_main_pub/planewaves.cpp`: mocked file reads and incomplete update helpers; standalone `semba_fdtd.cpp` has a duplicate simplified planewave path and an empty `applyPlaneWaveSource` | `src_main_pub/planewaves.F90`, `timestepping.F90` | required | Translate original planewave initialization and E/H corrections. Use strict byte tests for PEC/MUR/periodic boundaries. |
| Bulk-current observation | `src_main_pub/observation.cpp`: public-interface stub comments; standalone bulk-current implementation lives in `semba_fdtd.cpp` | `src_main_pub/observation.F90` | required | Migrate exact Fortran orientation and output naming. Keep the existing strict tests for probe bytes and add cases for each current orientation. |
| Classic wires | `src_wires_pub/wires.cpp`: `WarnErrReport`, `calc_wirehollandconstants`, `AdvanceWiresE`, `AdvanceWiresH`, `AdvanceWiresEcrank`, `StoreFieldsWires`, `DestroyWires`, `GetHwires`, `ReportWireJunctions` | `src_wires_pub/wires.F90` | required | Current standalone solver has a partial Holland path. Full migration still requires the original wire module. Coordinate with MTLN work only where wire/MTLN ownership overlaps. |
| Lumped elements | `src_main_pub/lumped.cpp`: current executable-gate update logic plus remaining restart/file-I/O gaps | `src_main_pub/lumped.F90` | required | Continue translating the full original module while keeping lumped tests strict against Fortran outputs. |
| CPML borders | `src_main_pub/borderscpml.cpp`: placeholder constants and incomplete routines | `src_main_pub/borderscpml.F90` | required | Needed for PML/CPML full-system tests and exact boundary behavior. Independent from SGBC if tests use MUR/PEC. |
| Other borders | `src_main_pub/bordersother.cpp`: minimal type stubs | `src_main_pub/bordersother.F90` | required | Translate periodic/PEC/PMC helpers and compare against current standalone boundary code. |
| MUR borders | `src_main_pub/bordersmur.cpp`: placeholder file-read hook | `src_main_pub/bordersmur.F90` | required | Mostly implemented in standalone solver; finish by migrating original serialization/restart hooks. |
| PML bodies | `src_main_pub/pml_bodies.cpp`: placeholder branch for material/body data | `src_main_pub/pml_bodies.F90` | required | Translate after CPML, because both touch absorbing regions. |
| Dispersive electric media | `src_main_pub/electricdispersive.cpp`: placeholder types and error handling | `src_main_pub/electricdispersive.F90` | required | Translate for dispersive-material cases. Independent from planewave/SGBC except when a case combines them. |
| Dispersive magnetic media | `src_main_pub/magneticdispersive.cpp`: structural stubs | `src_main_pub/magneticdispersive.F90` | required | Translate after electric dispersive media to reuse patterns. |
| Anisotropic media | `src_main_pub/anisotropic.cpp`: placeholder dependency definitions | `src_main_pub/anisotropic.F90` | required | Translate when anisotropic tests are enabled. |
| NFDE rotation | `src_main_pub/nfde_rotate.cpp`, `nfde_rotate.h`: simplified rotation/permutation and unsupported anisotropic rotation cases | `src_main_pub/nfde_rotate.F90` | required | Translate exact tensor rotation before anisotropic tests are made strict. |
| Geometry healing | `src_main_pub/healer.cpp`: file ends with placeholder implementation | `src_main_pub/healer.F90` | required | Large, shared preprocessing surface. Do this in a separate branch/package because it can affect every material assignment test. |
| Geometry preprocessing | `src_main_pub/preprocess_geom.cpp`: MPI barrier placeholder and visibly malformed generated sections | `src_main_pub/preprocess_geom.F90` | required | Required for full original executable path. Keep independent from standalone solver fixes. |
| Constants | `src_main_pub/calc_constants.cpp`: SGBC comments mirror the Fortran deferral to `InitSGBCs` | `src_main_pub/calc_constants.F90` | required | Revisit after SGBC translation to ensure material coefficients are not duplicated or stale. |
| Postprocess | `src_main_pub/postprocess.cpp`: simplified DTFT placeholder and generic I/O comments | `src_main_pub/postprocess.F90` | required | Translate after time-domain probes are exact. |
| Far field | `src_main_pub/farfield.cpp`: local `stoponerror`/`print11` empty stubs and placeholder for snippet-start code | `src_main_pub/farfield.F90` | required | Existing sphere/farfield tests use the standalone path, but full module migration remains required. |
| XDMF/observation facade | `src_main_pub/xdmf.cpp`: dummy `SGGFDTDINFO_t`, dummy output structures, and many no-op `createxdmf*`/observation calls | `src_main_pub/xdmf.F90` | artifact now, required later | Artifact for standalone compilation today. Replace with real translation before enabling the full `semba-outputs` path. |
| HDF5 XDMF writer | `src_main_pub/xdmf_h5.cpp`: comment says HDF5 calls are stubbed | `src_main_pub/xdmf_h5.F90` | required | HDF5 movie tests currently pass through `semba_fdtd.cpp`; migrate original writer for full output compatibility. |
| VTK writer | `src_main_pub/vtk.cpp`: implementation placeholders | `src_main_pub/vtk.F90` | required | Translate for map/VTK output tests. |
| Map VTK writer | `src_main_pub/mapvtk_writer.cpp`: no real stub found, mostly implemented | none direct | artifact | Keep as standalone support code unless it diverges from Fortran map output. |
| Error reporting core | `src_main_pub/errorreport_core.cpp`: empty `reportmedia`, help/default switch routines | `src_main_pub/errorreport.F90`, `interpresenta_switches.F90` | required | Low physics risk, but needed for CLI/report parity. |
| Error reporting full module | `src_main_pub/errorreport.cpp`: helper and MPI placeholders | `src_main_pub/errorreport.F90` | required | Translate or consolidate with `errorreport_core.cpp`; avoid two divergent reporting implementations. |
| CLI argument parsing | `src_main_pub/getargs.cpp`: placeholder argument string handling | `src_main_pub/getargs.F90` | required | Needed when C++ executable stops depending only on JSON wrapper calls. |
| Switch parsing | `src_main_pub/interpreta_switches.cpp`: external stubs | `src_main_pub/interpreta_switches.F90` | required | Translate with CLI parsing. |
| Resume/restart | `src_main_pub/resuming.cpp`: placeholder stream and stubbed environment assumptions | `src_main_pub/resuming.F90` | defer | Real migration, but independent unless restart tests are enabled. |
| MPI communication | `src_main_pub/mpicomm.cpp`: cut-off placeholder implementation and MPI shim comments | `src_main_pub/mpicomm.F90` | defer | MPI-specific package; keep separate from single-process full-system work. |
| Store geometry | `src_main_pub/storegeom.cpp`: placeholder type definitions | `src_main_pub/storegeom.F90` | required | Needed for report/output parity, not core timestepping. |
| Snapshot XDMF | `src_main_pub/snapxdmf.cpp`: placeholder constants/types | `src_main_pub/snapxdmf.F90` | required | Translate with XDMF/HDF output package. |
| Thin slot | `src_main_pub/dmma_thin_slot.cpp`: placeholder constants/types | `src_main_pub/dmma_thin_slot.F90` | required | Translate when DMMA/thin-slot tests are in scope. |
| Timestepping full module | `src_main_pub/timestepping.cpp`: forward declarations for methods not implemented in snippet | `src_main_pub/timestepping.F90` | required | Needed before `semba-main` can switch from the current executable-gate loop to the full translated Fortran time loop. |
| Conformal mapping | `src_conformal/conformal.cpp`: `addFaceRatio`, `addEdgeRatio`, `buildCellMap`, `buildSideMap` | `src_conformal/conformal.F90`, `cell_map.F90` | required | Independent package. Affects conformal/curved geometry tests. |
| JSON mesh helpers | `src_json_parser/mesh_m.h`: `allocateCoordinates`, `allocateElements` | JSON parser Fortran equivalents | artifact | Resolved as C++ `unordered_map::reserve()` calls. No Fortran translation needed because storage is map-backed, not array-backed. |
| MTLN utilities | `src_mtln/utils.cpp`, `src_mtln/preprocess.cpp`, MTLN headers with default empty vectors | `src_mtln/*.F90` | defer | Real migration, but keep as a separate MTLN package to avoid interference with non-MTLN work. |

## Current Checklist

- [x] Inventory obvious C++ stubs and placeholder comments under `src_cpp`.
- [x] Classify current standalone artifacts versus missing Fortran translations.
- [x] Remove dead standalone `applyPlaneWaveSource()` artifact; TF/SF updates are
      already handled by `advancePlaneWaveE/H`.
- [x] Resolve JSON parser allocation artifacts with map `reserve()` calls.
- [x] Translate SGBC scalar helpers in `maloney_nostoch.cpp`.
- [x] Add focused C++/Python checks for SGBC helper formulas against Fortran
      values.
- [x] Translate SGBC depth allocation and layer mapping.
- [ ] Attach the translated SGBC depth helper to `depth(...)` once the migrated
      C++ modules share one `SGGFDTDINFO_t`/`MediaData_t` definition instead of
      per-file generated stand-ins.
- [ ] Wire SGBC initialization into the standalone solver or switch the
      executable to the migrated component path once equivalent.
- [ ] Re-run `test_sgbc_shielding_effectiveness`,
      `test_sgbc_structured_resistance_single_wire`,
      `test_pec_overlapping_sgbcs`, and `test_sgbc_overlapping_sgbc`.
- [ ] Move to the planewave TF/SF work package after SGBC has strict tests.

## Verification Log

- `cpp_build_nomtln/bin/cpp_tests`: 149 passed, 1 skipped
  (`BordersMur.PulseMatchesFortranProbe`, missing golden probe).
- `pytest test/pyWrapper/test_full_system.py::test_sgbc_shielding_effectiveness`
  still fails. The C++ standalone solver currently behaves like the SGBC panel
  is transparent, so the next SGBC task is wiring the translated SGBC state and
  update routines into the executable path.
