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
| CPML borders | `src_main_pub/borderscpml.cpp`: placeholder constants and incomplete routines; active `semba_fdtd.cpp` has translated CPML border coefficient/state hooks and E/H updates, but the active region/sweep alignment is not byte-identical to Fortran yet | `src_main_pub/borderscpml.F90` | required | Finish matching Fortran `PMLc(...)`, `SINPML_fullsize(...)`, and full-size padded sweep behavior. PML must not fall back to MUR. |
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
| HDF5 XDMF writer | `src_main_pub/xdmf_h5.cpp`: real HDF5 C API (`openh5file` / `writeh5file` / `closeh5file`), linked in `semba-main` when `SEMBA_FDTD_ENABLE_HDF=ON` | `src_main_pub/xdmf_h5.F90` | required | Movie sampling in `semba_fdtd.cpp`; `test_movie_in_planewave_in_box` and `test_hdf_movie_dataset_shape_matches_fortran`. |
| VTK writer | `src_main_pub/vtk.cpp`: implementation placeholders | `src_main_pub/vtk.F90` | required | Translate for map/VTK output tests. |
| Map VTK writer | `src_main_pub/mapvtk_writer.cpp`: JSON-driven `-mapvtk` in active exe | `src_main_pub/vtk.F90` (`what==mapvtk`) | required | Byte parity: `test_mapvtk_vtk_file_matches_fortran_exact` on `pec_volume.fdtd.json`. |
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
| Conformal mapping | `src_cpp/src_conformal/conformal.cpp`: `addFaceRatio`, `addEdgeRatio`, `buildCellMap`, `buildSideMap` | `src_conformal/conformal.F90`, `cell_map.F90` | required | Independent package. Affects conformal/curved geometry tests. |
| JSON mesh helpers | `src_json_parser/mesh_m.h`: `allocateCoordinates`, `allocateElements` | JSON parser Fortran equivalents | artifact | Resolved as C++ `unordered_map::reserve()` calls. No Fortran translation needed because storage is map-backed, not array-backed. |
| MTLN utilities | `src_cpp/src_mtln/utils.cpp`, `src_cpp/src_mtln/preprocess.cpp`, MTLN headers with default empty vectors | `src_mtln/*.F90` | defer | Real migration, but keep as a separate MTLN package to avoid interference with non-MTLN work. |

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
- [x] Wire the migrated depth-zero SGBC/Maloney update into the standalone
      solver timestep (`advanceSgbcE` after wire E, `advanceSgbcH` before
      planewave H).
- [ ] Replace the standalone depth-zero SGBC node builder with the full
      translated `InitSGBCs`/`depth` path once the generated C++ modules share
      the executable's real field/media types.
- [ ] Translate and enable nonzero-depth SGBC, Crank-Nicolson SGBC, and
      dispersive SGBC internals from `maloney_nostoch.F90`.
- [x] Add hard failing C++ TODO tests for the remaining Maloney gaps:
      `MaloneyMissing.InitSgbcsBuildsFullFortranSurfaceState`,
      `MaloneyMissing.NonZeroDepthAdvancesInternalOneDimensionalFields`,
      `MaloneyMissing.CrankNicolsonSgbcMatchesFortranSheetSolve`, and
      `MaloneyMissing.DispersiveSgbcUpdatesPolarizationState`.
- [x] Add a hard failing full-system SGBC frequency sweep,
      `test_sgbc_surface_impedance_frequency_points_require_maloney_depth`,
      which samples a raw probe at several frequencies so depth-zero Maloney
      cannot pass by matching only the low-frequency response.
- [x] Add a hard failing byte-exact Fortran-vs-C++ probe comparison,
      `test_sgbc_surface_impedance_raw_probe_matches_fortran_exact`, for the
      same raw surface-impedance probe.
- [x] Stop aliasing JSON `"pml"` boundaries to Mur in the active standalone
      solver.
- [x] Add active timestep PML call sites in Fortran order:
      `advancePmlE()` after wire E, and `advancePmlBodyH()` plus
      `advanceMagneticCpml()` immediately after magnetic Yee H.
- [x] Add PML boundary classification and timestep-hook tests:
      `PmlBoundary.*`.
- [x] Add a hard failing byte-exact PML regression,
      `test_planewave_pml_probe_files_match_fortran_exact`, comparing Fortran
      and C++ probe files for a small planewave case with PML boundaries.
- [x] Replace the current PML no-op hooks with translated CPML border state
      initialization and E/H updates from `borderscpml.F90`.
- [ ] Match active CPML region/sweep bounds to Fortran `PMLc(...)` and
      `SINPML_fullsize(...)` exactly enough for byte-identical probe files.
- [ ] Translate PML body state initialization and H updates from
      `pml_bodies.F90`.
- [x] Re-run `test_sgbc_shielding_effectiveness`,
      `test_sgbc_structured_resistance_single_wire`,
      `test_pec_overlapping_sgbcs`, and `test_sgbc_overlapping_sgbc`.
- [ ] Move to the planewave TF/SF work package after SGBC has strict tests.

## Active Timestepping Coverage

The active C++ executable currently uses the standalone `timestepping()`
wrapper in `src_cpp/src_main_pub/semba_fdtd.cpp`, not the generated
`src_cpp/src_main_pub/timestepping.cpp`.  The wrapper keeps separate E and H
phase helpers ordered after `src_main_pub/timestepping.F90`.

- Present in the active loop: Yee E/H updates, MTLN E when enabled, Holland wire
  E, translated CPML border E/H hooks, PML-body H placeholder, depth-zero SGBC/Maloney E/H, lumped E, PEC E/H,
  planewave E/H, nodal E, magnetic PMC/periodic cloning, Mur H, and probe
  sampling.
- TODO: align active CPML bounds and padded-domain sweeps with Fortran exactly,
  then translate PML-body H if needed by cases with PML media bodies.
- TODO: migrate electric and magnetic dispersive material E/H phases.
- TODO: migrate Fortran multiport phases beyond the MTLN hook.
- TODO: migrate nodal H if cases require magnetic nodal source behavior.
- TODO: migrate full wire H/Crank paths from `src_wires_pub/wires.F90`; the
  active solver has only the current standalone Holland E path.
- TODO: replace depth-zero SGBC with full Maloney internals once the generated
  SGBC module shares the active solver field/media types.

## Verification Log

- `/tmp/semba_cpp_mtln_check/bin/cpp_tests`: 187 passed, 1 skipped
  (`BordersMur.PulseMatchesFortranProbe`, missing golden probe).
- After adding the hard failing Maloney TODO tests,
  `/tmp/semba_cpp_mtln_check/bin/cpp_tests '--gtest_filter=MaloneyMissing.*'`
  is expected to fail until full SGBC/Maloney translation is completed.
- The hard failing full-system SGBC frequency sweep
  `test_sgbc_surface_impedance_frequency_points_require_maloney_depth` is also
  expected to fail until the nonzero-depth/internal Maloney propagation is
  translated.  Current depth-zero output is nearly flat from 50 MHz to 900 MHz
  (`-45.55 dB` to `-45.83 dB`) while the analytical slab response falls to
  about `-73.0 dB`.
- The byte-exact migration regression
  `test_sgbc_surface_impedance_raw_probe_matches_fortran_exact` currently fails
  against C++ at the first differing raw probe sample, confirming the active
  time-domain SGBC output is not yet Fortran-identical.
- `/tmp/semba_cpp_mtln_check/bin/cpp_tests '--gtest_filter=PmlBoundary.*'`:
  4 passed.  These tests verify PML no longer enables Mur and that the active
  timestep reaches the PML E/H hook functions.
- `test_planewave_pml_probe_files_match_fortran_exact` is still a strict xfail.
  Active CPML border updates are translated, but the forced run first differs
  at line 21 of `pw-in-box.fdtd_before_Ex_3_3_1.dat`
  (`0.235246427E-033` Fortran vs `0.237195566E-033` C++).  The likely remaining
  work is exact Fortran `PMLc(...)`/padded-sweep alignment, not MUR fallback.
- Focused SGBC Python checks passed:
  `test_sgbc_shielding_effectiveness`,
  `test_sgbc_structured_resistance_single_wire`,
  `test_pec_overlapping_sgbcs`, `test_sgbc_overlapping_sgbc`,
  `test_can_assign_same_surface_impedance_to_multiple_geometries`, and
  `test_one_cell_SGBC_surface_Jprobe`.
- `pytest test/pyWrapper/test_integration.py -m sgbc`: 3 passed, 50 deselected.
- `pytest -q test/pyWrapper`: 93 passed, 37 skipped.
- Active executable path now initializes depth-zero SGBC nodes from
  `multilayeredSurface` material associations and calls the migrated
  `SGBC_nostoch_m::AdvanceSGBCE/H` depth-zero routines in the timestep loop.
  This is intentionally not the full Fortran Maloney translation yet:
  nonzero-depth, Crank-Nicolson, and dispersive SGBC remain open items.
