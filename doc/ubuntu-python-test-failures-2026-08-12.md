# Ubuntu Python Test Failure Report

## Scope

This report covers the latest completed Ubuntu GitHub Actions workflow that ran Python tests at the time of investigation:

- Workflow: [ubuntu #31588351969](https://github.com/OpenSEMBA/fdtd/actions/runs/31588351969)
- Commit: `ed9a1e708217b2fea1f8da66321782f6494e47e8`
- Branch: `feature/conformal_thin_panel_pec`
- Test command: `python -m pytest test/ -n auto --durations=20`

## Remediation Tracker

| Priority | Failure family | Status | Verification target |
| --- | --- | --- | --- |
| 1 | B: map media-type contract | Residual | All targeted map checks pass except `test_wires_collision`, which publishes two tag-`128` lines rather than four. |
| 2 | C: duplicate movie geometry tags | Verified | `test_three_surfaces_Jprobe` and `test_wires_collision_Jprobe` pass. |
| 3 | A: output discovery | Verified | Representative MTLN and wire probe tests pass. |
| 4 | D: conformal cylinder validation | Verified | `test_conformal_impedance_cylinder_unshielded` passes. |
| 5 | E: current-source amplitude | Verified | `test_current_generators_without_resistance` passes. |

### Remediation Verification

Completed on the current branch:

- MTLN representative: `test_movie_in_planewave_in_box` passes after accepting the deliberate `tagnumber` movie attribute.
- No-MTLN representative: `test_conformal_impedance_cylinder_unshielded`, `test_current_generators_without_resistance`, `test_three_surfaces`, `test_three_surfaces_Jprobe`, and `test_wires_collision_Jprobe` pass.
- Residual: `test_wires_collision` has the correct map file and PEC entries, but tag `128` appears twice rather than four times. This is the only failure in the six-test targeted no-MTLN set.

Six matrix configurations reached pytest and all failed.
Three Intel configurations failed while building, so they have no Python-test result:

- Intel, MPI=OFF, MTLN=ON, double=OFF
- Intel, MPI=ON, MTLN=ON, double=OFF
- Intel, MPI=ON, MTLN=OFF, double=ON

## Failure Families

### A. Probe-artifact discovery returns no files

Affected tests fail with an empty `getSolvedProbeFilenames(...)` result or by trying to open a non-existent expected output file.

The most likely cause is the output-artifact discovery refactor in `src_pyWrapper/pyWrapper.py`.
The discovery method filters filenames rather than obtaining a solver-produced manifest, and the run shows that it does not recognise a broad set of wire and MTLN outputs.
This is a wrapper/output naming regression, not evidence that each affected physics simulation failed.

Relevant code: `FDTD.getSolvedProbeFilenames()` in `src_pyWrapper/pyWrapper.py:386`.
The regression is consistent with the recent `FDTD | Fix | Discover probe output artifacts` change (`baf9238c`).

### B. VTK map media types are missing or changed

Affected tests raise `KeyError` for expected media keys such as `0.5`, `16.5`, `12`, or `128`.

The map output no longer publishes the legacy values expected by the tests for PEC, PMC, wire, and multiwire segments.
`mapVTKOutput_m` now derives edge types from the material selected at the sampled electric-field location.
At collisions and conformal boundaries this can resolve to a different material, and thus omits the expected type.

Relevant code: `edge_media_type()` in `src_output/mapVTKOutput.F90:399`.
The problem is most likely introduced by the map-output changes `50a910bf` and `06436b2f`.

### C. Movie geometry tags are duplicated

Affected current-movie tests report exactly double the expected geometry-tag counts:

- `24` rather than `12`
- `8` rather than `4`
- `4` rather than `2`

This points to duplicate geometry publication or aggregation, rather than an incorrect tag value.
The recent movie-geometry feature writes companion geometry artifacts and combines tag arrays across all matching `tagnumber` attributes.

Relevant code: `current_movie_geometry_all_tag_counts()` in `test/pyWrapper/test_integration.py:38` and movie/map output code under `src_output/`.
The likely regression range starts with `c4334ed4` (`FDTD | Feature | Publish movie geometry tags`).

### D. Conformal cylinder is rejected during preprocessing

`test_conformal_impedance_cylinder_unshielded` stops before time stepping with:

```
Invalid conformal volume PEC@layer1: combined surface patches traverse a shared edge in the same direction
```

This is a real input-validation/preprocessing failure.
It is independent of compiler, MPI, and MTLN mode because it occurs in every configuration that reaches the test.
The likely cause is the conformal thin-panel PEC validation introduced on this branch, which rejects the existing `conformal_impedance_cylinder_conformal.fdtd.json` geometry or applies the shared-edge orientation rule incorrectly.

Relevant implementation area: `src_main_pub/preprocess_geom.F90`.

### E. Current-source amplitude is 100 times too small

In non-MTLN configurations, `test_current_generators_without_resistance` receives a steady current of `0.01` while the test expects `1.0`.

This is a physics/output value regression, not a discovery error.
It suggests that source magnitude normalisation or the current-generator implementation applies a factor of `1/100`.
The MTLN configurations do not reach the value assertion because their probe-artifact discovery fails first.

Relevant test: `test/pyWrapper/test_full_system.py:1806` in the run revision.

## Failures By Compilation Mode

### GCC, MPI=ON, MTLN=ON, HDF=ON, double=OFF

Job: [94087314615](https://github.com/OpenSEMBA/fdtd/actions/runs/31588351969/job/94087314615)
Result: 32 failed, 69 passed, 25 skipped, 1 xfailed.

| Failed tests | Likely cause |
| --- | --- |
| `test_shieldedPair`, `test_shieldedPair_mpi`, `test_coated_antenna`, `test_holland`, `test_holland_short_terminals_match_open_terminals`, `test_holland_mtln_mpi`, `test_unshielded_multiwires`, `test_towelHanger`, `test_towelHanger_mpi`, `test_towel_rack_with_and_without_shorting_plane`, `test_current_generators_with_resistance`, `test_current_generators_without_resistance`, `test_voltage_generators` | A: probe-artifact discovery does not find expected wire/MTLN outputs. |
| `test_conductors_forming_y_on_panel_holland_vs_mtln` | A: expected generated probe file is absent; likely the same discovery/output naming regression. |
| `test_holland_case_checking_number_of_outputs_single_wire`, `test_holland_case_checking_number_of_outputs_wire`, `test_holland_case_checking_number_of_outputs_unshielded` | A: discovery returns zero outputs. |
| `test_three_surfaces`, `test_1_line`, `test_volume_and_surfaces`, `test_multiwire_z_collision_y`, `test_multiwire_z_collision_x`, `test_multiwire_y_collision_x`, `test_multiwire_y_collision_z`, `test_multiwire_x_collision_y`, `test_multiwire_x_collision_z`, `test_multiwire_x_long_collision_z` | B: map media-type keys are absent or differ from the legacy contract. |
| `test_three_surfaces_Jprobe`, `test_three_one_cell_surfaces_Jprobe`, `test_one_cell_PEC_surface_Jprobe`, `test_one_cell_SGBC_surface_Jprobe` | C: current-movie geometry tag counts are doubled. |
| `test_conformal_impedance_cylinder_unshielded` | D: conformal shared-edge orientation validation rejects the case. |

### GCC, MPI=OFF, MTLN=ON, HDF=ON, double=OFF

Job: [94087314632](https://github.com/OpenSEMBA/fdtd/actions/runs/31588351969/job/94087314632)
Result: 29 failed, 64 passed, 33 skipped, 1 xfailed.

| Failed tests | Likely cause |
| --- | --- |
| `test_shieldedPair`, `test_coated_antenna`, `test_holland`, `test_holland_short_terminals_match_open_terminals`, `test_unshielded_multiwires`, `test_towelHanger`, `test_towel_rack_with_and_without_shorting_plane`, `test_current_generators_with_resistance`, `test_current_generators_without_resistance`, `test_voltage_generators` | A: probe-artifact discovery/output naming regression. |
| `test_conductors_forming_y_on_panel_holland_vs_mtln` | A: expected generated probe file is absent. |
| `test_holland_case_checking_number_of_outputs_single_wire`, `test_holland_case_checking_number_of_outputs_wire`, `test_holland_case_checking_number_of_outputs_unshielded` | A: discovery returns zero outputs. |
| `test_three_surfaces`, `test_1_line`, `test_volume_and_surfaces`, `test_multiwire_z_collision_y`, `test_multiwire_z_collision_x`, `test_multiwire_y_collision_x`, `test_multiwire_y_collision_z`, `test_multiwire_x_collision_y`, `test_multiwire_x_collision_z`, `test_multiwire_x_long_collision_z` | B: map media-type regression. |
| `test_three_surfaces_Jprobe`, `test_three_one_cell_surfaces_Jprobe`, `test_one_cell_PEC_surface_Jprobe`, `test_one_cell_SGBC_surface_Jprobe` | C: duplicate current-movie geometry tags. |
| `test_conformal_impedance_cylinder_unshielded` | D: conformal shared-edge orientation validation. |

### GCC, MPI=ON, MTLN=OFF, HDF=ON, double=OFF

Job: [94087314670](https://github.com/OpenSEMBA/fdtd/actions/runs/31588351969/job/94087314670)
Result: 19 failed, 71 passed, 36 skipped, 1 xfailed.

| Failed tests | Likely cause |
| --- | --- |
| `test_current_generators_without_resistance` | E: current is `0.01` rather than `1.0`. |
| `test_three_surfaces`, `test_1_line`, `test_volume_and_surfaces`, `test_wires_collision`, `test_wire_x_collision_y`, `test_wire_x_collision_z`, `test_wire_x_long_collision_z`, `test_wire_y_collision_x`, `test_wire_y_collision_z`, `test_wire_z_collision_x`, `test_wire_z_collision_y` | B: map media-type regression. |
| `test_three_surfaces_Jprobe`, `test_three_one_cell_surfaces_Jprobe`, `test_one_cell_PEC_surface_Jprobe`, `test_one_cell_SGBC_surface_Jprobe`, `test_wires_collision_Jprobe`, `test_wire_x_collision_y_Jprobe` | C: duplicate or extra movie geometry tags. |
| `test_conformal_impedance_cylinder_unshielded` | D: conformal shared-edge orientation validation. |

### GCC, MPI=OFF, MTLN=OFF, HDF=ON, double=OFF

Job: [94087314708](https://github.com/OpenSEMBA/fdtd/actions/runs/31588351969/job/94087314708)
Result: 19 failed, 67 passed, 40 skipped, 1 xfailed.

The failing-test set and causes are identical to GCC, MPI=ON, MTLN=OFF, double=OFF:

- E: `test_current_generators_without_resistance`.
- B: `test_three_surfaces`, `test_1_line`, `test_volume_and_surfaces`, `test_wires_collision`, `test_wire_x_collision_y`, `test_wire_x_collision_z`, `test_wire_x_long_collision_z`, `test_wire_y_collision_x`, `test_wire_y_collision_z`, `test_wire_z_collision_x`, `test_wire_z_collision_y`.
- C: `test_three_surfaces_Jprobe`, `test_three_one_cell_surfaces_Jprobe`, `test_one_cell_PEC_surface_Jprobe`, `test_one_cell_SGBC_surface_Jprobe`, `test_wires_collision_Jprobe`, `test_wire_x_collision_y_Jprobe`.
- D: `test_conformal_impedance_cylinder_unshielded`.

### Intel, MPI=ON, MTLN=OFF, HDF=ON, double=OFF

Job: [94087314618](https://github.com/OpenSEMBA/fdtd/actions/runs/31588351969/job/94087314618)
Result: 4 failed, 86 passed, 36 skipped, 1 xfailed.

| Failed tests | Likely cause |
| --- | --- |
| `test_conformal_impedance_cylinder_unshielded` | D: conformal shared-edge orientation validation. |
| `test_current_generators_without_resistance` | E: current is `0.01` rather than `1.0`. |
| `test_fill_conformal_vtk_sphere`, `test_fill_conformal_vtk_corner` | B: conformal VTK media-type keys are not published as expected. |

### Intel, MPI=OFF, MTLN=OFF, HDF=ON, double=OFF

Job: [94087314705](https://github.com/OpenSEMBA/fdtd/actions/runs/31588351969/job/94087314705)
Result: 4 failed, 82 passed, 40 skipped, 1 xfailed.

The failing-test set and causes are identical to Intel, MPI=ON, MTLN=OFF, double=OFF:

- D: `test_conformal_impedance_cylinder_unshielded`.
- E: `test_current_generators_without_resistance`.
- B: `test_fill_conformal_vtk_sphere`, `test_fill_conformal_vtk_corner`.

## Prioritised Investigation Order

1. Restore the expected map media-type contract in `src_output/mapVTKOutput.F90`; it accounts for 46 direct failures across the six completed jobs.
2. Inspect movie geometry companion publication and ensure each geometry cell is emitted once; it accounts for 20 direct failures.
3. Repair output discovery or align solver-generated filenames with the wrapper contract; it blocks 31 wire/MTLN assertions and conceals any underlying solver-value regressions in those modes.
4. Decide whether the conformal-cylinder mesh should be reoriented or whether the new shared-edge validation is too strict.
5. Trace the current-source magnitude conversion responsible for the `0.01` result.
