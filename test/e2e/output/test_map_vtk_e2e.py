def test_geometry_map_publishes_payloads_without_metadata(run_output_case):
    process, output_root = run_output_case(
        "geometry-map",
        [],
        additional_arguments="-mapvtk",
    )

    assert process.returncode == 0, process.stdout + process.stderr
    geometry_paths = list(output_root.rglob("*.vtu"))
    assert len(geometry_paths) == 1
    geometry_path = geometry_paths[0]
    assert geometry_path.with_suffix(".txt").is_file()
    assert not geometry_path.with_suffix(".json").exists()
    assert not geometry_path.with_suffix(".h5").exists()
    assert not geometry_path.with_suffix(".xdmf").exists()
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()
