def test_geometry_map_publishes_payloads_without_metadata(run_output_case):
    process, output_root = run_output_case(
        "geometry-map",
        [],
        additional_arguments="-mapvtk",
    )

    assert process.returncode == 0, process.stdout + process.stderr
    map_folder = output_root / "common_geometry.fdtd__MAP"
    assert map_folder.is_dir()
    assert not [
        path
        for path in output_root.glob("common_geometry.fdtd__MAP_*")
        if path.is_dir()
    ]
    geometry_paths = list(map_folder.glob("*.vtu"))
    assert len(geometry_paths) == 1
    geometry_path = geometry_paths[0]
    assert not geometry_path.with_suffix(".txt").exists()
    assert not list(output_root.rglob("*.pvtu"))
    assert not geometry_path.with_suffix(".json").exists()
    assert not geometry_path.with_suffix(".h5").exists()
    assert not geometry_path.with_suffix(".xdmf").exists()
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()
