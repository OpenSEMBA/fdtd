def test_bulk_probe_publishes_only_flat_dat_output(run_output_case):
    process, output_root = run_output_case(
        "bulk",
        [{"name": "bulk_probe", "type": "bulkCurrent", "elementIds": [10], "direction": "x"}],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    outputs = sorted(output_root.glob("common_geometry.fdtd_bulk_probe*"))
    assert outputs
    assert all(path.is_file() for path in outputs)
    assert all(path.suffix == ".dat" for path in outputs)
    assert all(path.parent == output_root for path in outputs)
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()
