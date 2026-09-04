from test.utils.utils import mtln_skip, no_mtln_skip


@no_mtln_skip
def test_wire_probe_publishes_MTLN_data_without_metadata(run_output_case):
    process, output_root = run_output_case(
        "wire",
        [{"name": "wire_probe", "type": "wire", "field": "current", "elementIds": [1]}],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    outputs = sorted(output_root.glob("common_geometry.fdtd_wire_probe*.dat"))
    assert len(outputs) == 1
    assert outputs[0].is_file()
    assert outputs[0].parent == output_root
    assert not [
        path for path in output_root.rglob("*.json")
        if path.name != "common_geometry.fdtd.json"
    ]
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()


@mtln_skip
def test_wire_probe_publishes_only_flat_dat_output_without_mtln(run_output_case):
    process, output_root = run_output_case(
        "wire",
        [{"name": "wire_probe", "type": "wire", "field": "current", "elementIds": [1]}],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    outputs = sorted(output_root.glob("common_geometry.fdtd_wire_probe*"))
    assert len(outputs) == 1
    assert outputs[0].is_file()
    assert outputs[0].suffix == ".dat"
    assert outputs[0].parent == output_root
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()
