def test_line_probe_publishes_only_flat_dat_output(run_output_case):
    process, output_root = run_output_case(
        "line",
        [
            {
                "name": "line_probe",
                "type": "line",
                "elementIds": [2],
                "field": "electric",
                "domain": {"type": "time"},
            }
        ],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    outputs = sorted(output_root.glob("common_geometry.fdtd_line_probe*"))
    assert [path.name for path in outputs] == ["common_geometry.fdtd_line_probe_LI_tm.dat"]
    assert outputs[0].is_file()
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()
