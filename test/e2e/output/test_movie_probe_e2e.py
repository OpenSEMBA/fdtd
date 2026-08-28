import json


def test_movie_probe_publishes_complete_descriptor(run_output_case):
    process, output_root = run_output_case(
        "movie",
        [
            {
                "name": "movie_probe",
                "type": "movie",
                "field": "electric",
                "component": "x",
                "elementIds": [10],
                "domain": {"type": "time", "samplingPeriod": 1e-11},
            }
        ],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    descriptors = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in output_root.rglob("*.json")
        if path.name != "common_geometry.fdtd.json"
    ]
    assert any("lifecycle" in value for value in descriptors)


def test_mixed_scalar_and_movie_manifest_contains_only_coordinated_probe(run_output_case):
    process, output_root = run_output_case(
        "mixed",
        [
            {
                "name": "point_probe",
                "type": "point",
                "field": "electric",
                "elementIds": [1],
                "directions": ["x"],
                "domain": {"type": "time"},
            },
            {
                "name": "movie_probe",
                "type": "movie",
                "field": "electric",
                "component": "x",
                "elementIds": [10],
                "domain": {"type": "time", "samplingPeriod": 1e-11},
            },
        ],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    assert (output_root / "common_geometry.fdtd_point_probe_Ex_5_4_4_tm.dat").is_file()
    manifest = json.loads(
        (output_root / "common_geometry.fdtd_output_manifest.json").read_text(encoding="utf-8")
    )
    assert len(manifest["probes"]) == 1
    assert "movie_probe" in manifest["probes"][0]["probe_id"]
