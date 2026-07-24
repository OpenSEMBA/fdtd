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
