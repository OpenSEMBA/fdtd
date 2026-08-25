import json


def test_point_probe_publishes_complete_descriptor(run_output_case):
    process, output_root = run_output_case(
        "point",
        [
            {
                "name": "point_probe",
                "type": "point",
                "field": "electric",
                "elementIds": [1],
                "directions": ["x"],
                "domain": {"type": "time"},
            }
        ],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    descriptors = []
    for path in output_root.rglob("*.json"):
        if path.name == "common_geometry.fdtd.json":
            continue
        value = json.loads(path.read_text(encoding="utf-8"))
        if "lifecycle" in value:
            descriptors.append(value)
    assert descriptors
    assert len(descriptors) == 1
    descriptor = descriptors[0]
    assert descriptor["lifecycle"]["state"] == "complete"
    assert descriptor["artifacts"]
    assert descriptor["lower_bound"] == {"x": 5, "y": 4, "z": 4}
    assert descriptor["upper_bound"] == {"x": 5, "y": 4, "z": 4}
    assert descriptor["domain"] == "time"
    assert descriptor["ownership"] == {
        "participant_ranks": [0],
        "scalar_writer_rank": 0,
    }
