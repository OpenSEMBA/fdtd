import json


def test_far_field_probe_publishes_complete_descriptor(run_output_case):
    process, output_root = run_output_case(
        "far-field",
        [
            {
                "name": "far_field_probe",
                "type": "farField",
                "elementIds": [10],
                "theta": {"initial": 0, "final": 0, "step": 1},
                "phi": {"initial": 0, "final": 0, "step": 1},
                "domain": {
                    "type": "frequency",
                    "initialFrequency": 1e6,
                    "finalFrequency": 1e9,
                    "numberOfFrequencies": 30,
                    "frequencySpacing": "logarithmic"
                }
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
