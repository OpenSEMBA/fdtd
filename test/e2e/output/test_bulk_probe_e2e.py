import json


def test_bulk_probe_publishes_complete_descriptor(run_output_case):
    process, output_root = run_output_case(
        "bulk",
        [{"name": "bulk_probe", "type": "bulkCurrent", "elementIds": [10], "direction": "x"}],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    descriptors = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in output_root.rglob("*.json")
        if path.name != "common_geometry.fdtd.json"
    ]
    assert any("lifecycle" in value for value in descriptors)
