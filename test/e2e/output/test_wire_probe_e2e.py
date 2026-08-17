import json


def test_wire_probe_publishes_complete_descriptor(run_output_case):
    process, output_root = run_output_case(
        "wire",
        [{"name": "wire_probe", "type": "wire", "field": "current", "elementIds": [1]}],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    descriptors = [
        (path, json.loads(path.read_text(encoding="utf-8")))
        for path in output_root.rglob("*.json")
        if path.name != "common_geometry.fdtd.json"
    ]
    probe_descriptors = [(path, value) for path, value in descriptors if "lifecycle" in value]
    assert probe_descriptors
    descriptor_path, descriptor = probe_descriptors[0]
    assert descriptor["lifecycle"]["state"] == "complete"
    assert descriptor["quantity"] == "Wx"
    artifact = descriptor["artifacts"][0]
    assert artifact["kind"] == "text"
    assert (descriptor_path.parent / artifact["relative_path"]).is_file()
