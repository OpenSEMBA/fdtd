import json


def test_geometry_map_publishes_a_complete_descriptor(run_output_case):
    process, output_root = run_output_case(
        "geometry-map",
        [],
        additional_arguments="-mapvtk",
    )

    assert process.returncode == 0, process.stdout + process.stderr
    descriptors = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in output_root.rglob("*.json")
        if path.name != "common_geometry.fdtd.json"
    ]
    assert descriptors
