import json
from pathlib import Path


def test_geometry_map_publishes_a_complete_descriptor(run_output_case):
    process, output_root = run_output_case(
        "geometry-map",
        [],
        additional_arguments="-mapvtk",
    )

    assert process.returncode == 0, process.stdout + process.stderr
    descriptors = [
        (path, json.loads(path.read_text(encoding="utf-8")))
        for path in output_root.rglob("*.json")
        if path.name != "common_geometry.fdtd.json"
    ]
    descriptor_path, descriptor = next(
        (path, value) for path, value in descriptors if "lifecycle" in value
    )
    assert sum("lifecycle" in value for _, value in descriptors) == 1
    assert descriptor["lifecycle"]["state"] == "complete"
    assert descriptor["quantity"] == "geometry"
    assert descriptor["domain"] == "undefined"
    assert descriptor["ownership"] == {
        "participant_ranks": [0],
        "scalar_writer_rank": 0,
    }
    assert {artifact["kind"] for artifact in descriptor["artifacts"]} == {
        "geometry",
        "text",
    }
    for artifact in descriptor["artifacts"]:
        assert (descriptor_path.parent / artifact["relative_path"]).is_file()

    geometry_artifact = next(
        artifact for artifact in descriptor["artifacts"] if artifact["kind"] == "geometry"
    )
    geometry_path = descriptor_path.parent / geometry_artifact["relative_path"]
    assert not geometry_path.with_suffix(".h5").exists()
    assert not geometry_path.with_suffix(".xdmf").exists()

    manifest = next(value for _, value in descriptors if "probes" in value)
    manifest_probe = next(
        probe for probe in manifest["probes"] if probe["probe_id"] == descriptor["probe_id"]
    )
    expected_artifacts = {descriptor_path.name} | {
        artifact["relative_path"] for artifact in descriptor["artifacts"]
    }
    assert {
        Path(artifact["path"]).name for artifact in manifest_probe["artifacts"]
    } == expected_artifacts
