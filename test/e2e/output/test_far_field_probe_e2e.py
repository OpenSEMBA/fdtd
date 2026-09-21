import pytest


def test_far_field_probe_publishes_only_flat_dat_output(run_output_case):
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
                    "frequencySpacing": "logarithmic",
                },
            }
        ],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    outputs = sorted(output_root.glob("common_geometry.fdtd_far_field_probe*"))
    assert len(outputs) == 1
    assert outputs[0].is_file()
    assert outputs[0].suffix == ".dat"
    assert outputs[0].parent == output_root
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()

    records = []
    for line in outputs[0].read_text(encoding="utf-8").splitlines():
        try:
            records.append([float(value) for value in line.split()])
        except ValueError:
            continue

    assert len(records) == 30
    assert records[0][0] == pytest.approx(1e6)
    assert records[-1][0] == pytest.approx(1e9)
