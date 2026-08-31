def test_frequency_slice_probe_publishes_payloads_without_metadata(run_output_case):
    process, output_root = run_output_case(
        "frequency-slice",
        [
            {
                "name": "frequency_slice_probe",
                "type": "movie",
                "field": "electric",
                "component": "x",
                "elementIds": [10],
                "domain": {
                    "type": "frequency",
                    "initialFrequency": 1e9,
                    "finalFrequency": 1e9,
                    "numberOfFrequencies": 1,
                },
            }
        ],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    output_directory = next(output_root.glob("common_geometry.fdtd_frequency_slice_probe_*"))
    assert list(output_directory.glob("*.bin"))
    assert list(output_directory.glob("*.xdmf"))
    assert list(output_directory.glob("*.h5"))
    assert not list(output_directory.glob("*.json"))
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()
