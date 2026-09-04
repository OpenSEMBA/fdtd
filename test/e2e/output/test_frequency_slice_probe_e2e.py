from conftest import assert_static_point_attributes, read_xdmf_point_data_names


def test_frequency_slice_probe_publishes_xdmf_and_geometry_without_binary(run_output_case):
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
    output_name = output_directory.name
    assert not list(output_directory.glob("*.bin"))
    assert (output_directory / f"{output_name}.xdmf").is_file()
    assert (output_directory / f"{output_name}.h5").is_file()
    assert (output_directory / f"{output_name}_geometry.xdmf").is_file()
    assert (output_directory / f"{output_name}_geometry.h5").is_file()
    assert not list(output_directory.glob("*.json"))
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()
    frequency_xdmf = output_directory / f"{output_name}.xdmf"
    assert_static_point_attributes(frequency_xdmf, ("tagnumber", "mediatype"))
    assert {"tagnumber", "mediatype"} <= read_xdmf_point_data_names(frequency_xdmf)


def test_vector_frequency_slice_publishes_component_classification(run_output_case):
    process, output_root = run_output_case(
        "frequency-slice-vector",
        [
            {
                "name": "frequency_slice_probe",
                "type": "movie",
                "field": "electric",
                "component": "magnitude",
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
    output_name = output_directory.name
    assert not list(output_directory.glob("*.bin"))
    frequency_xdmf = output_directory / f"{output_name}.xdmf"
    assert_static_point_attributes(
        frequency_xdmf,
        (
            "tagnumber_x",
            "tagnumber_y",
            "tagnumber_z",
            "mediatype_x",
            "mediatype_y",
            "mediatype_z",
        ),
    )
