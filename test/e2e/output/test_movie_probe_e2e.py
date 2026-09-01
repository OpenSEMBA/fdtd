from conftest import assert_static_point_attributes, read_xdmf_point_data_names


def test_movie_probe_publishes_payloads_without_metadata(run_output_case):
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
    output_directory = next(output_root.glob("common_geometry.fdtd_movie_probe_*"))
    assert list(output_directory.glob("*.bin"))
    assert list(output_directory.glob("*.xdmf"))
    assert list(output_directory.glob("*.h5"))
    assert not list(output_directory.glob("*.json"))
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()
    movie_xdmf = next(
        path
        for path in output_directory.glob("*.xdmf")
        if not path.stem.endswith("_geometry")
    )
    assert_static_point_attributes(movie_xdmf, ("tagnumber", "mediatype"))
    assert {"tagnumber", "mediatype"} <= read_xdmf_point_data_names(movie_xdmf)


def test_mixed_scalar_and_movie_publishes_no_metadata(run_output_case):
    process, output_root = run_output_case(
        "mixed",
        [
            {
                "name": "point_probe",
                "type": "point",
                "field": "electric",
                "elementIds": [1],
                "directions": ["x"],
                "domain": {"type": "time"},
            },
            {
                "name": "movie_probe",
                "type": "movie",
                "field": "electric",
                "component": "x",
                "elementIds": [10],
                "domain": {"type": "time", "samplingPeriod": 1e-11},
            },
        ],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    assert (output_root / "common_geometry.fdtd_point_probe_Ex_5_4_4_tm.dat").is_file()
    assert not list(output_root.rglob("common_geometry.fdtd_*probe*.json"))
    assert not (output_root / "common_geometry.fdtd_output_manifest.json").exists()


def test_vector_movie_probe_publishes_component_classification(run_output_case):
    process, output_root = run_output_case(
        "movie-vector",
        [
            {
                "name": "movie_probe",
                "type": "movie",
                "field": "electric",
                "component": "magnitude",
                "elementIds": [10],
                "domain": {"type": "time", "samplingPeriod": 1e-11},
            }
        ],
    )

    assert process.returncode == 0, process.stdout + process.stderr
    output_directory = next(output_root.glob("common_geometry.fdtd_movie_probe_*"))
    movie_xdmf = next(
        path
        for path in output_directory.glob("*.xdmf")
        if not path.stem.endswith("_geometry")
    )
    assert_static_point_attributes(
        movie_xdmf,
        (
            "tagnumber_x",
            "tagnumber_y",
            "tagnumber_z",
            "mediatype_x",
            "mediatype_y",
            "mediatype_z",
        ),
    )
