import re

import h5py
import numpy as np

from conftest import assert_static_point_attributes, read_xdmf_point_data_names


def _movie_xdmf(output_root, case_root, probe_name):
    directory = next(output_root.glob(f"{case_root}_{probe_name}_*"))
    return next(
        path for path in directory.glob("*.xdmf") if not path.stem.endswith("_geometry")
    )


def _read_movie_attribute(xdmf_path, name):
    """Read one movie attribute by name from its XDMF/HDF5 payload."""
    contents = xdmf_path.read_text()
    match = re.search(
        rf'<Attribute Name="{name}".*?Format="HDF">\s*([^<\s]+)\s*</DataItem>',
        contents,
        re.DOTALL,
    )
    assert match is not None, f"Missing attribute {name} in {xdmf_path}"
    hdf_path = xdmf_path.with_suffix(".h5")
    resource = match.group(1)
    assert resource.startswith(hdf_path.name + ":"), resource
    dataset = resource.split(":", maxsplit=1)[1]
    with h5py.File(hdf_path, "r") as hdf_file:
        return hdf_file[dataset][()]


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


def test_movie_current_density_magnitude_reports_only_surface_components(run_output_case):
    """Magnitude movies must not leak non-surface current components.

    On a z-normal PEC surface only the in-plane (x, y) edges are part of the
    sheet. The z component is evaluated on an edge that is not a surface edge,
    so it must stay zero and match the per-component probe.
    """
    sampling = {"type": "time", "samplingPeriod": 1e-10}
    probes = [
        {
            "name": f"j_{component}",
            "type": "movie",
            "field": "currentDensity",
            "component": component,
            "elementIds": [12],
            "domain": sampling,
        }
        for component in ("magnitude", "x", "y", "z")
    ]

    process, output_root = run_output_case("pec_surface", probes)

    assert process.returncode == 0, process.stdout + process.stderr

    def attribute(probe_name, attribute_name):
        xdmf_path = _movie_xdmf(output_root, "pec_surface.fdtd", probe_name)
        return _read_movie_attribute(xdmf_path, attribute_name)

    magnitude_x = attribute("j_magnitude", "CurrenDensityX")
    magnitude_y = attribute("j_magnitude", "CurrenDensityY")
    magnitude_z = attribute("j_magnitude", "CurrenDensityZ")

    np.testing.assert_array_equal(magnitude_x, attribute("j_x", "CurrenDensityX"))
    np.testing.assert_array_equal(magnitude_y, attribute("j_y", "CurrenDensityY"))
    np.testing.assert_array_equal(attribute("j_z", "CurrenDensityZ"), 0.0)
    np.testing.assert_array_equal(magnitude_z, 0.0)

    # The magnitude movie must carry actual surface current, and the
    # classification must confirm that only the in-plane edges belong to the sheet.
    assert np.any(magnitude_x != 0.0)
    assert np.any(np.isclose(attribute("j_magnitude", "mediatype_x"), 0.5))
    assert not np.any(np.isclose(attribute("j_magnitude", "mediatype_z"), 0.5))
