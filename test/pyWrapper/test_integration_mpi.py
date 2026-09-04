import json
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import pyvista as pv

from test.utils.utils import *


def movie_binary_records(movie_probe):
    binary_path = movie_probe.getBinFile()
    assert binary_path is not None
    records = np.fromfile(binary_path, dtype="<f8")
    assert records.size % 7 == 0
    return records.reshape(-1, 7)


def movie_binary_coordinates(movie_probe):
    return movie_binary_records(movie_probe)[:, 1:4].astype(int)


def assert_static_point_attributes(xdmf_path, attribute_names):
    root = ET.parse(xdmf_path).getroot()
    for name in attribute_names:
        attributes = root.findall(f'.//Attribute[@Name="{name}"]')
        assert attributes
        for attribute in attributes:
            data_item = attribute.find("DataItem")
            assert data_item is not None
            assert data_item.attrib.get("ItemType") is None
            assert data_item.attrib.get("Format") == "HDF"


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.wires
def test_airplane_case_with_mpi(tmp_path):
    fn = CASES_FOLDER + "airplane/airplane.fdtd.json"
    solver = FDTD(
        fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
        mpi_command="mpirun -np 2",
    )
    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)
    assert Path(vtkmapfile).suffix == ".pvtu"

    map_directories = sorted(
        path for path in tmp_path.glob("airplane.fdtd__MAP_*") if path.is_dir()
    )
    assert [path.name for path in map_directories] == [
        "airplane.fdtd__MAP_0_0_0__49_49_24",
        "airplane.fdtd__MAP_0_0_25__49_49_49",
    ]
    assert all((path / f"{path.name}.vtu").is_file() for path in map_directories)
    assert not list(tmp_path.rglob("*MAP*.txt"))

    root = ET.parse(vtkmapfile).getroot()
    piece_sources = [piece.attrib["Source"] for piece in root.findall(".//Piece")]
    assert piece_sources == [
        f"{path.name}/{path.name}.vtu" for path in map_directories
    ]
    assert all((tmp_path / source).is_file() for source in piece_sources)

    parallel_map = pv.read(vtkmapfile)
    pieces = [pv.read(tmp_path / source) for source in piece_sources]
    assert parallel_map.n_cells == sum(piece.n_cells for piece in pieces)
    assert set(parallel_map.cell_data) == {"tagnumber", "mediatype"}


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.wires
@pytest.mark.probes
@pytest.mark.movie
def test_airplane_movie_is_published_canonically_with_mpi(tmp_path):
    fn = CASES_FOLDER + "airplane/airplane.fdtd.json"
    solver = FDTD(
        fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        mpi_command="mpirun -np 2",
    )
    solver.run()

    movie_probes = [Probe(path) for path in solver.getSolvedProbeFolders("Movie")]
    assert len(movie_probes) == 1
    probe = movie_probes[0]
    assert probe.getXDMFFile() is not None
    assert probe.getH5File() is not None
    assert_static_point_attributes(
        probe.getXDMFFile(),
        (
            "tagnumber_x",
            "tagnumber_y",
            "tagnumber_z",
            "mediatype_x",
            "mediatype_y",
            "mediatype_z",
        ),
    )

    movie_coordinates = np.unique(movie_binary_coordinates(probe), axis=0)
    assert np.all(movie_coordinates >= probe.cell_init)
    assert np.all(movie_coordinates <= probe.cell_end)

    assert np.all(movie_coordinates[:, 0] >= 40)
    assert np.all(movie_coordinates[:, 0] <= 49)
    assert np.all(movie_coordinates[:, 1] == 26)
    assert np.all(movie_coordinates[:, 2] >= 10)
    assert np.all(movie_coordinates[:, 2] <= 36)
    assert len(np.unique(movie_coordinates, axis=0)) == len(movie_coordinates)
    assert np.array_equal(np.unique(movie_coordinates[:, 2]), np.arange(10, 37))


    serial_folder = tmp_path / "serial"
    serial_folder.mkdir()
    serial_solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=serial_folder)
    serial_solver.run()
    serial_probe = Probe(serial_solver.getSolvedProbeFolders("Movie")[0])
    np.testing.assert_allclose(movie_binary_records(probe), movie_binary_records(serial_probe))


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.probes
@pytest.mark.hdf
def test_frequency_slice_is_published_canonically_with_mpi(tmp_path):
    fn = CASES_FOLDER + "planewave/pw-in-box-with-frequency-slice.fdtd.json"
    solver = FDTD(
        fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        mpi_command="mpirun -np 2",
    )
    solver.run()

    probe_folders = solver.getSolvedProbeFolders("electric_field_frequency_slice")
    assert len(probe_folders) == 1
    probe = Probe(probe_folders[0])
    h5_path = probe.getH5File()
    xdmf_path = probe.getXDMFFile()
    assert xdmf_path is not None
    assert h5_path is not None
    geometry_xdmf_path = Path(xdmf_path).with_name(
        Path(xdmf_path).stem + "_geometry.xdmf"
    )
    geometry_h5_path = geometry_xdmf_path.with_suffix(".h5")
    assert geometry_xdmf_path.is_file()
    assert geometry_h5_path.is_file()
    assert not list(Path(probe.folder).glob("*.bin"))
    assert_static_point_attributes(
        xdmf_path,
        ("tagnumber", "mediatype"),
    )

    with h5py.File(h5_path, "r") as h5_file:
        assert h5_file["attributes/a0001/values"].shape == (4, 120)
        assert np.max(np.abs(h5_file["attributes/a0001/values"][()])) > 0.0

    assert not list(Path(probe.folder).glob("*.json"))
    assert not (tmp_path / f"{solver.getCaseName()}_output_manifest.json").exists()

    serial_folder = tmp_path / "serial"
    serial_folder.mkdir()
    serial_solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=serial_folder)
    serial_solver.run()
    serial_probe = Probe(
        serial_solver.getSolvedProbeFolders("electric_field_frequency_slice")[0]
    )
    serial_xdmf_path = serial_probe.getXDMFFile()
    assert serial_xdmf_path is not None
    serial_geometry_xdmf_path = Path(serial_xdmf_path).with_name(
        Path(serial_xdmf_path).stem + "_geometry.xdmf"
    )
    serial_geometry_h5_path = serial_geometry_xdmf_path.with_suffix(".h5")
    assert serial_geometry_xdmf_path.is_file()
    assert serial_geometry_h5_path.is_file()
    assert not list(Path(serial_probe.folder).glob("*.bin"))
    with h5py.File(h5_path, "r") as mpi_h5, h5py.File(
        serial_probe.getH5File(), "r"
    ) as serial_h5:
        np.testing.assert_allclose(
            mpi_h5["attributes/a0001/values"][()],
            serial_h5["attributes/a0001/values"][()],
            atol=1e-14
        )


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.wires
def test_simple_cabin_initialization_with_mpi(tmp_path):
    fn = CASES_FOLDER + "simple_cabin/simple_cabin.fdtd.json"
    solver = FDTD(
        fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
        mpi_command="mpirun -np 2",
    )
    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)


@no_mpi_skip
@no_mtln_skip
@pytest.mark.mpi
@pytest.mark.mtln
@pytest.mark.probes
def test_mtln_non_root_writer_publishes_data_without_metadata(tmp_path):
    input_data = json.loads(
        (Path(CASES_FOLDER) / "mpi" / "bundles_for_mpi.fdtd.json").read_text()
    )
    input_data["mesh"]["elements"].append(
        {"id": 2, "type": "node", "coordinateIds": [2]}
    )
    input_data["probes"] = [
        {
            "name": "upper",
            "type": "wire",
            "field": "current",
            "elementIds": [2],
            "domain": {"type": "time"},
        }
    ]
    input_path = tmp_path / "mtln_non_root.fdtd.json"
    input_path.write_text(json.dumps(input_data))
    solver = FDTD(
        input_path,
        path_to_exe=SEMBA_EXE,
        mpi_command="mpirun -np 2",
    )

    solver.run()

    probe_path = Path(solver.getSolvedProbeFolders("upper")[0])
    assert probe_path.is_file()
    assert probe_path.suffix == ".dat"
