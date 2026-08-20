from test.utils.utils import *


def movie_binary_coordinates(movie_probe):
    binary_path = movie_probe.getBinFile()
    assert binary_path is not None
    records = np.fromfile(binary_path, dtype="<f8")
    assert records.size % 7 == 0
    return records.reshape(-1, 7)[:, 1:4].astype(int)


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

    map_directories = sorted(tmp_path.glob("airplane.fdtd__MAP_*"))
    assert [path.name for path in map_directories] == [
        "airplane.fdtd__MAP_0_0_0__49_49_24",
        "airplane.fdtd__MAP_0_0_25__49_49_49",
    ]
    assert all((path / f"{path.name}.vtu").is_file() for path in map_directories)


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.wires
@pytest.mark.probes
@pytest.mark.movie
def test_airplane_movie_coordinates_are_partitioned_with_mpi(tmp_path):
    fn = CASES_FOLDER + "airplane/airplane.fdtd.json"
    solver = FDTD(
        fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        mpi_command="mpirun -np 2",
    )
    solver.run()

    movie_probes = [Probe(path) for path in solver.getSolvedProbeFolders("Movie")]
    assert len(movie_probes) == 2
    assert all(probe.getXDMFFile() is not None for probe in movie_probes)
    assert all(probe.getH5File() is not None for probe in movie_probes)

    fragment_coordinates = [
        np.unique(movie_binary_coordinates(probe), axis=0) for probe in movie_probes
    ]
    for probe, coordinates in zip(movie_probes, fragment_coordinates):
        assert np.all(coordinates >= probe.cell_init)
        assert np.all(coordinates <= probe.cell_end)

    movie_coordinates = np.concatenate(fragment_coordinates)

    assert np.all(movie_coordinates[:, 0] >= 40)
    assert np.all(movie_coordinates[:, 0] <= 49)
    assert np.all(movie_coordinates[:, 1] == 26)
    assert np.all(movie_coordinates[:, 2] >= 10)
    assert np.all(movie_coordinates[:, 2] <= 36)
    assert len(np.unique(movie_coordinates, axis=0)) == len(movie_coordinates)
    assert np.array_equal(np.unique(movie_coordinates[:, 2]), np.arange(10, 37))


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
