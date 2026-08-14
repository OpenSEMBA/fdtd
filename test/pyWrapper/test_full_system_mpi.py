from pathlib import Path

from utils import *


def _get_solved_probe_folder(solver, probe_name, *, filename=None, contains=None) -> str:
    probe_files = solver.getSolvedProbeFolders(probe_name)
    if filename is not None:
        probe_files = [path for path in probe_files if Path(path).name == Path(filename).stem]
    if contains is not None:
        probe_files = [path for path in probe_files if contains in Path(path).name]
    assert len(probe_files) == 1, (
        f"Expected one artifact for probe {probe_name!r}, found {probe_files}"
    )
    return probe_files[0]


@no_mtln_skip
@no_mpi_skip
@pytest.mark.mtln
@pytest.mark.mpi
@pytest.mark.wires
@pytest.mark.multiwire
def test_bundles_mpi_n_ranks(tmp_path):
    fn = CASES_FOLDER + "mpi/bundles_for_mpi.fdtd.json"
    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        mpi_command="mpirun -np 2",
        run_in_folder=tmp_path,
    )
    solver.run()
    assert solver.hasFinishedSuccessfully()


@no_mtln_skip
@no_mpi_skip
@pytest.mark.mtln
@pytest.mark.mpi
@pytest.mark.wires
@pytest.mark.multiwire
def test_bundles_mpi_n_ranks_2(tmp_path):
    fn = CASES_FOLDER + "mpi/bundles_for_mpi_2.fdtd.json"
    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        mpi_command="mpirun -np 2",
        run_in_folder=tmp_path,
    )
    solver.run()
    assert solver.hasFinishedSuccessfully()


@no_mtln_skip
@no_mpi_skip
@pytest.mark.mtln
@pytest.mark.mpi
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.probes
def test_shieldedPair_mpi(tmp_path):
    fn = CASES_FOLDER + "shieldedPair/shieldedPair.fdtd.json"
    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        mpi_command="mpirun -np 2",
        run_in_folder=tmp_path,
    )
    solver.run()

    probe_files = [
        "shieldedPair.fdtd_wire_start_line_out_V_75_74_74.dat",
        "shieldedPair.fdtd_wire_start_line_out_I_75_74_74.dat",
        "shieldedPair.fdtd_wire_end_line_out_I_75_71_74.dat",
        "shieldedPair.fdtd_wire_end_line_out_V_75_71_74.dat",
    ]
    p_expected = [probe_from_fixture(tmp_path, filename) for filename in probe_files]
    p_solved = [
        Probe(
            _get_solved_probe_folder(
                solver,
                "wire_start" if "_wire_start_" in filename else "wire_end",
                filename=filename,
            )
        )
        for filename in probe_files
    ]

    for index in [0, 3]:
        for component in range(3):
            solved = np.interp(
                p_expected[index]["time"].to_numpy(),
                p_solved[index]["time"].to_numpy(),
                p_solved[index][f"voltage_{component}"].to_numpy(),
            )
            assert np.corrcoef(solved, p_expected[index][f"voltage_{component}"])[0, 1] > 0.999
    for index in [1, 2]:
        for component in range(3):
            solved = np.interp(
                p_expected[index]["time"].to_numpy(),
                p_solved[index]["time"].to_numpy(),
                p_solved[index][f"current_{component}"].to_numpy(),
            )
            assert np.corrcoef(solved, p_expected[index][f"current_{component}"])[0, 1] > 0.999


@no_mtln_skip
@no_mpi_skip
@pytest.mark.mtln
@pytest.mark.mpi
@pytest.mark.wires
@pytest.mark.probes
def test_holland_mtln_mpi(tmp_path):
    fn = CASES_FOLDER + "holland/holland1981_unshielded.fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    probe_mid_no_mpi = Probe(
        _get_solved_probe_folder(solver, "mid_point", contains="_I_")
    )

    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        mpi_command="mpirun -np 1",
        run_in_folder=tmp_path,
    )
    solver.cleanUp()
    solver.run()
    probe_mid_mpi_1 = Probe(
        _get_solved_probe_folder(solver, "mid_point", contains="_I_")
    )

    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        mpi_command="mpirun -np 2",
        run_in_folder=tmp_path,
    )
    solver.cleanUp()
    solver.run()
    probe_mid_mpi_2 = Probe(
        _get_solved_probe_folder(solver, "mid_point", contains="_I_")
    )

    expected_f = json.load(
        open(OUTPUTS_FOLDER + "holland1981_mid_point_expected_current.json")
    )
    expected_t, expected_i = np.array([]), np.array([])
    for data in expected_f["datasetColl"][0]["data"]:
        expected_t = np.append(expected_t, float(data["value"][0]))
        expected_i = np.append(expected_i, float(data["value"][1]))
    expected_i_interp = np.interp(
        probe_mid_no_mpi["time"] - 3.05 * 1e-9, expected_t, expected_i
    )

    for probe in [probe_mid_no_mpi, probe_mid_mpi_1, probe_mid_mpi_2]:
        assert np.allclose(
            expected_i_interp, probe["current_0"], rtol=1e-4, atol=5e-5
        )


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.wires
@pytest.mark.probes
def test_towelHanger_mpi(tmp_path):
    fn = CASES_FOLDER + "towelHanger/towelHanger_mpi.fdtd.json"
    for direction_index, direction in enumerate(["x", "y", "z"]):
        solver = FDTD(
            input_filename=fn,
            path_to_exe=SEMBA_EXE,
            run_in_folder=tmp_path,
            flags=["-mpidir " + direction],
            mpi_command="mpirun -np 1",
        )
        for coordinate in solver["mesh"]["coordinates"]:
            position = coordinate["relativePosition"]
            coordinate["relativePosition"] = [
                position[(axis + direction_index) % 3] for axis in range(3)
            ]

        element = solver["mesh"]["elements"][2]
        element["intervals"] = [
            [
                [endpoint[(axis + direction_index) % 3] for axis in range(3)]
                for endpoint in element["intervals"][0]
            ]
        ]
        solver.cleanUp()
        solver.run()

        p_solved = [
            Probe(_get_solved_probe_folder(solver, name))
            for name in ["wire_start", "wire_mid", "wire_end"]
        ]
        p_expected = [
            probe_from_fixture(tmp_path, filename)
            for filename in [
                "towelHanger.fdtd_wire_start_Wz_27_25_30_s1.dat",
                "towelHanger.fdtd_wire_mid_Wx_35_25_32_s5.dat",
                "towelHanger.fdtd_wire_end_Wz_43_25_30_s4.dat",
            ]
        ]
        for solved_probe, expected_probe in zip(p_solved, p_expected):
            solved = np.interp(
                expected_probe["time"].to_numpy(),
                solved_probe["time"].to_numpy(),
                solved_probe["current_0"].to_numpy(),
            )
            assert np.corrcoef(solved, expected_probe["current_0"])[0, 1] > 0.999
