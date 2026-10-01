import copy
from pathlib import Path

from test.utils.utils import *


def _get_solved_probe_folder(solver, probe_name, *, filename=None, contains=None) -> str:
    probe_files = solver.getSolvedProbeFolders(probe_name)
    if filename is not None:
        probe_files = [path for path in probe_files if Path(path).stem == Path(filename).stem]
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
@pytest.mark.codemodel
@pytest.mark.codemodel
def test_towelHanger_mpi(tmp_path):
    fn = CASES_FOLDER + "towelHanger/towelHanger_mpi.fdtd.json"
    setNgspice(tmp_path)
    print(SEMBA_EXE)
    for layers in [1, 2]:
        for direction_index, direction in enumerate(["x", "y", "z"]):
            solver = FDTD(
                input_filename=fn,
                path_to_exe=SEMBA_EXE,
                run_in_folder=tmp_path,
                flags=["-mpidir " + direction],
                mpi_command="mpirun -np " + str(layers),
            )
            for coordinate in solver["mesh"]["coordinates"]:
                position = coordinate["relativePosition"]
                coordinate["relativePosition"] = [
                    position[(axis - direction_index) % 3] for axis in range(3)
                ]

            element = solver["mesh"]["elements"][2]
            element["intervals"] = [
                [
                    [endpoint[(axis - direction_index) % 3] for axis in range(3)]
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


_POINT_PROBE_BASE_CASE = CASES_FOLDER + "output_e2e/common_geometry.fdtd.json"
_POINT_PROBE_CASE_NAME = "point_probe_interface.fdtd"
_POINT_PROBE_COMPONENTS = (
    ("electric", "x"),
    ("electric", "y"),
    ("electric", "z"),
    ("magnetic", "x"),
    ("magnetic", "y"),
    ("magnetic", "z"),
)


def _stage_point_probe_case(tmp_path: Path, planes) -> Path:
    """Copy the shared output case and add point probes at the given z planes."""
    case = copy.deepcopy(json.loads(Path(_POINT_PROBE_BASE_CASE).read_text()))
    for z in planes:
        coordinate_id = 100 + z
        element_id = 200 + z
        case["mesh"]["coordinates"].append(
            {"id": coordinate_id, "relativePosition": [5, 4, z]}
        )
        case["mesh"]["elements"].append(
            {"id": element_id, "type": "node", "coordinateIds": [coordinate_id]}
        )
        for field, direction in _POINT_PROBE_COMPONENTS:
            case["probes"].append(
                {
                    "name": f"p{z}{field[0]}{direction}",
                    "type": "point",
                    "field": field,
                    "elementIds": [element_id],
                    "directions": [direction],
                    "domain": {"type": "time"},
                }
            )

    case_dir = tmp_path / "case"
    case_dir.mkdir()
    shutil.copy(EXCITATIONS_FOLDER + "gauss.exc", case_dir)
    input_path = case_dir / (_POINT_PROBE_CASE_NAME + ".json")
    input_path.write_text(json.dumps(case, indent=2))
    return input_path


def _run_point_probe_case(input_path, run_dir, mpi_command=None, flags=None):
    run_dir.mkdir()
    solver = FDTD(
        input_path,
        path_to_exe=SEMBA_EXE,
        run_in_folder=run_dir,
        mpi_command=mpi_command,
        flags=flags or [],
    )
    solver.run()
    assert solver.hasFinishedSuccessfully()
    return solver


def _read_point_probe_outputs(run_dir: Path):
    outputs = {}
    for path in sorted(run_dir.glob(f"{_POINT_PROBE_CASE_NAME}_*_tm.dat")):
        outputs[path.name] = np.atleast_2d(np.loadtxt(path, skiprows=1))
    return outputs


def _assert_point_outputs_match_serial(serial, parallel):
    assert set(parallel) == set(serial), sorted(set(parallel) ^ set(serial))
    for name, expected in serial.items():
        got = parallel[name]
        assert got.shape == expected.shape, f"{name}: {got.shape} != {expected.shape}"
        times = got[:, 0]
        assert np.all(np.diff(times) > 0), f"{name}: duplicated or unordered samples"
        assert np.allclose(got[:, 1], expected[:, 1], rtol=1e-6, atol=1e-12), name


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.probes
@pytest.mark.parametrize(
    "ranks, flags, planes",
    [
        # -force pins the cut at z=5, so planes 4-6 lie in the Alloc overlap.
        (2, ["-force", "5"], [2, 4, 5, 6]),
        # A second cut position exercises interface ownership at z=3.
        (2, ["-force", "3"], [2, 3, 4]),
    ],
)
def test_point_probe_at_mpi_interface_is_written_once(tmp_path, ranks, flags, planes):
    input_path = _stage_point_probe_case(tmp_path, planes)

    serial_dir = tmp_path / "serial"
    _run_point_probe_case(input_path, serial_dir)
    serial = _read_point_probe_outputs(serial_dir)
    assert serial, "serial run produced no point probe outputs"

    parallel_dir = tmp_path / f"mpi_{ranks}"
    _run_point_probe_case(
        input_path,
        parallel_dir,
        mpi_command=f"mpirun -np {ranks}",
        flags=flags,
    )
    parallel = _read_point_probe_outputs(parallel_dir)

    _assert_point_outputs_match_serial(serial, parallel)
