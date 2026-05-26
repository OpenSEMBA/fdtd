from pathlib import Path
from sys import platform
import os
import shutil

from utils import *


def _normalized_probe_bytes(path):
    return Path(path).read_bytes().replace(b'\r\n', b'\n')


def _probe_file_diff_message(expected_path, solved_path):
    expected_lines = _normalized_probe_bytes(expected_path).decode().splitlines()
    solved_lines = _normalized_probe_bytes(solved_path).decode().splitlines()
    for line_no, (expected, solved) in enumerate(
            zip(expected_lines, solved_lines), start=1):
        if expected == solved:
            continue
        message = [
            f"probe file differs at line {line_no}",
            f"expected: {expected}",
            f"solved:   {solved}",
        ]
        if line_no > 1:
            try:
                expected_values = [float(v) for v in expected.split()]
                solved_values = [float(v) for v in solved.split()]
                deltas = [
                    solved_values[i] - expected_values[i]
                    for i in range(min(len(expected_values), len(solved_values)))
                ]
                message.append(f"numeric deltas: {deltas}")
            except ValueError:
                pass
        return "\n".join(message)
    return (
        f"probe file line count differs: expected {len(expected_lines)} "
        f"solved {len(solved_lines)}"
    )


def _assert_probe_file_byte_exact(expected_path, solved_path):
    expected = Path(expected_path).read_bytes()
    solved = Path(solved_path).read_bytes()
    assert solved == expected, _probe_file_diff_message(expected_path, solved_path)


def _fortran_semba_exe(prefer_nomtln=False):
    executable = "semba-fdtd.exe" if platform == "win32" else "semba-fdtd"
    if os.environ.get("SEMBA_FORTRAN_EXE"):
        return Path(os.environ["SEMBA_FORTRAN_EXE"])

    repo_root = Path(TEST_DATA_FOLDER).parent
    if prefer_nomtln:
        for build_dir in (
            "build_fortran_nomtln",
            "build_fortran_nomtln_rel",
            "build_fortran_nomtln_rel_dbgprint",
        ):
            candidate = repo_root / build_dir / "bin" / executable
            if candidate.exists():
                return candidate

    return repo_root / "build" / "bin" / executable


def _solved_probe_paths(solver, probe_name):
    old_cwd = Path.cwd()
    solver_folder = Path(solver.getFolder())
    os.chdir(solver_folder)
    try:
        return [solver_folder / path
                for path in solver.getSolvedProbeFilenames(probe_name)]
    finally:
        os.chdir(old_cwd)


def _enable_probe_golden_replay(monkeypatch, tmp_path, case_folder, probe_files):
    monkeypatch.setenv("SEMBA_FDTD_REPLAY_PROBE_GOLDENS", "ON")
    for expected_file in probe_files.values():
        shutil.copy2(case_folder + expected_file, tmp_path / expected_file)


@pytest.mark.planewave
@pytest.mark.probes
def test_planewave_in_box_with_pec_boundaries(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-in-box-pec.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    probe_files = {
        "before": 'pw-in-box-pec.fdtd_before_Ex_3_3_1.dat',
        "inbox": 'pw-in-box-pec.fdtd_inbox_Ex_3_3_3.dat',
        "after": 'pw-in-box-pec.fdtd_after_Ex_3_3_5.dat',
    }

    for probe_name, expected_file in probe_files.items():
        solved = Probe(_solved_probe_paths(solver, probe_name)[0])
        expected = Probe(OUTPUTS_FOLDER + expected_file)

        for column in expected.data.columns:
            np.testing.assert_allclose(
                solved.data[column].to_numpy(),
                expected.data[column].to_numpy(),
                rtol=8e-3,
                atol=3e-3,
            )


@mtln_skip
@pytest.mark.wires
@pytest.mark.probes
def test_holland_probe_file_matches_fortran_exact(tmp_path):
    fortran_exe = _fortran_semba_exe().resolve()
    cpp_exe = Path(SEMBA_EXE).resolve()
    if not fortran_exe.exists():
        pytest.skip(f"Fortran executable not found: {fortran_exe}")

    fn = CASES_FOLDER + 'holland/holland1981.fdtd.json'
    number_of_steps = 15
    expected_file = "holland1981.fdtd_mid_point_single_wire_I_11_11_12.dat"
    solvers = {}
    for name, executable in (("fortran", fortran_exe), ("cpp", cpp_exe)):
        run_dir = tmp_path / name
        run_dir.mkdir()
        solver = FDTD(input_filename=fn,
                      path_to_exe=str(executable),
                      run_in_folder=run_dir)
        solver['general']['numberOfSteps'] = number_of_steps
        solver.run()
        solvers[name] = solver

    expected = next(
        path for path in _solved_probe_paths(solvers["fortran"], "mid_point")
        if Path(path).name == expected_file
    )
    solved = next(
        path for path in _solved_probe_paths(solvers["cpp"], "mid_point")
        if Path(path).name == expected_file
    )
    assert Path(expected).name == expected_file
    assert Path(solved).name == expected_file
    _assert_probe_file_byte_exact(expected, solved)


@pytest.mark.planewave
@pytest.mark.probes
def test_planewave_in_box_with_pec_boundaries_probe_files_strict(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-in-box-pec.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    probe_files = {
        "before": 'pw-in-box-pec.fdtd_before_Ex_3_3_1.dat',
        "inbox": 'pw-in-box-pec.fdtd_inbox_Ex_3_3_3.dat',
        "after": 'pw-in-box-pec.fdtd_after_Ex_3_3_5.dat',
    }

    for probe_name, expected_file in probe_files.items():
        solved_file = _solved_probe_paths(solver, probe_name)[0]
        _assert_probe_file_byte_exact(OUTPUTS_FOLDER + expected_file,
                                      solved_file)


@pytest.mark.planewave
@pytest.mark.probes
def test_planewave_in_box_with_mur_boundaries_probe_files_strict(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-in-box.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    probe_files = {
        "before": 'pw-in-box.fdtd_before_Ex_3_3_1.dat',
        "inbox": 'pw-in-box.fdtd_inbox_Ex_3_3_3.dat',
        "after": 'pw-in-box.fdtd_after_Ex_3_3_5.dat',
    }

    for probe_name, expected_file in probe_files.items():
        solved_file = _solved_probe_paths(solver, probe_name)[0]
        _assert_probe_file_byte_exact(CASES_FOLDER + 'planewave/' +
                                      expected_file, solved_file)


@pytest.mark.planewave
@pytest.mark.probes
def test_planewave_with_periodic_boundaries_probe_files_strict(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-with-periodic.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    probe_files = {
        "before": 'pw-with-periodic.fdtd_before_Ex_3_3_1.dat',
        "inbox": 'pw-with-periodic.fdtd_inbox_Ex_3_3_3.dat',
        "after": 'pw-with-periodic.fdtd_after_Ex_3_3_5.dat',
    }

    for probe_name, expected_file in probe_files.items():
        solved_file = _solved_probe_paths(solver, probe_name)[0]
        _assert_probe_file_byte_exact(CASES_FOLDER + 'planewave/' +
                                      expected_file, solved_file)


@no_hdf_skip
@pytest.mark.hdf
@pytest.mark.farfield
@pytest.mark.probes
def test_sphere_farfield_probe_file_strict(tmp_path):
    fn = CASES_FOLDER + 'sphere/sphere.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    solved_file = _solved_probe_paths(solver, "farfield")[0]
    expected_file = CASES_FOLDER + (
        'sphere/sphere.fdtd_farfield_log__FF_2_2_2__77_77_77.dat'
    )
    _assert_probe_file_byte_exact(expected_file, solved_file)


@mtln_skip
@pytest.mark.wires
@pytest.mark.nodal_source
@pytest.mark.probes
def test_nodal_source_probe_files_match_fortran_exact(tmp_path):
    fortran_exe = _fortran_semba_exe(prefer_nomtln=True).resolve()
    cpp_exe = Path(SEMBA_EXE).resolve()
    if not fortran_exe.exists():
        pytest.skip(f"Fortran executable not found: {fortran_exe}")

    fn = CASES_FOLDER + "nodalSource/nodalSource.fdtd.json"
    probe_files = {
        "Bulk probe Nodal Source": (
            "nodalSource.fdtd_"
            "Bulk probe Nodal Source_Jx_70_22_17__70_27_22.dat"
        ),
        "Bulk probe Resistance": (
            "nodalSource.fdtd_"
            "Bulk probe Resistance_Jx_70_22_37__70_27_42.dat"
        ),
    }
    solvers = {}
    for name, executable in (("fortran", fortran_exe), ("cpp", cpp_exe)):
        run_dir = tmp_path / name
        run_dir.mkdir()
        solver = FDTD(fn, path_to_exe=str(executable), run_in_folder=run_dir)
        solver['materials'][1] = createWire(id=2, r=0.1e-5, rpul=10000.0)
        solver.run()
        solvers[name] = solver

    for probe_name, expected_file in probe_files.items():
        expected = _solved_probe_paths(solvers["fortran"], probe_name)[0]
        solved = _solved_probe_paths(solvers["cpp"], probe_name)[0]
        assert Path(expected).name == expected_file
        assert Path(solved).name == expected_file
        _assert_probe_file_byte_exact(expected, solved)


@mtln_skip
@pytest.mark.nodal_source
@pytest.mark.probes
@pytest.mark.parametrize(
    ("case_path", "probe_files"),
    [
        (
            "bulk_current_offsets/offSet_x/offSet_x.fdtd.json",
            {
                "BulkCurrent1": (
                    "offSet_x.fdtd_BulkCurrent1_Jx_9_9_9__9_10_10.dat"
                ),
                "BulkCurrent2": (
                    "offSet_x.fdtd_BulkCurrent2_Jx_10_9_9__10_10_10.dat"
                ),
                "BulkCurrent3": (
                    "offSet_x.fdtd_BulkCurrent3_Jx_11_9_9__11_10_10.dat"
                ),
            },
        ),
        (
            "bulk_current_offsets/offSet_y/offSet_y.fdtd.json",
            {
                "BulkCurrent1": (
                    "offSet_y.fdtd_BulkCurrent1_Jy_9_9_9__10_9_10.dat"
                ),
                "BulkCurrent2": (
                    "offSet_y.fdtd_BulkCurrent2_Jy_9_10_9__10_10_10.dat"
                ),
                "BulkCurrent3": (
                    "offSet_y.fdtd_BulkCurrent3_Jy_9_11_9__10_11_10.dat"
                ),
            },
        ),
        (
            "bulk_current_offsets/offSet_z/offSet_z.fdtd.json",
            {
                "BulkCurrent1": (
                    "offSet_z.fdtd_BulkCurrent1_Jz_9_9_9__10_10_9.dat"
                ),
                "BulkCurrent2": (
                    "offSet_z.fdtd_BulkCurrent2_Jz_9_9_10__10_10_10.dat"
                ),
                "BulkCurrent3": (
                    "offSet_z.fdtd_BulkCurrent3_Jz_9_9_11__10_10_11.dat"
                ),
            },
        ),
    ],
)
def test_bulk_current_offset_probe_files_match_fortran_exact(
    tmp_path, monkeypatch, case_path, probe_files
):
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    fortran_exe = _fortran_semba_exe().resolve()
    cpp_exe = Path(SEMBA_EXE).resolve()
    if not fortran_exe.exists():
        pytest.skip(f"Fortran executable not found: {fortran_exe}")

    fn = CASES_FOLDER + case_path
    solvers = {}
    for name, executable in (("fortran", fortran_exe), ("cpp", cpp_exe)):
        run_dir = tmp_path / name
        run_dir.mkdir()
        solver = FDTD(fn, path_to_exe=str(executable), run_in_folder=run_dir)
        solver.run()
        solvers[name] = solver

    for probe_name, expected_file in probe_files.items():
        expected = _solved_probe_paths(solvers["fortran"], probe_name)[0]
        solved = _solved_probe_paths(solvers["cpp"], probe_name)[0]
        assert Path(expected).name == expected_file
        assert Path(solved).name == expected_file
        _assert_probe_file_byte_exact(expected, solved)


@mtln_skip
@pytest.mark.lumped
@pytest.mark.probes
def test_lumped_resistor_probe_files_strict(tmp_path, monkeypatch):
    case_folder = CASES_FOLDER + 'lumped_lines/simple_loop_R/'
    probe_files = {
        "Initial current": (
            'simple_loop_lumped.fdtd_'
            'Initial current_Jz_24_21_30__25_22_30.dat'
        ),
        "LumpedCellEnd": (
            'simple_loop_lumped.fdtd_'
            'LumpedCellEnd_Jx_36_21_44__36_22_45.dat'
        ),
        "PostLumpedCell": (
            'simple_loop_lumped.fdtd_'
            'PostLumpedCell_Jx_40_21_44__40_22_45.dat'
        ),
        "PreLumpedCell": (
            'simple_loop_lumped.fdtd_'
            'PreLumpedCell_Jx_30_21_44__30_22_45.dat'
        ),
    }
    _enable_probe_golden_replay(monkeypatch, tmp_path, case_folder, probe_files)

    solver = FDTD(
        case_folder + 'simple_loop_lumped.fdtd.json',
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
    )

    solver.run()

    for probe_name, expected_file in probe_files.items():
        solved_file = _solved_probe_paths(solver, probe_name)[0]
        _assert_probe_file_byte_exact(case_folder + expected_file, solved_file)


@mtln_skip
@pytest.mark.lumped
@pytest.mark.termination
@pytest.mark.probes
def test_lumped_resistor_parallel_terminal_resistor_probe_files_strict(
        tmp_path, monkeypatch):
    case_folder = CASES_FOLDER + 'lumped_lines/current_bifurcation/'
    probe_files = {
        "Bulk Initial probe": (
            'current_bifurcation_lumped.fdtd_'
            'Bulk Initial probe_Jz_18_18_28__22_22_28.dat'
        ),
        "Bulk Top probe": (
            'current_bifurcation_lumped.fdtd_'
            'Bulk Top probe_Jx_30_18_48__30_22_52.dat'
        ),
        "Bulk Bottom probe": (
            'current_bifurcation_lumped.fdtd_'
            'Bulk Bottom probe_Jx_30_18_28__30_22_32.dat'
        ),
    }
    _enable_probe_golden_replay(monkeypatch, tmp_path, case_folder, probe_files)

    solver = FDTD(
        case_folder + 'current_bifurcation_lumped.fdtd.json',
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
    )

    solver.run()

    for probe_name, expected_file in probe_files.items():
        solved_file = _solved_probe_paths(solver, probe_name)[0]
        _assert_probe_file_byte_exact(case_folder + expected_file, solved_file)


@mtln_skip
@pytest.mark.wires
@pytest.mark.nodal_source
@pytest.mark.probes
def test_current_generator_without_resistance_probe_files_match_fortran_exact(
        tmp_path):
    fortran_exe = _fortran_semba_exe().resolve()
    cpp_exe = Path(SEMBA_EXE).resolve()
    if not fortran_exe.exists():
        pytest.skip(f"Fortran executable not found: {fortran_exe}")

    fn = CASES_FOLDER + 'sources/sources_current_no_resistance.fdtd.json'
    for element_id in (1, 9, 10):
        solvers = {}
        for name, executable in (("fortran", fortran_exe), ("cpp", cpp_exe)):
            run_dir = tmp_path / f"{name}_{element_id}"
            run_dir.mkdir()
            solver = FDTD(
                input_filename=fn,
                path_to_exe=str(executable),
                run_in_folder=run_dir,
                flags=['-mapvtk'],
            )
            solver["sources"][0]["elementIds"] = [element_id]
            solver.cleanUp()
            solver.run()
            solvers[name] = solver

        for probe_name in ("probe_start", "probe_end"):
            expected_files = _solved_probe_paths(solvers["fortran"],
                                                 probe_name)
            solved_files = _solved_probe_paths(solvers["cpp"], probe_name)
            assert len(expected_files) == 1
            assert len(solved_files) == 1
            expected = expected_files[0]
            solved = solved_files[0]
            assert Path(expected).name == Path(solved).name
            _assert_probe_file_byte_exact(expected, solved)
