from pathlib import Path
import json
import os
import shutil
import subprocess

import pytest

from utils import CASES_FOLDER, FDTD, SEMBA_EXE


pytestmark = [
    pytest.mark.mpi,
    pytest.mark.skipif(
        os.getenv("SEMBA_FDTD_ENABLE_MPI") != "ON",
        reason="requires a C++ build configured with SEMBA_FDTD_ENABLE_MPI=ON",
    ),
]


def _mpi_command(ranks: int) -> str:
    template = os.getenv("SEMBA_MPI_COMMAND", "mpirun -np {ranks}")
    return template.format(ranks=ranks)


def _run_case(input_file: str, run_dir: Path, path_to_exe: str, flags=None, mpi_command=None):
    run_dir.mkdir(parents=True, exist_ok=True)
    cwd = Path.cwd()
    try:
        solver = FDTD(
            input_filename=input_file,
            path_to_exe=path_to_exe,
            flags=flags or [],
            mpi_command=mpi_command,
            run_in_folder=run_dir,
        )
        solver.run()
        return solver.getCaseName()
    finally:
        os.chdir(cwd)


def _probe_files(run_dir: Path, case_name: str):
    return sorted(
        path
        for path in run_dir.glob(f"{case_name}_*.dat")
        if path.is_file() and "Energy" not in path.name
    )


def _run_raw_mpi_case(executable: str, run_dir: Path, case_file: str, ranks: int):
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    exe = str(Path(executable).resolve())
    command = [*_mpi_command(ranks).split(), exe, "-i", case_file, "-n", "2", "-mpidir", "z"]
    completed = subprocess.run(
        command,
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
    )
    return completed.returncode, completed.stdout + completed.stderr


def test_planewave_pec_z_mpi_probe_files_match_fortran_bytes(tmp_path):
    if shutil.which("mpirun") is None and "SEMBA_MPI_COMMAND" not in os.environ:
        pytest.skip("mpirun not available")
    fortran_exe = os.getenv("SEMBA_FORTRAN_EXE")
    if not fortran_exe:
        pytest.skip("SEMBA_FORTRAN_EXE is required for Fortran MPI byte parity")

    input_file = CASES_FOLDER + "planewave/pw-in-box-pec.fdtd.json"
    fortran_dir = tmp_path / "fortran_z"
    cpp_dir = tmp_path / "cpp_z"

    case_name = _run_case(
        input_file,
        fortran_dir,
        path_to_exe=fortran_exe,
        flags=["-mpidir", "z"],
        mpi_command=_mpi_command(2),
    )
    cpp_case_name = _run_case(
        input_file,
        cpp_dir,
        path_to_exe=SEMBA_EXE,
        flags=["-mpidir", "z"],
        mpi_command=_mpi_command(2),
    )
    assert cpp_case_name == case_name

    fortran_files = _probe_files(fortran_dir, case_name)
    cpp_files = _probe_files(cpp_dir, case_name)
    assert fortran_files
    assert [path.name for path in cpp_files] == [path.name for path in fortran_files]

    for fortran_path, cpp_path in zip(fortran_files, cpp_files):
        assert cpp_path.read_bytes() == fortran_path.read_bytes(), cpp_path.name


def test_towelhanger_mpi_pml_slice_failure_is_clean_and_matches_fortran_condition(tmp_path):
    if shutil.which("mpirun") is None and "SEMBA_MPI_COMMAND" not in os.environ:
        pytest.skip("mpirun not available")
    fortran_exe = os.getenv("SEMBA_FORTRAN_EXE")
    if not fortran_exe:
        pytest.skip("SEMBA_FORTRAN_EXE is required for Fortran MPI failure parity")

    case_dir = Path(CASES_FOLDER) / "towelHanger"
    base_case = case_dir / "towelHanger_mpi.fdtd.json"
    excitation = case_dir / "towelHanger.exc"
    base_data = json.loads(base_case.read_text())
    fail_data = json.loads(base_case.read_text())
    fail_data["boundary"]["all"]["layers"] = 30

    fortran_ok_dir = tmp_path / "fortran_pml_slice_ok"
    cpp_ok_dir = tmp_path / "cpp_pml_slice_ok"
    fortran_fail_dir = tmp_path / "fortran_pml_slice_fail"
    cpp_fail_dir = tmp_path / "cpp_pml_slice_fail"
    for path in (fortran_ok_dir, cpp_ok_dir, fortran_fail_dir, cpp_fail_dir):
        path.mkdir(parents=True, exist_ok=True)
        shutil.copy2(excitation, path / excitation.name)

    (fortran_ok_dir / "case.fdtd.json").write_text(json.dumps(base_data, indent=2) + "\n")
    (cpp_ok_dir / "case.fdtd.json").write_text(json.dumps(base_data, indent=2) + "\n")
    (fortran_fail_dir / "case.fdtd.json").write_text(json.dumps(fail_data, indent=2) + "\n")
    (cpp_fail_dir / "case.fdtd.json").write_text(json.dumps(fail_data, indent=2) + "\n")

    rc_fortran_ok, txt_fortran_ok = _run_raw_mpi_case(
        fortran_exe, fortran_ok_dir, "case.fdtd.json", ranks=3
    )
    rc_cpp_ok, txt_cpp_ok = _run_raw_mpi_case(
        SEMBA_EXE, cpp_ok_dir, "case.fdtd.json", ranks=3
    )
    assert rc_fortran_ok == 0, txt_fortran_ok
    assert rc_cpp_ok == 0, txt_cpp_ok

    rc_fortran_3, txt_fortran_3 = _run_raw_mpi_case(
        fortran_exe, fortran_fail_dir, "case.fdtd.json", ranks=3
    )
    rc_cpp_3, txt_cpp_3 = _run_raw_mpi_case(
        SEMBA_EXE, cpp_fail_dir, "case.fdtd.json", ranks=3
    )
    assert rc_fortran_3 != 0, txt_fortran_3
    assert rc_cpp_3 != 0, txt_cpp_3

    expected_fortran_message = (
        "Minimum slice sizes along MPI should be larger that PML number of layers"
    )
    assert expected_fortran_message in txt_fortran_3

    assert "Segmentation fault" not in txt_cpp_3
    assert "signal 11" not in txt_cpp_3.lower()
    txt_cpp_3_lc = txt_cpp_3.lower()
    assert "pml" in txt_cpp_3_lc
    assert ("slice" in txt_cpp_3_lc) or ("layer" in txt_cpp_3_lc)
