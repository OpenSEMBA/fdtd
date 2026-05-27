from pathlib import Path
import os
import shutil

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
