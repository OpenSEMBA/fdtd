from pathlib import Path
import os
import shlex
import shutil
import subprocess

import pytest

from utils import CASES_FOLDER


pytestmark = [pytest.mark.cpp_migration]


MOVIE_CASE_DIR = Path(CASES_FOLDER) / "planewave"
MOVIE_CASE_FILE = "pw-in-box-with-movie.fdtd.json"
MOVIE_CASE_SUPPORT = ("gauss_1GHz.exc",)
PAUL_CASE_DIR = Path(CASES_FOLDER) / "paul"
PAUL_CASE_FILE = "paul_8_6_square.fdtd.json"


def _required_executable(env_name: str) -> Path:
    value = os.getenv(env_name, "").strip()
    assert value, (
        f"Missing {env_name}. Set it to the executable path before running this "
        "test module."
    )
    exe = Path(value).expanduser().resolve()
    assert exe.exists(), f"{env_name} does not exist: {exe}"
    return exe


@pytest.fixture(scope="session")
def cpp_mtln_exe() -> Path:
    return _required_executable("SEMBA_CPP_MTLN_EXE")


@pytest.fixture(scope="session")
def cpp_mpi_exe() -> Path:
    return _required_executable("SEMBA_CPP_MPI_EXE")


@pytest.fixture(scope="session")
def fortran_mtln_exe() -> Path:
    return _required_executable("SEMBA_FORTRAN_MTLN_EXE")


@pytest.fixture(scope="session")
def fortran_mpi_exe() -> Path:
    return _required_executable("SEMBA_FORTRAN_MPI_EXE")


def _mpi_launcher(ranks: int) -> list[str]:
    template = os.getenv("SEMBA_MPI_COMMAND", "mpirun -np {ranks}")
    launcher = template.format(ranks=ranks)
    argv = shlex.split(launcher)
    assert argv, "SEMBA_MPI_COMMAND resolved to an empty command"
    assert shutil.which(argv[0]) is not None, f"MPI launcher not found: {argv[0]}"
    return argv


def _stage_movie_case(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(MOVIE_CASE_DIR / MOVIE_CASE_FILE, run_dir / MOVIE_CASE_FILE)
    for filename in MOVIE_CASE_SUPPORT:
        shutil.copy2(MOVIE_CASE_DIR / filename, run_dir / filename)


def _stage_paul_case(run_dir: Path) -> None:
    shutil.copytree(PAUL_CASE_DIR, run_dir)


def _run_case(exe: Path, args: list[str], run_dir: Path,
              mpi_ranks: int = 1) -> str:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    if mpi_ranks > 1:
        command = [*_mpi_launcher(mpi_ranks), str(exe), *args]
    else:
        command = [str(exe), *args]
    completed = subprocess.run(
        command,
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, (
        "CLI process failed\n"
        f"Command: {' '.join(command)}\n"
        f"stdout:\n{completed.stdout}\n"
        f"stderr:\n{completed.stderr}\n"
    )
    return (completed.stdout + completed.stderr).replace("\r\n", "\n").replace(
        "\r", "\n")


def _assert_markers_in_order(text: str, markers: list[str], label: str) -> None:
    cursor = -1
    for marker in markers:
        pos = text.find(marker)
        assert pos >= 0, f"{label}: missing marker: {marker}"
        assert pos > cursor, (
            f"{label}: marker order mismatch around '{marker}'.\n{text}"
        )
        cursor = pos


def _assert_common_cli_markers(text: str, case_file: str, label: str) -> None:
    _assert_markers_in_order(
        text,
        [
            "semba-fdtd",
            "Compilation date:",
            "Compiler Id:",
            "git commit:",
            "cmake build type:",
            "cmake compilation flags:",
            f"INIT interpreting geometrical data from {case_file}",
            "Switches ",
        ],
        label,
    )


def _assert_movie_outputs(run_dir: Path, require_h5bin: bool) -> None:
    names = [path.name for path in run_dir.iterdir() if path.is_file()]
    assert any("_electric_field_movie_" in name and name.endswith(".bin")
               for name in names), "Missing movie .bin output"
    assert any(name.endswith("_time.xdmf") for name in names), (
        "Missing movie _time.xdmf output"
    )
    if require_h5bin:
        assert any(name.endswith(".h5bin") for name in names), (
            "Missing movie .h5bin output"
        )
    else:
        assert any(name.endswith(".h5") for name in names), (
            "Missing movie .h5 output"
        )


@pytest.mark.mtln
@pytest.mark.movie
def test_cli_mtln_fdtd_movie_messages_match_fortran_markers(
    tmp_path: Path,
    cpp_mtln_exe: Path,
    fortran_mtln_exe: Path,
) -> None:
    fortran_dir = tmp_path / "fortran_mtln_movie"
    cpp_dir = tmp_path / "cpp_mtln_movie"
    _stage_movie_case(fortran_dir)
    _stage_movie_case(cpp_dir)

    args = ["-i", MOVIE_CASE_FILE, "-n", "1200"]
    fortran_text = _run_case(fortran_mtln_exe, args, fortran_dir)
    cpp_text = _run_case(cpp_mtln_exe, args, cpp_dir)

    _assert_common_cli_markers(fortran_text, MOVIE_CASE_FILE, "fortran-mtln")
    _assert_common_cli_markers(cpp_text, MOVIE_CASE_FILE, "cpp-mtln")

    assert "END PREPROCESSING. STARTING simulation." in cpp_text
    assert "Next info at step:" in cpp_text
    assert "Mcells/sec  :" in cpp_text
    assert "Running FDTD:" not in cpp_text
    assert "Step 0/" not in cpp_text

    _assert_movie_outputs(cpp_dir, require_h5bin=False)


@pytest.mark.mpi
@pytest.mark.movie
def test_cli_mpi_fdtd_movie_messages_match_fortran_markers(
    tmp_path: Path,
    cpp_mpi_exe: Path,
    fortran_mpi_exe: Path,
) -> None:
    fortran_dir = tmp_path / "fortran_mpi_movie"
    cpp_dir = tmp_path / "cpp_mpi_movie"
    _stage_movie_case(fortran_dir)
    _stage_movie_case(cpp_dir)

    args = ["-i", MOVIE_CASE_FILE, "-n", "1200", "-mpidir", "z"]
    fortran_text = _run_case(fortran_mpi_exe, args, fortran_dir, mpi_ranks=2)
    cpp_text = _run_case(cpp_mpi_exe, args, cpp_dir, mpi_ranks=2)

    _assert_common_cli_markers(fortran_text, MOVIE_CASE_FILE, "fortran-mpi")
    _assert_common_cli_markers(cpp_text, MOVIE_CASE_FILE, "cpp-mpi")

    assert "END PREPROCESSING. STARTING simulation." in cpp_text
    assert "Next info at step:" in cpp_text
    assert "Mcells/sec  :" in cpp_text
    assert "Running FDTD:" not in cpp_text
    assert "Step 0/" not in cpp_text

    # Root-only user-facing logging: markers should not duplicate per rank.
    assert cpp_text.count("Compilation date:") == 1
    assert cpp_text.count(
        f"INIT interpreting geometrical data from {MOVIE_CASE_FILE}"
    ) == 1

    _assert_movie_outputs(cpp_dir, require_h5bin=True)


@pytest.mark.mtln
@pytest.mark.mtln_standalone
def test_cli_mtln_standalone_messages_match_fortran_markers(
    tmp_path: Path,
    cpp_mtln_exe: Path,
    fortran_mtln_exe: Path,
) -> None:
    fortran_dir = tmp_path / "fortran_mtln_standalone"
    cpp_dir = tmp_path / "cpp_mtln_standalone"
    _stage_paul_case(fortran_dir)
    _stage_paul_case(cpp_dir)

    args = ["-i", PAUL_CASE_FILE]
    fortran_text = _run_case(fortran_mtln_exe, args, fortran_dir)
    cpp_text = _run_case(cpp_mtln_exe, args, cpp_dir)

    _assert_common_cli_markers(fortran_text, PAUL_CASE_FILE, "fortran-mtln-standalone")
    _assert_common_cli_markers(cpp_text, PAUL_CASE_FILE, "cpp-mtln-standalone")
    assert "MTLN simulation finished." in fortran_text
    assert "MTLN simulation finished." in cpp_text

    report_file = cpp_dir / "paul_8_6_square.fdtd_Report.txt"
    assert report_file.exists(), "MTLN standalone report file was not generated"
    report_text = report_file.read_text(errors="replace")
    _assert_common_cli_markers(
        report_text, PAUL_CASE_FILE, "cpp-mtln-standalone-report")
    assert "MTLN simulation finished." in report_text

