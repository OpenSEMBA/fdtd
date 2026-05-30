from dataclasses import dataclass
from pathlib import Path
import os
import re
import shutil
import subprocess

import pytest

from utils import CASES_FOLDER, SEMBA_EXE


pytestmark = [
    pytest.mark.cpp_migration,
    pytest.mark.planewave,
]


CASE_DIR = Path(CASES_FOLDER) / "planewave"
CASE_FILE = "pw-in-box-pec.fdtd.json"
CASE_SUPPORT_FILES = ("gauss_1GHz.exc",)

TEXT_EXTENSIONS = {".txt", ".dat", ".vtk", ".json", ".log", ".pl", ".xdmf", ".exc"}
VOLATILE_TEXT_PATTERNS = (
    re.compile(r"^Launched on\s+"),
    re.compile(r"^Compilation date:\s+"),
    re.compile(r"^Compiler Id:\s+"),
    re.compile(r"^git commit:\s+"),
    re.compile(r"^cmake build type:\s+"),
    re.compile(r"^cmake compilation flags:\s+"),
    re.compile(r"^Start Date/time\s+"),
    re.compile(r"^Date/time\s+"),
    re.compile(r"^BEGUN\s+"),
    re.compile(r"^ENDED\s+"),
)


SCENARIOS = (
    ("noargs", []),
    ("help", ["-h"]),
    ("input", ["-i", CASE_FILE]),
    ("input_missing", ["-i"]),
    ("unknown_flag", ["-zzz"]),
    ("mapvtk", ["-i", CASE_FILE, "-mapvtk"]),
    ("n2", ["-i", CASE_FILE, "-n", "2"]),
    ("prefix", ["-i", CASE_FILE, "-prefix", "ABC"]),
    ("mpidir_x", ["-i", CASE_FILE, "-mpidir", "x"]),
)


@dataclass(frozen=True)
class RunResult:
    returncode: int
    stdout: str
    stderr: str
    generated_files: tuple[str, ...]
    run_dir: Path
    exe_path: Path


def _resolve_fortran_exe() -> Path:
    env_path = os.environ.get("SEMBA_FORTRAN_EXE")
    if env_path:
        return Path(env_path)
    repo_root = Path(CASES_FOLDER).parent.parent
    exe_name = "semba-fdtd.exe" if os.name == "nt" else "semba-fdtd"
    for candidate_dir in (
        "build_fortran_nomtln",
        "build_fortran_nomtln_rel",
        "build_fortran_rel",
        "build_fortran",
        "build",
    ):
        candidate = repo_root / candidate_dir / "bin" / exe_name
        if candidate.exists():
            return candidate
    return repo_root / "build_fortran_nomtln" / "bin" / exe_name


@pytest.fixture(scope="session")
def fortran_exe() -> Path:
    exe = _resolve_fortran_exe().resolve()
    assert exe.exists(), (
        "Fortran executable is required for CLI parity tests. "
        f"Set SEMBA_FORTRAN_EXE or build Fortran binary. Missing: {exe}"
    )
    return exe


@pytest.fixture(scope="session")
def cpp_exe() -> Path:
    exe = Path(SEMBA_EXE).resolve()
    assert exe.exists(), f"C++ executable not found: {exe}"
    return exe


def _stage_case(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(CASE_DIR / CASE_FILE, run_dir / CASE_FILE)
    for filename in CASE_SUPPORT_FILES:
        shutil.copy2(CASE_DIR / filename, run_dir / filename)


def _list_files(run_dir: Path) -> set[str]:
    return {
        path.relative_to(run_dir).as_posix()
        for path in run_dir.rglob("*")
        if path.is_file()
    }


def _run_case(exe_path: Path, args: list[str], run_dir: Path) -> RunResult:
    _stage_case(run_dir)
    before = _list_files(run_dir)
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    completed = subprocess.run(
        [str(exe_path), *args],
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
    )
    after = _list_files(run_dir)
    generated = tuple(sorted(after - before))
    return RunResult(
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
        generated_files=generated,
        run_dir=run_dir,
        exe_path=exe_path,
    )


def _is_volatile_text_file(relative_path: str) -> bool:
    name = Path(relative_path).name
    return (
        name == "SEMBA_FDTD_temp.log"
        or name.endswith("_Report.txt")
        or name.endswith("_Warnings.txt")
    )


def _is_text_file(path: Path) -> bool:
    return path.suffix.lower() in TEXT_EXTENSIONS


def _normalize_text(
    text: str,
    run_dir: Path,
    exe_path: Path,
    drop_volatile_lines: bool,
) -> str:
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    normalized = normalized.replace(str(run_dir), "<RUN_DIR>")
    normalized = normalized.replace(str(exe_path), "<EXE>")
    lines: list[str] = []
    for raw_line in normalized.split("\n"):
        line = raw_line.rstrip()
        if drop_volatile_lines and any(pat.match(line) for pat in VOLATILE_TEXT_PATTERNS):
            continue
        lines.append(line)
    return "\n".join(lines).strip()


def _first_text_diff(expected: str, got: str) -> str:
    exp_lines = expected.splitlines()
    got_lines = got.splitlines()
    for line_no, (exp_line, got_line) in enumerate(zip(exp_lines, got_lines), start=1):
        if exp_line == got_line:
            continue
        return (
            f"line {line_no} differs\n"
            f"expected: {exp_line}\n"
            f"got:      {got_line}"
        )
    return (
        f"line count differs: expected {len(exp_lines)} lines "
        f"got {len(got_lines)} lines"
    )


def _compare_generated_files(fortran_run: RunResult, cpp_run: RunResult) -> None:
    assert cpp_run.generated_files == fortran_run.generated_files

    for rel_path in fortran_run.generated_files:
        fortran_file = fortran_run.run_dir / rel_path
        cpp_file = cpp_run.run_dir / rel_path
        assert cpp_file.exists(), f"missing C++ artifact: {rel_path}"

        if _is_volatile_text_file(rel_path):
            expected = _normalize_text(
                fortran_file.read_text(errors="replace"),
                fortran_run.run_dir,
                fortran_run.exe_path,
                drop_volatile_lines=True,
            )
            got = _normalize_text(
                cpp_file.read_text(errors="replace"),
                cpp_run.run_dir,
                cpp_run.exe_path,
                drop_volatile_lines=True,
            )
            assert got == expected, (
                f"normalized volatile text mismatch for {rel_path}\n"
                f"{_first_text_diff(expected, got)}"
            )
        elif _is_text_file(fortran_file):
            expected_bytes = fortran_file.read_bytes()
            got_bytes = cpp_file.read_bytes()
            assert got_bytes == expected_bytes, f"byte mismatch for {rel_path}"
        else:
            expected_bytes = fortran_file.read_bytes()
            got_bytes = cpp_file.read_bytes()
            assert got_bytes == expected_bytes, f"byte mismatch for {rel_path}"


@pytest.mark.parametrize(("scenario_name", "args"), SCENARIOS)
def test_cli_console_and_outputs_match_fortran(
    tmp_path: Path,
    scenario_name: str,
    args: list[str],
    fortran_exe: Path,
    cpp_exe: Path,
) -> None:
    fortran_run = _run_case(fortran_exe, args, tmp_path / f"fortran_{scenario_name}")
    cpp_run = _run_case(cpp_exe, args, tmp_path / f"cpp_{scenario_name}")

    assert cpp_run.returncode == fortran_run.returncode

    expected_stdout = _normalize_text(
        fortran_run.stdout,
        fortran_run.run_dir,
        fortran_run.exe_path,
        drop_volatile_lines=True,
    )
    got_stdout = _normalize_text(
        cpp_run.stdout,
        cpp_run.run_dir,
        cpp_run.exe_path,
        drop_volatile_lines=True,
    )
    assert got_stdout == expected_stdout, _first_text_diff(expected_stdout, got_stdout)

    expected_stderr = _normalize_text(
        fortran_run.stderr,
        fortran_run.run_dir,
        fortran_run.exe_path,
        drop_volatile_lines=True,
    )
    got_stderr = _normalize_text(
        cpp_run.stderr,
        cpp_run.run_dir,
        cpp_run.exe_path,
        drop_volatile_lines=True,
    )
    assert got_stderr == expected_stderr, _first_text_diff(expected_stderr, got_stderr)

    _compare_generated_files(fortran_run, cpp_run)


def _looks_like_simulation_output(rel_path: str) -> bool:
    lower = rel_path.lower()
    return lower.endswith((".dat", ".bin", ".vtk", ".h5", ".xdmf", ".pl"))


def test_cpp_does_not_default_to_simulation_without_input(
    tmp_path: Path,
    fortran_exe: Path,
    cpp_exe: Path,
) -> None:
    fortran_run = _run_case(fortran_exe, [], tmp_path / "fortran_no_input")
    cpp_run = _run_case(cpp_exe, [], tmp_path / "cpp_no_input")

    assert cpp_run.returncode == fortran_run.returncode

    cpp_console = _normalize_text(
        cpp_run.stdout + "\n" + cpp_run.stderr,
        cpp_run.run_dir,
        cpp_run.exe_path,
        drop_volatile_lines=True,
    )
    assert "Running FDTD:" not in cpp_console

    fortran_sim_outputs = sorted(
        rel for rel in fortran_run.generated_files if _looks_like_simulation_output(rel)
    )
    cpp_sim_outputs = sorted(
        rel for rel in cpp_run.generated_files if _looks_like_simulation_output(rel)
    )
    assert cpp_sim_outputs == fortran_sim_outputs
