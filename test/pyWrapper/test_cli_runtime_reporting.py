from pathlib import Path
import os
import pty
import re
import select
import shutil
import signal
import subprocess
import time

import pytest

from utils import CASES_FOLDER, SEMBA_EXE


pytestmark = [pytest.mark.cpp_migration]


PLANEWAVE_DIR = Path(CASES_FOLDER) / "planewave"
PLANEWAVE_CASE = "pw-in-box-pec.fdtd.json"
PLANEWAVE_SUPPORT = ("gauss_1GHz.exc",)
CONFORMAL_DIR = Path(CASES_FOLDER) / "conformal_impedance_cylinder"
CONFORMAL_CASE = "conformal_impedance_cylinder_conformal.fdtd.json"

METADATA_KEYS = (
    "Compilation date",
    "Compiler Id",
    "git commit",
    "cmake build type",
    "cmake compilation flags",
)
PLACEHOLDER_VALUES = {"cpp migration", "CXX", "cpp"}


def _cpp_exe() -> Path:
    exe = Path(SEMBA_EXE).resolve()
    assert exe.exists(), f"C++ executable not found: {exe}"
    return exe


def _normalize_text(text: str) -> str:
    return text.replace("\r\n", "\n").replace("\r", "\n")


def _extract_metadata(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in _normalize_text(text).split("\n"):
        if ":" not in line:
            continue
        left, right = line.split(":", 1)
        key = left.strip()
        if key in METADATA_KEYS and key not in values:
            values[key] = right.strip()
    return values


def _assert_real_metadata(values: dict[str, str]) -> None:
    missing = [key for key in METADATA_KEYS if key not in values]
    assert not missing, f"missing metadata keys: {missing}"
    for key in METADATA_KEYS:
        value = values[key]
        assert value, f"empty metadata value for {key}"
        assert value not in PLACEHOLDER_VALUES, (
            f"placeholder metadata for {key}: {value}"
        )
    assert values["cmake build type"].lower() != "cpp"
    assert values["Compiler Id"].lower() != "cxx"
    assert values["cmake compilation flags"].lower() != "cpp"


def _stage_planewave(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(PLANEWAVE_DIR / PLANEWAVE_CASE, run_dir / PLANEWAVE_CASE)
    for filename in PLANEWAVE_SUPPORT:
        shutil.copy2(PLANEWAVE_DIR / filename, run_dir / filename)


def _stage_conformal(run_dir: Path) -> None:
    shutil.copytree(CONFORMAL_DIR, run_dir)


def _run_with_pty(
    exe: Path,
    args: list[str],
    cwd: Path,
    read_timeout_s: float = 1.0,
) -> tuple[str, bool]:
    master_fd, slave_fd = pty.openpty()
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    proc = subprocess.Popen(
        [str(exe), *args],
        cwd=cwd,
        env=env,
        stdin=subprocess.DEVNULL,
        stdout=slave_fd,
        stderr=slave_fd,
        text=False,
        start_new_session=True,
        close_fds=True,
    )
    os.close(slave_fd)

    output = bytearray()
    deadline = time.time() + read_timeout_s
    while time.time() < deadline and proc.poll() is None:
        ready, _, _ = select.select([master_fd], [], [], 0.05)
        if not ready:
            continue
        try:
            chunk = os.read(master_fd, 8192)
        except OSError as exc:
            if exc.errno == 5:
                break
            raise
        if not chunk:
            break
        output.extend(chunk)

    still_running = proc.poll() is None
    if still_running:
        os.killpg(proc.pid, signal.SIGTERM)
        try:
            proc.wait(timeout=2.0)
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait(timeout=2.0)

    while True:
        ready, _, _ = select.select([master_fd], [], [], 0.01)
        if not ready:
            break
        try:
            chunk = os.read(master_fd, 8192)
        except OSError as exc:
            if exc.errno == 5:
                break
            raise
        if not chunk:
            break
        output.extend(chunk)
    os.close(master_fd)
    return output.decode(errors="replace"), still_running


@pytest.mark.planewave
def test_cpp_cli_uses_real_build_metadata_in_stdout_and_report(tmp_path):
    run_dir = tmp_path / "meta"
    _stage_planewave(run_dir)
    exe = _cpp_exe()
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    completed = subprocess.run(
        [str(exe), "-i", PLANEWAVE_CASE, "-n", "2"],
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr

    stdout_meta = _extract_metadata(completed.stdout)
    _assert_real_metadata(stdout_meta)

    report_path = run_dir / "pw-in-box-pec.fdtd_Report.txt"
    assert report_path.exists(), "expected _Report.txt was not generated"
    report_text = report_path.read_text(errors="replace")
    report_meta = _extract_metadata(report_text)
    _assert_real_metadata(report_meta)

    for key in METADATA_KEYS:
        assert report_meta[key] == stdout_meta[key], (
            f"{key} mismatch stdout/report: {stdout_meta[key]} vs {report_meta[key]}"
        )


@pytest.mark.conformal
def test_cpp_cli_streams_startup_before_long_run_finishes(tmp_path):
    run_dir = tmp_path / "live"
    _stage_conformal(run_dir)
    exe = _cpp_exe()
    output, still_running = _run_with_pty(
        exe,
        ["-i", CONFORMAL_CASE],
        run_dir,
        read_timeout_s=1.0,
    )
    normalized = _normalize_text(output)
    assert still_running, (
        "expected a long-running process to still be active after 1s; "
        "choose a heavier case if this becomes flaky"
    )
    assert normalized.strip(), "no streamed output captured while simulation was running"
    assert "semba-fdtd" in normalized
    assert "Compilation date:" in normalized
    assert (
        f"INIT interpreting geometrical data from {CONFORMAL_CASE}" in normalized
    )
    assert "Compilation date: cpp migration" not in normalized
    assert "Compiler Id: CXX" not in normalized
    assert "cmake build type: cpp" not in normalized


@pytest.mark.planewave
def test_cpp_cli_emits_runtime_status_and_speed_reports(tmp_path):
    run_dir = tmp_path / "status"
    _stage_planewave(run_dir)
    exe = _cpp_exe()
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    completed = subprocess.run(
        [str(exe), "-i", PLANEWAVE_CASE, "-n", "1200"],
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr

    stdout = _normalize_text(completed.stdout)
    assert re.search(r"Next info at step:\s+\d+", stdout), (
        "missing runtime status line 'Next info at step'"
    )
    assert "Mcells/sec  :" in stdout, "missing runtime speed report lines"
