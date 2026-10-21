"""End-to-end checks that sources with empty excitation files abort preprocessing."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pytest

from test.utils.build_resolver import build_feature_enabled, solver_executable


PROJECT_ROOT = Path(__file__).resolve().parents[3]
CASES = PROJECT_ROOT / "testData" / "cases"
SEMBA_EXE = solver_executable(PROJECT_ROOT)

no_mtln_skip = pytest.mark.skipif(
    not build_feature_enabled("SEMBA_FDTD_ENABLE_MTLN"),
    reason="MTLN is not available",
)


def _stage_case(tmp_path: Path, case_path: Path, excitations: list[str]) -> Path:
    """Copy a case and its excitation files into an isolated solver folder."""
    input_path = tmp_path / case_path.name
    shutil.copy2(case_path, input_path)
    for excitation in excitations:
        shutil.copy2(case_path.parent / excitation, tmp_path / excitation)
    return input_path


def _empty_file(path: Path) -> None:
    path.write_text("")


def _run_solver(cwd: Path, input_path: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(SEMBA_EXE), "-i", input_path.name],
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


def _assert_solver_reports_empty_file(process, excitation: str) -> None:
    assert process.returncode != 0, process.stdout + process.stderr
    output = (process.stdout + process.stderr).lower()
    assert excitation.lower() in output
    assert "empty" in output


def test_empty_excitation_file_fails_in_main_solver(tmp_path):
    """A nodal source whose magnitude file is empty must fail during preprocessing."""
    case = CASES / "nodalSource" / "nodalSource.fdtd.json"
    excitation = "predefinedExcitation.1.exc"
    input_path = _stage_case(tmp_path, case, [excitation])
    _empty_file(tmp_path / excitation)

    process = _run_solver(tmp_path, input_path)

    _assert_solver_reports_empty_file(process, excitation)


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
def test_empty_excitation_file_fails_in_mtln_terminal_generator(tmp_path):
    """An MTLN generator on a terminal node must fail when its file is empty."""
    case = CASES / "paul" / "paul_8_6_square.fdtd.json"
    excitation = "coaxial_line_paul_8_6_0.25_square.exc"
    input_path = _stage_case(tmp_path, case, [excitation])
    _empty_file(tmp_path / excitation)

    process = _run_solver(tmp_path, input_path)

    _assert_solver_reports_empty_file(process, excitation)


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
def test_empty_excitation_file_fails_in_mtln_wire_generator(tmp_path):
    """An MTLN generator on an interior wire node must fail when its file is empty."""
    case = CASES / "sources" / "sources_voltage.fdtd.json"
    excitation = "voltage_source_50V.exc"
    input_path = _stage_case(tmp_path, case, [excitation])
    _empty_file(tmp_path / excitation)

    process = _run_solver(tmp_path, input_path)

    _assert_solver_reports_empty_file(process, excitation)
