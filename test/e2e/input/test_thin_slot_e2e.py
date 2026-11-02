"""End-to-end checks for surface and zero-length thin slots."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from test.utils.build_resolver import solver_executable


PROJECT_ROOT = Path(__file__).resolve().parents[3]
EXAMPLES = PROJECT_ROOT / "testData" / "input_examples"
SEMBA_EXE = solver_executable(PROJECT_ROOT)


def _run_example(tmp_path: Path, filename: str) -> subprocess.CompletedProcess[str]:
    case = EXAMPLES / filename
    input_path = tmp_path / case.name
    shutil.copy2(case, input_path)
    return subprocess.run(
        [str(SEMBA_EXE), "-i", input_path.name],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )


def test_surface_thin_slot_is_accepted(tmp_path):
    """A surface-defined slot must expand during preprocessing and run."""
    process = _run_example(tmp_path, "thinSlot_surface.fdtd.json")

    assert process.returncode == 0, process.stdout + process.stderr


def test_negative_normal_surface_thin_slot_is_accepted(tmp_path):
    """A surface slot with a negative normal must expand during preprocessing."""
    process = _run_example(tmp_path, "thinSlot_surface_negative.fdtd.json")

    assert process.returncode == 0, process.stdout + process.stderr


def test_inconsistent_surface_thin_slot_normals_fail(tmp_path):
    """Adjacent coplanar faces with opposite normals must abort preprocessing."""
    process = _run_example(tmp_path, "thinSlot_surface_inconsistent.fdtd.json")

    assert process.returncode != 0, process.stdout + process.stderr
    output = (process.stdout + process.stderr).lower()
    assert "inconsistent normal orientation" in output


def test_zero_length_thin_slot_is_accepted(tmp_path):
    """A mesher-collapsed slot must parse and preprocess without components."""
    process = _run_example(tmp_path, "thinSlot_zeroLength.fdtd.json")

    assert process.returncode == 0, process.stdout + process.stderr
