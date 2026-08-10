"""Shared staging helpers for output end-to-end cases."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Callable

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_CASES = PROJECT_ROOT / "testData" / "cases" / "output_e2e"
EXCITATIONS = PROJECT_ROOT / "testData" / "excitations"


@pytest.fixture
def stage_output_case(tmp_path: Path) -> Callable[[str], Path]:
    """Copy an input and its source magnitudes into an isolated solver folder."""

    def stage(case_name: str) -> Path:
        source = OUTPUT_CASES / case_name
        destination = tmp_path / case_name
        shutil.copy2(source, destination)

        with destination.open(encoding="utf-8") as input_file:
            case = json.load(input_file)
        for source_definition in case.get("sources", []):
            magnitude = source_definition.get("magnitudeFile")
            if magnitude:
                shutil.copy2(EXCITATIONS / Path(magnitude).name, tmp_path / Path(magnitude).name)
        return destination

    return stage


@pytest.fixture
def output_root(tmp_path: Path) -> Path:
    """Return a nested output path containing a space for path tests."""

    return tmp_path / "nested output" / "results"


@pytest.fixture
def failed_output_root(tmp_path: Path) -> Path:
    """Return an isolated path used by publication failure tests."""

    return tmp_path / "failed output" / "results"


@pytest.fixture
def run_output_case(
    stage_output_case: Callable[[str], Path],
    tmp_path: Path,
) -> Callable[[str, list[dict], str], tuple[subprocess.CompletedProcess[str], Path]]:
    """Run a staged output case and return its process and output directory."""

    def run(
        case_name: str,
        probes: list[dict],
        additional_arguments: str = "",
    ) -> tuple[subprocess.CompletedProcess, Path]:
        input_path = stage_output_case("common_geometry.fdtd.json")
        with input_path.open(encoding="utf-8") as input_file:
            case = json.load(input_file)
        case["probes"] = probes
        if case_name == "wire":
            case["materials"] = [
                {
                    "id": 1,
                    "type": "wire",
                    "radius": 0.01,
                    "resistancePerMeter": 0.0,
                    "inductancePerMeter": 0.0,
                },
                {
                    "id": 2,
                    "type": "terminal",
                    "terminations": [{"type": "open"}],
                },
            ]
            case["materialAssociations"] = [
                {
                    "name": "e2e_wire",
                    "materialId": 1,
                    "initialTerminalId": 2,
                    "endTerminalId": 2,
                    "elementIds": [2],
                }
            ]
        if additional_arguments:
            case.setdefault("general", {})["additionalArguments"] = additional_arguments
        with input_path.open("w", encoding="utf-8") as input_file:
            json.dump(case, input_file)

        executable_name = "semba-fdtd.exe" if os.name == "nt" else "semba-fdtd"
        executable = PROJECT_ROOT / "build" / "bin" / executable_name
        process = subprocess.run(
            [str(executable), "-i", str(input_path)],
            cwd=tmp_path,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
            capture_output=True,
            text=True,
            check=False,
        )
        return process, tmp_path

    return run
