"""Shared staging helpers for output end-to-end cases."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_CASES = PROJECT_ROOT / "testData" / "cases" / "output_e2e"
EXCITATIONS = PROJECT_ROOT / "testData" / "excitations"


@pytest.fixture
def stage_output_case(tmp_path):
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
