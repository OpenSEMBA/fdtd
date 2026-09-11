"""Shared staging helpers for output end-to-end cases."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Callable

import h5py
import pyvista as pv
import pytest
from vtkmodules.vtkIOXdmf2 import vtkXdmfReader

from test.utils.build_resolver import solver_executable


PROJECT_ROOT = Path(__file__).resolve().parents[3]
OUTPUT_CASES = PROJECT_ROOT / "testData" / "cases" / "output_e2e"
EXCITATIONS = PROJECT_ROOT / "testData" / "excitations"


def assert_static_point_attributes(
    xdmf_path: Path,
    attribute_names: tuple[str, ...],
) -> None:
    """Require node attributes backed by one non-series HDF5 dataset."""
    root = ET.parse(xdmf_path).getroot()
    with h5py.File(xdmf_path.with_suffix(".h5"), "r") as hdf5:
        for name in attribute_names:
            dataset_paths = []
            for attribute in root.findall(f'.//Attribute[@Name="{name}"]'):
                data_item = attribute.find("DataItem")
                assert data_item is not None
                assert data_item.get("ItemType") != "HyperSlab", (
                    f"Node attribute {name} must be static"
                )
                assert data_item.text is not None
                dataset_paths.append(data_item.text.strip().split(":", maxsplit=1)[1])

            assert dataset_paths, f"Missing static node attribute {name}"
            assert len(set(dataset_paths)) == 1
            assert hdf5[dataset_paths[0]].ndim in (1, 3)


def read_xdmf_point_data_names(xdmf_path: Path) -> set[str]:
    """Read all leaf point arrays without assuming a time-series collection."""
    reader = vtkXdmfReader()
    reader.SetFileName(str(xdmf_path))
    reader.Update()
    pending = [pv.wrap(reader.GetOutputDataObject(0))]
    names = set()
    while pending:
        dataset = pending.pop()
        if isinstance(dataset, pv.MultiBlock):
            pending.extend(dataset)
        else:
            names.update(dataset.point_data)
    return names


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

        executable = solver_executable(PROJECT_ROOT)
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
