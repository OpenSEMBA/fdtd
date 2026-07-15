#!/usr/bin/env pvpython
"""Optional smoke test using ParaView's real XDMF reader."""

from __future__ import annotations

import sys
from pathlib import Path

from paraview import servermanager
from paraview.simple import Delete, XDMFReader


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def read_case(root: Path, name: str, expected_points: int | None = None) -> None:
    reader = XDMFReader(FileNames=[str(root / f"{name}.xdmf")])
    reader.UpdatePipeline()
    information = reader.GetDataInformation()
    require(information.GetNumberOfPoints() > 0, f"{name}: no points")
    if expected_points is not None:
        require(
            information.GetNumberOfPoints() == expected_points,
            f"{name}: expected {expected_points} points",
        )
    Delete(reader)


def first_leaf(data):
    if hasattr(data, "GetPointData"):
        return data
    for index in range(data.GetNumberOfBlocks()):
        block = data.GetBlock(index)
        if block is not None:
            leaf = first_leaf(block)
            if leaf is not None:
                return leaf
    return None


def read_volume(root: Path) -> None:
    reader = XDMFReader(FileNames=[str(root / "volume.xdmf")])
    reader.PointArrayStatus = ["electric-field-magnitude"]
    reader.UpdatePipeline()
    information = reader.GetDataInformation()
    require(information.GetNumberOfPoints() == 21**3, "volume: point count")
    data = first_leaf(servermanager.Fetch(reader))
    require(data is not None, "volume: missing grid")
    field = data.GetPointData().GetArray("electric-field-magnitude")
    require(field is not None, "volume: missing scalar field")
    require(field.GetRange()[0] >= 0.0, "volume: negative field magnitude")
    require(field.GetRange()[1] > 1.0, "volume: field maximum")
    Delete(reader)


def main() -> int:
    root = Path(sys.argv[1]).resolve()
    read_case(root, "uniform", 6)
    read_volume(root)
    read_case(root, "curvilinear", 6)
    read_case(root, "unstructured")
    read_case(root, "mixed", 8)

    reader = XDMFReader(FileNames=[str(root / "time-series.xdmf")])
    times = list(reader.TimestepValues)
    require(len(times) == 20, "time-series: timestep count")
    require(abs(times[0]) < 1.0e-6, "time-series: initial time")
    require(abs(times[-1] - 1.0) < 1.0e-6, "time-series: final time")
    reader.PointArrayStatus = ["electric-field-magnitude", "electric-field"]
    reader.UpdatePipeline(times[0])
    first_step = first_leaf(servermanager.Fetch(reader))
    require(first_step is not None, "time-series: missing first grid")
    first_magnitude = first_step.GetPointData().GetArray("electric-field-magnitude")
    first_electric = first_step.GetPointData().GetArray("electric-field")
    require(first_magnitude is not None, "time-series: missing first magnitude")
    require(first_electric is not None, "time-series: missing first electric field")
    require(first_magnitude.GetRange()[1] > 0.99, "time-series: first pulse maximum")
    require(first_electric.GetNumberOfComponents() == 3, "time-series: vector components")
    reader.UpdatePipeline(times[-1])
    information = reader.GetDataInformation()
    require(information.GetNumberOfPoints() == 25**3, "time-series: point count")
    require(
        set(reader.PointArrayStatus) == {"electric-field-magnitude", "electric-field"},
        "time-series: point arrays",
    )
    last_step = first_leaf(servermanager.Fetch(reader))
    require(last_step is not None, "time-series: missing last grid")
    last_magnitude = last_step.GetPointData().GetArray("electric-field-magnitude")
    require(last_magnitude is not None, "time-series: missing last magnitude")
    require(last_magnitude.GetRange()[1] > 0.99, "time-series: last pulse maximum")
    Delete(reader)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
