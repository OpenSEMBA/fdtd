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


def main() -> int:
    root = Path(sys.argv[1]).resolve()
    read_case(root, "uniform", 6)
    read_case(root, "curvilinear", 6)
    read_case(root, "unstructured")
    read_case(root, "mixed", 8)

    reader = XDMFReader(FileNames=[str(root / "time-series.xdmf")])
    require(list(reader.TimestepValues) == [0.1, 0.25], "time-series: timesteps")
    reader.PointArrayStatus = ["pressure", "velocity"]
    reader.UpdatePipeline(0.1)
    first_step = first_leaf(servermanager.Fetch(reader))
    require(first_step is not None, "time-series: missing first grid")
    first_pressure = first_step.GetPointData().GetArray("pressure")
    require(first_pressure is not None, "time-series: missing first pressure")
    require(first_pressure.GetRange() == (101.0, 124.0), "time-series: first values")
    reader.UpdatePipeline(0.25)
    information = reader.GetDataInformation()
    require(information.GetNumberOfPoints() == 24, "time-series: point count")
    require(
        set(reader.PointArrayStatus) == {"pressure", "velocity"},
        "time-series: point arrays",
    )
    second_step = first_leaf(servermanager.Fetch(reader))
    require(second_step is not None, "time-series: missing second grid")
    second_pressure = second_step.GetPointData().GetArray("pressure")
    require(second_pressure is not None, "time-series: missing second pressure")
    require(second_pressure.GetRange() == (201.0, 224.0), "time-series: second values")
    Delete(reader)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
