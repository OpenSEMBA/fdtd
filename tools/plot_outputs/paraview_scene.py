"""Render probe locations, optionally with the problem geometry, using ParaView.

Run by plot_outputs.py as:
  pvbatch paraview_scene.py probes.csv image.png [--geometry piece.vtu ...]
"""

import argparse

from paraview.simple import *  # noqa: F403


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("probes", help="CSV with x, y, z, name and peak columns")
    parser.add_argument("image", help="PNG file to write")
    parser.add_argument("--geometry", nargs="+", default=[], help="geometry .vtu pieces from a __MAP folder")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    view = GetActiveViewOrCreate("RenderView")
    view.OrientationAxesVisibility = 1

    if arguments.geometry:
        geometry = XMLUnstructuredGridReader(FileName=arguments.geometry)
        geometry_display = Show(geometry, view)
        geometry_display.SetRepresentationType("Surface With Edges")
        geometry_display.Opacity = 0.2
        geometry_display.ColorArrayName = [None, ""]

    table = CSVReader(FileName=[arguments.probes])
    points = TableToPoints(Input=table, XColumn="x", YColumn="y", ZColumn="z")
    points.KeepAllDataArrays = 1
    points.UpdatePipeline()
    bounds = points.GetDataInformation().GetBounds()
    extent = max(bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4])
    radius = 0.02 * extent or 0.01

    # Glyph defaults to sampling random points inside the bounds; one sphere per probe is required.
    probes = Glyph(
        Input=points,
        GlyphType="Sphere",
        GlyphMode="All Points",
        ScaleArray=["POINTS", "No scale array"],
        ScaleFactor=2 * radius,
    )
    probes_display = Show(probes, view)
    ColorBy(probes_display, ("POINTS", "peak"))
    probes_display.RescaleTransferFunctionToDataRange(True, False)
    probes_display.SetScalarBarVisibility(view, True)

    # The default camera looks straight down z and flattens the 3D layout.
    ResetCamera(view)
    focal = view.CameraFocalPoint
    distance = sum((p - f) ** 2 for p, f in zip(view.CameraPosition, focal)) ** 0.5
    direction = [1.0, -1.0, 0.8]
    norm = sum(d * d for d in direction) ** 0.5
    view.CameraPosition = [f + d / norm * distance for f, d in zip(focal, direction)]
    view.CameraViewUp = [0.0, 0.0, 1.0]
    ResetCamera(view)
    Render(view)
    SaveScreenshot(arguments.image, view, ImageResolution=[1600, 1000])


if __name__ == "__main__":
    main()
