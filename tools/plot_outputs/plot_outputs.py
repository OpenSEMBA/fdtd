#!/usr/bin/env python3
"""Plot semba-fdtd probe outputs.

Commands:
  probe      one PNG per probe .dat file
  group      one PNG per group of probes sharing a figure (see --group)
  3d-probes  probe locations in 3D, rendered with ParaView (pvbatch)
  3d-scene   probe locations plus the problem geometry, rendered with ParaView

Every command takes the solver output directory (or one probe file for
"probe") and writes PNG files to --output, which defaults to <source>/plots.
"""

import argparse
import csv
import json
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

PARAVIEW_SCRIPT = Path(__file__).with_name("paraview_scene.py")
# <case>.fdtd_<name>_<component>_<i>_<j>_<k>[__<i>_<j>_<k>]_<tm|fq>.dat
PROBE_FILE = re.compile(
    r"\.fdtd_(?P<name>.+)_(?P<component>[EHJ][xyz])_(?P<i>\d+)_(?P<j>\d+)_(?P<k>\d+)"
    r"(?:__\d+_\d+_\d+)?_(?P<domain>tm|fq)\.dat$"
)


@dataclass(frozen=True)
class Probe:
    path: Path
    name: str
    component: str
    index: tuple
    domain: str


@dataclass(frozen=True)
class Series:
    probe: Probe
    columns: list
    data: np.ndarray


def parse_probe(path):
    match = PROBE_FILE.search(path.name)
    if match is None:
        return None
    index = tuple(int(match[axis]) for axis in "ijk")
    return Probe(path, match["name"], match["component"], index, match["domain"])


def load_series(probe):
    lines = probe.path.read_text(encoding="utf-8").splitlines()
    rows = [line.split() for line in lines[1:] if line.strip()]
    if not rows:
        return None
    return Series(probe, lines[0].split(), np.array(rows, dtype=float))


def load_all(source):
    files = [source] if source.is_file() else sorted(source.glob("*.dat"))
    series = []
    for path in files:
        probe = parse_probe(path)
        if probe is None:
            continue
        loaded = load_series(probe)
        if loaded is None:
            print(f"skipped {path.name}: no samples")
        else:
            series.append(loaded)
    if not series:
        raise SystemExit(f"No supported probe files in {source}")
    return series


def traces(series):
    """Columns to draw: field/current over time, or magnitudes over frequency."""
    if series.probe.domain == "fq":
        return [(i, name) for i, name in enumerate(series.columns) if name.endswith("magnitude")]
    return [(i, name) for i, name in enumerate(series.columns) if i > 0]


def plot_series(series_list, title, domain, output):
    figure, axis = plt.subplots(figsize=(8, 4.5), layout="constrained")
    for series in series_list:
        x = series.data[:, 0]
        for i, column in traces(series):
            axis.plot(x, series.data[:, i], label=f"{series.probe.name} {series.probe.component} {column}")
    if domain == "tm":
        axis.set_xlabel("Time (s)")
        axis.set_ylabel("Value")
    else:
        axis.set_xlabel("Frequency (Hz)")
        axis.set_ylabel("Magnitude")
        axis.set_yscale("log")
    axis.set_title(title)
    axis.grid(alpha=0.3)
    axis.legend(fontsize=8)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=150)
    plt.close(figure)
    print(output)


def parse_group(text):
    title, _, names = text.partition("=")
    selected = [name.strip() for name in names.split(",") if name.strip()]
    if not title or not selected:
        raise argparse.ArgumentTypeError(f"expected TITLE=NAME[,NAME...], got '{text}'")
    return title, selected


def slug(text):
    return re.sub(r"\W+", "_", text).strip("_")


def output_directory(arguments):
    base = arguments.source if arguments.source.is_dir() else arguments.source.parent
    return arguments.output or base / "plots"


def command_probe(arguments):
    output = output_directory(arguments)
    for series in load_all(arguments.source):
        probe = series.probe
        plot_series([series], probe.path.stem, probe.domain, output / f"{probe.path.stem}.png")


def command_group(arguments):
    output = output_directory(arguments)
    by_name = defaultdict(list)
    for series in load_all(arguments.source):
        by_name[series.probe.name].append(series)
    groups = dict(arguments.group) if arguments.group else {name: [name] for name in by_name}

    for title, names in groups.items():
        unknown = [name for name in names if name not in by_name]
        if unknown:
            raise SystemExit(f"Unknown probe name(s) {', '.join(unknown)}; available: {', '.join(sorted(by_name))}")
        for domain in ("tm", "fq"):
            selected = [series for name in names for series in by_name[name] if series.probe.domain == domain]
            if selected:
                plot_series(selected, title, domain, output / f"{slug(title)}_{domain}.png")


def grid_edges(source_dir):
    """Cell-vertex coordinates along x, y and z, from the input's mesh.grid."""
    cases = sorted(source_dir.glob("*.fdtd.json"))
    if not cases:
        raise SystemExit(f"No .fdtd.json input in {source_dir}; it is needed to place the probes")
    grid = json.loads(cases[0].read_text(encoding="utf-8"))["mesh"]["grid"]
    origin = grid.get("origin", [0.0, 0.0, 0.0])
    edges = []
    for axis, count, start in zip("xyz", grid["numberOfCells"], origin):
        steps = grid["steps"][axis]
        if len(steps) == 1:
            steps = steps * count
        edges.append(start + np.concatenate(([0.0], np.cumsum(steps))))
    return edges


def probe_csv(source_dir, output):
    """Write one row per probe location: physical position and peak time-domain value."""
    edges = grid_edges(source_dir)
    peaks = {}
    for series in load_all(source_dir):
        if series.probe.domain != "tm":
            continue
        key = (series.probe.name, series.probe.index)
        peaks[key] = max(peaks.get(key, 0.0), float(np.abs(series.data[:, 1:]).max()))
    if not peaks:
        raise SystemExit(f"No time-domain probe files in {source_dir}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["x", "y", "z", "name", "peak"])
        for (name, index), peak in sorted(peaks.items()):
            writer.writerow([*(edges[axis][index[axis]] for axis in range(3)), name, peak])
    return output


def geometry_pieces(source_dir):
    pieces = sorted(piece for folder in source_dir.glob("*__MAP") for piece in folder.glob("*.vtu"))
    if not pieces:
        raise SystemExit(f"No __MAP geometry in {source_dir}; run the solver with -mapvtk")
    return pieces


def command_3d(arguments, with_geometry):
    source = arguments.source if arguments.source.is_dir() else arguments.source.parent
    output = output_directory(arguments)
    pvbatch = shutil.which("pvbatch")
    if pvbatch is None:
        raise SystemExit("pvbatch was not found in PATH; install ParaView to render 3D views")

    csv_path = probe_csv(source, output / "probes.csv")
    image = output / ("scene_3d.png" if with_geometry else "probes_3d.png")
    command = [pvbatch, "--force-offscreen-rendering", str(PARAVIEW_SCRIPT), str(csv_path), str(image)]
    if with_geometry:
        command += ["--geometry", *map(str, geometry_pieces(source))]
    subprocess.run(command, check=True)
    print(image)


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    def add_command(name, run, help_text):
        command = commands.add_parser(name, help=help_text, description=help_text)
        command.add_argument("source", type=Path, help="solver output directory, or one probe .dat file for 'probe'")
        command.add_argument("-o", "--output", type=Path, help="directory for the PNG files (default: <source>/plots)")
        command.set_defaults(run=run)
        return command

    add_command("probe", command_probe, "one plot per probe file")
    group = add_command("group", command_group, "probes sharing one figure")
    group.add_argument(
        "-g",
        "--group",
        action="append",
        type=parse_group,
        metavar="TITLE=NAME[,NAME...]",
        help="probe names drawn in one figure; repeat for more figures (default: one figure per probe name)",
    )
    add_command("3d-probes", lambda arguments: command_3d(arguments, False), "probe locations in 3D (ParaView)")
    add_command("3d-scene", lambda arguments: command_3d(arguments, True), "probes and geometry in 3D (ParaView)")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    arguments.run(arguments)


if __name__ == "__main__":
    main()
