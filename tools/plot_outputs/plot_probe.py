#!/usr/bin/env python3
"""Plot a two-column time-series probe written by semba-fdtd."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input",
        type=Path,
        help="Probe descriptor (.json) or two-column text output (.dat)",
    )
    parser.add_argument("--output", type=Path, help="Image file to write instead of displaying")
    parser.add_argument("--title", help="Plot title")
    return parser.parse_args()


def resolve_probe(path):
    if path.suffix != ".json":
        return path, path.stem

    descriptor = json.loads(path.read_text(encoding="utf-8"))
    artifact = next(
        (
            artifact
            for artifact in descriptor["artifacts"]
            if artifact["kind"] == "text" and artifact["role"] == "canonical"
        ),
        None,
    )
    if artifact is None:
        raise ValueError(f"{path} has no canonical text artifact")
    return path.parent / artifact["relative_path"], descriptor.get("quantity", path.stem)


def plot(data_path, quantity, title, output):
    values = np.loadtxt(data_path, ndmin=2)
    if values.shape[1] != 2:
        raise ValueError(f"Expected two columns (time and value), found {values.shape[1]} in {data_path}")

    figure, axis = plt.subplots(layout="constrained")
    axis.plot(values[:, 0], values[:, 1])
    axis.set_xlabel("Time (s)")
    axis.set_ylabel(quantity)
    axis.set_title(title or data_path.stem)
    axis.grid()

    if output:
        output.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output, dpi=150)
    else:
        plt.show()
    plt.close(figure)


def main():
    arguments = parse_arguments()
    data_path, quantity = resolve_probe(arguments.input)
    if not data_path.is_file():
        raise FileNotFoundError(f"Probe data does not exist: {data_path}")

    plot(data_path, quantity, arguments.title, arguments.output)


if __name__ == "__main__":
    main()
