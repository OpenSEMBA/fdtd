#!/usr/bin/env python3
"""Plot every canonical text-series probe descriptor below an output directory."""

import argparse
import json
from pathlib import Path

from plot_probe import plot, resolve_probe


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Solver output directory")
    parser.add_argument("output", type=Path, help="Directory for PNG plots")
    return parser.parse_args()


def main():
    arguments = parse_arguments()
    descriptors = sorted(arguments.input.rglob("*.json"))
    plotted = 0

    for descriptor in descriptors:
        try:
            data_path, quantity = resolve_probe(descriptor)
        except (json.JSONDecodeError, KeyError, ValueError):
            continue
        if not data_path.is_file():
            continue

        output_path = arguments.output / descriptor.relative_to(arguments.input).with_suffix(".png")
        plot(data_path, quantity, None, output_path)
        print(output_path)
        plotted += 1

    if not plotted:
        raise RuntimeError(f"No canonical text-series probe descriptors found in {arguments.input}")


if __name__ == "__main__":
    main()
