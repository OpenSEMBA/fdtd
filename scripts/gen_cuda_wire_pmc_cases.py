#!/usr/bin/env python3
"""Generate scaled Holland-wire / PMC / periodic nodal CUDA golden JSON cases."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "testData/cases/cuda"


def scaled(n: int, steps: int, boundary: str, wire: bool) -> dict:
    cx = cy = cz = n // 2
    z0, z1 = cz - 4, cz + 4
    px, py, pz = cx + 2, cy + 2, cz
    coords = [
        {"id": 1, "relativePosition": [px, py, pz]},
        {"id": 2, "relativePosition": [cx, cy, z0]},
        {"id": 3, "relativePosition": [cx, cy, z1]},
    ]
    if wire:
        coords.extend(
            [
                {"id": 4, "relativePosition": [cx, cy, cz]},
                {"id": 5, "relativePosition": [cx, cy, cz + 2]},
            ]
        )
    elements = [{"id": 1, "type": "node", "coordinateIds": [1]}]
    if wire:
        elements.extend(
            [
                {"id": 2, "type": "polyline", "coordinateIds": [2, 3]},
                {"id": 3, "type": "node", "coordinateIds": [2]},
                {"id": 4, "type": "node", "coordinateIds": [3]},
                {
                    "id": 5,
                    "type": "cell",
                    "intervals": [[[cx, cy, cz], [cx, cy, cz + 2]]],
                },
            ]
        )
        src_elements = [5]
    else:
        elements.append(
            {
                "id": 2,
                "type": "cell",
                "intervals": [[[cx, cy, z0], [cx, cy, z1]]],
            }
        )
        src_elements = [2]

    doc = {
        "format": "FDTD Input file",
        "__comments": f"{n}^3 CUDA golden ({boundary}" + (", Holland wire" if wire else "") + ").",
        "general": {"timeStep": 1e-11, "numberOfSteps": steps},
        "boundary": {"all": {"type": boundary}},
        "mesh": {
            "grid": {
                "numberOfCells": [n, n, n],
                "steps": {"x": [0.1], "y": [0.1], "z": [0.1]},
            },
            "coordinates": coords,
            "elements": elements,
        },
        "materials": [],
        "materialAssociations": [],
        "sources": [
            {
                "name": "nodal_soft",
                "type": "nodalSource",
                "hardness": "soft",
                "field": "current",
                "magnitudeFile": "gauss.exc",
                "elementIds": src_elements,
            }
        ],
        "probes": [
            {
                "name": "electric_field_point_probe",
                "type": "point",
                "elementIds": [1],
                "directions": ["x", "y", "z"],
                "domain": {"type": "time"},
            }
        ],
    }
    if wire:
        doc["materials"] = [
            {
                "name": "wireMaterial",
                "id": 1,
                "type": "wire",
                "radius": 0.001,
                "resistancePerMeter": 0.0229,
            },
            {
                "name": "openEnd",
                "id": 2,
                "type": "terminal",
                "terminations": [{"type": "open"}],
            },
        ]
        doc["materialAssociations"] = [
            {
                "name": "wirez",
                "elementIds": [2],
                "materialId": 1,
                "initialTerminalId": 2,
                "endTerminalId": 2,
            }
        ]
    return doc


def main() -> None:
    specs = [
        ("box_wire_holland_mur", "mur", True, {40: 200, 100: 300, 200: 300}),
        ("box_pmc_nodal", "pmc", False, {40: 200, 100: 300, 200: 300}),
        ("box_periodic_nodal", "periodic", False, {40: 200, 100: 300, 200: 300}),
    ]
    for prefix, bnd, wire, step_map in specs:
        for n, steps in step_map.items():
            path = OUT / f"{prefix}_{n}.fdtd.json"
            path.write_text(json.dumps(scaled(n, steps, bnd, wire), indent=2) + "\n")
            print("wrote", path.name)


if __name__ == "__main__":
    main()
