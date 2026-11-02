"""Generates the CPML reflection cases: same source/probe in a small and a large PML-terminated domain."""
import json
import numpy as np

DX = 0.01
DT = 1.5e-11
W0 = 1.5e-10
STEPS = 1200


def write_excitation(fn="gauss_cpml.exc"):
    t0 = 6 * W0
    t = np.arange(0, 2 * t0, DT)
    np.savetxt(fn, np.column_stack([t, np.exp(-((t - t0) / W0) ** 2)]))


def case(name, n):
    m = n // 2
    d = {
        "format": "FDTD Input file",
        "general": {"timeStep": DT, "numberOfSteps": STEPS},
        "boundary": {"all": {"type": "pml", "layers": 8}},
        "mesh": {
            "grid": {"numberOfCells": [n, n, n], "steps": {"x": [DX], "y": [DX], "z": [DX]}},
            "coordinates": [{"id": 1, "relativePosition": [m + 4, m, m]}],
            "elements": [
                {"id": 1, "type": "cell", "intervals": [[[m, m, m], [m, m, m + 1]]]},
                {"id": 2, "type": "node", "coordinateIds": [1]},
            ],
        },
        "materials": [],
        "sources": [{"type": "nodalSource", "magnitudeFile": "gauss_cpml.exc",
                     "elementIds": [1], "hardness": "soft"}],
        "materialAssociations": [],
        "probes": [{"name": "ez", "type": "point", "field": "electric",
                    "elementIds": [2], "directions": ["z"], "domain": {"type": "time"}}],
    }
    json.dump(d, open(name, "w"), indent=2)


if __name__ == "__main__":
    write_excitation()
    case("cpml_small.fdtd.json", 30)
    case("cpml_large.fdtd.json", 90)
