"""Generates a 1D-like cavity (PEC x walls) half filled with er=4 to check resonances across a dielectric interface."""
import json
import numpy as np

DX = 0.01
DT = 1.5e-11
W0 = 1.5e-10
N = [40, 2, 2]  # 0.40 m long, interface at x = 0.20 m


def write_excitation(fn="gauss_dielectric_cavity.exc"):
    t0 = 6 * W0
    t = np.arange(0, 2 * t0, DT)
    np.savetxt(fn, np.column_stack([t, np.exp(-((t - t0) / W0) ** 2)]))


def case(name, steps):
    d = {
        "format": "FDTD Input file",
        "general": {"timeStep": DT, "numberOfSteps": steps},
        "boundary": {"xLower": {"type": "pec"}, "xUpper": {"type": "pec"},
                     "yLower": {"type": "pmc"}, "yUpper": {"type": "pmc"},
                     "zLower": {"type": "pec"}, "zUpper": {"type": "pec"}},
        "mesh": {
            "grid": {"numberOfCells": N, "steps": {"x": [DX], "y": [DX], "z": [DX]}},
            "coordinates": [{"id": 1, "relativePosition": [14, 1, 1]}],
            "elements": [
                {"id": 1, "type": "cell", "intervals": [[[5, 1, 0], [5, 1, 2]]]},
                {"id": 2, "type": "node", "coordinateIds": [1]},
                {"id": 3, "type": "cell", "intervals": [[[20, 0, 0], [40, 2, 2]]]},
            ],
        },
        "materials": [{"id": 1, "name": "dielectric", "type": "isotropic", "relativePermittivity": 4.0}],
        "sources": [{"type": "nodalSource", "magnitudeFile": "gauss_dielectric_cavity.exc",
                     "elementIds": [1], "hardness": "soft"}],
        "materialAssociations": [{"materialId": 1, "elementIds": [3]}],
        "probes": [{"name": "ez", "type": "point", "field": "electric",
                    "elementIds": [2], "directions": ["z"], "domain": {"type": "time"}}],
    }
    json.dump(d, open(name, "w"), indent=2)


if __name__ == "__main__":
    write_excitation()
    case("dielectric_cavity.fdtd.json", 60000)
