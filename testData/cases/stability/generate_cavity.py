"""Generates the PEC cavity stability cases (lossless and lossy)."""
import json
import numpy as np

DX = 0.01
N = [20, 16, 12]          # 0.20 x 0.16 x 0.12 m
DT = 1.5e-11              # below the solver CFL correction (no auto dt change)
W0 = 1.5e-10


def write_excitation(fn="gauss_cavity.exc"):
    t0 = 6 * W0
    t = np.arange(0, 2 * t0, DT)
    np.savetxt(fn, np.column_stack([t, np.exp(-((t - t0) / W0) ** 2)]))


def case(name, steps, sigma=None):
    d = {
        "format": "FDTD Input file",
        "general": {"timeStep": DT, "numberOfSteps": steps},
        "boundary": {"all": {"type": "pec"}},
        "mesh": {
            "grid": {"numberOfCells": N, "steps": {"x": [DX], "y": [DX], "z": [DX]}},
            "coordinates": [
                {"id": 1, "relativePosition": [4, 5, 3]},
                {"id": 2, "relativePosition": [4, 5, 5]},
                {"id": 3, "relativePosition": [13, 9, 7]},
            ],
            "elements": [
                {"id": 1, "type": "cell", "intervals": [[[4, 5, 3], [4, 5, 5]]]},
                {"id": 2, "type": "node", "coordinateIds": [3]},
                {"id": 3, "type": "cell", "intervals": [[[0, 0, 0], [20, 16, 12]]]},
            ],
        },
        "materials": [],
        "sources": [{"type": "nodalSource", "magnitudeFile": "gauss_cavity.exc",
                     "elementIds": [1], "hardness": "soft"}],
        "materialAssociations": [],
        "probes": [{"name": "ez", "type": "point", "field": "electric",
                    "elementIds": [2], "directions": ["z"], "domain": {"type": "time"}}],
    }
    if sigma:
        d["materials"].append({"id": 1, "name": "lossy", "type": "isotropic",
                               "relativePermittivity": 1.0, "electricConductivity": sigma})
        d["materialAssociations"].append({"materialId": 1, "elementIds": [3]})
    json.dump(d, open(name, "w"), indent=2)


if __name__ == "__main__":
    write_excitation()
    case("pec_cavity.fdtd.json", 100000)
    case("lossy_cavity.fdtd.json", 8000, sigma=1e-3)
    case("courant_cavity.fdtd.json", 20000)
