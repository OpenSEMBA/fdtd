"""Generates a long CPML run with a zero-mean source placed right next to the absorbing layer (grazing waves)."""
import json
import numpy as np

DX = 0.01
DT = 1.5e-11
W0 = 1.5e-10
N = 40


def write_excitation(fn="dgauss_cpml.exc"):
    t0 = 6 * W0
    t = np.arange(0, 2 * t0, DT)
    np.savetxt(fn, np.column_stack([t, -2 * (t - t0) / W0 * np.exp(-((t - t0) / W0) ** 2)]))


def case(name, steps):
    d = {
        "format": "FDTD Input file",
        "general": {"timeStep": DT, "numberOfSteps": steps},
        "boundary": {"all": {"type": "pml", "layers": 8}},
        "mesh": {
            "grid": {"numberOfCells": [N, N, N], "steps": {"x": [DX], "y": [DX], "z": [DX]}},
            "coordinates": [{"id": 1, "relativePosition": [10, 30, 20]}],
            "elements": [
                {"id": 1, "type": "cell", "intervals": [[[10, 20, 20], [10, 20, 21]]]},
                {"id": 2, "type": "node", "coordinateIds": [1]},
            ],
        },
        "materials": [],
        "sources": [{"type": "nodalSource", "magnitudeFile": "dgauss_cpml.exc",
                     "elementIds": [1], "hardness": "soft"}],
        "materialAssociations": [],
        "probes": [{"name": "ez", "type": "point", "field": "electric",
                    "elementIds": [2], "directions": ["z"], "domain": {"type": "time"}}],
    }
    json.dump(d, open(name, "w"), indent=2)


if __name__ == "__main__":
    write_excitation()
    case("cpml_grazing.fdtd.json", 20000)
