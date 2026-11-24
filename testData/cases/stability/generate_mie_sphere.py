"""Staircase PEC sphere (r = 0.5 m, 20 cells of 2.5 cm) in a plane wave: bistatic RCS to be compared with the Mie series.

Reuses the plane wave box, far-field probe and boundary of the conformal sphere RCS case.
"""
import json
import numpy as np

BASE = "../conformal/conformal_sphere_rcs.fdtd.json"
REFINE = 2  # cells per base cell (base: 60^3 cells of 5 cm, r = 10 cells)
N, CENTER, RADIUS = 60 * REFINE, 30 * REFINE, 10 * REFINE  # in cells


def sphere_intervals():
    out = []
    for j in range(N):
        for k in range(N):
            dy, dz = j + 0.5 - CENTER, k + 0.5 - CENTER
            xs = [i for i in range(N) if (i + 0.5 - CENTER) ** 2 + dy**2 + dz**2 <= RADIUS**2]
            if xs:
                out.append([[xs[0], j, k], [xs[-1] + 1, j + 1, k + 1]])
    return out


if __name__ == "__main__":
    d = json.load(open(BASE))
    d["general"] = {"timeStep": 7e-11 / REFINE, "numberOfSteps": 3000 * REFINE}
    d["mesh"]["grid"]["numberOfCells"] = [N, N, N]
    d["mesh"]["grid"]["steps"] = {k: [0.05 / REFINE] for k in "xyz"}
    d["mesh"]["elements"] = [
        {"id": 1, "type": "cell", "intervals": sphere_intervals()},
        {"id": 2, "type": "cell", "intervals": [[[8 * REFINE] * 3, [52 * REFINE] * 3]]},
        {"id": 4, "type": "cell", "intervals": [[[4 * REFINE] * 3, [56 * REFINE] * 3]]},
    ]
    d["probes"] = [p for p in d["probes"] if p["name"] == "n2f"]
    d["sources"][0]["magnitudeFile"] = "gauss_mie.exc"
    json.dump(d, open("mie_sphere.fdtd.json", "w"), indent=1)
