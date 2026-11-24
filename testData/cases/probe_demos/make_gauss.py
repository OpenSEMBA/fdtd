"""Write the Gaussian excitation files used by the probe demo cases.

The pulse is exp(-((t - T0) / S)^2), the same shape as testData/excitations/gauss.exc,
sampled at each case's time step so the solver does not warn about undersampling:
  gauss.exc       30 ps step, 1000 steps (point, movie, line, bulk_current, wire)
  gauss_wide.exc  7.7 ps step, 3900 steps (farfield_small)
Run from any directory: python3 make_gauss.py
"""
import math
import os

T0 = 1.2e-9
S = 2.8284271247e-10
# wider pulse, same peak: more spectral energy below 0.3 GHz
WIDE_T0 = 3.5e-9
WIDE_S = 1.0e-9

# (file name, time step, number of steps) matching general.timeStep / numberOfSteps in each case
EXCITATIONS = [
    ("gauss.exc", 3e-11, 1000),
    ("gauss_wide.exc", 7.7e-12, 3900),
]

here = os.path.dirname(os.path.abspath(__file__))
for name, dt, n_steps in EXCITATIONS:
    t0, s = (WIDE_T0, WIDE_S) if name == "gauss_wide.exc" else (T0, S)
    with open(os.path.join(here, name), "w") as f:
        for i in range(n_steps + 1):
            t = i * dt
            f.write(f"{t:.10e} {math.exp(-((t - t0) / s) ** 2):.10e}\n")
