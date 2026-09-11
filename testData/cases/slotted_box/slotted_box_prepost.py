"""Compare DMMA FDTD and Robinson SEee references for the slotted PEC box.

This is the lossless, single-mode Robinson transmission-line model for a
centred aperture in a 300 mm x 300 mm x 120 mm PEC enclosure.  The script
runs the normal-incidence DMMA FDTD case and obtains the simulated electric
shielding effectiveness from the total and incident electric fields at the
cavity centre.  The FDTD model uses a normal-incidence plane wave and a
10 mm Yee-cell size, matching the Robinson reference assumption.

The formulation follows R. M. Robinson et al., "Analytical formulation for
the shielding effectiveness of enclosures with apertures", IEEE Transactions
on Electromagnetic Compatibility, 1998.
"""

# %% Imports and case paths
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.constants import c, epsilon_0, mu_0, pi


CASE_DIR = Path(__file__).resolve().parent
REPOSITORY_DIR = CASE_DIR.parents[2]
sys.path.insert(0, str(REPOSITORY_DIR / "src_pyWrapper"))

from pyWrapper import FDTD, Probe

def effective_slot_width(width, thickness):
    """Return Robinson's effective aperture width."""
    if thickness == 0.0:
        return width
    return width - 5.0 * thickness / (4.0 * pi) * (
        1.0 + np.log(4.0 * pi * width / thickness)
    )


def shielding_effectiveness(frequency):
    """Return Robinson electric and magnetic SE in dB at cavity centre."""
    frequency = np.asarray(frequency, dtype=float)
    k_0 = 2.0 * pi * frequency / c
    wavelength = c / frequency
    slot_width = effective_slot_width(SLOT_WIDTH, WALL_THICKNESS)

    fourth_root = np.power(1.0 - (slot_width / B) ** 2, 0.25)
    slot_impedance = 120.0 * pi**2 / np.log(
        2.0 * (1.0 + fourth_root) / (1.0 - fourth_root)
    )

    # Treat the below-cutoff TE10 propagation constant as complex.  The
    # ideal, lossless model has singular points at cavity resonances; those
    # are deliberately left unregularized and masked from the plot below.
    cutoff_factor = np.sqrt(1.0 - (wavelength / (2.0 * A)) ** 2 + 0.0j)
    guide_impedance = ETA_0 / cutoff_factor
    guide_wavenumber = k_0 * cutoff_factor

    aperture_impedance = 0.5j * (SLOT_LENGTH / A) * slot_impedance * np.tan(
        k_0 * SLOT_LENGTH / 2.0
    )
    v_1 = aperture_impedance / (ETA_0 + aperture_impedance)
    z_1 = ETA_0 * aperture_impedance / (ETA_0 + aperture_impedance)

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        propagation = np.cos(guide_wavenumber * PROBE_DEPTH) + 1j * (
            z_1 / guide_impedance
        ) * np.sin(guide_wavenumber * PROBE_DEPTH)
        v_2 = v_1 / propagation
        z_2 = (z_1 + 1j * guide_impedance * np.tan(
            guide_wavenumber * PROBE_DEPTH
        )) / (1.0 + 1j * (z_1 / guide_impedance) * np.tan(
            guide_wavenumber * PROBE_DEPTH
        ))
        z_3 = 1j * guide_impedance * np.tan(
            guide_wavenumber * (D - PROBE_DEPTH)
        )
        v_probe = v_2 * z_3 / (z_2 + z_3)
        i_probe = v_2 / (z_2 + z_3)
        se_electric = -20.0 * np.log10(np.abs(2.0 * v_probe))
        se_magnetic = -20.0 * np.log10(np.abs(2.0 * i_probe * ETA_0))

    return se_electric.real, se_magnetic.real



# %% Enclosure and aperture parameters
# Enclosure and aperture dimensions in metres.  The aperture is centred in
# the 300 mm x 120 mm front face and the observation point is at the cavity
# centre (p = d / 2).
A = 0.300
B = 0.120
D = 0.300
SLOT_LENGTH = 0.100
SLOT_WIDTH = 0.005
WALL_THICKNESS = 0.0
PROBE_DEPTH = D / 2.0
ETA_0 = np.sqrt(mu_0 / epsilon_0)


# %% Simulation
INPUT_FILE = CASE_DIR / "slotted_box_normal_incidence.fdtd.json"
SEMBA_EXE = REPOSITORY_DIR / "build" / "bin" / "semba-fdtd"

solver = FDTD(input_filename=INPUT_FILE, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

# %% Postproc
probes = {}
for filename in solver.getSolvedProbeFolders("Point probe"):
    probe = Probe(filename)
    if probe.field == "E":
        probes[probe.direction] = probe

if set(probes) != {"x", "y", "z"}:
    raise RuntimeError("Expected Ex, Ey, and Ez point-probe outputs.")
if any("incident" not in probe.data for probe in probes.values()):
    raise RuntimeError("Expected incident fields in the point-probe outputs.")

time = probes["x"]["time"].to_numpy()
for direction in ("y", "z"):
    if not np.array_equal(time, probes[direction]["time"].to_numpy()):
        raise RuntimeError("Point-probe components do not share a time grid.")

field = np.vstack([
    probes[direction]["field"].to_numpy() for direction in ("x", "y", "z")
])
incident = np.vstack([
    probes[direction]["incident"].to_numpy()
    for direction in ("x", "y", "z")
])

# %% Simulated-versus-analytical SEee plot
frequency = np.geomspace(300.0e6, 1.0e9, 101)
se_electric, _ = shielding_effectiveness(frequency)
valid_electric = np.isfinite(se_electric)


def dtft(signal, time, freqs):
    """Continuous-time DTFT approximation using trapezoidal integration."""
    signal = np.asarray(signal)
    time = np.asarray(time)
    res = np.empty_like(freqs, dtype=complex)
    trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
    for i, f in enumerate(freqs):
        res[i] = trapz(signal * np.exp(-1j * 2 * np.pi * f * time), time)
    return res


frequency_fdtd = frequency
total_spectrum = np.vstack([dtft(component, time, frequency_fdtd) for component in field])
incident_spectrum = np.vstack([dtft(component, time, frequency_fdtd) for component in incident])

with np.errstate(divide="ignore", invalid="ignore"):
    se_electric_fdtd = 20.0 * np.log10(
        np.linalg.norm(incident_spectrum, axis=0)
        / np.linalg.norm(total_spectrum, axis=0)
    )


# %%
plt.figure()
plt.semilogx(frequency[valid_electric], se_electric[valid_electric],
             label=r"$SE_{ee}$ Robinson reference (normal incidence)")
plt.semilogx(frequency_fdtd, se_electric_fdtd, ".-",
             label=r"$SE_{ee}$ DMMA FDTD (45 deg, p-polarized)")
plt.xlabel("Frequency (Hz)")
plt.ylabel("Shielding effectiveness (dB)")
plt.title("Electric shielding effectiveness of a slotted PEC enclosure")
plt.grid(which="both")
plt.legend()
plt.tight_layout()
plt.savefig(CASE_DIR /"slotted_box_seee_comparison.png", dpi=200)
plt.show()

# %%
