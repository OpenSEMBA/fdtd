"""RCS convergence tests for the sphere PEC/PMC symmetry cases.

Methodology (per case):
1. Launch the solver with a high numberOfSteps cap.
2. Monitor the "watchdog" point probe (a node outside the sphere) which
   records the time-domain electric field and the incident field.
3. When the scattered transient has decayed (trailing 5 ns peak below 1% of
   the global peak, never before 40 ns), place a void 'stopflushing' file in
   the run folder: the solver then ends gracefully and flushes all
   observation data, including the far-field transform.
4. Compare the backscatter RCS (theta=pi/2, phi=3*pi/2, the direction back
   towards the +y propagating incident plane wave) against the analytical
   Mie series of the 0.5 m radius PEC sphere.

All angles in these cases (planewave direction/polarization and farField
probe theta/phi) are expressed in radians, the unit the solver core expects.
"""

import threading
import time
from pathlib import Path

import numpy as np
import pytest
from scipy.special import spherical_jn, spherical_yn

from test.utils.utils import *

SYMMETRY_FOLDER = "sphere/symmetry/"

SPHERE_SYMMETRY_CASES = [
    "sphere_full_polx",
    "sphere_xpec_polx",
    "sphere_zpmc_polx",
    "sphere_xpec_zpmc_polx",
    "sphere_full_polz",
    "sphere_xpmc_polz",
    "sphere_zpec_polz",
    "sphere_xpmc_zpec_polz",
]

SPHERE_RADIUS_M = 0.5
C0 = 299792458.0

# Convergence criterion for the watchdog scattered field.
T_MIN_S = 40e-9        # never stop before the pulse has fully passed
T_CAP_S = 150e-9       # stop regardless past this simulated time
DECAY_WINDOW_S = 5e-9  # trailing window for the decay check
DECAY_RATIO = 1e-2     # trailing peak below 1% of the global peak
POLL_S = 10.0          # watchdog polling period
WALL_CAP_S = 30 * 60   # per-case wall-clock safety cap

# Acceptance tolerances vs Mie (measured: rms <= 0.89 dB, max <= 3.07 dB).
RMS_TOL_DB = 1.5
MAX_TOL_DB = 4.0


def _mie_backscatter_rcs(frequency_hz, radius_m=SPHERE_RADIUS_M):
    """Monostatic RCS [m^2] of a PEC sphere via the Mie series."""
    frequency_hz = np.atleast_1d(np.asarray(frequency_hz, dtype=float))
    result = np.empty_like(frequency_hz)
    for i, f in enumerate(frequency_hz):
        x = 2.0 * np.pi * f * radius_m / C0
        n = np.arange(1, int(np.ceil(x + 4.05 * x ** (1.0 / 3.0) + 10.0)) + 1)
        jn = spherical_jn(n, x)
        yn = spherical_yn(n, x)
        psi = x * jn
        xi = x * (jn + 1j * yn)
        psi_d = jn + x * spherical_jn(n, x, derivative=True)
        xi_d = psi_d + 1j * (yn + x * spherical_yn(n, x, derivative=True))
        s = np.sum((2 * n + 1) * (-1.0) ** n * (psi / xi - psi_d / xi_d))
        result[i] = np.pi / (x / radius_m) ** 2 * abs(s) ** 2
    return result


def _read_watchdog(folder):
    """Return (t, max |E_scattered|) from the watchdog point probe files."""
    series = []
    for f in sorted(Path(folder).glob("*_watchdog_E*_tm.dat")):
        try:
            data = np.loadtxt(f, skiprows=1)
        except Exception:
            return None
        if data.ndim != 2 or data.shape[0] < 10:
            return None
        series.append(data)  # columns: t, field, incident
    if not series:
        return None
    t = series[0][:, 0]
    scat = np.zeros_like(t)
    for s in series:
        scat = np.maximum(scat, np.abs(s[:, 1] - s[:, 2]))
    return t, scat


def _converged(t, scat):
    if t[-1] < T_MIN_S:
        return False
    peak = scat.max()
    if peak <= 0.0:
        return False
    return scat[t >= t[-1] - DECAY_WINDOW_S].max() < DECAY_RATIO * peak


def _read_backscatter_rcs(folder):
    logs = sorted(Path(folder).glob("*_farfield_log_*.dat"))
    assert len(logs) == 1, f"Expected one far-field log in {folder}, found {len(logs)}"
    data = np.loadtxt(logs[0], skiprows=1)
    if data.ndim == 1:
        data = data[None, :]
    mask = np.isclose(data[:, 1], np.pi / 2, atol=1e-4) & np.isclose(
        data[:, 2], 3 * np.pi / 2, atol=1e-4
    )
    sel = data[mask]
    assert len(sel) > 0, "No backscatter direction found in far-field log"
    order = np.argsort(sel[:, 0])
    return sel[order, 0], sel[order, 8]  # freq, RCS(GEOM)


@pytest.mark.farfield
@pytest.mark.planewave
@pytest.mark.probes
@pytest.mark.parametrize("case", SPHERE_SYMMETRY_CASES)
def test_sphere_symmetry_rcs_converges_to_mie(tmp_path, case):
    """Run to convergence via stopflushing; backscatter RCS must match Mie."""
    fn = CASES_FOLDER + SYMMETRY_FOLDER + case + ".fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 20000  # hard cap; stopflushing ends earlier

    run_error = []

    def _run():
        try:
            solver.run()
        except AssertionError as exc:
            run_error.append(exc)

    thread = threading.Thread(target=_run, daemon=True)
    thread.start()

    folder = solver.getFolder()
    stop_file = Path(folder) / "stopflushing"
    t0 = time.time()
    stopped = False
    while thread.is_alive():
        time.sleep(POLL_S)
        wd = _read_watchdog(folder)
        if wd is not None and not stopped:
            t, scat = wd
            if _converged(t, scat) or t[-1] >= T_CAP_S:
                stop_file.touch()
                stopped = True
        if time.time() - t0 > WALL_CAP_S and not stopped:
            stop_file.touch()
            stopped = True
    thread.join(timeout=600)

    assert not run_error, f"Solver run failed: {run_error[0]}"
    assert solver.hasFinishedSuccessfully()
    assert stopped, "Convergence criterion was never met; stopflushing never placed"

    freq, rcs = _read_backscatter_rcs(folder)
    mie = _mie_backscatter_rcs(freq)
    err_db = 10 * np.log10(np.maximum(rcs, 1e-30)) - 10 * np.log10(mie)
    rms_db = float(np.sqrt(np.mean(err_db**2)))
    max_db = float(np.max(np.abs(err_db)))
    print(f"\n{case}: RCS vs Mie  rms={rms_db:.3f} dB  max|err|={max_db:.3f} dB")
    assert np.all(np.isfinite(rcs)) and np.all(rcs > 0.0)
    assert rms_db < RMS_TOL_DB, f"RMS error vs Mie {rms_db:.2f} dB exceeds {RMS_TOL_DB} dB"
    assert max_db < MAX_TOL_DB, f"Max error vs Mie {max_db:.2f} dB exceeds {MAX_TOL_DB} dB"
