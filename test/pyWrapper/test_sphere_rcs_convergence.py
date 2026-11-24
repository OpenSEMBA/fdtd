"""Compare converged backscatter RCS across eight sphere symmetry cases.

All pairwise Semba RMS differences must pass before each case is checked
against the analytical Mie solution up to 300 MHz.
"""

import threading
import time
from itertools import combinations
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

# RMS limit for pairwise backscatter-RCS differences, measured in dB.
RMS_TOL_DB = 1.5
MIE_RMS_TOL_DB = 1.5
MIE_MAX_FREQUENCY_HZ = 300e6
SPHERE_RADIUS_M = 0.5
C0 = 299792458.0

T_MIN_S = 40e-9
T_CAP_S = 150e-9
DECAY_WINDOW_S = 5e-9
DECAY_RATIO = 1e-2
POLL_S = 10.0
WALL_CAP_S = 30 * 60


def _read_watchdog(folder):
    """Return (t, max |E_scattered|) from the watchdog point probe files."""
    series = []
    for path in sorted(Path(folder).glob("*_watchdog_E*_tm.dat")):
        try:
            data = np.loadtxt(path, skiprows=1)
        except Exception:
            return None
        if data.ndim != 2 or data.shape[0] < 10:
            return None
        series.append(data)
    if not series:
        return None
    t = series[0][:, 0]
    scattered = np.zeros_like(t)
    for data in series:
        scattered = np.maximum(scattered, np.abs(data[:, 1] - data[:, 2]))
    return t, scattered


def _converged(t, scattered):
    if t[-1] < T_MIN_S:
        return False
    peak = scattered.max()
    if peak <= 0.0:
        return False
    return scattered[t >= t[-1] - DECAY_WINDOW_S].max() < DECAY_RATIO * peak


def _read_backscatter_rcs(folder):
    logs = sorted(Path(folder).glob("*_farfield_log_*.dat"))
    assert len(logs) == 1, f"Expected one far-field log in {folder}, found {len(logs)}"
    data = np.loadtxt(logs[0], skiprows=1)
    if data.ndim == 1:
        data = data[None, :]
    mask = np.isclose(data[:, 1], np.pi / 2, atol=1e-4) & np.isclose(
        data[:, 2], 3 * np.pi / 2, atol=1e-4
    )
    selected = data[mask]
    assert len(selected) > 0, "No backscatter direction found in far-field log"
    selected = selected[np.argsort(selected[:, 0])]
    return selected[:, 0], selected[:, 8]


def _read_backscatter_rcs_csvs(folder):
    """Load the eight precomputed case CSVs without launching SEMBA."""
    folder = Path(folder)
    if not folder.is_dir():
        raise NotADirectoryError(f"Sphere RCS CSV directory does not exist: {folder}")

    results = {}
    for case in SPHERE_SYMMETRY_CASES:
        csv_name = case[len("sphere_") :] + "_rcs.csv"
        path = folder / csv_name
        if not path.is_file():
            raise FileNotFoundError(f"Missing RCS results for {case}: {path}")
        data = np.loadtxt(path, delimiter=",", skiprows=1)
        if data.ndim == 1:
            data = data[None, :]
        if data.ndim != 2 or data.shape[1] < 2:
            raise ValueError(f"Unexpected RCS CSV format in {path}")
        results[case] = (data[:, 0], data[:, 1])
    return results


def _mie_backscatter_rcs(frequency_hz, radius_m=SPHERE_RADIUS_M):
    """Return analytical PEC-sphere monostatic RCS [m^2] via the Mie series."""
    frequency_hz = np.atleast_1d(np.asarray(frequency_hz, dtype=float))
    result = np.empty_like(frequency_hz)
    for index, frequency in enumerate(frequency_hz):
        x = 2.0 * np.pi * frequency * radius_m / C0
        n = np.arange(1, int(np.ceil(x + 4.05 * x ** (1.0 / 3.0) + 10.0)) + 1)
        jn = spherical_jn(n, x)
        yn = spherical_yn(n, x)
        psi = x * jn
        xi = x * (jn + 1j * yn)
        psi_d = jn + x * spherical_jn(n, x, derivative=True)
        xi_d = psi_d + 1j * (yn + x * spherical_yn(n, x, derivative=True))
        series = np.sum((2 * n + 1) * (-1.0) ** n * (psi / xi - psi_d / xi_d))
        result[index] = np.pi / (x / radius_m) ** 2 * abs(series) ** 2
    return result


def _rms_db_vs_mie(freq, rcs, mie_rcs, max_frequency_hz):
    if freq.shape != rcs.shape or freq.shape != mie_rcs.shape:
        raise ValueError("Frequency and RCS arrays must have matching shapes")
    in_band = freq <= max_frequency_hz
    if not np.any(in_band):
        raise ValueError(f"No frequency samples at or below {max_frequency_hz:g} Hz")
    difference_db = 10.0 * (np.log10(rcs[in_band]) - np.log10(mie_rcs[in_band]))
    return float(np.sqrt(np.mean(difference_db**2)))


def _pairwise_rms_db(results):
    """Calculate RMS dBsm differences for every pair of case RCS curves."""
    pairwise_rms = {}
    first_case = next(iter(results))
    reference_freq, _ = results[first_case]
    for case, (freq, rcs) in results.items():
        if not np.array_equal(reference_freq, freq):
            raise ValueError(f"{case} has a different frequency grid from {first_case}")
        if freq.shape != rcs.shape:
            raise ValueError(f"{case} has different numbers of frequencies and RCS samples")
        if not np.all(np.isfinite(rcs)) or np.any(rcs <= 0.0):
            raise ValueError(f"{case} contains non-finite or non-positive RCS values")

    for left, right in combinations(results, 2):
        left_rcs = results[left][1]
        right_rcs = results[right][1]
        difference_db = 10.0 * (np.log10(left_rcs) - np.log10(right_rcs))
        pairwise_rms[(left, right)] = float(np.sqrt(np.mean(difference_db**2)))
    return pairwise_rms


def _rms_failures(rms_by_comparison, tolerance_db=RMS_TOL_DB):
    return {
        comparison: rms_db
        for comparison, rms_db in rms_by_comparison.items()
        if rms_db >= tolerance_db
    }


def test_mie_rms_uses_only_frequencies_through_300_mhz():
    frequencies = np.array([100e6, 300e6, 300.1e6])
    mie_rcs = np.ones(3)
    fdtd_rcs = mie_rcs * 10 ** (np.array([0.5, -0.5, 20.0]) / 10)

    rms_db = _rms_db_vs_mie(frequencies, fdtd_rcs, mie_rcs, MIE_MAX_FREQUENCY_HZ)
    over_limit_rcs = mie_rcs * 10 ** (np.array([1.6, 1.6, 20.0]) / 10)
    over_limit_rms_db = _rms_db_vs_mie(
        frequencies, over_limit_rcs, mie_rcs, MIE_MAX_FREQUENCY_HZ
    )

    assert rms_db == pytest.approx(0.5)
    assert _rms_failures({"case": rms_db}, MIE_RMS_TOL_DB) == {}
    assert _rms_failures({"case": over_limit_rms_db}, MIE_RMS_TOL_DB) == {
        "case": pytest.approx(1.6)
    }


def test_pairwise_rms_db_and_pass_fail_threshold():
    frequencies = np.array([1.0, 2.0])
    reference = np.array([1.0, 2.0])
    results = {
        "reference": (frequencies, reference),
        "within_limit": (frequencies, reference * 10 ** (1.4 / 10)),
        "over_limit": (frequencies, reference * 10 ** (1.6 / 10)),
    }

    pairwise_rms = _pairwise_rms_db(results)

    assert pairwise_rms[("reference", "within_limit")] == pytest.approx(1.4)
    assert pairwise_rms[("reference", "over_limit")] == pytest.approx(1.6)
    assert _rms_failures(pairwise_rms) == {("reference", "over_limit"): pytest.approx(1.6)}


@pytest.mark.farfield
@pytest.mark.planewave
@pytest.mark.probes
def test_sphere_symmetry_rcs_agrees_across_cases(tmp_path, request):
    """Compare eight cases from SEMBA runs or existing CSVs."""
    csv_dir = request.config.getoption("--sphere-rcs-csv-dir")
    if csv_dir is not None:
        results = _read_backscatter_rcs_csvs(csv_dir)
        print(f"Using precomputed RCS CSVs from {csv_dir}; SEMBA was not launched.")
    else:
        results = {}
        for case in SPHERE_SYMMETRY_CASES:
            fn = CASES_FOLDER + SYMMETRY_FOLDER + case + ".fdtd.json"
            run_folder = tmp_path / case
            run_folder.mkdir()
            solver = FDTD(
                input_filename=fn,
                path_to_exe=SEMBA_EXE,
                run_in_folder=run_folder,
            )
            solver["general"]["numberOfSteps"] = 20000

            run_error = []

            def _run():
                try:
                    solver.run()
                except Exception as exc:
                    run_error.append(exc)

            thread = threading.Thread(target=_run, daemon=True)
            thread.start()

            folder = solver.getFolder()
            stop_file = Path(folder) / "stopflushing"
            t0 = time.time()
            stop_reason = None
            while thread.is_alive():
                time.sleep(POLL_S)
                watchdog = _read_watchdog(folder)
                if watchdog is not None:
                    t, scattered = watchdog
                    if _converged(t, scattered):
                        stop_reason = "convergence"
                    elif t[-1] >= T_CAP_S:
                        stop_reason = "simulation-time cap"
                    if stop_reason is not None:
                        stop_file.touch()
                        break
                if time.time() - t0 > WALL_CAP_S:
                    stop_file.touch()
                    stop_reason = "wall-clock cap"
                    break
            thread.join(timeout=600)

            assert not thread.is_alive(), f"{case}: solver did not stop after stopflushing"
            assert not run_error, f"{case}: solver run failed: {run_error[0]}"
            assert solver.hasFinishedSuccessfully(), f"{case}: solver did not finish successfully"
            assert stop_reason in {"convergence", "simulation-time cap"}, (
                f"{case}: watchdog did not converge before the wall-clock cap"
            )

            freq, rcs = _read_backscatter_rcs(folder)
            assert np.all(np.isfinite(rcs)) and np.all(rcs > 0.0), (
                f"{case}: RCS contains non-finite or non-positive values"
            )
            results[case] = (freq, rcs)
            print(f"\n{case}: completed ({stop_reason})")

    pairwise_rms = _pairwise_rms_db(results)
    failures = _rms_failures(pairwise_rms)
    print(f"\nPairwise FDTD backscatter-RCS RMS (limit: {RMS_TOL_DB:.2f} dB):")
    for (left, right), rms_db in pairwise_rms.items():
        status = "FAIL" if rms_db >= RMS_TOL_DB else "PASS"
        print(f"{status}: {left} vs {right}: {rms_db:.3f} dB")
    print(f"Pairwise Semba result: {'FAIL' if failures else 'PASS'}")
    assert not failures, (
        f"Pairwise RMS must be below {RMS_TOL_DB:.2f} dB; "
        + ", ".join(
            f"{left} vs {right}: {rms_db:.3f} dB"
            for (left, right), rms_db in failures.items()
        )
    )

    print(
        f"\nAll pairwise Semba comparisons passed; checking Mie through "
        f"{MIE_MAX_FREQUENCY_HZ / 1e6:.0f} MHz:"
    )
    mie_rms = {}
    for case, (freq, rcs) in results.items():
        in_band = freq <= MIE_MAX_FREQUENCY_HZ
        mie_freq = freq[in_band]
        mie_rcs = _mie_backscatter_rcs(mie_freq)
        mie_rms[case] = _rms_db_vs_mie(
            mie_freq, rcs[in_band], mie_rcs, MIE_MAX_FREQUENCY_HZ
        )
    mie_failures = _rms_failures(mie_rms, MIE_RMS_TOL_DB)
    for case, rms_db in mie_rms.items():
        status = "FAIL" if rms_db >= MIE_RMS_TOL_DB else "PASS"
        print(f"{status}: {case} vs Mie: {rms_db:.3f} dB")
    print(f"Overall result: {'FAIL' if mie_failures else 'PASS'}")
    assert not mie_failures, (
        f"Mie RMS through {MIE_MAX_FREQUENCY_HZ / 1e6:.0f} MHz must be below "
        f"{MIE_RMS_TOL_DB:.2f} dB; "
        + ", ".join(f"{case}: {rms_db:.3f} dB" for case, rms_db in mie_failures.items())
    )
