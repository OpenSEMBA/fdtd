"""Compare converged backscatter RCS across all pairs of eight symmetry cases.

Each pair's RMS difference is calculated in dBsm over the shared frequency
grid; every pair must be below the 1.5 dB acceptance limit.
"""

import threading
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import pytest

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


def _rms_failures(pairwise_rms):
    return {
        pair: rms_db
        for pair, rms_db in pairwise_rms.items()
        if rms_db >= RMS_TOL_DB
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
def test_sphere_symmetry_rcs_agrees_across_cases(tmp_path):
    """Run all eight JSON cases and pass only if every pairwise RMS is below 1.5 dB."""
    results = {}
    for case in SPHERE_SYMMETRY_CASES:
        fn = CASES_FOLDER + SYMMETRY_FOLDER + case + ".fdtd.json"
        solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path / case)
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
    print(f"Overall result: {'FAIL' if failures else 'PASS'}")
    assert not failures, (
        f"Pairwise RMS must be below {RMS_TOL_DB:.2f} dB; "
        + ", ".join(
            f"{left} vs {right}: {rms_db:.3f} dB"
            for (left, right), rms_db in failures.items()
        )
    )
