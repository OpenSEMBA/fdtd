"""Plot full-domain versus PEC/PMC-symmetry sphere RCS from solver logs."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
RUNS = ROOT / "build-rls-mpi" / "sphere_rcs_comparison" / "mpi4_visible"


def find_log(folder):
    logs = sorted(folder.glob("sphere.fdtd_farfield_log_*.dat"))
    if len(logs) != 1:
        raise RuntimeError(f"Expected one far-field log in {folder}, found {len(logs)}")
    return logs[0]


def read_backscatter(path):
    data = np.loadtxt(path, skiprows=1)
    if data.ndim != 2 or data.shape[1] < 9:
        raise ValueError(f"Unexpected far-field columns in {path}")
    mask = np.isclose(data[:, 1], 90.0) & np.isclose(data[:, 2], 180.0)
    samples = data[mask]
    if samples.size == 0 or not np.all(np.isfinite(samples[:, [0, 8]])):
        raise ValueError(f"No finite backscatter RCS samples in {path}")
    return samples[:, 0], samples[:, 8]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-log", type=Path, default=None)
    parser.add_argument("--symmetry-log", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=RUNS)
    args = parser.parse_args()

    full_log = args.full_log or find_log(RUNS / "full")
    symmetry_log = args.symmetry_log or find_log(RUNS / "pec_pmc_symmetry")
    full_freq, full_rcs = read_backscatter(full_log)
    symmetry_freq, symmetry_rcs = read_backscatter(symmetry_log)
    if not np.array_equal(full_freq, symmetry_freq):
        raise ValueError("Full and symmetry logs contain different frequency samples")

    floor = np.finfo(float).tiny
    full_dbsm = 10.0 * np.log10(np.maximum(full_rcs, floor))
    symmetry_dbsm = 10.0 * np.log10(np.maximum(symmetry_rcs, floor))
    difference_db = symmetry_dbsm - full_dbsm

    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "sphere_rcs_comparison.csv"
    np.savetxt(
        csv_path,
        np.column_stack((full_freq, full_rcs, symmetry_rcs, difference_db)),
        delimiter=",",
        header="frequency_hz,full_rcs_geom_m2,symmetry_rcs_geom_m2,symmetry_minus_full_db",
        comments="",
    )

    fig, (rcs_ax, diff_ax) = plt.subplots(2, 1, sharex=True, constrained_layout=True)
    freq_mhz = full_freq / 1e6
    rcs_ax.semilogx(freq_mhz, full_dbsm, label="Full domain")
    rcs_ax.semilogx(freq_mhz, symmetry_dbsm, "--", label="PEC/PMC symmetry")
    rcs_ax.set_ylabel("Backscatter RCS [dBsm]")
    rcs_ax.grid(True, which="both", alpha=0.3)
    rcs_ax.legend()
    diff_ax.semilogx(freq_mhz, difference_db, color="tab:green")
    diff_ax.axhline(0.0, color="black", linewidth=0.8)
    diff_ax.set_xlabel("Frequency [MHz]")
    diff_ax.set_ylabel("Symmetry - full [dB]")
    diff_ax.grid(True, which="both", alpha=0.3)
    png_path = args.output_dir / "sphere_rcs_comparison.png"
    fig.savefig(png_path, dpi=180)
    plt.close(fig)

    worst = int(np.argmax(np.abs(difference_db)))
    print(f"Figure: {png_path}")
    print(f"Data:   {csv_path}")
    print(
        "Maximum absolute RCS difference: "
        f"{abs(difference_db[worst]):.3f} dB at {freq_mhz[worst]:.3f} MHz"
    )


if __name__ == "__main__":
    main()
