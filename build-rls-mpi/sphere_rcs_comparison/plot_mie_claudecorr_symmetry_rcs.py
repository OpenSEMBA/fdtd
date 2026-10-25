"""Compare Mie, full-domain and PEC/PMC quarter-domain backscatter RCS."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
MIE = HERE / "full_sphere_dx10mm_pml15further_1rank_20threads_precise" / "sphere_rcs_mie_radius_1m.csv"
FULL = HERE / "sphere_claudecorr_fast_mpi1_omp20_6000_20261025_1342"
SYMMETRY = HERE / "sphere_claudecorr_symmetries_pecY_pmcZ_mpi1_omp20_6000_20261025_1349"
OUTPUT = HERE / "sphere_rcs_mie_claudecorr_symmetry_comparison.png"


def load_farfield(folder):
    files = sorted(folder.glob("*_farfield_log_*.dat"))
    if len(files) != 1:
        raise RuntimeError(f"Expected one far-field file in {folder}, found {len(files)}")
    data = np.loadtxt(files[0], skiprows=1)
    order = np.argsort(data[:, 0])
    return data[order, 0], data[order, 8]


def to_dbsm(rcs):
    return 10.0 * np.log10(np.maximum(rcs, np.finfo(float).tiny))


def main():
    mie = np.loadtxt(MIE, delimiter=",", skiprows=1)
    mie_frequency, mie_rcs = mie[:, 0], mie[:, 1]
    full_frequency, full_rcs = load_farfield(FULL)
    symmetry_frequency, symmetry_rcs = load_farfield(SYMMETRY)

    mie_full = np.interp(full_frequency, mie_frequency, mie_rcs)
    mie_symmetry = np.interp(symmetry_frequency, mie_frequency, mie_rcs)
    full_error = to_dbsm(full_rcs) - to_dbsm(mie_full)
    symmetry_error = to_dbsm(symmetry_rcs) - to_dbsm(mie_symmetry)

    figure, (rcs_axis, error_axis) = plt.subplots(
        2, 1, figsize=(10, 8), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [2, 1]},
    )
    rcs_axis.semilogx(mie_frequency / 1e6, to_dbsm(mie_rcs), "k-", linewidth=2, label="Mie (analytical)")
    rcs_axis.semilogx(full_frequency / 1e6, to_dbsm(full_rcs), "-", color="tab:blue", label="FDTD full domain (claudecorr)")
    rcs_axis.semilogx(symmetry_frequency / 1e6, to_dbsm(symmetry_rcs), "--", color="tab:red",
                      label="FDTD quarter domain, PEC y / PMC z symmetry")
    rcs_axis.set_ylabel("Backscatter RCS [dBsm]")
    rcs_axis.set_title("PEC sphere, radius 1 m, dx = 10 mm: backscatter (θ = 90°, φ = 180°)")
    rcs_axis.grid(True, which="both", alpha=0.3)
    rcs_axis.legend()

    error_axis.semilogx(full_frequency / 1e6, full_error, color="tab:blue", label="Full domain - Mie")
    error_axis.semilogx(symmetry_frequency / 1e6, symmetry_error, "--", color="tab:red", label="Symmetry - Mie")
    error_axis.axhline(0.0, color="black", linewidth=0.8)
    error_axis.set_xlabel("Frequency [MHz]")
    error_axis.set_ylabel("Error vs Mie [dB]")
    error_axis.grid(True, which="both", alpha=0.3)
    error_axis.legend()

    figure.savefig(OUTPUT, dpi=180)
    plt.close(figure)

    print(f"Figure: {OUTPUT}")
    print(f"Full domain  max |error| vs Mie: {np.max(np.abs(full_error)):.3f} dB, "
          f"RMS: {np.sqrt(np.mean(full_error**2)):.3f} dB")
    print(f"Symmetry     max |error| vs Mie: {np.max(np.abs(symmetry_error)):.3f} dB, "
          f"RMS: {np.sqrt(np.mean(symmetry_error**2)):.3f} dB")


if __name__ == "__main__":
    main()
