"""Compare full-domain and PEC/PMC-symmetry RCS for a PEC sphere."""

import argparse
import copy
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
CASE = ROOT / "testData" / "cases" / "sphere" / "sphere.fdtd.json"
DEFAULT_EXE = ROOT / "build-rls" / "bin" / "semba-fdtd.exe"
OUTPUT = ROOT / "build-rls" / "sphere_rcs_comparison"
SCALE = 2
SYMMETRY_INDEX = 80
FREQUENCY_RANGE = (10e6, 900e6)
NUMBER_OF_FREQUENCIES = 80
NUMBER_OF_STEPS = 3000

sys.path.insert(0, str(ROOT / "src_pyWrapper"))
from pyWrapper import FDTD, Probe


def configure_case(solver, symmetric):
    grid = solver["mesh"]["grid"]
    grid["steps"] = {"x": [0.025], "y": [0.025], "z": [0.025]}

    if symmetric:
        grid["numberOfCells"] = [160, 80, 80]
        grid["origin"] = [0.0, 2.0, 2.0]
    else:
        grid["numberOfCells"] = [160, 160, 160]
        grid["origin"] = [0.0, 0.0, 0.0]

    for element in solver["mesh"]["elements"]:
        transformed_intervals = []
        for interval in element.get("intervals", []):
            points = np.asarray(interval, dtype=np.int64) * SCALE
            low = points.min(axis=0)
            high = points.max(axis=0)

            if symmetric:
                if high[1] < SYMMETRY_INDEX or high[2] < SYMMETRY_INDEX:
                    continue
                low = np.maximum(low, [0, SYMMETRY_INDEX, SYMMETRY_INDEX])
                high = np.minimum(high, [160, 160, 160])
                low[1:] -= SYMMETRY_INDEX
                high[1:] -= SYMMETRY_INDEX

            transformed_intervals.append([low.tolist(), high.tolist()])
        element["intervals"] = transformed_intervals

    pml = {"type": "pml", "layers": 10, "order": 2, "reflection": 0.001}
    if symmetric:
        boundary = {
            "xLower": copy.deepcopy(pml),
            "xUpper": copy.deepcopy(pml),
            "yLower": {"type": "pec"},
            "yUpper": copy.deepcopy(pml),
            "zLower": {"type": "pmc"},
            "zUpper": copy.deepcopy(pml),
        }
    else:
        boundary = {"all": pml}
    solver["boundary"].clear()
    solver["boundary"].update(boundary)

    solver["general"]["numberOfSteps"] = NUMBER_OF_STEPS
    solver["probes"][:] = [
        probe for probe in solver["probes"] if probe["type"] == "farField"
    ]
    farfield = solver["probes"][0]
    farfield["theta"] = {"initial": 90.0, "final": 90.0, "step": 1.0}
    farfield["phi"] = {"initial": 180.0, "final": 180.0, "step": 1.0}
    farfield["domain"] = {
        "type": "frequency",
        "initialFrequency": FREQUENCY_RANGE[0],
        "finalFrequency": FREQUENCY_RANGE[1],
        "numberOfFrequencies": NUMBER_OF_FREQUENCIES,
        "frequencySpacing": "logarithmic",
    }


def run_case(executable, output_folder, symmetric):
    output_folder.mkdir(parents=True, exist_ok=True)
    solver = FDTD(
        input_filename=CASE,
        path_to_exe=str(executable),
        run_in_folder=str(output_folder),
    )
    configure_case(solver, symmetric)
    solver.cleanUp()
    solver.run()

    probe_folders = solver.getSolvedProbeFolders("farfield")
    if len(probe_folders) != 1:
        raise RuntimeError(f"Expected one far-field output, found {len(probe_folders)}")

    data = Probe(probe_folders[0]).data
    mask = np.isclose(data["Theta"], 90.0) & np.isclose(data["Phi"], 180.0)
    samples = data.loc[mask, ["freq", "rcs_geom"]].sort_values("freq")
    if len(samples) != NUMBER_OF_FREQUENCIES:
        raise RuntimeError(f"Expected {NUMBER_OF_FREQUENCIES} RCS samples, found {len(samples)}")
    return samples["freq"].to_numpy(), samples["rcs_geom"].to_numpy()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exe", type=Path, default=DEFAULT_EXE, help="Path to semba-fdtd.exe")
    args = parser.parse_args()
    if not args.exe.is_file():
        raise FileNotFoundError(f"Solver executable not found: {args.exe}")

    full_frequency, full_rcs = run_case(args.exe, OUTPUT / "full", symmetric=False)
    symmetry_frequency, symmetry_rcs = run_case(
        args.exe, OUTPUT / "pec_pmc_symmetry", symmetric=True
    )
    if not np.array_equal(full_frequency, symmetry_frequency):
        raise RuntimeError("The full and symmetry cases produced different frequency samples")
    if not np.all(np.isfinite(full_rcs)) or not np.all(np.isfinite(symmetry_rcs)):
        raise RuntimeError("RCS output contains non-finite values")

    floor = np.finfo(float).tiny
    full_dbsm = 10.0 * np.log10(np.maximum(full_rcs, floor))
    symmetry_dbsm = 10.0 * np.log10(np.maximum(symmetry_rcs, floor))
    difference_db = symmetry_dbsm - full_dbsm

    OUTPUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT / "sphere_rcs_comparison.csv"
    np.savetxt(
        csv_path,
        np.column_stack((full_frequency, full_rcs, symmetry_rcs, difference_db)),
        delimiter=",",
        header="frequency_hz,full_rcs_m2,symmetry_rcs_m2,symmetry_minus_full_db",
        comments="",
    )

    figure, (rcs_axis, difference_axis) = plt.subplots(
        2, 1, figsize=(9, 7), sharex=True, constrained_layout=True
    )
    frequency_mhz = full_frequency / 1e6
    rcs_axis.semilogx(frequency_mhz, full_dbsm, label="Full domain")
    rcs_axis.semilogx(frequency_mhz, symmetry_dbsm, "--", label="PEC/PMC quarter domain")
    rcs_axis.set_ylabel("Backscatter RCS [dBsm]")
    rcs_axis.set_title("1 m radius PEC sphere, Gaussian plane wave, PML boundaries")
    rcs_axis.grid(True, which="both", alpha=0.3)
    rcs_axis.legend()
    difference_axis.semilogx(frequency_mhz, difference_db, color="tab:green")
    difference_axis.axhline(0.0, color="black", linewidth=0.8)
    difference_axis.set_xlabel("Frequency [MHz]")
    difference_axis.set_ylabel("Symmetry - full [dB]")
    difference_axis.grid(True, which="both", alpha=0.3)

    figure_path = OUTPUT / "sphere_rcs_comparison.png"
    figure.savefig(figure_path, dpi=180)
    plt.close(figure)

    print(f"Figure: {figure_path}")
    print(f"Data:   {csv_path}")
    print(f"Maximum absolute RCS difference: {np.max(np.abs(difference_db)):.3f} dB")


if __name__ == "__main__":
    main()
