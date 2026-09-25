#!/usr/bin/env python3
"""Run all CUDA golden cases: CPU OMP=1, CPU full OMP, and CUDA; tabulate accuracy."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
EXE = Path(os.environ.get("SEMBA_EXE", ROOT / "build-rls-cuda/bin/semba-fdtd"))
CASES = ROOT / "testData/cases/cuda"
CASE_NAMES = [
    "box_pml_smallpw_50",
    "box_pml_smallpw_100",
    "box_pml_smallpw_256",
    "box_pml_smallpw_300",
    "box_mur_smallpw_100",
    "box_nodal_soft_pml_40",
    "box_nodal_hard_pml_40",
    "box_wire_holland_mur_40",
    "box_wire_holland_mur_100",
    "box_wire_holland_mur_200",
    "box_pmc_nodal_40",
    "box_pmc_nodal_100",
    "box_pmc_nodal_200",
    "box_periodic_nodal_40",
    "box_periodic_nodal_100",
    "box_periodic_nodal_200",
]
WORKDIR = Path(os.environ.get("TMPDIR", "/tmp")) / "fdtd-cuda-full-bench"
NPROC = int(os.environ.get("OMP_FULL", str(os.cpu_count() or 1)))


def load_probe_matrix(path: Path) -> np.ndarray:
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line[0] in "#%" or line[0].isalpha():
            continue
        try:
            rows.append([float(x) for x in line.split()])
        except ValueError:
            continue
    return np.asarray(rows, dtype=float)


def max_reldiff_vs_ref(ref_dir: Path, other_dir: Path) -> tuple[float, str]:
    """Max over probes of max|Δ|/max|ref| (same logic as test_cuda.py)."""
    worst = 0.0
    worst_name = ""
    probes = sorted(ref_dir.glob("*_tm.dat"))
    if not probes:
        raise RuntimeError(f"no probes in {ref_dir}")
    loaded = []
    for p in probes:
        o = other_dir / p.name
        if not o.is_file():
            raise RuntimeError(f"missing {o.name} in {other_dir}")
        a = load_probe_matrix(p)
        b = load_probe_matrix(o)
        if a.shape != b.shape:
            raise RuntimeError(f"{p.name}: shape {a.shape} vs {b.shape}")
        loaded.append((p.name, a[:, 1:], b[:, 1:]))

    dom = max(float(np.max(np.abs(a))) for _, a, _ in loaded)
    null_floor = 1e-4 * dom

    for name, a_fields, b_fields in loaded:
        scale = float(np.max(np.abs(a_fields)))
        maxdiff = float(np.max(np.abs(a_fields - b_fields)))
        if scale < null_floor:
            reldiff = 0.0 if maxdiff <= 1e-6 else float("inf")
        else:
            reldiff = maxdiff / scale
        if reldiff > worst:
            worst = reldiff
            worst_name = name
    return worst, worst_name


def run_case(case_json: Path, run_dir: Path, omp: int, cuda: bool) -> float:
    run_dir.mkdir(parents=True, exist_ok=True)
    case_dir = case_json.parent
    bn = case_json.name
    if not (run_dir / bn).exists():
        shutil.copy2(case_json, run_dir / bn)
    exc = case_dir / "gauss.exc"
    if exc.is_file() and not (run_dir / "gauss.exc").exists():
        shutil.copy2(exc, run_dir / "gauss.exc")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp)
    if cuda:
        env["SEMBA_FDTD_DEVICE"] = "cuda"
    else:
        env.pop("SEMBA_FDTD_DEVICE", None)

    t0 = time.perf_counter()
    r = subprocess.run(
        [str(EXE), "-i", bn],
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - t0
    if r.returncode != 0:
        print(r.stdout, file=sys.stderr)
        print(r.stderr, file=sys.stderr)
        raise RuntimeError(f"run failed {case_json.name} omp={omp} cuda={cuda}")
    return elapsed


def main() -> int:
    if not EXE.is_file():
        print(f"Missing binary: {EXE}", file=sys.stderr)
        return 1

    rows = []
    for name in CASE_NAMES:
        case_json = CASES / f"{name}.fdtd.json"
        base = WORKDIR / name
        t1 = run_case(case_json, base / "cpu_omp1", 1, False)
        t_full = run_case(case_json, base / "cpu_full", NPROC, False)
        t_gpu = run_case(case_json, base / "cuda", 1, True)

        ref_omp1 = base / "cpu_omp1"
        ref_full = base / "cpu_full"
        re_full, _ = max_reldiff_vs_ref(ref_omp1, ref_full)
        re_cuda_vs_omp1, _ = max_reldiff_vs_ref(ref_omp1, base / "cuda")
        re_cuda_vs_full, _ = max_reldiff_vs_ref(ref_full, base / "cuda")
        rows.append(
            {
                "case": name,
                "t_omp1_s": t1,
                "t_full_s": t_full,
                "t_cuda_s": t_gpu,
                "reldiff_full_vs_omp1": re_full,
                "reldiff_cuda_vs_full": re_cuda_vs_full,
                "reldiff_cuda_vs_omp1": re_cuda_vs_omp1,
                "speedup_full": t1 / t_full if t_full > 0 else float("nan"),
                "speedup_cuda": t1 / t_gpu if t_gpu > 0 else float("nan"),
            }
        )

    print(f"Binary: {EXE}")
    print(f"CPU full: OMP_NUM_THREADS={NPROC}")
    print(f"Workdir: {WORKDIR}")
    print()
    hdr = (
        "| Case | t OMP=1 (s) | t CPU full (s) | t CUDA (s) | "
        "max rel err CUDA vs CPU full | max rel err CPU full vs OMP=1 | "
        "speedup full vs OMP=1 | speedup CUDA vs OMP=1 |"
    )
    sep = "|" + "|".join(["---"] * 8) + "|"
    print(hdr)
    print(sep)
    for r in rows:
        print(
            f"| {r['case']} | {r['t_omp1_s']:.3f} | {r['t_full_s']:.3f} | {r['t_cuda_s']:.3f} | "
            f"{r['reldiff_cuda_vs_full']:.3e} | {r['reldiff_full_vs_omp1']:.3e} | "
            f"{r['speedup_full']:.2f}× | {r['speedup_cuda']:.2f}× |"
        )

    gtest = subprocess.run(
        [str(EXE.parent / "fdtd_tests"), "--gtest_filter=cuda.*"],
        capture_output=True,
        text=True,
    )
    print()
    if gtest.returncode == 0:
        print("GoogleTest `cuda.*`: PASSED (exit 0)")
    else:
        print("GoogleTest `cuda.*`: FAILED", file=sys.stderr)
        print(gtest.stdout[-2000:], file=sys.stderr)
        print(gtest.stderr[-2000:], file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
