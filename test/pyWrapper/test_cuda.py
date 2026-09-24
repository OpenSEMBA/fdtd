"""CUDA-only case smoke and CPU-golden probe checks.

Collected/run only when the selected build has SEMBA_FDTD_ENABLE_CUDA=ON
(same gating pattern as MPI: set SEMBA_FDTD_ENABLE_CUDA=ON and SEMBA_EXE to
the rls-cuda binary). Otherwise tests are skipped and never launched.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from test.utils.utils import (
    CASES_FOLDER,
    SEMBA_EXE,
    no_cuda_skip,
    FDTD,
)

CUDA_CASES = Path(CASES_FOLDER) / "cuda"


def _load_probe_matrix(path: Path) -> np.ndarray:
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


def _assert_probes_close(cpu_dir: Path, gpu_dir: Path, rtol=1e-3):
    """Compare CPU vs CUDA probe series (max-norm relative).

    Uses max|Δ| / max|CPU| rather than elementwise np.allclose: mid-ramp
    samples with tiny absolute values make atol=1e-6*scale too tight even when
    the global relative error is well under rtol.

    rtol=1e-3 covers FP32 Yee/CPML/PW on device for the dominant field
    components (see 07-cuda-validation.md). Near-null polarization components
    (scale ≪ dominant) use an absolute floor so relative noise is not a false
    fail — still far below "garbage" on the real signal.
    """
    probes = sorted(cpu_dir.glob("*_tm.dat"))
    assert probes, f"no *_tm.dat probes under {cpu_dir}"
    loaded = []
    for cpu_probe in probes:
        gpu_probe = gpu_dir / cpu_probe.name
        assert gpu_probe.is_file(), f"missing CUDA probe {gpu_probe.name}"
        a = _load_probe_matrix(cpu_probe)
        b = _load_probe_matrix(gpu_probe)
        assert a.shape == b.shape, f"{cpu_probe.name}: shape {a.shape} vs {b.shape}"
        assert a.ndim == 2 and a.shape[1] >= 2, f"{cpu_probe.name}: unexpected layout"
        assert np.allclose(a[:, 0], b[:, 0], rtol=0.0, atol=1e-18), (
            f"{cpu_probe.name}: time axis mismatch"
        )
        loaded.append((cpu_probe.name, a[:, 1:], b[:, 1:]))

    dom = max(float(np.max(np.abs(a))) for _, a, _ in loaded)
    assert dom > 1e-12, "all CPU probes ~empty; not a valid golden"
    null_floor = 1e-4 * dom

    for name, a_fields, b_fields in loaded:
        scale = float(np.max(np.abs(a_fields)))
        maxdiff = float(np.max(np.abs(a_fields - b_fields)))
        if scale < null_floor:
            assert maxdiff <= 1e-6, (
                f"{name}: near-null abs maxdiff={maxdiff:.3e} "
                f"(scale={scale:.3e}, dom={dom:.3e})"
            )
            continue
        reldiff = maxdiff / scale
        assert reldiff <= rtol, (
            f"{name}: maxdiff={maxdiff:.3e} reldiff={reldiff:.3e} "
            f"scale={scale:.3e} (limit rtol={rtol})"
        )


def _run_cpu_vs_cuda_golden(tmp_path, monkeypatch, case_name: str):
    """OMP=1 CPU golden vs CUDA; probes must match (not just 'finished')."""
    input_filename = str(CUDA_CASES / f"{case_name}.fdtd.json")
    assert Path(input_filename).is_file(), input_filename
    cpu_dir = tmp_path / "cpu"
    gpu_dir = tmp_path / "gpu"
    cpu_dir.mkdir()
    gpu_dir.mkdir()

    monkeypatch.delenv("SEMBA_FDTD_DEVICE", raising=False)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    cpu = FDTD(
        input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=cpu_dir,
    )
    cpu.run()
    assert cpu.hasFinishedSuccessfully()

    monkeypatch.setenv("SEMBA_FDTD_DEVICE", "cuda")
    gpu = FDTD(
        input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=gpu_dir,
    )
    gpu.run()
    assert gpu.hasFinishedSuccessfully()

    _assert_probes_close(cpu_dir, gpu_dir)


@no_cuda_skip
@pytest.mark.cuda
@pytest.mark.planewave
@pytest.mark.parametrize(
    "case_name",
    [
        "box_pml_smallpw_50",
        "box_pml_smallpw_100",
        "box_pml_smallpw_256",
        "box_pml_smallpw_300",
        "box_mur_smallpw_100",
    ],
)
def test_cuda_matches_cpu_golden_pml(tmp_path, case_name, monkeypatch):
    """Correctness gate: OMP=1 CPU golden vs CUDA probes.

    rtol 1e-3 on dominant probes (device Yee+CPML/Mur+PW). Includes 256^3 / 300^3
    so large-case speed claims are backed by matching non-trivial probe signals.
    Mur uses the same device-resident first-order ABC as the host path.
    """
    _run_cpu_vs_cuda_golden(tmp_path, monkeypatch, case_name)
