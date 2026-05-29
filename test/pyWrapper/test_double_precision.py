"""Tests for the SEMBA_FDTD_ENABLE_DOUBLE_PRECISION (CompileWithReal8) flag.

These tests build off two parallel build trees:
  - SP build: ``build/`` (Fortran) and ``cpp_build_nomtln/`` (C++) compiled
    without ``-DSEMBA_FDTD_ENABLE_DOUBLE_PRECISION=ON``.
  - DP build: ``build_dp/`` (Fortran) and ``cpp_build_dp/`` (C++) compiled
    with the flag enabled.

You can override the discovery via env vars:
  - SEMBA_FORTRAN_EXE      (SP Fortran)
  - SEMBA_EXE              (SP C++)
  - SEMBA_FORTRAN_EXE_DP   (DP Fortran)
  - SEMBA_CPP_EXE_DP       (DP C++)

If a build is missing the corresponding test is skipped, never failed.
"""
from pathlib import Path
import os
import shutil

import pytest

from utils import CASES_FOLDER, FDTD, SEMBA_EXE, TEST_DATA_FOLDER


REPO_ROOT = Path(TEST_DATA_FOLDER).parent


def _exe(env_var: str, default_rel: str):
    val = os.environ.get(env_var) or str(REPO_ROOT / default_rel)
    p = Path(val)
    if not p.exists():
        pytest.skip(f"{env_var} not found at {p}")
    return p.resolve()


def _fortran_sp():
    return _exe("SEMBA_FORTRAN_EXE", "build/bin/semba-fdtd")


def _fortran_dp():
    return _exe("SEMBA_FORTRAN_EXE_DP", "build_dp/bin/semba-fdtd")


def _cpp_sp():
    val = os.environ.get("SEMBA_EXE") or SEMBA_EXE
    p = Path(val)
    if not p.exists():
        pytest.skip(f"SEMBA_EXE not found at {p}")
    return p.resolve()


def _cpp_dp():
    return _exe("SEMBA_CPP_EXE_DP", "cpp_build_dp/bin/semba-fdtd-cpp")


def _solved_probe_paths(solver, probe_name):
    """Resolve probe filenames to absolute paths in the solver's run dir."""
    cwd = os.getcwd()
    folder = Path(solver.getFolder())
    os.chdir(folder)
    try:
        return [folder / p for p in solver.getSolvedProbeFilenames(probe_name)]
    finally:
        os.chdir(cwd)


def _run_short_pwbox(executable, run_dir: Path, num_steps: int = 10):
    """Run the planewave-in-box-pec case for a few steps and return the solver."""
    fn = CASES_FOLDER + 'planewave/pw-in-box-pec.fdtd.json'
    run_dir.mkdir(parents=True, exist_ok=True)
    solver = FDTD(input_filename=fn,
                  path_to_exe=str(executable),
                  run_in_folder=run_dir)
    solver['general']['numberOfSteps'] = num_steps
    solver.run()
    return solver


def _read_first_data_line(probe_path: Path):
    with open(probe_path) as f:
        f.readline()  # header
        return f.readline().rstrip("\n")


def _columns(line: str):
    """Split a Fortran fixed-width probe line into its columns.

    The first column is always 27 chars wide (time, ``e27.17e3``).
    Remaining columns are either all 27 (DP) or all 19 (SP).
    """
    time_col = line[:27]
    rest = line[27:]
    if not rest:
        return [time_col]
    field_width = 27 if (len(rest) % 27 == 0) else 19
    fields = [rest[i:i + field_width] for i in range(0, len(rest), field_width)]
    return [time_col] + fields


def _mantissa_digits_after_point(token: str) -> int:
    """Count digits after the decimal point in a Fortran ``e``-format token.

    For ``0.294470937E-006`` this returns 9; for ``0.29447082901330123E-006``
    it returns 17.
    """
    stripped = token.strip()
    point = stripped.find('.')
    if point < 0:
        return 0
    after = stripped[point + 1:]
    digits = 0
    for ch in after:
        if ch in ('e', 'E'):
            break
        if ch.isdigit():
            digits += 1
    return digits


# ---------------------------------------------------------------------------
# Probe output format must match the precision the C++ build was compiled in
# ---------------------------------------------------------------------------
def test_cpp_sp_build_uses_9_digit_field_format(tmp_path):
    """SP C++ build writes fields with the SP Fortran format ``e19.9e3``."""
    cpp = _cpp_sp()
    solver = _run_short_pwbox(cpp, tmp_path / "cpp_sp", num_steps=5)
    probe = _solved_probe_paths(solver, "inbox")[0]
    line = _read_first_data_line(probe)
    cols = _columns(line)
    assert len(cols) >= 2, f"unexpected probe row: {line!r}"
    # SP fields are 19 chars wide and have exactly 9 mantissa digits.
    assert len(cols[1]) == 19, (
        f"expected SP field width 19, got {len(cols[1])} for {cols[1]!r}")
    assert _mantissa_digits_after_point(cols[1]) == 9, (
        f"expected 9 mantissa digits in SP field, got {cols[1]!r}")


def test_cpp_dp_build_uses_17_digit_field_format(tmp_path):
    """DP C++ build writes fields with the DP Fortran format ``e27.17e3``."""
    cpp = _cpp_dp()
    solver = _run_short_pwbox(cpp, tmp_path / "cpp_dp", num_steps=5)
    probe = _solved_probe_paths(solver, "inbox")[0]
    line = _read_first_data_line(probe)
    cols = _columns(line)
    assert len(cols) >= 2, f"unexpected probe row: {line!r}"
    # DP fields are 27 chars wide and have 17 mantissa digits.
    assert len(cols[1]) == 27, (
        f"expected DP field width 27 (CompileWithReal8), got {len(cols[1])} "
        f"for {cols[1]!r}; the C++ build is ignoring the double precision flag"
    )
    assert _mantissa_digits_after_point(cols[1]) == 17, (
        f"expected 17 mantissa digits in DP field, got {cols[1]!r}")


def test_cpp_dp_format_matches_fortran_dp_format(tmp_path):
    """Top column widths must match between Fortran DP and C++ DP."""
    fortran = _fortran_dp()
    cpp = _cpp_dp()

    f_solver = _run_short_pwbox(fortran, tmp_path / "f_dp", num_steps=5)
    c_solver = _run_short_pwbox(cpp, tmp_path / "c_dp", num_steps=5)

    f_line = _read_first_data_line(_solved_probe_paths(f_solver, "inbox")[0])
    c_line = _read_first_data_line(_solved_probe_paths(c_solver, "inbox")[0])
    assert len(f_line) == len(c_line), (
        f"DP probe line lengths differ: Fortran={len(f_line)} C++={len(c_line)}\n"
        f"F: {f_line!r}\nC: {c_line!r}")
    f_cols = _columns(f_line)
    c_cols = _columns(c_line)
    assert [len(c) for c in f_cols] == [len(c) for c in c_cols], (
        f"DP column widths differ: F={[len(c) for c in f_cols]} "
        f"C={[len(c) for c in c_cols]}")


# ---------------------------------------------------------------------------
# The flag must actually change the output (existence check)
# ---------------------------------------------------------------------------
def test_cpp_double_precision_flag_changes_output(tmp_path):
    """SP and DP C++ builds must produce *different* probe outputs."""
    cpp_sp = _cpp_sp()
    cpp_dp = _cpp_dp()
    sp_solver = _run_short_pwbox(cpp_sp, tmp_path / "cpp_sp")
    dp_solver = _run_short_pwbox(cpp_dp, tmp_path / "cpp_dp")
    sp_probe = Path(_solved_probe_paths(sp_solver, "inbox")[0])
    dp_probe = Path(_solved_probe_paths(dp_solver, "inbox")[0])
    assert sp_probe.read_bytes() != dp_probe.read_bytes(), (
        "SP and DP C++ builds produced identical probe output; the "
        "double-precision flag is being ignored at runtime.")


def test_fortran_double_precision_flag_changes_output(tmp_path):
    """Sanity: Fortran SP vs DP also produce different outputs."""
    f_sp = _fortran_sp()
    f_dp = _fortran_dp()
    sp_solver = _run_short_pwbox(f_sp, tmp_path / "f_sp")
    dp_solver = _run_short_pwbox(f_dp, tmp_path / "f_dp")
    sp_probe = Path(_solved_probe_paths(sp_solver, "inbox")[0])
    dp_probe = Path(_solved_probe_paths(dp_solver, "inbox")[0])
    assert sp_probe.read_bytes() != dp_probe.read_bytes()


# ---------------------------------------------------------------------------
# Cross-precision F vs C++ comparison
# ---------------------------------------------------------------------------
def _max_relative_error(path_a, path_b):
    """Compute the max relative error between two probe data files."""
    a_lines = Path(path_a).read_text().splitlines()
    b_lines = Path(path_b).read_text().splitlines()
    assert len(a_lines) == len(b_lines), (
        f"probe line counts differ: {len(a_lines)} vs {len(b_lines)}")
    max_rel = 0.0
    max_abs = 0.0
    for a, b in zip(a_lines[1:], b_lines[1:]):  # skip header
        a_vals = [float(t) for t in a.split()]
        b_vals = [float(t) for t in b.split()]
        for x, y in zip(a_vals, b_vals):
            d = abs(x - y)
            if d > max_abs:
                max_abs = d
            denom = max(abs(x), abs(y))
            if denom > 0:
                rel = d / denom
                if rel > max_rel:
                    max_rel = rel
    return max_rel, max_abs


@pytest.mark.parametrize("probe_name", ["before", "inbox", "after"])
def test_pwbox_dp_fortran_vs_cpp_close(tmp_path, probe_name):
    """At DP, F and C++ probe traces must agree closely.

    A small residual remains because Fortran builds the cell-step array
    by accumulation while C++ takes the JSON value directly, which yields
    ``dxmin`` values that differ by a few parts per billion. After 80
    steps this propagates to a relative error in the field around 1e-5.
    """
    fortran = _fortran_dp()
    cpp = _cpp_dp()
    f_solver = _run_short_pwbox(fortran, tmp_path / "f_dp", num_steps=80)
    c_solver = _run_short_pwbox(cpp, tmp_path / "c_dp", num_steps=80)
    f_probe = _solved_probe_paths(f_solver, probe_name)[0]
    c_probe = _solved_probe_paths(c_solver, probe_name)[0]
    rel, _abs = _max_relative_error(f_probe, c_probe)
    assert rel < 5e-5, (
        f"DP F vs C++ relative error too large for {probe_name}: rel={rel:.3e}")


@pytest.mark.parametrize("probe_name", ["before", "inbox", "after"])
def test_pwbox_sp_fortran_vs_cpp_close(tmp_path, probe_name):
    """At SP, F and C++ probe traces must agree to within float resolution."""
    fortran = _fortran_sp()
    cpp = _cpp_sp()
    f_solver = _run_short_pwbox(fortran, tmp_path / "f_sp", num_steps=80)
    c_solver = _run_short_pwbox(cpp, tmp_path / "c_sp", num_steps=80)
    f_probe = _solved_probe_paths(f_solver, probe_name)[0]
    c_probe = _solved_probe_paths(c_solver, probe_name)[0]
    rel, _abs = _max_relative_error(f_probe, c_probe)
    assert rel < 1e-3, (
        f"SP F vs C++ relative error too large for {probe_name}: rel={rel:.3e}")


def test_dp_more_accurate_than_sp(tmp_path):
    """SP and DP must give close-but-different results; DP error must be smaller.

    This is the main correctness test for the double-precision flag: enabling
    it should *materially* improve agreement between the C++ and Fortran
    builds versus a high-accuracy reference (the Fortran DP run).
    """
    sp_cpp = _cpp_sp()
    dp_cpp = _cpp_dp()
    f_dp = _fortran_dp()

    sp_solver = _run_short_pwbox(sp_cpp, tmp_path / "cpp_sp", num_steps=80)
    dp_solver = _run_short_pwbox(dp_cpp, tmp_path / "cpp_dp", num_steps=80)
    f_solver = _run_short_pwbox(f_dp, tmp_path / "f_dp_ref", num_steps=80)

    sp_probe = _solved_probe_paths(sp_solver, "inbox")[0]
    dp_probe = _solved_probe_paths(dp_solver, "inbox")[0]
    f_probe = _solved_probe_paths(f_solver, "inbox")[0]

    rel_sp, _ = _max_relative_error(sp_probe, f_probe)
    rel_dp, _ = _max_relative_error(dp_probe, f_probe)
    assert rel_dp < rel_sp, (
        f"DP build is no better than SP build vs Fortran DP reference: "
        f"sp_rel={rel_sp:.3e} dp_rel={rel_dp:.3e}")
