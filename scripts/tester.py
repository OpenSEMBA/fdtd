#!/usr/bin/env python3
"""semba-fdtd test runner — builds, runs unit + integration tests, produces JSON report + terminal summary."""

import subprocess
import sys
import os
import json
import re
import time
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
BUILD_DIR = ROOT / "build"
TMP_DIR = ROOT / "tmp"
BIN = BUILD_DIR / "bin" / "semba-fdtd"
UNIT_TESTS = BUILD_DIR / "bin" / "fdtd_tests"


def run(cmd, cwd=None, env=None, capture=True):
    """Run a command, return (exit_code, stdout, stderr)."""
    full_env = os.environ.copy()
    if env:
        full_env.update(env)
    result = subprocess.run(
        cmd,
        cwd=cwd or str(ROOT),
        env=full_env,
        capture_output=capture,
        text=True,
    )
    return result.returncode, result.stdout, result.stderr


def detect_cmake_flags():
    """Read CMakeCache.txt to detect which features were enabled."""
    cache = BUILD_DIR / "CMakeCache.txt"
    flags = {}
    if not cache.exists():
        return {"MPI": "OFF", "MTLN": "OFF", "HDF": "OFF", "SMBJSON": "OFF"}
    content = cache.read_text()
    for key, short in [
        ("SEMBA_FDTD_ENABLE_MPI", "MPI"),
        ("SEMBA_FDTD_ENABLE_MTLN", "MTLN"),
        ("SEMBA_FDTD_ENABLE_HDF", "HDF"),
        ("SEMBA_FDTD_ENABLE_SMBJSON", "SMBJSON"),
    ]:
        for line in content.splitlines():
            if line.startswith(key + ":"):
                val = line.split("=", 1)[1].strip().strip('"')
                flags[short] = val
                break
        else:
            flags[short] = "OFF"
    return flags


def parse_git_commit():
    rc, out, _ = run(["git", "rev-parse", "--short", "HEAD"], cwd=str(ROOT))
    if rc == 0:
        return out.strip()
    return "unknown"


def count_pattern(text, pattern):
    return len(re.findall(pattern, text))


def run_build(cmake_flags=None):
    """Configure and build the project. Returns (success, duration_sec)."""
    if cmake_flags is None:
        cmake_flags = {}

    # Detect existing flags from cache, override with user-provided
    existing = detect_cmake_flags() if BUILD_DIR.exists() else {}
    combined = {**existing, **cmake_flags}

    cmake_cmd = [
        "cmake", "-S", ".", "-B", "build",
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    flag_map = {
        "MPI": "SEMBA_FDTD_ENABLE_MPI",
        "MTLN": "SEMBA_FDTD_ENABLE_MTLN",
        "HDF": "SEMBA_FDTD_ENABLE_HDF",
        "SMBJSON": "SEMBA_FDTD_ENABLE_SMBJSON",
    }
    for short, cmake_var in flag_map.items():
        if short in combined:
            cmake_cmd.append(f"-D{cmake_var}={combined[short]}")

    print("=" * 60)
    print("BUILD")
    print("=" * 60)
    print(f"  Command: {' '.join(cmake_cmd)}")

    t0 = time.time()
    rc, out, err = run(cmake_cmd)
    if rc != 0:
        print(f"  cmake configure: FAIL (exit {rc})")
        return False, time.time() - t0
    print(f"  cmake configure: OK")

    build_cmd = ["cmake", "--build", "build", "-j"]
    rc, out, err = run(build_cmd)
    duration = time.time() - t0

    if rc != 0:
        print(f"  cmake build: FAIL (exit {rc})")
        return False, duration

    # Count compiled targets from build output
    lines = out.splitlines()
    targets = count_pattern(out, r"Building [CF]+ object")
    print(f"  cmake build: OK ({targets} targets, {duration:.1f}s)")
    return True, duration, combined


def run_unit_tests():
    """Run GoogleTest suite. Returns (success, duration_sec, output_lines)."""
    print()
    print("=" * 60)
    print("UNIT TESTS (GoogleTest)")
    print("=" * 60)

    if not UNIT_TESTS.exists():
        print(f"  SKIP: {UNIT_TESTS} not found (build first)")
        return False, 0, 0

    t0 = time.time()
    rc, out, err = run([str(UNIT_TESTS)])
    duration = time.time() - t0
    lines = out.splitlines() if out else []

    # Parse GoogleTest summary
    passed = count_pattern(out, r"^\[  PASSED  \] [0-9]+ test")
    failed = count_pattern(out, r"^\[  FAILED  \] [0-9]+ test")
    total = passed + failed

    # Also count individual test results
    passed_tests = count_pattern(out, r"\[  PASSED  \]")
    failed_tests = count_pattern(out, r"\[  FAILED  \]")

    print(f"  Binary: {UNIT_TESTS}")
    print(f"  Duration: {duration:.1f}s")
    print(f"  Output lines: {len(lines)}")
    print(f"  Result: {'PASS' if rc == 0 else 'FAIL'} (exit {rc})")

    return rc == 0, duration, len(lines)


def run_integration_tests(env_vars=None):
    """Run pytest integration tests. Returns (success, duration_sec, output_lines)."""
    print()
    print("=" * 60)
    print("INTEGRATION TESTS (pytest)")
    print("=" * 60)

    pytest_env = os.environ.copy()
    if env_vars:
        pytest_env.update(env_vars)

    t0 = time.time()
    rc, out, err = run(
        [sys.executable or "python3", "-m", "pytest", "test/", "--durations=20", "-v"],
        env=pytest_env,
    )
    duration = time.time() - t0
    lines = out.splitlines() if out else []

    # Parse pytest summary
    passed = count_pattern(out, r" passed")
    failed = count_pattern(out, r" failed")
    skipped = count_pattern(out, r" skipped")

    print(f"  Duration: {duration:.1f}s")
    print(f"  Output lines: {len(lines)}")
    print(f"  Result: {'PASS' if rc == 0 else 'FAIL'} (exit {rc})")

    return rc == 0, duration, len(lines)


def print_summary(build_ok, build_dur, unit_ok, unit_dur, unit_lines,
                  integ_ok, integ_dur, integ_lines, cmake_flags):
    """Print a formatted terminal summary."""
    print()
    print("=" * 60)
    print("TEST REPORT SUMMARY")
    print("=" * 60)
    print(f"  {'Stage':<25} {'Status':<8} {'Duration':<12} {'Details'}")
    print(f"  {'─' * 25} {'─' * 8} {'─' * 12} {'─' * 20}")

    build_status = "PASS" if build_ok else "FAIL"
    print(f"  {'Build':<25} {build_status:<8} {build_dur:.1f}s{'':<6}")

    unit_status = "PASS" if unit_ok else "FAIL"
    print(f"  {'Unit Tests':<25} {unit_status:<8} {unit_dur:.1f}s  {unit_lines} lines")

    integ_status = "PASS" if integ_ok else "FAIL"
    print(f"  {'Integration Tests':<25} {integ_status:<8} {integ_dur:.1f}s  {integ_lines} lines")

    print(f"  {'─' * 25} {'─' * 8} {'─' * 12} {'─' * 20}")

    overall = "PASS" if (build_ok and unit_ok and integ_ok) else "FAIL"
    print(f"  {'OVERALL':<25} {overall:<8}")
    print("=" * 60)


def write_report(build_ok, build_dur, unit_ok, unit_dur, unit_lines,
                 integ_ok, integ_dur, integ_lines, cmake_flags):
    """Write JSON report to tmp/test_report_<timestamp>.json."""
    TMP_DIR.mkdir(exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    report = {
        "timestamp": datetime.now().isoformat(),
        "commit": parse_git_commit(),
        "build_type": "Release",
        "cmake_flags": cmake_flags,
        "build": {"success": build_ok, "duration_sec": round(build_dur, 2)},
        "unit_tests": {"success": unit_ok, "duration_sec": round(unit_dur, 2),
                       "output_lines": unit_lines},
        "integration_tests": {"success": integ_ok, "duration_sec": round(integ_dur, 2),
                              "output_lines": integ_lines},
        "overall": "pass" if (build_ok and unit_ok and integ_ok) else "fail",
    }
    report_path = TMP_DIR / f"test_report_{ts}.json"
    report_path.write_text(json.dumps(report, indent=2))
    print(f"\n  JSON report: {report_path}")
    return report


def main():
    # Accept optional cmake flags as arguments: e.g., --MPI=OFF --MTLN=OFF
    cmake_flags = {}
    i = 1
    while i < len(sys.argv):
        arg = sys.argv[i]
        if arg.startswith("--"):
            key_val = arg[2:]
            if "=" in key_val:
                k, v = key_val.split("=", 1)
                cmake_flags[k] = v
        i += 1

    # Build
    build_result = run_build(cmake_flags)
    if len(build_result) == 3:
        build_ok, build_dur, cmake_flags = build_result
    else:
        build_ok, build_dur = build_result[0], build_result[1]
        cmake_flags = detect_cmake_flags()

    if not build_ok:
        print("\nBuild failed. Skipping tests.")
        unit_ok, unit_dur, unit_lines = False, 0, 0
        integ_ok, integ_dur, integ_lines = False, 0, 0
    else:
        # Unit tests
        unit_result = run_unit_tests()
        unit_ok, unit_dur, unit_lines = unit_result

        # Integration tests
        env_vars = {f"SEMBA_FDTD_ENABLE_{k}": v for k, v in cmake_flags.items()
                    if k in ("MPI", "MTLN", "HDF")}
        integ_result = run_integration_tests(env_vars)
        integ_ok, integ_dur, integ_lines = integ_result

    # Summary + report
    print_summary(build_ok, build_dur, unit_ok, unit_dur, unit_lines,
                  integ_ok, integ_dur, integ_lines, cmake_flags)
    write_report(build_ok, build_dur, unit_ok, unit_dur, unit_lines,
                 integ_ok, integ_dur, integ_lines, cmake_flags)

    # Exit code
    all_ok = build_ok and unit_ok and integ_ok
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
