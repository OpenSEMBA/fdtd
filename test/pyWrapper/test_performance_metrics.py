import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
CASES_ROOT = REPO_ROOT / "testData" / "cases"
DEFAULT_MAX_STEPS = 2000
DEFAULT_CASES = {
    "nodal_source": CASES_ROOT / "nodalSource" / "nodalSource.fdtd.json",
}


pytestmark = pytest.mark.performance


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    input_file: Path
    steps: int


@dataclass(frozen=True)
class BenchmarkExecutable:
    name: str
    path: Path


@dataclass
class BenchmarkResult:
    case: str
    executable: str
    exe_path: str
    build_type: str
    threads: int
    repeat: int
    requested_steps: int
    reported_steps: int
    elapsed_s: float
    steps_per_s: float
    returncode: int
    run_dir: str
    stdout_tail: str
    stderr_tail: str


def _executable_name() -> str:
    return "semba-fdtd.exe" if sys.platform == "win32" else "semba-fdtd"


def _default_fortran_exe() -> Path:
    for env_name in ("SEMBA_PERF_FORTRAN_EXE", "SEMBA_FORTRAN_EXE"):
        if os.environ.get(env_name):
            return Path(os.environ[env_name])

    exe = _executable_name()
    for build_dir in (
        "build_fortran_rel",
        "build_fortran_nomtln_rel",
        "build_fortran",
        "build",
    ):
        candidate = REPO_ROOT / build_dir / "bin" / exe
        if candidate.exists():
            return candidate
    return REPO_ROOT / "build_fortran" / "bin" / exe


def _default_cpp_exe() -> Path:
    for env_name in ("SEMBA_PERF_CPP_EXE", "SEMBA_EXE"):
        if os.environ.get(env_name):
            return Path(os.environ[env_name])

    exe = "semba-fdtd-cpp.exe" if sys.platform == "win32" else "semba-fdtd-cpp"
    for build_dir in (
        "cpp_build_nomtln_rel",
        "cpp_build_release",
        "cpp_build_nomtln",
        "cpp_build",
    ):
        candidate = REPO_ROOT / build_dir / "bin" / exe
        if candidate.exists():
            return candidate
    return REPO_ROOT / "cpp_build_nomtln" / "bin" / exe


def _build_type_for_exe(exe_path: Path) -> str:
    build_dir = exe_path.resolve().parent.parent
    cache = build_dir / "CMakeCache.txt"
    if not cache.exists():
        return "unknown"

    for line in cache.read_text(errors="replace").splitlines():
        if line.startswith("CMAKE_BUILD_TYPE:"):
            return line.split("=", 1)[1] or "unknown"
    return "unknown"


def _used_files(input_data: dict) -> list[str]:
    result = []

    for source in input_data.get("sources", []):
        if "magnitudeFile" in source:
            result.append(source["magnitudeFile"])

    for probe in input_data.get("probes", []):
        if "magnitudeFile" in probe:
            result.append(probe["magnitudeFile"])

    for material in input_data.get("materials", []):
        for termination in material.get("terminations", []):
            if "file" in termination and termination["file"] not in result:
                result.append(termination["file"])

    return result


def _prepare_case_run(input_file: Path, run_dir: Path, max_steps: int) -> BenchmarkCase:
    run_dir.mkdir(parents=True, exist_ok=True)
    case_dir = input_file.parent
    input_data = json.loads(input_file.read_text())

    original_steps = int(input_data.get("general", {}).get("numberOfSteps", max_steps))
    steps = min(original_steps, max_steps)
    input_data.setdefault("general", {})["numberOfSteps"] = steps

    for used_file in _used_files(input_data):
        source = case_dir / used_file
        if source.exists():
            shutil.copy2(source, run_dir / Path(used_file).name)

    target_input = run_dir / input_file.name
    target_input.write_text(json.dumps(input_data, indent=4) + "\n")
    return BenchmarkCase(name=input_file.stem.replace(".fdtd", ""), input_file=target_input, steps=steps)


def _parse_reported_steps(stdout: str, fallback: int) -> int:
    match = re.search(r"Running FDTD:\s+(\d+)\s+steps", stdout)
    if match:
        return int(match.group(1))
    return fallback


def _tail(text: str, max_chars: int = 2000) -> str:
    if len(text) <= max_chars:
        return text
    return text[-max_chars:]


def _run_one(
    case_name: str,
    input_file: Path,
    requested_steps: int,
    executable: BenchmarkExecutable,
    threads: int,
    repeat: int,
    timeout_s: float,
    run_dir: Path,
) -> BenchmarkResult:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env.setdefault("OMP_PROC_BIND", "false")

    command = [str(executable.path.resolve()), "-i", input_file.name]
    start = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=run_dir,
        env=env,
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    elapsed = time.perf_counter() - start
    reported_steps = _parse_reported_steps(completed.stdout, requested_steps)

    if completed.returncode != 0:
        steps_per_s = 0.0
    else:
        steps_per_s = reported_steps / elapsed if elapsed > 0.0 else 0.0

    return BenchmarkResult(
        case=case_name,
        executable=executable.name,
        exe_path=str(executable.path.resolve()),
        build_type=_build_type_for_exe(executable.path),
        threads=threads,
        repeat=repeat,
        requested_steps=requested_steps,
        reported_steps=reported_steps,
        elapsed_s=elapsed,
        steps_per_s=steps_per_s,
        returncode=completed.returncode,
        run_dir=str(run_dir),
        stdout_tail=_tail(completed.stdout),
        stderr_tail=_tail(completed.stderr),
    )


def run_benchmark_matrix(
    cases: dict[str, Path],
    executables: list[BenchmarkExecutable],
    threads: list[int],
    max_steps: int = DEFAULT_MAX_STEPS,
    repeats: int = 1,
    timeout_s: float = 600.0,
    work_dir: Path | None = None,
    keep_workdirs: bool = False,
) -> list[BenchmarkResult]:
    if work_dir is None:
        manager = tempfile.TemporaryDirectory(prefix="semba_perf_")
        root = Path(manager.name)
    else:
        manager = None
        root = work_dir
        root.mkdir(parents=True, exist_ok=True)

    try:
        results = []
        for case_name, case_input in cases.items():
            if not case_input.exists():
                raise FileNotFoundError(f"benchmark case not found: {case_input}")

            for executable in executables:
                if not executable.path.exists():
                    raise FileNotFoundError(
                        f"{executable.name} executable not found: {executable.path}"
                    )

                for thread_count in threads:
                    for repeat in range(1, repeats + 1):
                        run_dir = root / case_name / executable.name / f"t{thread_count}" / f"r{repeat}"
                        prepared = _prepare_case_run(case_input, run_dir, max_steps)
                        result = _run_one(
                            case_name=case_name,
                            input_file=prepared.input_file,
                            requested_steps=prepared.steps,
                            executable=executable,
                            threads=thread_count,
                            repeat=repeat,
                            timeout_s=timeout_s,
                            run_dir=run_dir,
                        )
                        results.append(result)
                        if result.returncode != 0:
                            raise RuntimeError(
                                f"{result.executable} failed for {result.case} "
                                f"with OMP_NUM_THREADS={result.threads}; "
                                f"see {result.run_dir}"
                            )

        if keep_workdirs:
            manager = None
        return results
    finally:
        if manager is not None:
            manager.cleanup()


def _default_threads() -> list[int]:
    if os.environ.get("SEMBA_PERF_THREADS"):
        return _parse_threads(os.environ["SEMBA_PERF_THREADS"])

    cpu_count = os.cpu_count() or 1
    if cpu_count == 1:
        return [1]
    return [1, cpu_count]


def _parse_threads(value: str) -> list[int]:
    threads = []
    for item in value.split(","):
        item = item.strip()
        if not item:
            continue
        threads.append(int(item))
    if not threads:
        raise ValueError("at least one thread count is required")
    return threads


def _parse_case_args(case_args: list[str]) -> dict[str, Path]:
    if not case_args:
        return DEFAULT_CASES

    cases = {}
    for item in case_args:
        if "=" in item:
            name, path = item.split("=", 1)
        else:
            path = item
            name = Path(path).stem.replace(".fdtd", "")

        input_path = Path(path)
        if not input_path.is_absolute():
            input_path = REPO_ROOT / input_path
        cases[name] = input_path
    return cases


def _results_as_dict(results: list[BenchmarkResult]) -> list[dict]:
    return [asdict(result) for result in results]


def _format_markdown(results: list[BenchmarkResult]) -> str:
    lines = [
        "| case | executable | build | threads | repeat | steps | elapsed_s | steps_per_s |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    for result in results:
        lines.append(
            f"| {result.case} | {result.executable} | {result.build_type} | "
            f"{result.threads} | {result.repeat} | {result.reported_steps} | "
            f"{result.elapsed_s:.6f} | {result.steps_per_s:.3f} |"
        )
    return "\n".join(lines) + "\n"


def _print_release_warnings(results: list[BenchmarkResult]) -> None:
    for result in results:
        if result.build_type.lower() != "release":
            print(
                f"warning: {result.executable} build type is {result.build_type}; "
                "use Release builds for performance numbers",
                file=sys.stderr,
            )


@pytest.mark.skipif(
    os.environ.get("SEMBA_RUN_PERFORMANCE_BENCHMARK") != "1",
    reason="set SEMBA_RUN_PERFORMANCE_BENCHMARK=1 to run performance benchmarks",
)
def test_fortran_cpp_performance_metrics(tmp_path):
    results = run_benchmark_matrix(
        cases=DEFAULT_CASES,
        executables=[
            BenchmarkExecutable("fortran", _default_fortran_exe()),
            BenchmarkExecutable("cpp", _default_cpp_exe()),
        ],
        threads=_default_threads(),
        max_steps=int(os.environ.get("SEMBA_PERF_MAX_STEPS", DEFAULT_MAX_STEPS)),
        repeats=int(os.environ.get("SEMBA_PERF_REPEATS", "1")),
        timeout_s=float(os.environ.get("SEMBA_PERF_TIMEOUT", "600")),
        work_dir=tmp_path,
        keep_workdirs=True,
    )
    _print_release_warnings(results)
    print(_format_markdown(results))
    assert all(result.returncode == 0 for result in results)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare Fortran and C++ SEMBA-FDTD runtime metrics."
    )
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        help=(
            "Benchmark case as relative/path.fdtd.json or name=relative/path.fdtd.json. "
            "Defaults to nodalSource/nodalSource.fdtd.json."
        ),
    )
    parser.add_argument("--fortran-exe", default=str(_default_fortran_exe()))
    parser.add_argument("--cpp-exe", default=str(_default_cpp_exe()))
    parser.add_argument(
        "--threads",
        default=",".join(str(v) for v in _default_threads()),
        help="Comma-separated OMP_NUM_THREADS values.",
    )
    parser.add_argument("--max-steps", type=int, default=DEFAULT_MAX_STEPS)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--work-dir", type=Path, default=None)
    parser.add_argument("--keep-workdirs", action="store_true")
    parser.add_argument("--json", type=Path, default=None, help="Optional JSON output file.")
    parser.add_argument(
        "--markdown",
        type=Path,
        default=None,
        help="Optional Markdown table output file.",
    )
    args = parser.parse_args(argv)

    results = run_benchmark_matrix(
        cases=_parse_case_args(args.case),
        executables=[
            BenchmarkExecutable("fortran", Path(args.fortran_exe)),
            BenchmarkExecutable("cpp", Path(args.cpp_exe)),
        ],
        threads=_parse_threads(args.threads),
        max_steps=args.max_steps,
        repeats=args.repeats,
        timeout_s=args.timeout,
        work_dir=args.work_dir,
        keep_workdirs=args.keep_workdirs,
    )

    _print_release_warnings(results)
    markdown = _format_markdown(results)
    print(markdown, end="")

    if args.json is not None:
        args.json.write_text(json.dumps(_results_as_dict(results), indent=2) + "\n")

    if args.markdown is not None:
        args.markdown.write_text(markdown)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
