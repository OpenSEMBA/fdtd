import argparse
import json
import os
import re
import shutil
import shlex
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
DEFAULT_MPI_CASES = {
    "holland": CASES_ROOT / "holland_mpi" / "holland1981.fdtd.json",
    "towel": CASES_ROOT / "towelHanger" / "towelHanger_mpi.fdtd.json",
}
DEFAULT_CYBONERA_CASE = Path(
    os.environ.get(
        "SEMBA_PERF_CYBONERA_CASE",
        "/home/luis/cybonera/analysis/cybonera_w4_corrected/cybonera_w4_corrected.fdtd.json",
    )
)


pytestmark = pytest.mark.performance


@dataclass(frozen=True)
class BenchmarkCase:
    name: str
    input_file: Path
    steps: int
    mcells_per_step: float


@dataclass(frozen=True)
class BenchmarkExecutable:
    name: str
    path: Path
    mode: str = "default"


@dataclass
class BenchmarkResult:
    case: str
    executable: str
    mode: str
    exe_path: str
    build_type: str
    ranks: int
    threads: int
    total_workers: int
    mpi_command: str
    repeat: int
    requested_steps: int
    reported_steps: int
    mcells_per_step: float
    elapsed_s: float
    steps_per_s: float
    mcells_per_s: float
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


def _cmake_cache_value(exe_path: Path, key: str) -> str | None:
    build_dir = exe_path.resolve().parent.parent
    cache = build_dir / "CMakeCache.txt"
    if not cache.exists():
        return None
    prefix = f"{key}:"
    for line in cache.read_text(errors="replace").splitlines():
        if line.startswith(prefix):
            return line.split("=", 1)[1].strip()
    return None


def _cmake_bool(exe_path: Path, key: str) -> bool | None:
    raw = _cmake_cache_value(exe_path, key)
    if raw is None:
        return None
    value = raw.strip().upper()
    if value in {"ON", "TRUE", "1", "YES"}:
        return True
    if value in {"OFF", "FALSE", "0", "NO"}:
        return False
    return None


def _solver_signature(exe_path: Path) -> dict[str, bool | None]:
    return {
        "mpi": _cmake_bool(exe_path, "SEMBA_FDTD_ENABLE_MPI"),
        "mtln": _cmake_bool(exe_path, "SEMBA_FDTD_ENABLE_MTLN"),
        "double_precision": _cmake_bool(exe_path, "SEMBA_FDTD_ENABLE_DOUBLE_PRECISION"),
        "strict_rounding": _cmake_bool(exe_path, "SEMBA_FDTD_ENABLE_STRICT_FORTRAN_ROUNDING"),
    }


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


def _extract_mcells_per_step(input_data: dict) -> float:
    try:
        number_of_cells = input_data["mesh"]["grid"]["numberOfCells"]
        if len(number_of_cells) != 3:
            return 0.0
        total_cells = (
            float(number_of_cells[0]) *
            float(number_of_cells[1]) *
            float(number_of_cells[2])
        )
        return total_cells / 1.0e6
    except Exception:
        return 0.0


def _sanitize_case_for_benchmark(input_data: dict, strip_case_flags: bool) -> None:
    if not strip_case_flags:
        return
    general = input_data.setdefault("general", {})
    if isinstance(general, dict):
        general.pop("additionalArguments", None)


def _prepare_case_run(
    input_file: Path,
    run_dir: Path,
    max_steps: int,
    strip_case_flags: bool,
) -> BenchmarkCase:
    run_dir.mkdir(parents=True, exist_ok=True)
    case_dir = input_file.parent
    input_data = json.loads(input_file.read_text())
    _sanitize_case_for_benchmark(input_data, strip_case_flags)

    original_steps = int(input_data.get("general", {}).get("numberOfSteps", max_steps))
    steps = min(original_steps, max_steps)
    input_data.setdefault("general", {})["numberOfSteps"] = steps

    for used_file in _used_files(input_data):
        source = case_dir / used_file
        if source.exists():
            shutil.copy2(source, run_dir / Path(used_file).name)

    target_input = run_dir / input_file.name
    target_input.write_text(json.dumps(input_data, indent=4) + "\n")
    return BenchmarkCase(
        name=input_file.stem.replace(".fdtd", ""),
        input_file=target_input,
        steps=steps,
        mcells_per_step=_extract_mcells_per_step(input_data),
    )


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
    mcells_per_step: float,
    executable: BenchmarkExecutable,
    ranks: int,
    threads: int,
    repeat: int,
    timeout_s: float,
    run_dir: Path,
    mpi_command_template: str | None,
    force_launcher_for_rank1: bool,
) -> BenchmarkResult:
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)
    env.setdefault("OMP_PROC_BIND", "false")

    use_launcher = ranks > 1 or mpi_command_template is not None
    if ranks == 1 and force_launcher_for_rank1:
        use_launcher = True
    if use_launcher:
        template = mpi_command_template or "mpirun -np {ranks}"
        mpi_command = template.format(ranks=ranks)
        command = shlex.split(mpi_command) + [
            str(executable.path.resolve()), "-i", input_file.name
        ]
    else:
        mpi_command = ""
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
        mcells_per_s = 0.0
    else:
        steps_per_s = reported_steps / elapsed if elapsed > 0.0 else 0.0
        mcells_per_s = (
            (mcells_per_step * float(reported_steps)) / elapsed
            if elapsed > 0.0 else 0.0
        )

    return BenchmarkResult(
        case=case_name,
        executable=executable.name,
        mode=executable.mode,
        exe_path=str(executable.path.resolve()),
        build_type=_build_type_for_exe(executable.path),
        ranks=ranks,
        threads=threads,
        total_workers=ranks * threads,
        mpi_command=mpi_command,
        repeat=repeat,
        requested_steps=requested_steps,
        reported_steps=reported_steps,
        mcells_per_step=mcells_per_step,
        elapsed_s=elapsed,
        steps_per_s=steps_per_s,
        mcells_per_s=mcells_per_s,
        returncode=completed.returncode,
        run_dir=str(run_dir),
        stdout_tail=_tail(completed.stdout),
        stderr_tail=_tail(completed.stderr),
    )


def run_benchmark_matrix(
    cases: dict[str, Path],
    executables: list[BenchmarkExecutable],
    threads: list[int],
    ranks: list[int] | None = None,
    max_steps: int = DEFAULT_MAX_STEPS,
    repeats: int = 1,
    timeout_s: float = 600.0,
    work_dir: Path | None = None,
    keep_workdirs: bool = False,
    mpi_command: str | None = None,
    strip_case_flags: bool = True,
) -> list[BenchmarkResult]:
    if ranks is None:
        ranks = [1]
    if work_dir is None:
        manager = tempfile.TemporaryDirectory(prefix="semba_perf_")
        root = Path(manager.name)
    else:
        manager = None
        root = work_dir
        root.mkdir(parents=True, exist_ok=True)

    force_launcher_for_rank1 = bool(mpi_command) or any(rank > 1 for rank in ranks)
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

                for rank_count in ranks:
                    for thread_count in threads:
                        for repeat in range(1, repeats + 1):
                            run_dir = (
                                root / case_name / executable.name /
                                f"np{rank_count}_t{thread_count}" / f"r{repeat}"
                            )
                            prepared = _prepare_case_run(
                                case_input,
                                run_dir,
                                max_steps,
                                strip_case_flags=strip_case_flags,
                            )
                            result = _run_one(
                                case_name=case_name,
                                input_file=prepared.input_file,
                                requested_steps=prepared.steps,
                                mcells_per_step=prepared.mcells_per_step,
                                executable=executable,
                                ranks=rank_count,
                                threads=thread_count,
                                repeat=repeat,
                                timeout_s=timeout_s,
                                run_dir=run_dir,
                                mpi_command_template=mpi_command,
                                force_launcher_for_rank1=force_launcher_for_rank1,
                            )
                            results.append(result)
                            if result.returncode != 0:
                                raise RuntimeError(
                                    f"{result.executable} failed for {result.case} "
                                    f"with ranks={result.ranks} and "
                                    f"OMP_NUM_THREADS={result.threads}; "
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


def _default_ranks() -> list[int]:
    if os.environ.get("SEMBA_PERF_RANKS"):
        return _parse_threads(os.environ["SEMBA_PERF_RANKS"])
    return [1]


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


def _default_mpi_cases() -> dict[str, Path]:
    cases = dict(DEFAULT_MPI_CASES)
    if DEFAULT_CYBONERA_CASE.exists():
        cases["cybonera"] = DEFAULT_CYBONERA_CASE
    return cases


def _resolve_case_preset(case_preset: str) -> dict[str, Path]:
    if case_preset == "default":
        return DEFAULT_CASES
    if case_preset == "mpi":
        return _default_mpi_cases()
    raise ValueError(f"unsupported case preset: {case_preset}")


def _reference_mcells(results: list[BenchmarkResult]) -> dict[tuple[str, int, int, int], float]:
    refs: dict[tuple[str, int, int, int], float] = {}
    for result in results:
        if result.mode != "reference" and result.executable != "fortran":
            continue
        refs[(result.case, result.ranks, result.threads, result.repeat)] = result.mcells_per_s
    return refs


def _results_as_dict(results: list[BenchmarkResult]) -> list[dict]:
    return [asdict(result) for result in results]


def _format_markdown(results: list[BenchmarkResult]) -> str:
    refs = _reference_mcells(results)
    lines = [
        "| case | executable | mode | build | ranks | threads | workers | repeat | steps | elapsed_s | steps_per_s | mcells_per_s | vs_fortran |",
        "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for result in results:
        ref = refs.get((result.case, result.ranks, result.threads, result.repeat), 0.0)
        ratio = ""
        if ref > 0.0 and result.returncode == 0:
            ratio = f"{(result.mcells_per_s / ref):.3f}x"
        lines.append(
            f"| {result.case} | {result.executable} | {result.mode} | {result.build_type} | "
            f"{result.ranks} | {result.threads} | {result.total_workers} | "
            f"{result.repeat} | {result.reported_steps} | {result.elapsed_s:.6f} | "
            f"{result.steps_per_s:.3f} | {result.mcells_per_s:.3f} | {ratio} |"
        )
    return "\n".join(lines) + "\n"


def _print_release_warnings(results: list[BenchmarkResult]) -> None:
    warned_builds: set[tuple[str, str]] = set()
    for result in results:
        key = (result.executable, result.exe_path)
        if key in warned_builds:
            continue
        warned_builds.add(key)
        if result.build_type.lower() != "release":
            print(
                f"warning: {result.executable} build type is {result.build_type}; "
                "use Release builds for performance numbers",
                file=sys.stderr,
            )

    reference: BenchmarkResult | None = next(
        (
            result for result in results
            if result.executable == "fortran" or result.mode == "reference"
        ),
        None,
    )
    if reference is None:
        return
    reference_signature = _solver_signature(Path(reference.exe_path))

    warned_signature: set[tuple[str, str]] = set()
    for result in results:
        key = (result.executable, result.exe_path)
        if key in warned_signature:
            continue
        warned_signature.add(key)
        if result.exe_path == reference.exe_path:
            continue
        signature = _solver_signature(Path(result.exe_path))
        mismatches = []
        for field, ref_value in reference_signature.items():
            value = signature.get(field)
            if ref_value is None or value is None:
                continue
            if value != ref_value:
                mismatches.append(f"{field}: ref={ref_value} candidate={value}")
        if mismatches:
            print(
                f"warning: solver configuration mismatch for {result.executable} "
                f"({result.exe_path}) vs reference ({reference.exe_path}): "
                + ", ".join(mismatches),
                file=sys.stderr,
            )


@pytest.mark.skipif(
    os.environ.get("SEMBA_RUN_PERFORMANCE_BENCHMARK") != "1",
    reason="set SEMBA_RUN_PERFORMANCE_BENCHMARK=1 to run performance benchmarks",
)
def test_fortran_cpp_performance_metrics(tmp_path):
    strip_case_flags_env = os.environ.get("SEMBA_PERF_KEEP_CASE_FLAGS", "").strip().lower()
    strip_case_flags = strip_case_flags_env not in {"1", "true", "yes", "on"}
    results = run_benchmark_matrix(
        cases=DEFAULT_CASES,
        executables=[
            BenchmarkExecutable("fortran", _default_fortran_exe(), mode="reference"),
            BenchmarkExecutable("cpp", _default_cpp_exe(), mode="strict"),
        ],
        threads=_default_threads(),
        ranks=_default_ranks(),
        max_steps=int(os.environ.get("SEMBA_PERF_MAX_STEPS", DEFAULT_MAX_STEPS)),
        repeats=int(os.environ.get("SEMBA_PERF_REPEATS", "1")),
        timeout_s=float(os.environ.get("SEMBA_PERF_TIMEOUT", "600")),
        work_dir=tmp_path,
        keep_workdirs=True,
        mpi_command=os.environ.get("SEMBA_PERF_MPI_COMMAND"),
        strip_case_flags=strip_case_flags,
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
    parser.add_argument(
        "--case-preset",
        choices=("default", "mpi"),
        default="default",
        help=(
            "Case preset used when --case is not provided. "
            "`mpi` includes holland/towel and cybonera if available."
        ),
    )
    parser.add_argument("--fortran-exe", default=str(_default_fortran_exe()))
    parser.add_argument("--cpp-exe", default=str(_default_cpp_exe()))
    parser.add_argument(
        "--cpp-perf-exe",
        default=None,
        help="Optional second C++ executable built in performance mode.",
    )
    parser.add_argument(
        "--threads",
        default=",".join(str(v) for v in _default_threads()),
        help="Comma-separated OMP_NUM_THREADS values.",
    )
    parser.add_argument(
        "--ranks",
        default=",".join(str(v) for v in _default_ranks()),
        help="Comma-separated MPI rank counts.",
    )
    parser.add_argument(
        "--mpi-command",
        default=os.environ.get("SEMBA_PERF_MPI_COMMAND"),
        help=(
            "MPI launcher template, for example 'mpirun -np {ranks}'. "
            "When omitted, rank>1 uses 'mpirun -np {ranks}'. "
            "If any requested rank is >1, rank=1 points also use MPI launcher for consistency."
        ),
    )
    parser.add_argument("--max-steps", type=int, default=DEFAULT_MAX_STEPS)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=600.0)
    parser.add_argument("--work-dir", type=Path, default=None)
    parser.add_argument("--keep-workdirs", action="store_true")
    parser.add_argument(
        "--keep-case-flags",
        action="store_true",
        help="Do not strip case-level additionalArguments when preparing benchmark inputs.",
    )
    parser.add_argument("--json", type=Path, default=None, help="Optional JSON output file.")
    parser.add_argument(
        "--markdown",
        type=Path,
        default=None,
        help="Optional Markdown table output file.",
    )
    args = parser.parse_args(argv)

    executables = [
        BenchmarkExecutable("fortran", Path(args.fortran_exe), mode="reference"),
        BenchmarkExecutable("cpp", Path(args.cpp_exe), mode="strict"),
    ]
    if args.cpp_perf_exe is not None:
        executables.append(
            BenchmarkExecutable("cpp_perf", Path(args.cpp_perf_exe), mode="perf")
        )

    cases = _parse_case_args(args.case) if args.case else _resolve_case_preset(args.case_preset)
    results = run_benchmark_matrix(
        cases=cases,
        executables=executables,
        threads=_parse_threads(args.threads),
        ranks=_parse_threads(args.ranks),
        max_steps=args.max_steps,
        repeats=args.repeats,
        timeout_s=args.timeout,
        work_dir=args.work_dir,
        keep_workdirs=args.keep_workdirs,
        mpi_command=args.mpi_command,
        strip_case_flags=not args.keep_case_flags,
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
