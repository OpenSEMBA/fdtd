"""Resolve one compatible solver build for Python tests."""

from __future__ import annotations

import os
from pathlib import Path
from sys import platform


PROJECT_ROOT = Path(__file__).resolve().parents[2]
BUILD_DIRECTORIES = (
    "build-rls",
    "build-dbg",
    "build-rls-mpi",
    "build-dbg-mpi",
    "build-rls-nomtln",
    "build-dbg-nomtln",
    "build-intel-rls",
    "build-intel-dbg",
    "build-intel-rls-nomtln",
    "build",
)
BUILD_FEATURES = (
    "SEMBA_FDTD_ENABLE_MPI",
    "SEMBA_FDTD_ENABLE_MTLN",
)


def _cache_values(build_directory: Path) -> dict[str, str]:
    values = {}
    try:
        with (build_directory / "CMakeCache.txt").open(encoding="utf-8") as cache:
            for line in cache:
                for feature in BUILD_FEATURES:
                    if line.startswith(feature + ":"):
                        values[feature] = line.rstrip().rsplit("=", 1)[-1]
    except OSError:
        pass
    return values


def _matches_requested_features(build_directory: Path) -> bool:
    requested = {
        feature: os.environ[feature]
        for feature in BUILD_FEATURES
        if feature in os.environ
    }
    if not requested:
        return True
    cache_values = _cache_values(build_directory)
    return all(
        cache_values.get(feature) == value for feature, value in requested.items()
    )


def resolve_build(project_root: Path | None = None) -> tuple[Path, Path]:
    """Return a compatible build directory and its solver executable."""
    configured_executable = os.environ.get("SEMBA_EXE")
    if configured_executable:
        executable = Path(configured_executable).resolve()
        return executable.parent.parent, executable

    root = (project_root or PROJECT_ROOT).resolve()
    executable_name = "semba-fdtd.exe" if platform == "win32" else "semba-fdtd"
    for directory_name in BUILD_DIRECTORIES:
        build_directory = root / directory_name
        executable = build_directory / "bin" / executable_name
        if executable.is_file() and _matches_requested_features(build_directory):
            return build_directory, executable

    build_directory = root / "build"
    return build_directory, build_directory / "bin" / executable_name


def solver_executable(project_root: Path | None = None) -> Path:
    return resolve_build(project_root)[1]


def build_feature_enabled(feature: str, project_root: Path | None = None) -> bool:
    configured_value = os.environ.get(feature)
    if configured_value is not None:
        return configured_value == "ON"
    build_directory, _ = resolve_build(project_root)
    return _cache_values(build_directory).get(feature) == "ON"
