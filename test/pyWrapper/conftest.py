"""Pytest support for retaining solver artifacts after a test run."""

import os
import re
import shutil
import sys
from pathlib import Path

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def pytest_addoption(parser):
    parser.addoption(
        "--keep-artifacts",
        action="store_true",
        help="Copy each tmp_path test directory to test-artifacts for inspection.",
    )


def _keep_artifacts(config):
    value = os.getenv("SEMBA_FDTD_KEEP_ARTIFACTS", "").lower()
    return config.getoption("--keep-artifacts") or value in {"1", "true", "on", "yes"}


@pytest.fixture(autouse=True)
def retain_test_artifacts(request, tmp_path):
    yield

    if not _keep_artifacts(request.config):
        return

    destination = Path(
        os.getenv("SEMBA_FDTD_ARTIFACT_DIR", request.config.rootpath / "test-artifacts")
    )
    test_name = re.sub(r"[^A-Za-z0-9_.-]+", "_", request.node.nodeid)
    target = destination / test_name
    shutil.copytree(tmp_path, target, dirs_exist_ok=True)
    print(f"Retained test artifacts: {target.resolve()}")
