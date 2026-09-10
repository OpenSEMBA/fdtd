"""Shared pytest configuration."""

import os

from test.utils.build_resolver import build_feature_enabled


def pytest_configure():
    if build_feature_enabled("SEMBA_FDTD_ENABLE_MPI"):
        os.environ.setdefault("OMP_NUM_THREADS", "1")
