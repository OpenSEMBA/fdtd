"""Tests for PEC/PMC symmetry-plane truncation of the sphere RCS case.

The cases under testData/cases/sphere/symmetry/ contain the same PEC sphere
(radius 0.5 m) illuminated by a plane wave propagating towards +y, with the
domain truncated at the sphere-center symmetry planes:

- pol x (E along x): xLower = PEC (E normal to the plane), zLower = PMC (H normal)
- pol z (E along z): xLower = PMC (H normal), zLower = PEC (E normal)

Each case also carries a time-domain point probe ("watchdog") outside the
sphere, which can be used to monitor convergence and stop the run early with
a void 'stopflushing' file. All cases reproduce the same analytical Mie
backscatter RCS as the full-domain case when run to convergence.
"""

import numpy as np
import pytest

from test.utils.utils import *

SYMMETRY_FOLDER = "sphere/symmetry/"

SPHERE_SYMMETRY_CASES = [
    "sphere_full_polx",
    "sphere_xpec_polx",
    "sphere_zpmc_polx",
    "sphere_xpec_zpmc_polx",
    "sphere_full_polz",
    "sphere_xpmc_polz",
    "sphere_zpec_polz",
    "sphere_xpmc_zpec_polz",
]

NUMBER_OF_FREQUENCIES = 60


@pytest.mark.farfield
@pytest.mark.planewave
@pytest.mark.probes
@pytest.mark.parametrize("case", SPHERE_SYMMETRY_CASES)
def test_sphere_symmetry_case_produces_farfield_and_watchdog(tmp_path, case):
    """Run a symmetry case briefly and check far-field and point probe output."""
    fn = CASES_FOLDER + SYMMETRY_FOLDER + case + ".fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 30

    solver.run()

    farfield_folders = solver.getSolvedProbeFolders("farfield")
    assert len(farfield_folders) == 1
    farfield = Probe(farfield_folders[0])
    assert farfield.type == "farField"
    farfield_file = farfield.getTextFile()
    assert farfield_file is not None
    data = np.loadtxt(farfield_file, skiprows=1)
    # the solver writes one row less than the requested frequency samples
    assert data.shape[0] == NUMBER_OF_FREQUENCIES - 1
    rcs_geom = data[:, 8]
    assert np.all(np.isfinite(rcs_geom))
    assert np.all(rcs_geom > 0.0)

    watchdog_folders = solver.getSolvedProbeFolders("watchdog")
    # one output per requested direction (x, y, z)
    assert len(watchdog_folders) == 3
    for folder in watchdog_folders:
        watchdog = Probe(folder)
        assert watchdog.type == "point"
        watchdog_file = watchdog.getDatFile()
        assert watchdog_file is not None
        watchdog_data = np.loadtxt(watchdog_file, skiprows=1)
        assert watchdog_data.shape[0] > 10
        assert np.all(np.isfinite(watchdog_data))
