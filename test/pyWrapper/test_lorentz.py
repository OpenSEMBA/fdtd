"""End-to-end checks of the JSON `lorentz` material against stored references.

Runs the PEC-backed 2 m Lorentz slabs of ``testData/cases/lorentz`` and compares
their S11 against the committed dgtd and scikit-rf vectors in
``testData/cases/lorentz/references``.

An all-vacuum run with the same mesh, source and probe supplies the incident
field, so S11 = DFT(E_total - E_incident) / DFT(E_incident) at the probe plane.
The two cases cover the underdamped Lorentz resonance (``lorentz_res2``) and the
cold-plasma/Drude limit (``lorentz_plasma``). Both slabs are PEC-backed, so S21
is identically zero and S11 fully characterizes them.
"""

import numpy as np
import pytest

from test.utils.utils import CASES_FOLDER, FDTD, Probe, SEMBA_EXE

CASE_FOLDER = CASES_FOLDER + "lorentz/"
REFERENCE_FOLDER = CASE_FOLDER + "references/"
F_MIN = 1.0e6
F_MAX = 400.0e6
MAX_ERROR_DB = 0.5
MAX_ERROR_DEG = 12.0


def _run_case(case_name, folder):
    solver = FDTD(
        input_filename=CASE_FOLDER + case_name + ".fdtd.json",
        path_to_exe=SEMBA_EXE,
        run_in_folder=folder,
    )
    solver.run()
    assert solver.hasFinishedSuccessfully()
    return solver


def _reflection_probe(solver):
    folders = solver.getSolvedProbeFolders("reflection")
    assert len(folders) == 1, f"Expected one reflection probe, found {folders}"
    return Probe(folders[0])


def _dft(time, values, frequencies):
    return np.exp(-2j * np.pi * np.outer(frequencies, time)) @ values


def _read_reference(case_name):
    return np.loadtxt(
        REFERENCE_FOLDER + case_name + "_s11.csv", delimiter=",", skiprows=1
    )


@pytest.fixture(scope="module")
def incident(tmp_path_factory):
    """Incident field at the probe plane from an all-vacuum run."""
    solver = _run_case("vacuum_reference", tmp_path_factory.mktemp("lorentz-vacuum"))
    probe = _reflection_probe(solver)
    return probe["time"].to_numpy(), probe["field"].to_numpy()


@pytest.mark.dielectric
@pytest.mark.probes
@pytest.mark.parametrize("case_name", ["lorentz_res2", "lorentz_plasma"])
def test_lorentz_slab_s11(case_name, incident, tmp_path):
    """Verify a Lorentz slab S11 matches the dgtd and scikit-rf vectors."""
    solver = _run_case(case_name, tmp_path)
    probe = _reflection_probe(solver)
    time = probe["time"].to_numpy()
    field = probe["field"].to_numpy()

    reference = _read_reference(case_name)
    frequencies = reference[:, 0]
    assert np.all((frequencies >= F_MIN) & (frequencies <= F_MAX))

    time_incident, field_incident = incident
    assert len(time) == len(time_incident)
    assert np.allclose(time, time_incident)

    s11 = _dft(time, field - field_incident, frequencies) / _dft(
        time, field_incident, frequencies
    )
    magnitude = 20.0 * np.log10(np.abs(s11))
    phase = np.rad2deg(np.angle(s11))

    for db, deg in ((reference[:, 1], reference[:, 2]), (reference[:, 3], reference[:, 4])):
        assert np.all(np.abs(magnitude - db) < MAX_ERROR_DB)
        phase_error = (phase - deg + 180.0) % 360.0 - 180.0
        assert np.all(np.abs(phase_error) < MAX_ERROR_DEG)
