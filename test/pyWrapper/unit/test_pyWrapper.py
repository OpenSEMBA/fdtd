from utils import *
import utils
from pathlib import Path

import numpy as np
import pytest

from src_pyWrapper.pyWrapper import FDTD, Probe


@pytest.mark.probes
@pytest.mark.wires
def test_read_wire_probe():
    p = Probe(
        "testData/outputs/fakeCurrentProbe.fdtd_mid_point_Wz_11_11_11_s2.dat"
    )

    assert p.case_name == "fakeCurrentProbe"
    assert p.name == "mid_point"
    assert p.type == "wire"
    assert p.domainType == "time"
    assert np.all(p.cell == np.array([11, 11, 11]))
    assert p.segment == 2

    assert len(p["time"]) == 3
    assert p["time"][0] == 0.0
    assert np.isclose(p["time"].iat[-1], 0.59999998025528356e-010)

    assert len(p["current"]) == 3
    assert p["current"][0] == 0.0
    assert p["current"].iat[-1] == -0.000000000e000


@pytest.mark.probes
@pytest.mark.wires
def test_read_probe_from_NFDE():
    p = Probe(
        "testData/outputs/fakeCurrentProbe.fdtd_mid_point_Wz_11_11_11_s2.dat"
    )

    assert p.type == "wire"


@pytest.mark.probes
@pytest.mark.wires
def test_read_frequency_probe_from_NFDE():
    p = Probe("testData/outputs/edelcadfixZ_COR2_log__Wz_21_21_28_s10_df.dat")

    assert p.type == "wire"
    assert p.domainType == "frequency"
    assert np.all(p.cell == np.array([21, 21, 28]))
    assert p.segment == 10


@pytest.mark.probes
def test_read_point_probe():
    p = Probe("testData/outputs/shieldingEffectiveness.fdtd_front_Ex_1_1_1.dat")

    assert p.case_name == "shieldingEffectiveness"
    assert p.name == "front"
    assert p.type == "point"
    assert p.domainType == "time"
    assert p.direction == "x"
    assert p.field == "E"
    assert np.all(p.cell == np.array([1, 1, 1]))

    assert len(p["time"]) == 5193
    assert p["time"][0] == 0.0
    assert np.isclose(p["time"].iat[-1], 0.19997851853637005e-007)

    assert len(p["field"]) == 5193
    assert p["field"][0] == 0.0
    assert np.isclose(p["field"].iat[-1], 0.120145380e000)

    assert len(p["incident"]) == 5193
    assert np.isclose(p["incident"][0], 0.134010895e-005)
    assert p["incident"].iat[-1] == 0.0


@pytest.mark.probes
def test_read_point_probe_without_planewave():
    p = Probe("testData/outputs/twoWires.fdtd_ProbeEnd_Ey_25_13_5.dat")

    assert p.case_name == "twoWires"
    assert p.name == "ProbeEnd"
    assert p.type == "point"
    assert p.domainType == "time"
    assert p.direction == "y"
    assert p.field == "E"
    assert np.all(p.cell == np.array([25, 13, 5]))


@pytest.mark.probes
def test_read_bulk_current_probe():
    p = Probe("testData/outputs/twoWires.fdtd_Bulk probe_Jx_15_11_13__15_13_17.dat")

    assert p.case_name == "twoWires"
    assert p.name == "Bulk probe"
    assert p.type == "bulkCurrent"
    assert p.domainType == "time"
    assert p.direction == "x"


@pytest.mark.probes
@pytest.mark.parametrize(
    ("domain_marker", "expected_domain"),
    [("_tm", "time"), ("_fq", "frequency")],
)
def test_read_point_probe_with_output_domain_marker(
    tmp_path, domain_marker, expected_domain
):
    probe_folder = tmp_path / "case.fdtd_sample_Ex_1_2_3"
    probe_folder.mkdir()
    probe_path = probe_folder / (f"case.fdtd_sample_Ex_1_2_3{domain_marker}.dat")
    probe_path.write_text("value_a value_b\n0.0 1.0\n")

    probe = Probe(probe_path)

    assert probe.name == "sample"
    assert probe.domainType == expected_domain
    assert np.all(probe.cell == np.array([1, 2, 3]))


@pytest.mark.probes
@pytest.mark.farfield
def test_read_extensionless_far_field_probe(tmp_path):
    probe_folder = tmp_path / "case.fdtd_farfield_FF_1_1_1__2_2_2"
    probe_folder.mkdir()
    probe_path = probe_folder / "case.fdtd_farfield_FF_1_1_1__2_2_2"
    probe_path.write_text(
        "frequency Theta Phi Etheta_mod Etheta_phase Ephi_mod "
        "Ephi_phase RCS_arithmetic RCS_geometric\n"
        "1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0\n"
    )

    probe = Probe(probe_path)

    assert probe.name == "farfield"
    assert probe.type == "farField"
    assert np.all(probe.cell_init == np.array([1, 1, 1]))
    assert np.all(probe.cell_end == np.array([2, 2, 2]))


@pytest.mark.probes
def test_probe_discovery_is_scoped_to_solver_folder(tmp_path, monkeypatch):
    run_folder = tmp_path / "own_run"
    other_folder = tmp_path / "other_run"
    run_folder.mkdir()
    other_folder.mkdir()

    input_path = run_folder / "case.fdtd.json"
    input_path.write_text(
        '{"probes": [{"name": "sample"}, {"name": "probe[0]"}, {"name": "farfield"}]}'
    )
    monkeypatch.chdir(tmp_path)
    solver = FDTD(Path("own_run/case.fdtd.json"), path_to_exe=input_path)

    probe_folder = run_folder / "case.fdtd_sample_Ex_1_2_3"
    probe_folder.mkdir()
    nested_text = probe_folder / "case.fdtd_sample_Ex_1_2_3_tm.dat"
    nested_binary = probe_folder / "case.fdtd_sample_Ex_1_2_3_tm.bin"
    descriptor = probe_folder / "case.fdtd_sample_Ex_1_2_3.json"
    flat_text = run_folder / "case.fdtd_sample_line_I_1_2_3.dat"
    unrelated_extensionless = run_folder / "case.fdtd_sample_notes"
    literal_name = run_folder / "case.fdtd_probe[0]_Ex_4_5_6.dat"
    regex_like_name = run_folder / "case.fdtd_probe0_Ex_4_5_6.dat"
    wrong_run_text = other_folder / "case.fdtd_sample_Ex_1_2_3.dat"
    for path in (
        nested_text,
        nested_binary,
        descriptor,
        flat_text,
        unrelated_extensionless,
        literal_name,
        regex_like_name,
        wrong_run_text,
    ):
        path.touch()

    far_field = run_folder / "case.fdtd_farfield_FF_1_1_1__2_2_2"
    far_field.touch()
    monkeypatch.chdir(other_folder)

    assert solver.getSolvedProbeFilenames("sample") == sorted(
        [
            str(flat_text.resolve()),
            str(nested_text.resolve()),
        ]
    )
    assert solver.getSolvedProbeFilenames("sample", include_binary=True) == sorted(
        [
            str(flat_text.resolve()),
            str(nested_binary.resolve()),
            str(nested_text.resolve()),
        ]
    )
    assert solver.getSolvedProbeFilenames("farfield") == [str(far_field.resolve())]
    assert solver.getSolvedProbeFilenames("farfield_FF_1_1_1__2_2_2") == [
        str(far_field.resolve())
    ]
    assert solver.getSolvedProbeFilenames("probe[0]") == [str(literal_name.resolve())]
    assert solver.getSolvedProbeFilenames("sample_Ex_1_2_3_tm") == [
        str(nested_text.resolve())
    ]


@pytest.mark.spice
@pytest.mark.mtln
def test_fdtd_get_used_files():
    solver = FDTD(
        "testData/cases/multilines_opamp/multilines_opamp.fdtd.json",
        path_to_exe="build/bin/semba-fdtd",
    )

    used_files = solver.getUsedFiles()

    assert len(used_files) == 2
    assert used_files[0] == 'spice_4port_pulse_start_75.exc'
    assert used_files[1] == 'opamp.model'


def test_default_semba_exe_prefers_environment_override(tmp_path, monkeypatch):
    configured_exe = tmp_path / 'configured' / 'semba-fdtd'
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('SEMBA_EXE', str(configured_exe))

    assert utils._default_semba_exe() == str(configured_exe)


@pytest.mark.parametrize(
    ('environment', 'build_dir'),
    [
        ({}, 'build-rls'),
        ({'SEMBA_FDTD_ENABLE_MPI': 'ON'}, 'build-rls-mpi'),
        ({'SEMBA_FDTD_ENABLE_MTLN': 'OFF'}, 'build-rls-nomtln'),
    ],
)
def test_default_semba_exe_selects_compatible_preset(
    tmp_path, monkeypatch, environment, build_dir
):
    executable = tmp_path / build_dir / 'bin' / 'semba-fdtd'
    executable.parent.mkdir(parents=True)
    executable.touch()
    legacy_executable = tmp_path / 'build' / 'bin' / 'semba-fdtd'
    legacy_executable.parent.mkdir(parents=True)
    legacy_executable.touch()
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv('SEMBA_EXE', raising=False)
    monkeypatch.delenv('SEMBA_FDTD_ENABLE_MPI', raising=False)
    monkeypatch.delenv('SEMBA_FDTD_ENABLE_MTLN', raising=False)
    for name, value in environment.items():
        monkeypatch.setenv(name, value)

    assert utils._default_semba_exe() == str(executable)

