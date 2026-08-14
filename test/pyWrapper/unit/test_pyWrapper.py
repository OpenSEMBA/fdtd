from utils import *
import utils
from pathlib import Path

import numpy as np
import pytest

from src_pyWrapper.pyWrapper import FDTD, Probe


OUTPUTS_PATH = Path(OUTPUTS_FOLDER)


PROBE_TYPES: dict = {
    "wire_time": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "wire_probe", "expected": "wire_probe"},
        "type": {"code": "_Wx_", "expected": "wire"},
        "region": {"code": "11_12_13", "expected": [11, 12, 13]},
        "segment": {"code": "_s2", "expected": 2},
        "domain": {"code": "_tm", "expected": "time"},
        "field": "W",
        "direction": "x",
        "expected_extensions": [".dat", ".bin"],
        "expected_dat_columns": [
            "t",
            "current",
            "delta_voltage",
            "plus_voltage",
            "minus_voltage",
            "voltage_difference",
        ],
    },
    "wire_frequency": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "wire_probe", "expected": "wire_probe"},
        "type": {"code": "_Wz_", "expected": "wire"},
        "region": {"code": "21_22_23", "expected": [21, 22, 23]},
        "segment": {"code": "_s10", "expected": 10},
        "domain": {"code": "_fq", "expected": "frequency"},
        "field": "W",
        "direction": "z",
        "expected_extensions": [".dat", ".bin"],
        "expected_dat_columns": ["frequency", "magnitude", "phase"],
    },
    "bulk_current": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "bulk_probe", "expected": "bulk_probe"},
        "type": {"code": "_Jy_", "expected": "bulkCurrent"},
        "region": {"code": "1_2_3__4_5_6", "expected": ([1, 2, 3], [4, 5, 6])},
        "domain": {"code": "_tm", "expected": "time"},
        "field": "J",
        "direction": "y",
        "expected_extensions": [".dat", ".bin"],
        "expected_dat_columns": ["t", "current"],
    },
    "line_integral": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "line_probe", "expected": "line_probe"},
        "type": {"code": "_LI_", "expected": "lineIntegral"},
        "region": {"code": "4_5_6", "expected": [4, 5, 6]},
        "domain": {"code": "_tm", "expected": "time"},
        "field": "L",
        "direction": "I",
        "expected_extensions": [".dat", ".bin"],
        "expected_dat_columns": ["t", "lineIntegral"],
    },
    "point_time": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "point_probe", "expected": "point_probe"},
        "type": {"code": "_Ex_", "expected": "point"},
        "region": {"code": "7_8_9", "expected": [7, 8, 9]},
        "domain": {"code": "_tm", "expected": "time"},
        "field": "E",
        "direction": "x",
        "expected_extensions": [".dat", ".bin"],
        "expected_dat_columns": ["t", "field"],
    },
    "point_frequency": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "point_probe", "expected": "point_probe"},
        "type": {"code": "_Hz_", "expected": "point"},
        "region": {"code": "10_11_12", "expected": [10, 11, 12]},
        "domain": {"code": "_fq", "expected": "frequency"},
        "field": "H",
        "direction": "z",
        "expected_extensions": [".dat", ".bin"],
        "expected_dat_columns": ["frequency", "real", "imaginary"],
    },
    "far_field": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "far_field", "expected": "far_field"},
        "type": {"code": "_FF_", "expected": "farField"},
        "region": {"code": "1_2_3__4_5_6", "expected": ([1, 2, 3], [4, 5, 6])},
        "domain": {"code": "_fq", "expected": "frequency"},
        "expected_extensions": ["", ".bin"],
        "expected_dat_columns": [
            "frequency",
            "Theta",
            "Phi",
            "Etheta_mod",
            "Etheta_phase",
            "Ephi_mod",
            "Ephi_phase",
            "RCS_arithmetic",
            "RCS_geometric",
        ],
    },
    "movie": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "movie_probe", "expected": "movie_probe"},
        "type": {"code": "_ExC_", "expected": "movie"},
        "region": {"code": "1_2_3__4_5_6", "expected": ([1, 2, 3], [4, 5, 6])},
        "domain": {"code": "_tm", "expected": "time"},
        "expected_extensions": [".bin", ".xdmf", ".h5"],
        "extra_outputs": ["_geometry.xdmf", "_geometry.h5"],
        "expected_dat_columns": None,
    },
    "mtln_time": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "mtln_probe", "expected": "mtln_probe"},
        "type": {"code": "_V_", "expected": "mtln"},
        "region": {"code": "13_14_15", "expected": [13, 14, 15]},
        "domain": {"code": "_tm", "expected": "time"},
        "expected_extensions": [".dat"],
        "expected_dat_columns": ["t", "voltage_0"],
    },
    "mtln_frequency": {
        "case": {"code": "case_name", "expected": "case_name"},
        "name": {"code": "mtln_probe", "expected": "mtln_probe"},
        "type": {"code": "_I_", "expected": "mtln"},
        "region": {"code": "16_17_18", "expected": [16, 17, 18]},
        "domain": {"code": "_fq", "expected": "frequency"},
        "expected_extensions": [".dat"],
        "expected_dat_columns": ["frequency", "current_0"],
    },
}


@pytest.mark.probes
@pytest.mark.parametrize("probe_type", PROBE_TYPES)
def test_read_probe(tmp_path, probe_type):
    probe_case = PROBE_TYPES[probe_type]
    probe_folder = "{case}.fdtd_{name}{type}{region}{segment}{domain}/".format(
        case=probe_case["case"]["code"],
        name=probe_case["name"]["code"],
        type=probe_case["type"]["code"],
        region=probe_case["region"]["code"],
        segment=probe_case.get("segment", {}).get("code", ""),
        domain=probe_case["domain"]["code"],
    )
    probe_path = tmp_path / "results" / probe_folder
    probe_path.parent.mkdir()
    probe_path.mkdir()

    probe = Probe(probe_path)

    assert probe.case_name == probe_case["case"]["expected"]
    assert probe.name == probe_case["name"]["expected"]
    assert probe.type == probe_case["type"]["expected"]
    assert probe.domainType == probe_case["domain"]["expected"]
    if "field" in probe_case:
        assert probe.field == probe_case["field"]
        assert probe.direction == probe_case["direction"]
    if probe.type in ("bulkCurrent", "farField", "movie"):
        expected_cell_init, expected_cell_end = probe_case["region"]["expected"]
        assert np.all(probe.cell_init == expected_cell_init)
        assert np.all(probe.cell_end == expected_cell_end)
    else:
        assert np.all(probe.cell == probe_case["region"]["expected"])
    if "segment" in probe_case:
        assert probe.segment == probe_case["segment"]["expected"]


@pytest.mark.probes
@pytest.mark.parametrize("probe_type", PROBE_TYPES)
def test_getExpectedOutputs_returns_complete_set(tmp_path, probe_type):
    probe_case = PROBE_TYPES[probe_type]
    probe_folder = tmp_path / "{case}.fdtd_{name}{type}{region}{segment}{domain}".format(
        case=probe_case["case"]["code"],
        name=probe_case["name"]["code"],
        type=probe_case["type"]["code"],
        region=probe_case["region"]["code"],
        segment=probe_case.get("segment", {}).get("code", ""),
        domain=probe_case["domain"]["code"],
    )
    probe_folder.mkdir()
    expected_outputs = [
        probe_folder / f"{probe_folder.name}{extension}"
        for extension in probe_case["expected_extensions"]
    ]
    expected_outputs.extend(
        probe_folder / f"{probe_folder.name}{suffix}"
        for suffix in probe_case.get("extra_outputs", [])
    )
    for output in expected_outputs:
        output.touch()
    (probe_folder / "unrelated.dat").touch()

    probe = Probe(probe_folder)

    assert probe.getExpectedOutputs() == sorted(
        str(output.resolve()) for output in expected_outputs
    )


@pytest.mark.probes
@pytest.mark.parametrize("probe_type", PROBE_TYPES)
def test_get_expected_columns_match_schematic_text_output(tmp_path, probe_type):
    probe_case = PROBE_TYPES[probe_type]
    expected_dat_columns = probe_case["expected_dat_columns"]
    if expected_dat_columns is None:
        pytest.skip("Movie probes do not produce a text .dat output")

    probe_stem = "{case}.fdtd_{name}{type}{region}{segment}{domain}".format(
        case=probe_case["case"]["code"],
        name=probe_case["name"]["code"],
        type=probe_case["type"]["code"],
        region=probe_case["region"]["code"],
        segment=probe_case.get("segment", {}).get("code", ""),
        domain=probe_case["domain"]["code"],
    )
    probe_folder = tmp_path / probe_stem
    probe_folder.mkdir()
    extension = ".dat" if ".dat" in probe_case["expected_extensions"] else ""
    probe_path = probe_folder / f"{probe_stem}{extension}"
    probe_path.write_text(
        " ".join(expected_dat_columns) + "\n" + " ".join("0" for _ in expected_dat_columns) + "\n"
    )

    probe = Probe(probe_path)

    assert probe.getExpectedColumns() == expected_dat_columns
    assert probe_path.read_text().splitlines()[0].split() == expected_dat_columns


@pytest.mark.probes
@pytest.mark.wires
def test_read_wire_probe():
    probe = Probe(
        OUTPUTS_PATH / "fakeCurrentProbe.fdtd_mid_point_Wz_11_11_11_s2.dat"
    )

    assert probe.case_name == "fakeCurrentProbe"
    assert probe.name == "mid_point"
    assert probe.type == "wire"
    assert probe.domainType == "time"
    assert np.array_equal(probe.cell, [11, 11, 11])
    assert probe.segment == 2

    assert len(probe["time"]) == 3
    assert probe["time"].iat[0] == 0.0
    assert np.isclose(probe["time"].iat[-1], 0.59999998025528356e-10)

    assert len(probe["current"]) == 3
    assert probe["current"].iat[0] == 0.0
    assert probe["current"].iat[-1] == 0.0

@pytest.mark.probes
@pytest.mark.wires
def test_read_frequency_probe():
    probe = Probe(
        OUTPUTS_PATH / "edelcadfixZ_COR2_log__Wz_21_21_28_s10_df.dat"
    )

    assert probe.type == "wire"
    assert probe.domainType == "frequency"
    assert np.array_equal(probe.cell, [21, 21, 28])
    assert probe.segment == 10
    assert list(probe.data.columns[:3]) == ["frequency", "magnitude", "phase"]


@pytest.mark.probes
def test_read_point_probe():
    probe = Probe(
        OUTPUTS_PATH / "shieldingEffectiveness.fdtd_front_Ex_1_1_1.dat"
    )

    assert probe.case_name == "shieldingEffectiveness"
    assert probe.name == "front"
    assert probe.type == "point"
    assert probe.domainType == "time"
    assert probe.direction == "x"
    assert probe.field == "E"
    assert np.array_equal(probe.cell, [1, 1, 1])

    assert len(probe["time"]) == 5193
    assert probe["time"].iat[0] == 0.0
    assert np.isclose(probe["time"].iat[-1], 0.19997851853637005e-7)

    assert len(probe["field"]) == 5193
    assert probe["field"].iat[0] == 0.0
    assert np.isclose(probe["field"].iat[-1], 0.120145380)

    assert len(probe["incident"]) == 5193
    assert np.isclose(probe["incident"].iat[0], 0.134010895e-5)
    assert probe["incident"].iat[-1] == 0.0


@pytest.mark.probes
def test_read_point_probe_without_planewave():
    probe = Probe(OUTPUTS_PATH / "twoWires.fdtd_ProbeEnd_Ey_25_13_5.dat")

    assert probe.case_name == "twoWires"
    assert probe.name == "ProbeEnd"
    assert probe.type == "point"
    assert probe.domainType == "time"
    assert probe.direction == "y"
    assert probe.field == "E"
    assert np.array_equal(probe.cell, [25, 13, 5])
    assert "incident" not in probe.data.columns




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
    assert np.array_equal(probe.cell_init, [1, 1, 1])
    assert np.array_equal(probe.cell_end, [2, 2, 2])


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
    mtln_text = run_folder / "case.fdtd_start_voltage_wire_V_5_5_1.dat"
    wire_text = run_folder / "case.fdtd_wire_start_Wz_27_25_30_s3_tm.dat"
    unrelated_extensionless = run_folder / "case.fdtd_sample_notes"
    literal_name = run_folder / "case.fdtd_probe[0]_Ex_4_5_6.dat"
    regex_like_name = run_folder / "case.fdtd_probe0_Ex_4_5_6.dat"
    wrong_run_text = other_folder / "case.fdtd_sample_Ex_1_2_3.dat"
    for path in (
        nested_text,
        nested_binary,
        descriptor,
        flat_text,
        mtln_text,
        wire_text,
        unrelated_extensionless,
        literal_name,
        regex_like_name,
        wrong_run_text,
    ):
        path.touch()

    far_field = run_folder / "case.fdtd_farfield_FF_1_1_1__2_2_2"
    far_field.touch()
    monkeypatch.chdir(other_folder)

    expected_results = [
        ("sample", False, [flat_text, nested_text]),
        ("sample", True, [flat_text, nested_binary, nested_text]),
        ("farfield", False, [far_field]),
        ("farfield_FF_1_1_1__2_2_2", False, [far_field]),
        ("probe[0]", False, [literal_name]),
        ("sample_Ex_1_2_3_tm", False, [nested_text]),
        ("start_voltage", False, [mtln_text]),
        ("wire_start", False, [wire_text]),
    ]

    for probe_name, include_binary, expected_paths in expected_results:
        expected = sorted(str(path.resolve()) for path in expected_paths)
        assert solver.getSolvedProbeFilenames(probe_name, include_binary) == expected


@pytest.mark.spice
@pytest.mark.mtln
def test_fdtd_get_used_files():
    input_path = (
        Path(CASES_FOLDER) / "multilines_opamp" / "multilines_opamp.fdtd.json"
    )
    solver = FDTD(input_path, path_to_exe=input_path)

    assert solver.getUsedFiles() == ["spice_4port_pulse_start_75.exc", "opamp.model"]


def test_default_semba_exe_prefers_environment_override(tmp_path, monkeypatch):
    configured_exe = tmp_path / "configured" / "semba-fdtd"
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv("SEMBA_EXE", str(configured_exe))

    assert utils._default_semba_exe() == str(configured_exe)


@pytest.mark.parametrize(
    ("environment", "build_dir"),
    [
        ({}, "build-rls"),
        ({"SEMBA_FDTD_ENABLE_MPI": "ON"}, "build-rls-mpi"),
        ({"SEMBA_FDTD_ENABLE_MTLN": "OFF"}, "build-rls-nomtln"),
    ],
)
def test_default_semba_exe_selects_compatible_preset(
    tmp_path, monkeypatch, environment, build_dir
):
    executable = tmp_path / build_dir / "bin" / "semba-fdtd"
    executable.parent.mkdir(parents=True)
    executable.touch()
    legacy_executable = tmp_path / "build" / "bin" / "semba-fdtd"
    legacy_executable.parent.mkdir(parents=True)
    legacy_executable.touch()
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("SEMBA_EXE", raising=False)
    monkeypatch.delenv("SEMBA_FDTD_ENABLE_MPI", raising=False)
    monkeypatch.delenv("SEMBA_FDTD_ENABLE_MTLN", raising=False)
    for name, value in environment.items():
        monkeypatch.setenv(name, value)

    assert utils._default_semba_exe() == str(executable)
