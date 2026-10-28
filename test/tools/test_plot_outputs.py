"""Unit tests for tools/plot_outputs/plot_outputs.py on synthetic probe files."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "tools" / "plot_outputs" / "plot_outputs.py"


@pytest.fixture(scope="module")
def plot_outputs():
    spec = importlib.util.spec_from_file_location("plot_outputs", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parse_probe_reads_name_component_index_and_domain(plot_outputs):
    probe = plot_outputs.parse_probe(Path("pw-in-box.fdtd_corner_000_Ex_0_0_0_tm.dat"))
    assert (probe.name, probe.component, probe.index, probe.domain) == ("corner_000", "Ex", (0, 0, 0), "tm")


def test_parse_probe_accepts_bulk_current_range(plot_outputs):
    probe = plot_outputs.parse_probe(Path("case.fdtd_bulk_edge_x1_Jx_1_0_0__1_19_19_tm.dat"))
    assert (probe.name, probe.component, probe.index, probe.domain) == ("bulk_edge_x1", "Jx", (1, 0, 0), "tm")


@pytest.mark.parametrize(
    ("name", "content"),
    [
        ("case.fdtd_Energy.dat", "7.6E-009 3142.2\n"),
        ("case.fdtd_far_field_FF_1_1_1__18_18_18.dat", "f_at_0.2E-007 Theta Phi\n"),
    ],
)
def test_parse_probe_rejects_non_probe_outputs(plot_outputs, tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content, encoding="utf-8")
    assert plot_outputs.parse_probe(path) is None


def write_case(directory, steps_x, origin):
    grid = {
        "numberOfCells": [len(steps_x), 2, 2],
        "steps": {"x": steps_x, "y": [0.5], "z": [0.5]},
        "origin": origin,
    }
    (directory / "case.fdtd.json").write_text(json.dumps({"mesh": {"grid": grid}}), encoding="utf-8")


def test_probe_csv_places_rectilinear_grid_and_keeps_peak_per_location(plot_outputs, tmp_path):
    write_case(tmp_path, steps_x=[0.1, 0.2], origin=[1.0, 0.0, 0.0])
    (tmp_path / "case.fdtd_p_Ex_1_1_1_tm.dat").write_text("t field\n0.0 1.0\n1.0 -3.0\n", encoding="utf-8")
    (tmp_path / "case.fdtd_p_Ey_1_1_1_tm.dat").write_text("t field\n0.0 2.0\n", encoding="utf-8")
    (tmp_path / "case.fdtd_p_Ex_1_1_1_fq.dat").write_text("frequency magnitude\n1e9 50.0\n", encoding="utf-8")

    output = plot_outputs.probe_csv(tmp_path, tmp_path / "plots" / "probes.csv")

    rows = output.read_text(encoding="utf-8").splitlines()
    assert rows[0] == "x,y,z,name,peak"
    x, y, z, name, peak = rows[1].split(",")
    assert float(x) == pytest.approx(1.0 + 0.1 + 0.0)  # lower vertex of cell 1 along x
    assert (float(y), float(z), name) == (0.5, 0.5, "p")
    assert float(peak) == pytest.approx(3.0)  # frequency-domain file is not used for 3D


def test_traces_skip_phase_columns_in_frequency_domain(plot_outputs):
    probe = plot_outputs.Probe(Path("x"), "p", "Ex", (0, 0, 0), "fq")
    columns = ["frequency", "magnitude", "phase", "incident_magnitude", "incident_phase"]
    series = plot_outputs.Series(probe, columns, np.zeros((1, 5)))
    assert [name for _, name in plot_outputs.traces(series)] == ["magnitude", "incident_magnitude"]


@pytest.mark.parametrize(
    ("file_name", "expected"),
    [
        ("case.fdtd_line_voltage_LI_tm.dat", ("line_voltage", "LI", None, "tm")),
        ("case.fdtd_wire_probe_Wx_5_4_4_s2_tm.dat", ("wire_probe", "Wx", (5, 4, 4), "tm")),
        ("case.fdtd_wire_charge_Qz_3_3_3_s1_fq.dat", ("wire_charge", "Qz", (3, 3, 3), "fq")),
        ("case.fdtd_bulk_mag_Mx_1_1_1__1_8_8_tm.dat", ("bulk_mag", "Mx", (1, 1, 1), "tm")),
    ],
)
def test_parse_probe_accepts_line_wire_and_circulation_outputs(plot_outputs, file_name, expected):
    probe = plot_outputs.parse_probe(Path(file_name))
    assert (probe.name, probe.component, probe.index, probe.domain) == expected


def test_load_all_reads_mtln_wire_output_by_its_time_header(plot_outputs, tmp_path):
    (tmp_path / "case.fdtd_wire_probe_bundle_I_5_4_4.dat").write_text(
        "time conductor_1 conductor_2\n0.0 1.0 2.0\n1.0 -3.0 0.5\n", encoding="utf-8"
    )
    (tmp_path / "case.fdtd_Energy.dat").write_text("7.6E-009 3142.2\n", encoding="utf-8")

    (series,) = plot_outputs.load_all(tmp_path)

    assert (series.probe.name, series.probe.index, series.probe.domain) == ("wire_probe_bundle_I_5_4_4", None, "tm")
    assert series.columns == ["time", "conductor_1", "conductor_2"]
