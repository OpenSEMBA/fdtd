from test.utils.utils import *
from typing import Dict
from pathlib import Path
import os
from sys import platform
from scipy import signal
import matplotlib.pyplot as plt


def _get_solved_probe_folder(
    solver,
    probe_name,
    *,
    filename=None,
    contains=None,
    suffix=None,
    include_binary=False,
) -> str:
    probe_files = solver.getSolvedProbeFolders(probe_name)
    if filename is not None:
        probe_files = [path for path in probe_files if Path(path).stem == Path(filename).stem]
    if contains is not None:
        probe_files = [path for path in probe_files if contains in Path(path).name]
    if suffix is not None:
        artifact_getter = {
            ".bin": Probe.getBinFile,
            ".dat": Probe.getDatFile,
            ".h5": Probe.getH5File,
            ".xdmf": Probe.getXDMFFile,
        }[suffix]
        artifacts = [artifact_getter(Probe(path)) for path in probe_files]
        artifacts = [path for path in artifacts if path is not None]
        assert len(artifacts) == 1, (
            f"Expected one {suffix} artifact for probe {probe_name!r}, found {artifacts}"
        )
        return artifacts[0]
    assert len(probe_files) == 1, (
        f"Expected one artifact for probe {probe_name!r}, found {probe_files}"
    )
    return probe_files[0]


def _get_source_magnitude_file(solver, source_index=0) -> str:
    return solver.resolveInputPath(solver["sources"][source_index]["magnitudeFile"])


def is_debugging() -> bool:
    """Return whether opt-in debug diagnostics should be generated."""
    return os.getenv("SEMBA_FDTD_GENERATE_DEBUG_DATA", "").lower() in {
        "1",
        "true",
        "on",
        "yes",
    }


# compiled without mtln uses classic wires
# compiled with mtln, wire is treated as an unshielded multiwire
@pytest.mark.skip(reason='newOutput module has not solved this problem yet')
def test_lineIntegralProbe(tmp_path):
    """Verify a line-integral probe matches its reference waveform."""
    def generate_debug_data():
        _, axes = plt.subplots(2, sharey=True, figsize=(10, 6))
        axes[0].plot(expected_time, expected_value, label="Expected")
        axes[0].plot(solved_time, solved_value, "--", label="li_probe")
        axes[0].set_title("Line integral probe")
        axes[0].set_xlabel("Time [s]")
        axes[0].set_ylabel("Line integral")
        axes[0].legend()
        axes[0].grid()
        axes[1].plot(expected_time[window], expected_value[window], label="Expected")
        axes[1].plot(solved_time[window], solved_value[window], "--", label="li_probe")
        axes[1].set_title("Line integral probe near peak")
        axes[1].set_xlabel("Time [s]")
        axes[1].set_ylabel("Line integral")
        axes[1].legend()
        axes[1].grid()
        plt.tight_layout()
        plt.savefig(tmp_path / "lineIntegralProbe_comparison.png")
        plt.close()

    fn = CASES_FOLDER + "lineIntegralProbe/lineIntegralProbe_plates.fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["materials"][0] = createWire(id=1, r=0.001)
    solver.run()

    pf = "lineIntegralProbe_plates.fdtd_vprobe_LI_20_20_10.dat"
    li_probe = Probe(_get_solved_probe_folder(solver, "vprobe"))
    expected = probe_from_fixture(tmp_path, pf)

    solved_time = li_probe["time"].to_numpy()
    solved_value = li_probe["lineIntegral"].to_numpy()
    expected_time = expected["time"].to_numpy()
    expected_value = expected["lineIntegral"].to_numpy()
    peak = np.argmax(np.abs(expected_value))
    window = slice(max(0, peak - 100), min(len(expected_value), peak + 101))

    if is_debugging():
        generate_debug_data()

    assert np.allclose(
        solved_value,
        expected_value,
        rtol=0.01,
        atol=0.01,
    )


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.probes
def test_shieldedPair(tmp_path):
    """Verify shielded-pair voltage and current probes match reference signals."""
    fn = CASES_FOLDER + "shieldedPair/shieldedPair.fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    probe_files = [
        "shieldedPair.fdtd_wire_start_line_out_V_75_74_74.dat",
        "shieldedPair.fdtd_wire_start_line_out_I_75_74_74.dat",
        "shieldedPair.fdtd_wire_end_line_out_I_75_71_74.dat",
        "shieldedPair.fdtd_wire_end_line_out_V_75_71_74.dat",
    ]

    p_expected = []
    for pf in probe_files:
        p_expected.append(probe_from_fixture(tmp_path, pf))

    p_solved = []
    for pf in probe_files:
        probe_name = "wire_start" if "_wire_start_" in pf else "wire_end"
        p_solved.append(
            Probe(_get_solved_probe_folder(solver, probe_name, filename=pf))
        )

    for i in [0, 3]:
        solved = np.interp(
            p_expected[i]["time"].to_numpy(),
            p_solved[i]["time"].to_numpy(),
            p_solved[i]["voltage_0"].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[i]["voltage_0"])[0, 1] > 0.999
        solved = np.interp(
            p_expected[i]["time"].to_numpy(),
            p_solved[i]["time"].to_numpy(),
            p_solved[i]["voltage_1"].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[i]["voltage_1"])[0, 1] > 0.999
        solved = np.interp(
            p_expected[i]["time"].to_numpy(),
            p_solved[i]["time"].to_numpy(),
            p_solved[i]["voltage_2"].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[i]["voltage_2"])[0, 1] > 0.999
    for i in [1, 2]:
        solved = np.interp(
            p_expected[i]["time"].to_numpy(),
            p_solved[i]["time"].to_numpy(),
            p_solved[i]["current_0"].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[i]["current_0"])[0, 1] > 0.999
        solved = np.interp(
            p_expected[i]["time"].to_numpy(),
            p_solved[i]["time"].to_numpy(),
            p_solved[i]["current_1"].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[i]["current_1"])[0, 1] > 0.999
        solved = np.interp(
            p_expected[i]["time"].to_numpy(),
            p_solved[i]["time"].to_numpy(),
            p_solved[i]["current_2"].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[i]["current_2"])[0, 1] > 0.999


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.dielectric
@pytest.mark.probes
def test_coated_antenna(tmp_path):
    """Test for a coated antenna with MTLN wires reproducing Fig. 2 in:
    A. Rubio Bretones, R. Gomez Martin, A. Salinas and I. Sanchez,
    "Time domain analysis of dielectric coated wire scatterers and antennas,"
    Proceedings of MELECON '94. Mediterranean Electrotechnical Conference,
    Antalya, Turkey, 1994, pp. 1174-1176 vol.3, doi: 10.1109/MELCON.1994.380859.
    """
    fn = CASES_FOLDER + "coated_antenna/coated_antenna.fdtd.json"

    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    probe_files = solver.getSolvedProbeFolders("mid_point")
    assert {os.path.basename(probe_file) for probe_file in probe_files} == {
        "coated_antenna.fdtd_mid_point_half_1_I_11_11_12",
        "coated_antenna.fdtd_mid_point_half_2_I_11_11_12",
    }

    p_expected = probe_from_fixture(
        tmp_path, "coated_antenna.fdtd_mid_point_half_1_I_11_11_12.dat"
    )

    p_solved = Probe(
        _get_solved_probe_folder(
            solver,
            "mid_point",
            filename="coated_antenna.fdtd_mid_point_half_1_I_11_11_12.dat",
        )
    )

    solved = np.interp(
        p_expected["time"].to_numpy(),
        p_solved["time"].to_numpy(),
        p_solved["current_0"].to_numpy(),
    )
    assert np.corrcoef(solved, p_expected["current_0"])[0, 1] > 0.999


# compiled without mtln uses classic wires
# compiled with mtln, wire is treated as an unshielded multiwire
@pytest.mark.wires
@pytest.mark.probes
def test_holland(tmp_path):
    """Verify the Holland wire current agrees with the reference solution."""
    fn = CASES_FOLDER + "holland/holland1981.fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    p = Probe(_get_solved_probe_folder(solver, "mid_point"))

    expected_f = json.load(
        open(OUTPUTS_FOLDER + "holland1981_mid_point_expected_current.json")
    )
    expected_t, expected_i = np.array([]), np.array([])
    for data in expected_f["datasetColl"][0]["data"]:
        expected_t = np.append(expected_t, float(data["value"][0]))
        expected_i = np.append(expected_i, float(data["value"][1]))

    expected_i_interp = np.interp(p["time"] - 3.05 * 1e-9, expected_t, expected_i)
    assert np.allclose(expected_i_interp, p["current"], rtol=1e-4, atol=5e-5)


@pytest.mark.wires
@pytest.mark.termination
@pytest.mark.probes
def test_holland_short_terminals_match_open_terminals(tmp_path):
    """Verify short terminal definitions preserve this Holland-case response."""
    fn = CASES_FOLDER + "holland/holland1981.fdtd.json"
    number_of_steps = 1000

    folder_open = os.path.join(tmp_path, "open_terminals")
    os.makedirs(folder_open)
    solver_open = FDTD(
        input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=folder_open
    )
    solver_open["general"]["numberOfSteps"] = number_of_steps
    solver_open.run()

    probe_open = Probe(_get_solved_probe_folder(solver_open, "mid_point"))

    folder_short = os.path.join(tmp_path, "short_terminals")
    os.makedirs(folder_short)
    solver_short = FDTD(
        input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=folder_short
    )
    solver_short["general"]["numberOfSteps"] = number_of_steps
    solver_short["materials"][1]["terminations"][0]["type"] = "short"
    solver_short.run()

    probe_short = Probe(_get_solved_probe_folder(solver_short, "mid_point"))

    assert np.allclose(
        probe_open["time"].to_numpy(),
        probe_short["time"].to_numpy(),
        rtol=0.0,
        atol=0.0,
    )
    assert np.allclose(
        probe_open["current"].to_numpy(),
        probe_short["current"].to_numpy(),
        rtol=1e-4,
        atol=5e-5,
    )


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.probes
def test_unshielded_multiwires(tmp_path):
    """Verify unshielded-multiwire current outputs match reference signals."""
    fn = CASES_FOLDER + "unshielded_multiwires/unshielded_multiwires_berenger.fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    p_solved = Probe(_get_solved_probe_folder(solver, "mid_point", contains="_I_"))
    p_expected = probe_from_fixture(
        tmp_path,
        "unshielded_multiwires_berenger.fdtd_mid_point_unshielded_two_wire_I_2_11_14.dat",
    )

    for current in ("current_0", "current_1"):
        solved = np.interp(
            p_expected["time"].to_numpy(),
            p_solved["time"].to_numpy(),
            p_solved[current].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[current])[0, 1] > 0.999


@pytest.mark.wires
@pytest.mark.probes
def test_towelHanger(tmp_path):
    """Verify towel-hanger wire currents match the stored reference probes."""
    fn = CASES_FOLDER + "towelHanger/towelHanger.fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p_solved = [
        Probe(_get_solved_probe_folder(solver, "wire_start")),
        Probe(_get_solved_probe_folder(solver, "wire_mid")),
        Probe(_get_solved_probe_folder(solver, "wire_end")),
    ]

    p_expected = [
        probe_from_fixture(tmp_path, "towelHanger.fdtd_wire_start_Wz_27_25_30_s1.dat"),
        probe_from_fixture(tmp_path, "towelHanger.fdtd_wire_mid_Wx_35_25_32_s5.dat"),
        probe_from_fixture(tmp_path, "towelHanger.fdtd_wire_end_Wz_43_25_30_s4.dat"),
    ]

    for i in range(3):
        solved = np.interp(
            p_expected[i]["time"].to_numpy(),
            p_solved[i]["time"].to_numpy(),
            p_solved[i]["current_0"].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[i]["current_0"])[0, 1] > 0.999


@pytest.mark.wires
@pytest.mark.termination
@pytest.mark.probes
def test_towel_rack_with_and_without_shorting_plane(tmp_path):
    """Verify a shorting plane leaves low-frequency input impedance unchanged."""
    def generate_debug_data():
        plt.figure()
        plt.semilogx(freqs, 20 * np.log10(np.abs(Z_in_w)), label="With shorting plane")
        plt.semilogx(freqs, 20 * np.log10(np.abs(Z_in_wo)), label="Without shorting plane")
        plt.xlabel("Frequency [Hz]")
        plt.ylabel("|Z(j2πf)| [dB]")
        plt.grid(which="both")
        plt.legend()
        plt.tight_layout()
        plt.savefig(tmp_path / "input_impedance_comparison.png")
        plt.close()

    fn = (
        CASES_FOLDER
        + "towel_rack_with_shorting_plane/towel_rack_with_shorting_plane.fdtd.json"
    )

    # --- excitation ---
    dt = 1e-12
    w0 = 0.1e-8
    t0 = 10 * w0
    t = np.arange(0, t0 + 20 * w0, dt)
    excitation = np.exp(-np.power(t - t0, 2) / w0**2)
    exc_file = os.path.join(tmp_path, "gauss.exc")
    np.savetxt(exc_file, np.column_stack([t, excitation]))

    freqs = np.geomspace(1e3, 1e9, 61)

    # --- with shorting plane ---
    folder_with = os.path.join(tmp_path, "with_shorting_plane")
    os.makedirs(folder_with)
    solver_w = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=folder_with)
    solver_w.run()

    probe_w = Probe(_get_solved_probe_folder(solver_w, "Wire probe"))
    time_I_w = probe_w["time"].to_numpy()
    current_w = probe_w["current_0"].to_numpy()
    V_interp_w = np.interp(time_I_w, t, excitation)
    Z_in_w = dtft(V_interp_w, time_I_w, freqs) / dtft(current_w, time_I_w, freqs)

    # --- without shorting plane ---
    folder_without = os.path.join(tmp_path, "without_shorting_plane")
    os.makedirs(folder_without)
    solver_wo = FDTD(
        input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=folder_without
    )
    solver_wo["materialAssociations"][0]["elementIds"] = [1]
    solver_wo.run()

    probe_wo = Probe(_get_solved_probe_folder(solver_wo, "Wire probe"))
    time_I_wo = probe_wo["time"].to_numpy()
    current_wo = probe_wo["current_0"].to_numpy()
    V_interp_wo = np.interp(time_I_wo, t, excitation)
    Z_in_wo = dtft(V_interp_wo, time_I_wo, freqs) / dtft(current_wo, time_I_wo, freqs)

    if is_debugging():
        generate_debug_data()

    # Expect the shorting plane not changing the impedance at low frequencies.
    assert np.allclose(
        np.abs(Z_in_w[freqs < 1e6]), np.abs(Z_in_wo[freqs < 1e6]), rtol=0.1
    )


@pytest.mark.hdf
@pytest.mark.farfield
@pytest.mark.movie
@pytest.mark.probes
def test_sphere(tmp_path):
    """Verify far-field and movie probes produce the expected HDF artifacts."""
    fn = CASES_FOLDER + "sphere/sphere.fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 20
    solver["probes"][0]["domain"]["initialFrequency"] = 1e8
    solver["probes"][0]["domain"]["finalFrequency"] = 1e9

    solver.run()

    p = Probe(_get_solved_probe_folder(solver, "farfield"))
    assert p.type == "farField"

    electric_field_movie_folders = solver.getSolvedProbeFolders("electric_field_movie")
    assert len(electric_field_movie_folders) == 1
    p = Probe(electric_field_movie_folders[0])
    assert p.getH5File() is not None
    assert p.getXDMFFile() is not None
    assert p.type == "movie"


@pytest.mark.hdf
@pytest.mark.planewave
@pytest.mark.movie
def test_movie_in_planewave_in_box(tmp_path):
    """Verify time-domain movie HDF and XDMF metadata and field data."""
    import h5py
    import xml.etree.ElementTree as ET

    fn = CASES_FOLDER + "planewave/pw-in-box-with-movie.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    h5file = _get_solved_probe_folder(solver, "electric_field_movie", suffix=".h5")
    xdmffile = _get_solved_probe_folder(
        solver, "electric_field_movie", suffix=".xdmf"
    )
    geometry_xdmffile = str(
        Path(xdmffile).with_name(Path(xdmffile).stem + "_geometry.xdmf")
    )
    geometry_h5file = str(Path(geometry_xdmffile).with_suffix(".h5"))

    root = ET.parse(xdmffile).getroot()
    h5_name = os.path.basename(h5file)
    temporal = root.find("./Domain/Grid")
    assert temporal.attrib == {
        "Name": "Time Series",
        "GridType": "Collection",
        "CollectionType": "Temporal",
    }
    steps = temporal.findall("./Grid")
    time_ds = np.array([float(step.find("./Time").attrib["Value"]) for step in steps])
    assert len(time_ds) > 1
    assert np.all(np.diff(time_ds) > 0.0)

    with h5py.File(h5file, "r") as f:
        required_datasets = {
            "series/values",
            "grids/g0001/axis_x",
            "grids/g0001/axis_y",
            "grids/g0001/axis_z",
        }
        assert all(dataset in f for dataset in required_datasets)
        np.testing.assert_allclose(f["series/values"][()], time_ds)

        attributes = {
            attribute.attrib["Name"]: attribute
            for attribute in steps[0].findall("./Grid/Attribute")
        }
        assert set(attributes) == {
            "ElectricFieldX",
            "ElectricFieldY",
            "ElectricFieldZ",
            "tagnumber_x",
            "tagnumber_y",
            "tagnumber_z",
            "mediatype_x",
            "mediatype_y",
            "mediatype_z",
        }

        def dataset_for(attribute):
            source = attribute.find("./DataItem").findall("./DataItem")[1]
            source_name, dataset = source.text.strip().split(":/", 1)
            assert source_name == h5_name
            assert dataset in f
            return dataset

        for component in ("ElectricFieldX", "ElectricFieldY", "ElectricFieldZ"):
            field_ds = f[dataset_for(attributes[component])][()]
            assert field_ds.ndim == 4
            assert field_ds.shape[0] == len(time_ds)
            assert np.all(np.isfinite(field_ds))

        electric_field = f[dataset_for(attributes["ElectricFieldY"])][()]
        assert electric_field.shape[1:] == (10, 30, 30)
        assert np.max(np.abs(electric_field[1000])) > 1e-2

        # The dielectric slows the x-propagating pulse so it remains in the
        # movie volume through timestep 1,000.
        early_profile = np.mean(np.abs(electric_field[400]), axis=(0, 1))
        late_profile = np.mean(np.abs(electric_field[1000]), axis=(0, 1))
        assert np.argmax(late_profile[3:]) > np.argmax(early_profile[3:])

        steps = temporal.findall("./Grid")
        assert len(steps) == len(time_ds)

        for index, step in enumerate(steps):
            assert step.attrib["Name"] == f"Step {index + 1}"
            assert np.isclose(
                float(step.find("./Time").attrib["Value"]), time_ds[index]
            )
            grid = step.find("./Grid")
            assert grid.attrib == {"Name": "movieProbe", "GridType": "Uniform"}

            for attribute in grid.findall("./Attribute"):
                hyperslab = attribute.find("./DataItem")
                if hyperslab.attrib.get("ItemType") is None:
                    source_name, dataset = hyperslab.text.strip().split(":/", 1)
                    assert source_name == h5_name
                    assert dataset in f
                    assert f[dataset].shape == (10, 30, 30)
                    continue
                selection, source = hyperslab.findall("./DataItem")
                assert hyperslab.attrib["ItemType"] == "HyperSlab"
                dataset = dataset_for(attribute)
                assert (
                    tuple(
                        int(value) for value in hyperslab.attrib["Dimensions"].split()
                    )
                    == f[dataset].shape[1:]
                )
                assert selection.attrib == {"Dimensions": "3 4", "Format": "XML"}
                offset, stride, count = [
                    tuple(int(value) for value in line.split())
                    for line in selection.text.splitlines()
                    if line.strip()
                ]
                assert offset == (index, 0, 0, 0)
                assert stride == (1, 1, 1, 1)
                assert count == (1, *f[dataset].shape[1:])
                assert (
                    tuple(int(value) for value in source.attrib["Dimensions"].split())
                    == f[dataset].shape
                )

    geometry_root = ET.parse(geometry_xdmffile).getroot()
    point_data_items = geometry_root.findall(".//Geometry/DataItem")
    assert point_data_items
    geometry_points = []
    with h5py.File(geometry_h5file, "r") as f:
        for point_data_item in point_data_items:
            source_name, dataset = point_data_item.text.strip().split(":/", 1)
            assert source_name == Path(geometry_h5file).name
            points = f[dataset][()]
            if points.shape[0] == 3:
                points = points.T
            assert points.ndim == 2
            assert points.shape[1] == 3
            geometry_points.append(points)

    points = np.concatenate(geometry_points)
    assert np.all(points >= -1e-12)
    np.testing.assert_allclose(points / 0.01, np.rint(points / 0.01))
    assert np.max(points[:, 0]) <= 0.30 + 1e-12
    assert np.max(points[:, 1]) <= 0.30 + 1e-12
    assert np.max(points[:, 2]) <= 0.10 + 1e-12


@pytest.mark.hdf
def test_frequency_slice_in_planewave_in_box(tmp_path):
    """Verify frequency-slice HDF and XDMF metadata and field data."""
    import h5py
    import xml.etree.ElementTree as ET

    fn = CASES_FOLDER + "planewave/pw-in-box-with-frequency-slice.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    h5file = _get_solved_probe_folder(
        solver, "electric_field_frequency_slice", suffix=".h5"
    )
    xdmffile = _get_solved_probe_folder(
        solver, "electric_field_frequency_slice", suffix=".xdmf"
    )

    with h5py.File(h5file, "r") as f:
        root = ET.parse(xdmffile).getroot()
        series = root.find("./Domain/Grid")
        assert series.attrib == {
            "Name": "Parameter Series",
            "GridType": "Collection",
            "CollectionType": "Spatial",
        }

        steps = series.findall("./Grid")
        assert len(steps) == 4
        frequencies = [
            float(step.find("./Information").attrib["Value"]) for step in steps
        ]
        np.testing.assert_allclose(
            frequencies, [5e8, 5e8 + 1e9 / 3, 5e8 + 2e9 / 3, 1.5e9]
        )

        for index, step in enumerate(steps):
            assert step.find("./Information").attrib["Name"] == "Frequency"
            grid = step.find("./Grid")
            topology = grid.find("./Topology")
            assert topology.attrib == {
                "TopologyType": "Polyvertex",
                "NumberOfElements": "120",
                "NodesPerElement": "1",
            }
            geometry = grid.find("./Geometry/DataItem")
            assert geometry.attrib["Dimensions"] == "120 3"
            _, geometry_dataset = geometry.text.strip().split(":/", 1)
            assert f[geometry_dataset].shape == (120, 3)

            attributes = grid.findall("./Attribute")
            assert {attribute.attrib["Name"] for attribute in attributes} == {
                "xMagnitude",
                "yMagnitude",
                "zMagnitude",
                "xPhase",
                "yPhase",
                "zPhase",
                "tagnumber",
                "mediatype",
            }
            for attribute in attributes:
                hyperslab = attribute.find("./DataItem")
                if hyperslab.attrib.get("ItemType") is None:
                    _, dataset = hyperslab.text.strip().split(":/", 1)
                    assert f[dataset].shape == (120,)
                    assert np.all(np.isfinite(f[dataset][()]))
                    continue
                selection, source = hyperslab.findall("./DataItem")
                offset, stride, count = [
                    tuple(int(value) for value in line.split())
                    for line in selection.text.splitlines()
                    if line.strip()
                ]
                assert offset == (index, 0)
                assert stride == (1, 1)
                assert count == (1, 120)
                _, dataset = source.text.strip().split(":/", 1)
                assert (
                    tuple(int(value) for value in source.attrib["Dimensions"].split())
                    == f[dataset].shape
                )
                assert np.all(np.isfinite(f[dataset][()]))

        assert np.max(np.abs(f["attributes/a0001/values"][()])) > 0.0


@pytest.mark.hdf
def test_central_dipole_frequency_slice(tmp_path):
    """Verify dipole frequency slices exhibit the expected equatorial field."""
    import h5py
    import xml.etree.ElementTree as ET

    fn = CASES_FOLDER + "antenna_frequency/central_dipole_frequency_slice.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    h5file = _get_solved_probe_folder(solver, "dipole_frequency_slice", suffix=".h5")
    xdmffile = _get_solved_probe_folder(
        solver, "dipole_frequency_slice", suffix=".xdmf"
    )

    with h5py.File(h5file, "r") as f:
        root = ET.parse(xdmffile).getroot()
        steps = root.findall("./Domain/Grid/Grid")
        frequencies = [
            float(step.find("./Information").attrib["Value"]) for step in steps
        ]
        np.testing.assert_allclose(frequencies, [7.5e8, 1e9, 1.25e9])

        def component(step, name):
            attribute = next(
                item
                for item in step.findall("./Grid/Attribute")
                if item.attrib["Name"] == name
            )
            source = attribute.findall("./DataItem")[0].findall("./DataItem")[1]
            _, dataset = source.text.strip().split(":/", 1)
            return f[dataset][1]

        electric_z = component(steps[1], "zMagnitude") * np.exp(
            1j * component(steps[1], "zPhase")
        )
        electric_z = electric_z.reshape((48, 48, 48), order="F")
        amplitude = np.abs(electric_z)
        assert electric_z.shape == (48, 48, 48)
        assert np.max(amplitude) > 1e-9

        # The z-oriented dipole has an azimuthally symmetric equatorial field
        # whose magnitude decays and phase changes with radial distance.
        equatorial_near = electric_z[24, 24, 32]
        equatorial_far = electric_z[24, 24, 36]
        equatorial_far_y = electric_z[24, 36, 24]
        assert np.isclose(np.abs(equatorial_far), np.abs(equatorial_far_y), rtol=1e-3)
        assert np.abs(equatorial_far) < np.abs(equatorial_near)
        assert np.abs(np.angle(equatorial_far / equatorial_near)) > 0.1


@pytest.mark.planewave
@pytest.mark.probes
def test_planewave_in_box(tmp_path):
    """Verify a plane wave is incident only within the simulation box."""
    fn = CASES_FOLDER + "planewave/pw-in-box.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    before = Probe(_get_solved_probe_folder(solver, "before"))
    inbox = Probe(_get_solved_probe_folder(solver, "inbox"))
    after = Probe(_get_solved_probe_folder(solver, "after"))

    assert (
        np.corrcoef(inbox.data["field"].to_numpy(), inbox.data["incident"].to_numpy())[
            0, 1
        ]
        > 0.999
    )
    zeros = np.zeros_like(before.data["field"])
    assert np.allclose(before.data["field"].to_numpy(), zeros, atol=5e-4)
    assert np.allclose(after.data["field"].to_numpy(), zeros, atol=5e-4)


@pytest.mark.planewave
@pytest.mark.probes
def test_planewave_with_periodic_boundaries(tmp_path):
    """Verify a plane wave remains confined with periodic boundaries."""
    fn = CASES_FOLDER + "planewave/pw-with-periodic.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    before = Probe(_get_solved_probe_folder(solver, "before"))
    inbox = Probe(_get_solved_probe_folder(solver, "inbox"))
    after = Probe(_get_solved_probe_folder(solver, "after"))

    assert (
        np.corrcoef(inbox.data["field"].to_numpy(), inbox.data["incident"].to_numpy())[
            0, 1
        ]
        > 0.999
    )
    zeros = np.zeros_like(before.data["field"])
    assert np.allclose(before.data["field"].to_numpy(), zeros, atol=1.5e-3)
    assert np.allclose(after.data["field"].to_numpy(), zeros, atol=1.5e-3)


@pytest.mark.sgbc
@pytest.mark.probes
def test_sgbc_shielding_effectiveness(tmp_path):
    """Verify SGBC transmission agrees with the analytical slab response."""
    def generate_debug_data():
        plt.figure()
        plt.plot(f, 20 * np.log10(np.abs(fdtd_s21)), ".-", label="FDTD S21")
        plt.plot(f, 20 * np.log10(np.abs(slab.s[:, 0, 1])), ".-", label="Analytical S21")
        plt.grid(which="both")
        plt.xlim(f[0], f[-1])
        plt.xscale("log")
        plt.legend()
        plt.savefig(tmp_path / "sgbc_s21_comparison.png")
        plt.close()

    fn = CASES_FOLDER + "sgbcShieldingEffectiveness/shieldingEffectiveness.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    # FDTD results
    back = Probe(_get_solved_probe_folder(solver, "back"))

    t = back.data["time"]
    dt = t[1] - t[0]
    fq = fftfreq(len(t)) / dt
    INC = fft(back.data["incident"])
    BACK = fft(back.data["field"])
    S21 = BACK / INC

    fmin = 8e6
    fmax = 1e9
    idx_min = (np.abs(fq - fmin)).argmin()
    idx_max = (np.abs(fq - fmax)).argmin()
    f = fq[idx_min:idx_max]
    fdtd_s21 = S21[idx_min:idx_max]

    # Analytical results
    from skrf.media import Freespace
    from skrf.frequency import Frequency
    import scipy.constants

    freq = Frequency.from_f(f, unit="Hz")
    air = Freespace(freq)

    sigma = 100
    width = 10e-3
    mat_ep_r = 1 + sigma / (1j * freq.w * scipy.constants.epsilon_0)
    conductive_material = Freespace(freq, ep_r=mat_ep_r)

    slab = air.thru() ** conductive_material.line(width, unit="m") ** air.thru()

    if is_debugging():
        generate_debug_data()

    fdtd_s21_db = 20 * np.log10(np.abs(fdtd_s21))
    anal_s21_db = 20 * np.log10(np.abs(slab.s[:, 0, 1]))

    assert np.allclose(fdtd_s21_db, anal_s21_db, rtol=0.05)


@pytest.mark.sgbc
@pytest.mark.probes
def test_current_orientation(tmp_path):
    """Verify bulk-current sign follows source rather than mesh orientation."""
    fn = CASES_FOLDER + "current_orientation/currentOrientation.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver["mesh"]["elements"][12]["coordinateIds"] = [1, 2]
    solver["sources"][0]["elementIds"] = [1]
    solver.cleanUp()
    solver.run()
    i = Probe(_get_solved_probe_folder(solver, "Bulk_probe")).data["current"]
    assert np.all(i >= 0)

    solver["mesh"]["elements"][12]["coordinateIds"] = [2, 1]
    solver["sources"][0]["elementIds"] = [1]
    solver.cleanUp()
    solver.run()
    i = Probe(_get_solved_probe_folder(solver, "Bulk_probe")).data["current"]
    assert np.all(i >= 0)

    solver["mesh"]["elements"][12]["coordinateIds"] = [1, 2]
    solver["sources"][0]["elementIds"] = [3]
    solver.cleanUp()
    solver.run()
    i = Probe(_get_solved_probe_folder(solver, "Bulk_probe")).data["current"]
    assert np.all(i <= 0)

    solver["mesh"]["elements"][12]["coordinateIds"] = [2, 1]
    solver["sources"][0]["elementIds"] = [3]
    solver.cleanUp()
    solver.run()
    i = Probe(_get_solved_probe_folder(solver, "Bulk_probe")).data["current"]
    assert np.all(i <= 0)


# compiled without mtln uses classic wires
# compiled with mtln, wire is treated as an unshielded multiwire
@pytest.mark.sgbc
@pytest.mark.wires
@pytest.mark.probes
def test_sgbc_structured_resistance_single_wire(tmp_path):
    """Verify structured SGBC resistance produces the expected wire current."""
    fn = CASES_FOLDER + "sgbcResistance/sgbcResistance.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver["materials"][2] = createWire(id=3, r=1e-4)
    solver.run()

    i = Probe(_get_solved_probe_folder(solver, "Bulk_probe")).data["current"]
    assert np.allclose(i.array[-101:-1], np.ones(100) * i.array[-100], rtol=1e-3)
    assert np.allclose(1 / i.array[-101:-1], np.ones(100) * (50 + 45), rtol=0.05)


# compiled without mtln uses classic wires
# compiled with mtln, wire is treated as an unshielded multiwire
@pytest.mark.sgbc
@pytest.mark.probes
def test_pec_overlapping_sgbcs(tmp_path):
    """Test that PEC surfaces overlapping SGBC surfaces prioritize PEC."""
    def generate_debug_data():
        plt.figure()
        plt.plot(t, iSGBC, ".-", label="SGBC case")
        plt.plot(t, iPEC, ".-", label="PEC overlapping")
        plt.grid(which="both")
        plt.legend()
        plt.savefig(tmp_path / "pec_sgbc_overlap_comparison.png")
        plt.close()

    fn = CASES_FOLDER + "sgbcOverlapping/sgbcOverlapping.fdtd.json"
    # Runs case without overlap.
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p = Probe(_get_solved_probe_folder(solver, "Bulk probe"))
    t = p["time"].to_numpy()
    iSGBC = p["current"].to_numpy()

    # Adds current SGBC elements as PEC. Now both are defined over same surface.
    sgbcElementIds = solver["materialAssociations"][1]["elementIds"]
    solver["materialAssociations"][0]["elementIds"].extend(sgbcElementIds)
    solver.cleanUp()
    solver.run()
    iPEC = Probe(_get_solved_probe_folder(solver, "Bulk probe"))["current"].to_numpy()

    if is_debugging():
        generate_debug_data()

    # Checks values are different due to PEC prioritization.
    assert np.all(np.greater(np.abs(iPEC[1000:]), np.abs(iSGBC[1000:])))


# compiled without mtln uses classic wires
# compiled with mtln, wire is treated as an unshielded multiwire
@pytest.mark.sgbc
@pytest.mark.probes
def test_sgbc_overlapping_sgbc(tmp_path):
    """Test that SGBC surfaces overlapping SGBC surfaces prioritize first in MatAss."""
    def generate_debug_data():
        plt.figure()
        plt.plot(t, iSGBC_top, ".-", label="SGBC sigma = 40 S/m, top")
        plt.plot(t, iSGBC_bottom, ".-", label="SGBC sigma = 20 S/m, bottom")
        plt.grid(which="both")
        plt.legend()
        plt.savefig(tmp_path / "sgbc_overlap_comparison.png")
        plt.close()

    fn = CASES_FOLDER + "sgbcOverlapping/sgbcOverlapping.fdtd.json"

    # Runs case without overlap.
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    # Changes materialId in first SGBC in MatAss to material with larger conductivity.
    solver["materialAssociations"][1]["materialId"] = 6
    solver.cleanUp()
    solver.run()
    p = Probe(_get_solved_probe_folder(solver, "Bulk probe"))

    t = p["time"].to_numpy()
    iSGBC_top = p["current"].to_numpy()

    # Changes materialId in second SGBC in MatAss to material with larger conductivity.
    solver["materialAssociations"][1]["materialId"] = 2
    solver["materialAssociations"][2]["materialId"] = 6
    solver.cleanUp()
    solver.run()
    iSGBC_bottom = Probe(_get_solved_probe_folder(solver, "Bulk probe"))[
        "current"
    ].to_numpy()

    if is_debugging():
        generate_debug_data()

    # Checks values are different due to prioritization of first written.
    assert np.all(np.greater(np.abs(iSGBC_top[1000:]), np.abs(iSGBC_bottom[1000:])))


@pytest.mark.dielectric
@pytest.mark.probes
def test_dielectric_transmission(tmp_path):
    """Verify dielectric reflection, transmission, and propagation delay."""
    _FIELD_TOLERANCE = 0.05

    def getPointProbe(probeName: str):
        dat_filename = _get_solved_probe_folder(solver, probeName, suffix=".dat")
        samples = np.loadtxt(dat_filename, skiprows=1)
        return samples[:, 0], samples[:, 1]

    def getIncidentField(time: np.ndarray, field: np.ndarray) -> Dict:
        idx = field.argmin()
        time = time[idx]
        value = field[idx]
        return {"time": time, "value": value}

    def getReflectedField(time: np.ndarray, field: np.ndarray) -> Dict:
        idx = field.argmax()
        time = time[idx]
        value = field[idx]
        return {"time": time, "value": value}

    def getTransmittedField(time: np.ndarray, field: np.ndarray) -> Dict:
        idx = field.argmin()
        time = time[idx]
        value = field[idx]
        return {"time": time, "value": value}

    def getReflectedDelay(incidentTime: float, reflectedTime: float):
        timeToSurface: float = ((reflectedTime - incidentTime) / 2) + incidentTime
        reflectedDelay: float = reflectedTime - timeToSurface
        return reflectedDelay

    def getTransmittedDelay(
        incidentTime: float, reflectedTime: float, transmittedTime: float
    ):
        timeToSurface: float = ((reflectedTime - incidentTime) / 2) + incidentTime
        transmitedDelay = transmittedTime - timeToSurface
        return transmitedDelay

    fn = CASES_FOLDER + "dielectric/dielectricTransmission.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    relativePermittivity = solver.getMaterialProperties("DielectricMaterial")[
        "relativePermittivity"
    ]
    materialRelativeImpedance = np.sqrt(1 / relativePermittivity)

    expectedReflectedCoeff = (materialRelativeImpedance - 1) / (
        materialRelativeImpedance + 1
    )
    expectedtransmittedCoeff = 1 + expectedReflectedCoeff
    expectedDelayRatio = 1 / np.sqrt(relativePermittivity)

    outsideTime, outsideField = getPointProbe("outside")
    insideTime, insideField = getPointProbe("inside")

    incidentField = getIncidentField(outsideTime, outsideField)
    reflectedField = getReflectedField(outsideTime, outsideField)
    transmittedField = getTransmittedField(insideTime, insideField)

    assert (
        incidentField["value"] - transmittedField["value"] + reflectedField["value"]
    ) < _FIELD_TOLERANCE
    assert np.allclose(
        reflectedField["value"] / incidentField["value"],
        expectedReflectedCoeff,
        rtol=_FIELD_TOLERANCE,
    )
    assert np.allclose(
        transmittedField["value"] / incidentField["value"],
        expectedtransmittedCoeff,
        rtol=_FIELD_TOLERANCE,
    )

    reflectedDelay: float = getReflectedDelay(
        incidentField["time"], reflectedField["time"]
    )
    transmitedDelay: float = getTransmittedDelay(
        incidentField["time"], reflectedField["time"], transmittedField["time"]
    )

    assert np.allclose(
        reflectedDelay / transmitedDelay, expectedDelayRatio, rtol=_FIELD_TOLERANCE
    )


@pytest.mark.conformal
@pytest.mark.probes
def test_rectilinear_mode(tmp_path):
    """Verify rectilinear conformal mode preserves pulse amplitude and timing."""
    _FIELD_TOLERANCE = 4
    _TIME_TOLERANCE = 4

    def getPeakPulse(probe: Probe) -> Dict:
        idx = probe["field"].argmax()
        time = probe["time"][idx]
        value = probe["field"][idx]
        return {"time": time, "value": value}

    rectilinearModeFile = CASES_FOLDER + "rectilinear_mode/rectilinearMode.fdtd.json"
    noRectilinearModeFile = (
        CASES_FOLDER + "rectilinear_mode/noRectilinearMode.fdtd.json"
    )

    rectilinearModeFolder = os.path.join(tmp_path, "rectilinear")
    noRectilinearModeFolder = os.path.join(tmp_path, "noRectilinear")

    os.mkdir(rectilinearModeFolder)
    os.mkdir(noRectilinearModeFolder)

    solverRectilinear = FDTD(
        rectilinearModeFile, path_to_exe=SEMBA_EXE, run_in_folder=rectilinearModeFolder
    )
    solverRectilinear.run()
    rectilinearFrontProbe = Probe(
        _get_solved_probe_folder(solverRectilinear, "Front probe", contains="_Ex_")
    )
    rectilinearVertexProbe = Probe(
        _get_solved_probe_folder(solverRectilinear, "Vertex probe", contains="_Ex_")
    )

    solverNoRectilinear = FDTD(
        noRectilinearModeFile,
        path_to_exe=SEMBA_EXE,
        run_in_folder=noRectilinearModeFolder,
    )
    solverNoRectilinear.run()
    noRectilinearFrontProbe = Probe(
        _get_solved_probe_folder(solverNoRectilinear, "Front probe", contains="_Ex_")
    )
    noRectilinearVertexProbe = Probe(
        _get_solved_probe_folder(solverNoRectilinear, "Vertex probe", contains="_Ex_")
    )

    np.testing.assert_almost_equal(
        getPeakPulse(rectilinearFrontProbe)["value"],
        getPeakPulse(noRectilinearFrontProbe)["value"],
        decimal=_FIELD_TOLERANCE,
    )
    np.testing.assert_almost_equal(
        getPeakPulse(rectilinearFrontProbe)["time"],
        getPeakPulse(noRectilinearFrontProbe)["time"],
        decimal=_TIME_TOLERANCE,
    )
    np.testing.assert_almost_equal(
        getPeakPulse(rectilinearVertexProbe)["value"],
        getPeakPulse(noRectilinearVertexProbe)["value"],
        decimal=_FIELD_TOLERANCE,
    )
    np.testing.assert_almost_equal(
        getPeakPulse(rectilinearVertexProbe)["time"],
        getPeakPulse(noRectilinearVertexProbe)["time"],
        decimal=_TIME_TOLERANCE,
    )


@pytest.mark.dielectric
@pytest.mark.probes
@pytest.mark.vtk
def test_can_execute_fdtd_from_folder_with_spaces_and_can_process_additional_arguments(
    tmp_path,
):
    """Verify execution from paths with spaces and VTK argument handling."""
    folderWithSpaces: str = os.path.join(tmp_path, "spaced bin")
    os.mkdir(folderWithSpaces)
    if platform == "win32":
        shutil.copy2(NGSPICE_DLL, folderWithSpaces)

    sembaExecutable = os.path.basename(SEMBA_EXE)
    pathToExe: str = os.path.join(folderWithSpaces, sembaExecutable)
    shutil.copy2(SEMBA_EXE, pathToExe)

    fn = CASES_FOLDER + "dielectric/dielectricTransmission.fdtd.json"
    solver = FDTD(fn, path_to_exe=pathToExe, run_in_folder=tmp_path, flags=["-mapvtk"])
    solver.run()
    assert Probe(_get_solved_probe_folder(solver, "outside")) is not None
    vtk_map_path = solver.getVTKMap()
    assert vtk_map_path is not None and os.path.isfile(vtk_map_path)


# compiled without mtln uses classic wires
# compiled with mtln, wire is treated as an unshielded multiwire
@pytest.mark.wires
@pytest.mark.nodal_source
@pytest.mark.probes
def test_nodal_source(tmp_path):
    """Verify a nodal source matches its excitation and resistive-wire current."""
    def generate_debug_data():
        plt.figure()
        plt.plot(resistanceBulkProbe["time"], resistanceBulkProbe["current"], label="Resistance")
        plt.plot(excitation.data["time"], excitation.data["value"], label="Excitation")
        plt.plot(nodalBulkProbe["time"], -nodalBulkProbe["current"], label="Nodal source")
        plt.legend()
        plt.savefig(tmp_path / "nodal_source_current_comparison.png")
        plt.close()

    fn = CASES_FOLDER + "nodalSource/nodalSource.fdtd.json"
    assert os.path.isfile(fn)
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["materials"][1] = createWire(id=2, r=0.1e-5, rpul=10000.0)
    solver.run()

    resistanceBulkProbe = Probe(
        _get_solved_probe_folder(solver, "Bulk probe Resistance")
    )
    nodalBulkProbe = Probe(
        _get_solved_probe_folder(solver, "Bulk probe Nodal Source")
    )
    excitation = ExcitationFile(
        excitation_filename=solver.getExcitationFile("predefinedExcitation")[0]
    )

    if is_debugging():
        generate_debug_data()

    exc = np.interp(
        nodalBulkProbe["time"].to_numpy(),
        excitation.data["time"].to_numpy(),
        excitation.data["value"].to_numpy(),
    )
    assert np.corrcoef(exc, -nodalBulkProbe["current"])[0, 1] > 0.999
    assert (
        np.corrcoef(-nodalBulkProbe["current"], resistanceBulkProbe["current"])[0, 1]
        > 0.998
    )


@pytest.mark.wires
@pytest.mark.nodal_source
@pytest.mark.probes
def test_nodal_source_with_total_resistance(tmp_path):
    """Verify that totalResistance in materialAssociation overrides resistancePerMeter from material.

    The nodalSource wire spans 10 cells of 0.001 m each (total 0.01 m).
    totalResistance = 100.0 Ohm  <=>  resistancePerMeter = 10000.0 Ohm/m.
    The material's resistancePerMeter is set to zero and the total resistance is
    supplied through the materialAssociation instead.
    """
    fn = CASES_FOLDER + "nodalSource/nodalSource.fdtd.json"
    assert os.path.isfile(fn)
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["materials"][1]["resistancePerMeter"] = 0.0
    solver["materialAssociations"][1]["totalResistance"] = 100.0
    solver.run()

    resistanceBulkProbe = Probe(
        _get_solved_probe_folder(solver, "Bulk probe Resistance")
    )
    nodalBulkProbe = Probe(
        _get_solved_probe_folder(solver, "Bulk probe Nodal Source")
    )
    excitation = ExcitationFile(
        excitation_filename=solver.getExcitationFile("predefinedExcitation")[0]
    )

    exc = np.interp(
        nodalBulkProbe["time"].to_numpy(),
        excitation.data["time"].to_numpy(),
        excitation.data["value"].to_numpy(),
    )
    assert np.corrcoef(exc, -nodalBulkProbe["current"])[0, 1] > 0.999
    assert (
        np.corrcoef(-nodalBulkProbe["current"], resistanceBulkProbe["current"])[0, 1]
        > 0.998
    )


@pytest.mark.sgbc
def test_can_assign_same_surface_impedance_to_multiple_geometries(tmp_path):
    """Verify one surface-impedance material can serve multiple geometries."""
    fn = CASES_FOLDER + "multipleAssigments/multipleSurfaceImpedance.fdtd.json"

    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    assert Probe(_get_solved_probe_folder(solver, "BulkProbeEntry")) is not None


@pytest.mark.dielectric
def test_can_assign_same_dielectric_material_to_multiple_geometries(tmp_path):
    """Verify one dielectric material can serve multiple geometries."""
    fn = CASES_FOLDER + "multipleAssigments/multipleDielectricMaterial.fdtd.json"

    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    assert Probe(_get_solved_probe_folder(solver, "BulkProbeEntry")) is not None


@mtln_skip
@pytest.mark.lumped
@pytest.mark.probes
def test_lumped_resistor(tmp_path):
    """Verify a lumped resistor matches terminal and theoretical responses."""
    # This test validates the behavior of lumped resistor materials in a simplified circuit.
    # The circuit consists of a 40mm x 40mm simple loop with a lumped resistor line inserted along one edge.
    # Due to the geometry, the circuit naturally exhibits a parasitic inductance of approximately 1.65e-7 H.
    #
    # Current probes are placed at the start and end of the lumped line, as well as a few cells before and after it.
    # These measurements are used to evaluate the accuracy of the lumped material model.
    #
    # For validation, the results are compared against two reference cases:
    # 1. A simple loop circuit with the same dimensions where a terminal  is inserted in place of the lumped line,
    #    using the same resistance of the lumped line.
    # 2. Theoretical current response calculated using Laplace transforms from the initial pulse excitation.
    #
    # For better interaction with the case, the user can go to the file: testData/cases/lumped_lines/simple_loop_R/simple_loop_prepost.py

    fn_lumped = CASES_FOLDER + "lumped_lines/simple_loop_R/simple_loop_lumped.fdtd.json"
    fn_terminal = (
        CASES_FOLDER + "lumped_lines/simple_loop_R/simple_loop_terminal.fdtd.json"
    )

    solver_lumped = FDTD(fn_lumped, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver_terminal = FDTD(fn_terminal, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver_lumped.run()
    solver_terminal.run()

    StartTerminalProbe = Probe(
        _get_solved_probe_folder(solver_terminal, "TerminalCellStart")
    )
    StartLumpedProbe = Probe(
        _get_solved_probe_folder(solver_lumped, "LumpedCellStart")
    )

    EndTerminalProbe = Probe(
        _get_solved_probe_folder(solver_terminal, "TerminalCellEnd")
    )
    EndLumpedProbe = Probe(_get_solved_probe_folder(solver_lumped, "LumpedCellEnd"))

    AdjacentPostLumpedProbe = Probe(
        _get_solved_probe_folder(solver_lumped, "PostLumpedCell")
    )
    AdjacentPostTerminalProbe = Probe(
        _get_solved_probe_folder(solver_terminal, "PostTerminalCell")
    )

    AdjacentPreLumpedProbe = Probe(
        _get_solved_probe_folder(solver_lumped, "PreLumpedCell")
    )
    AdjacentPreTerminalProbe = Probe(
        _get_solved_probe_folder(solver_terminal, "PreTerminalCell")
    )

    assert (
        np.corrcoef(
            StartLumpedProbe["current"].to_numpy(),
            StartTerminalProbe["current"].to_numpy(),
        )[0, 1]
        > 0.999
    )
    assert (
        np.corrcoef(
            EndLumpedProbe["current"].to_numpy(), EndTerminalProbe["current"].to_numpy()
        )[0, 1]
        > 0.999
    )
    assert (
        np.corrcoef(
            AdjacentPostLumpedProbe["current"].to_numpy(),
            AdjacentPostTerminalProbe["current"].to_numpy(),
        )[0, 1]
        > 0.999
    )
    assert (
        np.corrcoef(
            AdjacentPreLumpedProbe["current"].to_numpy(),
            AdjacentPreTerminalProbe["current"].to_numpy(),
        )[0, 1]
        > 0.999
    )

    R = solver_lumped.getMaterialProperties("lumped_line")["resistance"]
    L = 1.65e-7  # parasitic inductance mentioned above

    num = [1]
    den = [L, R]
    system = signal.TransferFunction(num, den)
    tout, I_out, _ = signal.lsim(
        system,
        U=np.loadtxt(_get_source_magnitude_file(solver_lumped), usecols=1),
        T=np.loadtxt(_get_source_magnitude_file(solver_lumped), usecols=0),
    )

    I_theo = np.interp(AdjacentPreLumpedProbe["time"], tout, I_out)

    assert (
        np.corrcoef(AdjacentPostLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999
    )
    assert (
        np.corrcoef(AdjacentPreLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999
    )
    assert np.corrcoef(StartLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(EndLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999


@mtln_skip
@pytest.mark.lumped
@pytest.mark.probes
def test_lumped_capacitor(tmp_path):
    """Verify a lumped capacitor matches the theoretical circuit response."""
    # This test validates the behavior of lumped capacitor materials in a simplified circuit. The lumped capacitor
    # can be modeled as a capacitor in parallel with a resistor.
    # The circuit consists of a 40mm x 40mm simple loop with a lumped capacitor line inserted along one edge.
    # Due to the geometry, the circuit naturally exhibits a parasitic inductance of approximately 1.65e-7 H.
    #
    # Current probes are placed at the start and end of the lumped line, as well as a few cells before and after it.
    # These measurements are used to evaluate the accuracy of the lumped material model.
    #
    # For validation, the results are compared against only one reference case:
    # 1. Theoretical current response calculated using Laplace transforms from the initial pulse excitation.
    #
    # For better interaction with the case, the user can go to the file: testData/cases/lumped_lines/simple_loop_RC/simple_loop_prepost.py

    fn_lumped = (
        CASES_FOLDER + "lumped_lines/simple_loop_RC/simple_loop_lumped.fdtd.json"
    )
    solver_lumped = FDTD(fn_lumped, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver_lumped.run()

    StartLumpedProbe = Probe(
        _get_solved_probe_folder(solver_lumped, "LumpedCellStart")
    )
    EndLumpedProbe = Probe(_get_solved_probe_folder(solver_lumped, "LumpedCellEnd"))
    AdjacentPostLumpedProbe = Probe(
        _get_solved_probe_folder(solver_lumped, "PostLumpedCell")
    )
    AdjacentPreLumpedProbe = Probe(
        _get_solved_probe_folder(solver_lumped, "PreLumpedCell")
    )

    R = solver_lumped.getMaterialProperties("lumped_RC")["resistance"]
    C = solver_lumped.getMaterialProperties("lumped_RC")["capacitance"]
    L = 1.65e-7  # parasitic inductance mentioned above

    num = [R * C, 1]
    den = [L * R * C, L, R]
    system = signal.TransferFunction(num, den)
    tout, I_out, _ = signal.lsim(
        system,
        U=np.loadtxt(_get_source_magnitude_file(solver_lumped), usecols=1),
        T=np.loadtxt(_get_source_magnitude_file(solver_lumped), usecols=0),
    )

    I_theo = np.interp(AdjacentPreLumpedProbe["time"], tout, I_out)

    assert (
        np.corrcoef(AdjacentPostLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999
    )
    assert (
        np.corrcoef(AdjacentPreLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999
    )
    assert np.corrcoef(StartLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(EndLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999


@mtln_skip
@pytest.mark.lumped
@pytest.mark.probes
def test_lumped_inductor(tmp_path):
    """Verify a lumped inductor matches terminal and theoretical responses."""
    # This test validates the behavior of lumped inductor materials in a simplified circuit. The lumped inductor
    # can be modeled as an inductor in series with a resistor.
    # The circuit consists of a 40mm x 40mm simple loop with a lumped resistor line inserted along one edge.
    # Due to the geometry, the circuit naturally exhibits a parasitic inductance of approximately 1.65e-7 H.
    #
    # Current probes are placed at the start and end of the lumped line, as well as a few cells before and after it.
    # These measurements are used to evaluate the accuracy of the lumped material model.
    #
    # For validation, the results are compared against two reference cases:
    # 1. A simple loop circuit with the same dimensions where a terminal  is inserted in place of the lumped line,
    #    using the same resistance and inductance of the lumped line.
    # 2. Theoretical current response calculated using Laplace transforms from the initial pulse excitation.
    #
    # For better interaction with the case, the user can go to the file: testData/cases/lumped_lines/simple_loop_RL/simple_loop_prepost.py

    fn_lumped = (
        CASES_FOLDER + "lumped_lines/simple_loop_RL/simple_loop_lumped.fdtd.json"
    )
    fn_terminal = (
        CASES_FOLDER + "lumped_lines/simple_loop_RL/simple_loop_terminal.fdtd.json"
    )

    solver_lumped = FDTD(fn_lumped, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver_terminal = FDTD(fn_terminal, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver_lumped.run()
    solver_terminal.run()

    StartTerminalProbe = Probe(
        _get_solved_probe_folder(solver_terminal, "TerminalCellStart")
    )
    StartLumpedProbe = Probe(
        _get_solved_probe_folder(solver_lumped, "LumpedCellStart")
    )

    EndTerminalProbe = Probe(
        _get_solved_probe_folder(solver_terminal, "TerminalCellEnd")
    )
    EndLumpedProbe = Probe(_get_solved_probe_folder(solver_lumped, "LumpedCellEnd"))

    AdjacentPostLumpedProbe = Probe(
        _get_solved_probe_folder(solver_lumped, "PostLumpedCell")
    )
    AdjacentPostTerminalProbe = Probe(
        _get_solved_probe_folder(solver_terminal, "PostTerminalCell")
    )

    AdjacentPreLumpedProbe = Probe(
        _get_solved_probe_folder(solver_lumped, "PreLumpedCell")
    )
    AdjacentPreTerminalProbe = Probe(
        _get_solved_probe_folder(solver_terminal, "PreTerminalCell")
    )

    assert (
        np.corrcoef(
            StartLumpedProbe["current"].to_numpy(),
            StartTerminalProbe["current"].to_numpy(),
        )[0, 1]
        > 0.999
    )
    assert (
        np.corrcoef(
            EndLumpedProbe["current"].to_numpy(), EndTerminalProbe["current"].to_numpy()
        )[0, 1]
        > 0.999
    )
    assert (
        np.corrcoef(
            AdjacentPostLumpedProbe["current"].to_numpy(),
            AdjacentPostTerminalProbe["current"].to_numpy(),
        )[0, 1]
        > 0.999
    )
    assert (
        np.corrcoef(
            AdjacentPreLumpedProbe["current"].to_numpy(),
            AdjacentPreTerminalProbe["current"].to_numpy(),
        )[0, 1]
        > 0.999
    )

    R = solver_lumped.getMaterialProperties("lumped_RL")["resistance"]
    L = (
        solver_lumped.getMaterialProperties("lumped_RL")["inductance"] + 1.65e-7
    )  # inductance of the lumped line + parasitic inductance mentioned above

    num = [1]
    den = [L, R]
    system = signal.TransferFunction(num, den)
    tout, I_out, _ = signal.lsim(
        system,
        U=np.loadtxt(_get_source_magnitude_file(solver_lumped), usecols=1),
        T=np.loadtxt(_get_source_magnitude_file(solver_lumped), usecols=0),
    )

    I_theo = np.interp(AdjacentPreLumpedProbe["time"], tout, I_out)

    assert (
        np.corrcoef(AdjacentPostLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999
    )
    assert (
        np.corrcoef(AdjacentPreLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999
    )
    assert np.corrcoef(StartLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(EndLumpedProbe["current"].to_numpy(), I_theo)[0, 1] > 0.999


@mtln_skip
@pytest.mark.lumped
@pytest.mark.termination
@pytest.mark.probes
def test_lumped_resistor_parallel_terminal_resistor(tmp_path):
    """Verify current conservation in parallel lumped and terminal resistors."""
    # This test verifies current splitting behavior in a parallel resistive configuration.
    # The setup consists of a 40mm x 40mm circuit with two parallel elements:
    # - A lumped resistor line
    # - A terminal with the same resistance
    #
    # Current measurements are taken at three key locations:
    # 1. At the source (input) (Bulk Initial probe)
    # 2. At the terminal branch (Bulk Top probe)
    # 3. At the lumped resistor line branch (Bulk Bottom probe)
    #
    # The goal is to validate that the current divides between the two parallel resistive paths, i.e., the
    # current at the terminal branch plus the lumped resistor line branch should be equal to the current at the source.
    # The test also compares the current at the source with the theoretical current response calculated using Laplace transforms.
    #
    # For better interaction with the case, the user can go to the file: testData/cases/lumped_lines/current_bifurcation/current_bifurcation_lumped_prepost.py

    fn = (
        CASES_FOLDER
        + "lumped_lines/current_bifurcation/current_bifurcation_lumped.fdtd.json"
    )
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    InitialBulk_probe = Probe(_get_solved_probe_folder(solver, "Bulk Initial probe"))
    TopBulk_probe = Probe(_get_solved_probe_folder(solver, "Bulk Top probe"))
    BottomBulk_probe = Probe(_get_solved_probe_folder(solver, "Bulk Bottom probe"))

    R_lumped = solver.getMaterialProperties("lumped_resistor")["resistance"]
    R_terminal = solver.getMaterialProperties("Terminal_R")["terminations"][0][
        "resistance"
    ]
    R = 1 / (1 / R_lumped + 1 / R_terminal)
    L = 1.65e-7  # parasitic inductance mentioned above

    num = [1]
    den = [L, R]
    system = signal.TransferFunction(num, den)
    tout, I_out, _ = signal.lsim(
        system,
        U=np.loadtxt(_get_source_magnitude_file(solver), usecols=1),
        T=np.loadtxt(_get_source_magnitude_file(solver), usecols=0),
    )

    I_theo = np.interp(InitialBulk_probe["time"], tout, I_out)

    assert (
        np.corrcoef(
            TopBulk_probe["current"].to_numpy()
            + BottomBulk_probe["current"].to_numpy(),
            I_theo,
        )[0, 1]
        > 0.999
    )
    assert np.corrcoef(InitialBulk_probe["current"].to_numpy(), I_theo)[0, 1] > 0.999


@pytest.mark.nodal_source
@pytest.mark.probes
def test_bulk_current_offset_normal_in_x(tmp_path):
    """Verify bulk-current probes apply the normal x-direction offset."""
    # This test verifies the positive offset along the normal vector (in x-direction) respect to the bulk plane
    # used to measure the current values of the system.
    # The setup consists in a polyline with points [(0 mm,0 mm,0 mm), (20 mm,0 mm,0 mm), (20 mm,-20 mm,0 mm)]
    # as nodal source and three bulk planes defined at x=18 mm, x=20 mm and x=22 mm.
    # The test checks that only the plane defined at x=18 mm has non-zero current values.

    fn = CASES_FOLDER + "bulk_current_offsets/offSet_x/offSet_x.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    I_in = np.loadtxt(_get_source_magnitude_file(solver), usecols=1)
    time = np.loadtxt(_get_source_magnitude_file(solver), usecols=0)

    probe_at_x_18 = Probe(_get_solved_probe_folder(solver, "BulkCurrent1"))
    probe_at_x_20 = Probe(_get_solved_probe_folder(solver, "BulkCurrent2"))
    probe_at_x_22 = Probe(_get_solved_probe_folder(solver, "BulkCurrent3"))

    I_interp = np.interp(probe_at_x_18["time"].to_numpy(), time, I_in)

    assert np.corrcoef(probe_at_x_18["current"].to_numpy(), I_interp)[0, 1] > 0.999
    assert np.allclose(probe_at_x_20["current"].to_numpy(), 0.0, atol=1.5e-3)
    assert np.allclose(probe_at_x_22["current"].to_numpy(), 0.0, atol=1.5e-3)


@pytest.mark.nodal_source
@pytest.mark.probes
def test_bulk_current_offset_normal_in_y(tmp_path):
    """Verify bulk-current probes apply the normal y-direction offset."""
    # This test verifies the positive offset along the normal vector (in y-direction) respect to the bulk plane
    # used to measure the current values of the system.
    # The setup consists in a polyline with points [(0 mm,0 mm,0 mm), (0 mm,0 mm,20 mm), (0 mm,-20 mm,20 mm)]
    # as nodal source and three bulk planes defined at y=-2 mm, y=0 mm and y=2 mm.
    # The test checks that only the plane defined at y=-2 mm has non-zero current values.

    fn = CASES_FOLDER + "bulk_current_offsets/offSet_y/offSet_y.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    I_in = np.loadtxt(_get_source_magnitude_file(solver), usecols=1)
    time = np.loadtxt(_get_source_magnitude_file(solver), usecols=0)

    probe_at_y_m2 = Probe(_get_solved_probe_folder(solver, "BulkCurrent1"))
    probe_at_y_0 = Probe(_get_solved_probe_folder(solver, "BulkCurrent2"))
    probe_at_y_2 = Probe(_get_solved_probe_folder(solver, "BulkCurrent3"))

    I_interp = np.interp(probe_at_y_m2["time"].to_numpy(), time, I_in)

    assert (
        np.corrcoef(np.abs(probe_at_y_m2["current"].to_numpy()), I_interp)[0, 1] > 0.999
    )
    assert np.allclose(probe_at_y_0["current"].to_numpy(), 0.0, atol=1.5e-3)
    assert np.allclose(probe_at_y_2["current"].to_numpy(), 0.0, atol=1.5e-3)


@pytest.mark.nodal_source
@pytest.mark.probes
def test_bulk_current_offset_normal_in_z(tmp_path):
    """Verify bulk-current probes apply the normal z-direction offset."""
    # This test verifies the positive offset along the normal vector (in z-direction) respect to the bulk plane
    # used to measure the current values of the system.
    # The setup consists in a polyline with points [(0 mm,0 mm,0 mm), (0 mm,0 mm,20 mm), (0 mm,-20 mm,20 mm)]
    # as nodal source and three bulk planes defined at z=18 mm, z=20 mm and z=22 mm.
    # The test checks that only the plane defined at z=18 mm has non-zero current values.

    fn = CASES_FOLDER + "bulk_current_offsets/offSet_z/offSet_z.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    I_in = np.loadtxt(_get_source_magnitude_file(solver), usecols=1)
    time = np.loadtxt(_get_source_magnitude_file(solver), usecols=0)

    probe_at_z_18 = Probe(_get_solved_probe_folder(solver, "BulkCurrent1"))
    probe_at_z_20 = Probe(_get_solved_probe_folder(solver, "BulkCurrent2"))
    probe_at_z_22 = Probe(_get_solved_probe_folder(solver, "BulkCurrent3"))

    I_interp = np.interp(probe_at_z_18["time"].to_numpy(), time, I_in)

    assert np.corrcoef(probe_at_z_18["current"].to_numpy(), I_interp)[0, 1] > 0.999
    assert np.allclose(probe_at_z_20["current"].to_numpy(), 0.0, atol=1.5e-3)
    assert np.allclose(probe_at_z_22["current"].to_numpy(), 0.0, atol=1.5e-3)


@pytest.mark.nodal_source
@pytest.mark.probes
def test_bulk_current_offset_perpendicular_in_x(tmp_path):
    """Verify x-normal bulk-current probes apply perpendicular offsets."""
    # This test verifies the negative offset presented in the y and z directions when the bulk plane is defined
    # with a normal vector in the x-direction.
    # The setup consists in three nodal sources defined as follow:
    #   - Line [(0 mm,18 mm,20 mm), (50 mm,18 mm,20 mm)]. Gaussian pulse with 1 A amplitude.
    #   - Line [(0 mm,20 mm,20 mm), (50 mm,20 mm,20 mm)]. Gaussian pulse with 2 A amplitude.
    #   - Line [(0 mm,18 mm,18 mm), (50 mm,18 mm,18 mm)]. Gaussian pulse with 3 A amplitude.
    #
    # The nodal lines are separeted by one cell in y and z directions. We define four bulk planes:
    #   - Plane at x=30 mm, (16 mm, 20 mm) and (20 mm, 24 mm) in y and z directions. For nodal source 1.
    #   - Plane at x=30 mm, (20 mm, 24 mm) and (20 mm, 24 mm) in y and z directions. For nodal source 2.
    #   - Plane at x=30 mm, (16 mm, 20 mm) and (16 mm, 20 mm) in y and z directions. For nodal source 3.
    #   - Plane at x=30 mm, (16 mm, 24 mm) and (16 mm, 24 mm) in y and z directions. For the total current.
    #
    # The test checks that the bulk planes only measure the current values of the respective nodal sources
    # and if we move one negative cell in the y and z directions, the current values are zero for the bulk planes
    # 1 and 3 respectively; similar behavior if we move one positive cell in the y and z directions for
    # the bulk planes 1 and 2. This proves that the bulk have a negative offset in the y and z directions.

    fn = (
        CASES_FOLDER
        + "bulk_current_offsets/threeLines_offSet_x_Perpendicular/threeLines.fdtd.json"
    )
    fn_negative = (
        CASES_FOLDER
        + "bulk_current_offsets/threeLines_offSet_x_Perpendicular/threeLinesNegative.fdtd.json"
    )
    fn_positive = (
        CASES_FOLDER
        + "bulk_current_offsets/threeLines_offSet_x_Perpendicular/threeLinesPositive.fdtd.json"
    )

    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver_negative = FDTD(
        input_filename=fn_negative, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path
    )
    solver_positive = FDTD(
        input_filename=fn_positive, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path
    )

    solver.run()
    solver_negative.run()
    solver_positive.run()

    probe1 = Probe(_get_solved_probe_folder(solver, "Bulk probe1"))
    probe2 = Probe(_get_solved_probe_folder(solver, "Bulk probe2"))
    probe3 = Probe(_get_solved_probe_folder(solver, "Bulk probe3"))
    probeTotal = Probe(_get_solved_probe_folder(solver, "Bulk probe total"))

    assert (
        np.corrcoef(
            probe1["current"].to_numpy()
            + probe2["current"].to_numpy()
            + probe3["current"].to_numpy(),
            probeTotal["current"].to_numpy(),
        )[0, 1]
        > 0.999
    )

    probe1_negative = Probe(_get_solved_probe_folder(solver_negative, "Bulk probe1"))
    probe2_negative = Probe(_get_solved_probe_folder(solver_negative, "Bulk probe2"))
    probe3_negative = Probe(_get_solved_probe_folder(solver_negative, "Bulk probe3"))

    I_in_2 = np.loadtxt(_get_source_magnitude_file(solver, 1), usecols=1)
    time_2 = np.loadtxt(_get_source_magnitude_file(solver, 1), usecols=0)

    I_2_interp = np.interp(probe2_negative["time"].to_numpy(), time_2, I_in_2)

    assert np.corrcoef(probe2_negative["current"].to_numpy(), I_2_interp)[0, 1] > 0.999
    assert np.allclose(probe1_negative["current"].to_numpy(), 0.0, atol=1.5e-2)
    assert np.allclose(probe3_negative["current"].to_numpy(), 0.0, atol=1.5e-2)

    probe1_positive = Probe(_get_solved_probe_folder(solver_positive, "Bulk probe1"))
    probe2_positive = Probe(_get_solved_probe_folder(solver_positive, "Bulk probe2"))
    probe3_positive = Probe(_get_solved_probe_folder(solver_positive, "Bulk probe3"))

    I_in_3 = np.loadtxt(_get_source_magnitude_file(solver, 2), usecols=1)
    time_3 = np.loadtxt(_get_source_magnitude_file(solver, 2), usecols=0)

    I_3_interp = np.interp(probe3_positive["time"].to_numpy(), time_3, I_in_3)

    assert np.corrcoef(probe3_positive["current"].to_numpy(), I_3_interp)[0, 1] > 0.999
    assert np.allclose(probe1_positive["current"].to_numpy(), 0.0, atol=1.5e-2)
    assert np.allclose(probe2_positive["current"].to_numpy(), 0.0, atol=1.5e-2)


@pytest.mark.nodal_source
@pytest.mark.probes
def test_bulk_current_negative_offset_in_x(tmp_path):
    """Verify y-normal bulk-current probes apply a negative x offset."""
    # Following the previous test, we have seen that the bulk surfaces has negative offsets in the directions
    # perpendicular to the normal vector of the bulk surface. The previous test checks the negative offset in the
    # y and z directions when the normal vector is in the x-direction. Now we check the negative offset in the
    # x-direction when the normal vector is in the y-direction.
    # The setup consists in a nodal source placed in a line defined by [(0 mm,0 mm,0 mm), (0 mm,50 mm,0 mm)].
    # The nodal source is a Gaussian pulse with 1 A amplitude. And three bulk planes defined at:
    #  - Plane at y=36 mm, (-4 mm, -2 mm) and (0 mm, 2 mm) in x and z directions.
    #  - Plane at y=38 mm, (0 mm, -2 mm) and (4 mm, 2 mm) in x and z directions.
    #  - Plane at y=40 mm, (-2 mm, -2 mm) and (2 mm, 2 mm) in x and z directions.
    #
    # All the bulk planes measure the current values of the nodal source, however, the first only captures the
    # current values at the right edge of the plane, similarly, the second captures the current values at the
    # left edge of the plane. This test checks that the second and third bulk planes measure correctly the current values
    # while the first bulk plane has zero current values. This proves that the bulk have a negative offset in the x-direction.

    fn = (
        CASES_FOLDER
        + "bulk_current_offsets/negative_offSet_x/offSet_negative_x.fdtd.json"
    )
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    I_in = np.loadtxt(_get_source_magnitude_file(solver), usecols=1)
    time = np.loadtxt(_get_source_magnitude_file(solver), usecols=0)

    probeR = Probe(_get_solved_probe_folder(solver, "Bulk_right"))
    probeL = Probe(_get_solved_probe_folder(solver, "Bulk_left"))
    probeTotal = Probe(_get_solved_probe_folder(solver, "BulkTotal"))

    I_interp = np.interp(probeTotal["time"].to_numpy(), time, I_in)

    assert np.corrcoef(probeTotal["current"].to_numpy(), I_interp)[0, 1] > 0.999
    assert np.corrcoef(probeL["current"].to_numpy(), I_interp)[0, 1] > 0.999
    assert np.allclose(probeR["current"].to_numpy(), 0.0, atol=3e-3)


def _run_four_probes(tmp_path, json_filename):
    """Helper: run a four-probe bulk-current case and return the four Probe objects
    together with the interpolated excitation evaluated at the BC_LL time grid."""
    solver = FDTD(
        input_filename=CASES_FOLDER + "bulk_current_four_probes/" + json_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
    )
    solver.run()

    exc_time = np.loadtxt(_get_source_magnitude_file(solver), usecols=0)
    exc_val = np.loadtxt(_get_source_magnitude_file(solver), usecols=1)

    probe_LL = Probe(_get_solved_probe_folder(solver, "BC_LL"))
    probe_LU = Probe(_get_solved_probe_folder(solver, "BC_LU"))
    probe_UU = Probe(_get_solved_probe_folder(solver, "BC_UU"))
    probe_UL = Probe(_get_solved_probe_folder(solver, "BC_UL"))

    probe_time = probe_LL["time"].to_numpy()
    exc_interp = np.interp(probe_time, exc_time, exc_val)

    return probe_LL, probe_LU, probe_UU, probe_UL, exc_interp


@pytest.mark.nodal_source
@pytest.mark.probes
def test_bulk_current_four_probes_X_oriented(tmp_path):
    """Verify x-oriented current is captured only by its intersecting probe."""
    # A nodal current source runs along X through cell (23,23).
    # Four bulk-current probes are arranged in the YZ plane at x=4:
    #   BC_LL and BC_LU and BC_UL lie outside the wire path -> near zero.
    #   BC_UU contains the wire -> correlates with the excitation.
    probe_LL, probe_LU, probe_UU, probe_UL, exc_interp = _run_four_probes(
        tmp_path, "bulk_currents_X_oriented.fdtd.json"
    )

    assert np.corrcoef(exc_interp, probe_UU["current"].to_numpy())[0, 1] > 0.9999
    assert np.allclose(probe_LL["current"].to_numpy(), 0.0, atol=2e-3)
    assert np.allclose(probe_LU["current"].to_numpy(), 0.0, atol=2e-3)
    assert np.allclose(probe_UL["current"].to_numpy(), 0.0, atol=2e-3)


@pytest.mark.nodal_source
@pytest.mark.probes
def test_bulk_current_four_probes_Y_oriented(tmp_path):
    """Verify y-oriented current is captured only by its intersecting probe."""
    # A nodal current source runs along Y through cell (23,23).
    # Four bulk-current probes are arranged in the XZ plane at y=4:
    #   BC_LL and BC_LU and BC_UL lie outside the wire path -> near zero.
    #   BC_UU contains the wire -> correlates with the excitation.
    probe_LL, probe_LU, probe_UU, probe_UL, exc_interp = _run_four_probes(
        tmp_path, "bulk_currents_Y_oriented.fdtd.json"
    )

    assert np.corrcoef(exc_interp, probe_UU["current"].to_numpy())[0, 1] > 0.9999
    assert np.allclose(probe_LL["current"].to_numpy(), 0.0, atol=2e-3)
    assert np.allclose(probe_LU["current"].to_numpy(), 0.0, atol=2e-3)
    assert np.allclose(probe_UL["current"].to_numpy(), 0.0, atol=2e-3)


@pytest.mark.nodal_source
@pytest.mark.probes
def test_bulk_current_four_probes_Z_oriented(tmp_path):
    """Verify z-oriented current is captured only by its intersecting probe."""
    # A nodal current source runs along Z through cell (23,23).
    # Four bulk-current probes are arranged in the XY plane at z=4:
    #   BC_LL and BC_LU and BC_UL lie outside the wire path -> near zero.
    #   BC_UU contains the wire -> correlates with the excitation.
    probe_LL, probe_LU, probe_UU, probe_UL, exc_interp = _run_four_probes(
        tmp_path, "bulk_currents_Z_oriented.fdtd.json"
    )

    assert np.corrcoef(exc_interp, probe_UU["current"].to_numpy())[0, 1] > 0.9999
    assert np.allclose(probe_LL["current"].to_numpy(), 0.0, atol=2e-3)
    assert np.allclose(probe_LU["current"].to_numpy(), 0.0, atol=2e-3)
    assert np.allclose(probe_UL["current"].to_numpy(), 0.0, atol=2e-3)


# compiled without mtln uses classic wires
# compiled with mtln, wire is treated as an unshielded multiwire
@pytest.mark.conformal
@pytest.mark.wires
@pytest.mark.probes
def test_conformal_impedance_cylinder_unshielded(tmp_path):
    """Verify conformal-cylinder impedance matches the reference spectrum."""
    case_name = "conformal_impedance_cylinder_conformal"
    solver = FDTD(
        input_filename=TEST_DATA_FOLDER
        + "cases/conformal_impedance_cylinder/"
        + case_name
        + ".fdtd.json",
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
    )
    solver.cleanUp()
    solver.run()
    assert solver.hasFinishedSuccessfully()
    bulk_conf = Probe(_get_solved_probe_folder(solver, "BulkProbe"))

    # discrete fourier transforms
    exc_file = solver.getExcitationFile("predefinedExcitation")[0]
    exc = pd.read_csv(exc_file, sep="\\s+")
    exc = exc.rename(columns={exc.columns[0]: "time", exc.columns[1]: "V"})
    new_freqs = np.geomspace(1e3, 1e6, num=100)
    Vexc = exc["V"].to_numpy()
    texc = exc["time"].to_numpy()
    dt_exc = texc[1] - texc[0]
    Vfexc = dt_exc * np.array(
        [np.sum(Vexc * np.exp(-1j * 2 * np.pi * f * texc)) for f in new_freqs]
    )

    Ibulk_conf = bulk_conf["current"].to_numpy()
    tbulk_conf = bulk_conf["time"].to_numpy()
    dt_bulk_conf = tbulk_conf[1] - tbulk_conf[0]
    Ifbulk_conf = dt_bulk_conf * np.array(
        [
            np.sum(Ibulk_conf * np.exp(-1j * 2 * np.pi * f * tbulk_conf))
            for f in new_freqs
        ]
    )

    data = pd.read_csv(
        OUTPUTS_FOLDER + "conformal_cylinder_impedance_output.dat", sep=" ", header=0
    )
    data.columns = ["f", "z"]

    assert np.corrcoef(data["z"], np.abs(Vfexc / Ifbulk_conf))[0, 1] > 0.999


@pytest.mark.conformal
@pytest.mark.farfield
@pytest.mark.probes
def test_conformal_sphere_rcs(tmp_path):
    """Verify conformal-sphere RCS agrees with the analytical solution."""
    case_name = "conformal_sphere_rcs"
    solver = FDTD(
        input_filename=TEST_DATA_FOLDER + "cases/conformal/" + case_name + ".fdtd.json",
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
    )
    solver.run()
    assert solver.hasFinishedSuccessfully()

    far = Probe(_get_solved_probe_folder(solver, "n2f"))
    ra = far.data.loc[
        (far.data["Theta"] == 90.0) & (far.data["Phi"] == 0.0), "rcs_arit"
    ]
    rg = far.data.loc[
        (far.data["Theta"] == 90.0) & (far.data["Phi"] == 0.0), "rcs_geom"
    ]
    ffar = far.data.loc[(far.data["Theta"] == 90.0) & (far.data["Phi"] == 0.0), "freq"]

    # analytical RCS
    f = np.linspace(1e7, 7e8, 200)
    r = 0.5  # in meters
    rcs = RCS(fspace=f, radius=r)
    # simulated, interpolated to analytical frequencies
    rcs_interp = np.interp(f, ffar, rg)

    assert np.corrcoef(rcs[5:150], rcs_interp[5:150])[0, 1] > 0.98


@pytest.mark.conformal
@pytest.mark.probes
def test_conformal_delay(tmp_path):
    """Verify conformal geometry produces the expected propagation delay."""
    fn = CASES_FOLDER + "conformal/conformal.fdtd.json"
    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )

    solver["materialAssociations"][0]["elementIds"] = [4]
    solver["mesh"]["elements"][3]["intervals"] = [[[0, 0, 4], [2, 2, 4]]]
    solver.cleanUp()
    solver.run()
    front = Probe(_get_solved_probe_folder(solver, "front"))
    t = front["time"]
    t4 = t[front["field"].argmin()]

    solver["materialAssociations"][0]["elementIds"] = [5]
    n = 4
    for i in range(1, n):
        solver["mesh"]["coordinates"][6]["relativePosition"][2] = 4.0 + i * 1.0 / n
        solver["mesh"]["coordinates"][7]["relativePosition"][2] = 4.0 + i * 1.0 / n
        solver["mesh"]["coordinates"][8]["relativePosition"][2] = 4.0 + i * 1.0 / n
        solver["mesh"]["coordinates"][9]["relativePosition"][2] = 4.0 + i * 1.0 / n
        solver["mesh"]["coordinates"][15]["relativePosition"][2] = 4.0 + i * 1.0 / n
        solver["mesh"]["coordinates"][16]["relativePosition"][2] = 4.0 + i * 1.0 / n
        solver["mesh"]["coordinates"][17]["relativePosition"][2] = 4.0 + i * 1.0 / n
        solver["mesh"]["coordinates"][18]["relativePosition"][2] = 4.0 + i * 1.0 / n
        solver["mesh"]["coordinates"][19]["relativePosition"][2] = 4.0 + i * 1.0 / n

        solver.cleanUp()
        solver.run()
        front = Probe(_get_solved_probe_folder(solver, "front"))
        t = front["time"]
        delay = t[front["field"].argmin()]
        tdelta = t4 + 2 * (i * 1.0 / n) * 0.02 / 3e8
        assert np.abs(delay - tdelta) / tdelta < 0.01


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.nodal_source
@pytest.mark.probes
def test_current_generators_with_resistance(tmp_path):
    """Verify resistive-wire current-source voltage and current distribution."""
    # Checks current and voltage of probes at the extremes of a wire
    # with a current generator in the middle of the wire

    fn = CASES_FOLDER + "sources/sources_current.fdtd.json"
    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["sources"][0]["elementIds"] = [1]
    solver.cleanUp()
    solver.run()
    Iend = Probe(_get_solved_probe_folder(solver, "probe_end", contains="_I_"))
    Istart = Probe(_get_solved_probe_folder(solver, "probe_start", contains="_I_"))
    Vend = Probe(_get_solved_probe_folder(solver, "probe_end", contains="_V_"))
    Vstart = Probe(_get_solved_probe_folder(solver, "probe_start", contains="_V_"))

    assert np.allclose(Iend["current_0"][-100:-1], 1.0 / 3.0, rtol=0.005)
    assert np.allclose(Istart["current_0"][-100:-1], 1.0 / 3.0, rtol=0.005)
    assert np.allclose(Vend["voltage_0"][-100:-1], 16.666, rtol=0.005)
    assert np.allclose(Vstart["voltage_0"][-100:-1], -16.666, rtol=0.005)


@pytest.mark.wires
@pytest.mark.nodal_source
@pytest.mark.probes
def test_current_generators_without_resistance(tmp_path):
    """Verify ideal-wire current-source sign and magnitude at each position."""
    # Checks current probes at the extremes of a wire
    # with a current generator in the middle of the wire and on the extremes of the wire

    fn = CASES_FOLDER + "sources/sources_current_no_resistance.fdtd.json"
    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["sources"][0]["elementIds"] = [1]  # wire center
    solver.cleanUp()
    solver.run()
    Iend = Probe(_get_solved_probe_folder(solver, "probe_end"))
    Istart = Probe(_get_solved_probe_folder(solver, "probe_start"))
    assert np.allclose(Iend["current_0"][-100:-1], 1.0, rtol=0.005)
    assert np.allclose(Istart["current_0"][-100:-1], 1.0, rtol=0.005)

    solver["sources"][0]["elementIds"] = [9]  # wire start
    solver.cleanUp()
    solver.run()
    Iend = Probe(_get_solved_probe_folder(solver, "probe_end"))
    Istart = Probe(_get_solved_probe_folder(solver, "probe_start"))

    assert np.allclose(Iend["current_0"][-100:-1], 1.0, rtol=0.005)
    assert np.allclose(Istart["current_0"][-100:-1], 1.0, rtol=0.005)

    solver["sources"][0]["elementIds"] = [10]  # wire end
    solver.cleanUp()
    solver.run()
    Iend = Probe(_get_solved_probe_folder(solver, "probe_end"))
    Istart = Probe(_get_solved_probe_folder(solver, "probe_start"))

    assert np.allclose(Iend["current_0"][-100:-1], -1.0, rtol=0.005)
    assert np.allclose(Istart["current_0"][-100:-1], -1.0, rtol=0.005)


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.nodal_source
@pytest.mark.probes
def test_voltage_generators(tmp_path):
    """Verify MTLN voltage-source currents and voltages at wire endpoints."""
    # Checks current and voltage of probes at the extremes of bundle (1 conductor + 1 shield)
    # with a voltage generator in the middle of the inner conductor
    fn = CASES_FOLDER + "sources/sources_voltage.fdtd.json"
    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["sources"][0]["elementIds"] = [1]
    solver.cleanUp()
    solver.run()

    Iend = Probe(_get_solved_probe_folder(solver, "probe_end", contains="_I_"))
    Istart = Probe(_get_solved_probe_folder(solver, "probe_start", contains="_I_"))
    Vend = Probe(_get_solved_probe_folder(solver, "probe_end", contains="_V_"))
    Vstart = Probe(_get_solved_probe_folder(solver, "probe_start", contains="_V_"))

    assert np.allclose(Iend["current_0"][-100:-1], 0.0, rtol=0.005)
    assert np.allclose(Istart["current_0"][-100:-1], 0.0, rtol=0.005)
    assert np.allclose(Vend["voltage_0"][-100:-1], 0.0, rtol=0.005)
    assert np.allclose(Vstart["voltage_0"][-100:-1], 0.0, rtol=0.005)

    assert np.allclose(Iend["current_1"][-100:-1], 1.0 / 3.0, rtol=0.005)
    assert np.allclose(Istart["current_1"][-100:-1], 1.0 / 3.0, rtol=0.005)
    assert np.allclose(Vend["voltage_1"][-100:-1], -16.666, rtol=0.005)
    assert np.allclose(Vstart["voltage_1"][-100:-1], -16.666, rtol=0.005)


@pytest.mark.probes
def test_bulk_current_outputs(tmp_path):
    """Verify bulk-current output folders and probe directions for all shapes."""
    # This test uses bulk_probe_cases_over_nodal_source.fdtd from input_examples as input.
    # Verifies all kind of bulk probes are recognised and setted properly by checking outputFile format.
    fn = PROBES_INPUT_EXAMPLE + "bulk_probe_cases_over_nodal_source.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    bulkXPlaneFiles = solver.getSolvedProbeFolders("BulkXPlane")
    bulkYPlaneFiles = solver.getSolvedProbeFolders("BulkYPlane")
    bulkZPlaneFiles = solver.getSolvedProbeFolders("BulkZPlane")
    bulkYPointFiles = solver.getSolvedProbeFolders("BulkYPoint")
    bulkZVolumeFiles = solver.getSolvedProbeFolders("BulkZVolume")

    assert len(bulkXPlaneFiles) == 1
    assert len(bulkYPlaneFiles) == 1
    assert len(bulkZPlaneFiles) == 1
    assert len(bulkYPointFiles) == 1
    assert len(bulkZVolumeFiles) == 10

    probeBulkXPlane = Probe(_get_solved_probe_folder(solver, "BulkXPlane"))
    probeBulkYPlane = Probe(_get_solved_probe_folder(solver, "BulkYPlane"))
    probeBulkZPlane = Probe(_get_solved_probe_folder(solver, "BulkZPlane"))
    probeBulkYPoint = Probe(_get_solved_probe_folder(solver, "BulkYPoint"))
    probeBulkZVolumes = [Probe(path) for path in bulkZVolumeFiles]

    assert probeBulkXPlane.direction == "x"
    assert probeBulkYPlane.direction == "y"
    assert probeBulkZPlane.direction == "z"
    assert probeBulkYPoint.direction == "y"
    assert all(probe.direction == "z" for probe in probeBulkZVolumes)


@no_mtln_skip
@pytest.mark.mtln
def test_conductors_forming_y_on_panel_holland_vs_mtln(tmp_path):
    """Verify Holland and MTLN models agree for conductors forming a Y."""
    holland_run = tmp_path / "holland"
    mtln_run = tmp_path / "mtln"
    holland_run.mkdir()
    mtln_run.mkdir()

    holland_solver = FDTD(
        input_filename=CASES_FOLDER
        + "conductors_forming_y_on_panel/Conductors_Holland_50ohm_terminals.fdtd.json",
        path_to_exe=SEMBA_EXE,
        run_in_folder=holland_run,
    )
    holland_solver.run()

    mtln_solver = FDTD(
        input_filename=CASES_FOLDER
        + "conductors_forming_y_on_panel/Conductors_MTLN_50ohm_terminals.fdtd.json",
        path_to_exe=SEMBA_EXE,
        run_in_folder=mtln_run,
    )
    mtln_solver.run()

    holland_yplus = Probe(_get_solved_probe_folder(holland_solver, "curr_yplus"))
    holland_yminus = Probe(_get_solved_probe_folder(holland_solver, "curr_yminus"))
    holland_joined_1 = Probe(
        _get_solved_probe_folder(holland_solver, "curr_joined_1")
    )
    holland_joined_2 = Probe(
        _get_solved_probe_folder(holland_solver, "curr_joined_2")
    )

    mtln_yplus = Probe(_get_solved_probe_folder(mtln_solver, "curr_yplus"))
    mtln_yminus = Probe(_get_solved_probe_folder(mtln_solver, "curr_yminus"))
    mtln_joined = Probe(_get_solved_probe_folder(mtln_solver, "curr_joined"))

    corr_yplus = corrcoef_on_common_time(
        holland_yplus["time"].to_numpy(),
        holland_yplus["current"].to_numpy(),
        mtln_yplus["time"].to_numpy(),
        mtln_yplus["current"].to_numpy(),
    )
    corr_yminus = corrcoef_on_common_time(
        holland_yminus["time"].to_numpy(),
        holland_yminus["current"].to_numpy(),
        mtln_yminus["time"].to_numpy(),
        mtln_yminus["current"].to_numpy(),
    )

    # MTLN joined probe contains two current columns for the two conductors.
    corr_j1 = corrcoef_on_common_time(
        holland_joined_1["time"].to_numpy(),
        holland_joined_1["current"].to_numpy(),
        mtln_joined["time"].to_numpy(),
        mtln_joined["current_0"].to_numpy(),
    )
    corr_j2 = corrcoef_on_common_time(
        holland_joined_2["time"].to_numpy(),
        holland_joined_2["current"].to_numpy(),
        mtln_joined["time"].to_numpy(),
        mtln_joined["current_1"].to_numpy(),
    )

    assert corr_yplus > 0.999
    assert corr_yminus > 0.999
    assert corr_j1 > 0.999
    assert corr_j2 > 0.999
