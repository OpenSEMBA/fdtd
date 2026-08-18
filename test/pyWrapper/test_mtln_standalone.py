from test.utils.utils import *


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.probes
def test_paul_8_6_square(tmp_path):
    fn = CASES_FOLDER + "paul/paul_8_6_square.fdtd.json"

    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p_expected = probe_from_fixture(
        tmp_path, "paul_8_6_square.fdtd_start_voltage_wire_V_5_5_1.dat"
    )

    probe_voltage = solver.getSolvedProbeFolders("start_voltage")[0]
    probe_current = solver.getSolvedProbeFolders("end_current")[0]
    probe_files = [probe_voltage, probe_current]
    p_solved = Probe(probe_files[0])

    solved = np.interp(
        p_expected["time"].to_numpy(),
        p_solved["time"].to_numpy(),
        p_solved["voltage_0"].to_numpy(),
    )
    assert np.corrcoef(solved, p_expected["voltage_0"])[0, 1] > 0.999


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.probes
def test_paul_8_6_triangle(tmp_path):
    fn = CASES_FOLDER + "paul/paul_8_6_triangle.fdtd.json"

    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p_expected = probe_from_fixture(
        tmp_path, "paul_8_6_triangle.fdtd_start_voltage_wire_V_5_5_1.dat"
    )

    probe_voltage = solver.getSolvedProbeFolders("start_voltage")[0]
    probe_current = solver.getSolvedProbeFolders("end_current")[0]
    probe_files = [probe_voltage, probe_current]
    p_solved = Probe(probe_files[0])

    solved = np.interp(
        p_expected["time"].to_numpy(),
        p_solved["time"].to_numpy(),
        p_solved["voltage_0"].to_numpy(),
    )
    assert np.corrcoef(solved, p_expected["voltage_0"])[0, 1] > 0.999


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.probes
def test_paul_9_6(tmp_path):
    fn = CASES_FOLDER + "paul/paul_9_6.fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p_expected = [
        probe_from_fixture(tmp_path, "paul_9_6.fdtd_start_voltage_two_wires_V_5_5_1.dat"),
        probe_from_fixture(tmp_path, "paul_9_6.fdtd_end_voltage_two_wires_V_5_5_795.dat"),
    ]

    probe_voltage_left = solver.getSolvedProbeFolders("start_voltage")[0]
    probe_voltage_right = solver.getSolvedProbeFolders("end_voltage")[0]
    probe_files = [probe_voltage_left, probe_voltage_right]

    p_solved = [Probe(probe_files[0]), Probe(probe_files[1])]

    for i in range(2):
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


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.spice
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.probes
def test_spice_multilines_opamp(tmp_path):
    fn = CASES_FOLDER + "multilines_opamp/multilines_opamp.fdtd.json"

    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p_expected = [
        probe_from_fixture(tmp_path, "multilines_opamp.fdtd_line_end_s2_V_5_5_102.dat")
    ]

    probe_files = [solver.getSolvedProbeFolders("line_end")[0]]

    p_solved = [Probe(probe_files[0]), Probe(probe_files[0])]

    solved = np.interp(
        p_expected[0]["time"].to_numpy(),
        p_solved[0]["time"].to_numpy(),
        p_solved[0]["voltage_0"].to_numpy(),
    )
    assert np.corrcoef(solved, p_expected[0]["voltage_0"])[0, 1] > 0.999


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.spice
@pytest.mark.wires
@pytest.mark.probes
def test_spice_connectors_diode(tmp_path):
    fn = CASES_FOLDER + "spice_connectors/spice_connectors.fdtd.json"

    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p_expected = [
        probe_from_fixture(tmp_path, "spice_connectors.fdtd_start_voltage_wire_V_10_10_8.dat"),
        probe_from_fixture(tmp_path, "spice_connectors.fdtd_end_voltage_wire_V_10_10_12.dat"),
    ]

    probe_voltage_left = solver.getSolvedProbeFolders("start_voltage_wire")[0]
    probe_voltage_right = solver.getSolvedProbeFolders("end_voltage_wire")[0]
    probe_files = [probe_voltage_left, probe_voltage_right]

    p_solved = [Probe(probe_files[0]), Probe(probe_files[1])]

    for i in range(2):
        t_exp = p_expected[i].data["time"].to_numpy()[:-1]
        t_sol = p_solved[i].data["time"].to_numpy()[:-1]
        v_exp = p_expected[i].data["voltage_0"].to_numpy()[:-1]
        v_sol = p_solved[i].data["voltage_0"].to_numpy()[:-1]
        v_sol_interp = np.interp(t_exp, t_sol, v_sol)
        assert np.corrcoef(v_exp, v_sol_interp)[0, 1] > 0.99999


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.probes
def test_line_multiline_junction(tmp_path):
    fn = CASES_FOLDER + "line_multiline_junction/line_multiline_junction.fdtd.json"
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    p_expected = [
        probe_from_fixture(tmp_path, "line_multiline_junction.fdtd_s4_end_s4_V_5_5_159.dat"),
        probe_from_fixture(tmp_path, "line_multiline_junction.fdtd_s5_end_s5_V_5_5_159.dat"),
        probe_from_fixture(tmp_path, "line_multiline_junction.fdtd_s2_start_s2_V_5_5_2.dat"),
    ]

    probe_s2 = solver.getSolvedProbeFolders("s2_start")[0]
    probe_s4 = solver.getSolvedProbeFolders("s4_end")[0]
    probe_s5 = solver.getSolvedProbeFolders("s5_end")[0]
    p_solved = [Probe(probe_s4), Probe(probe_s5), Probe(probe_s2)]

    for i in range(3):
        solved = np.interp(
            p_expected[i]["time"].to_numpy(),
            p_solved[i]["time"].to_numpy(),
            p_solved[i]["voltage_0"].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[i]["voltage_0"])[0, 1] > 0.999


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.codemodel
@pytest.mark.spice
@pytest.mark.wires
@pytest.mark.probes
def test_spice_opamp_saturation(tmp_path):
    fn = CASES_FOLDER + "opamp_saturation/opamp_saturation.fdtd.json"
    setNgspice(tmp_path)

    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p_expected = probe_from_fixture(
        tmp_path, "opamp_saturation.fdtd_opamp_voltage_wire1_V_10_10_7.dat"
    )
    p_solved = Probe(solver.getSolvedProbeFolders("opamp_voltage_wire1")[0])

    t_exp = p_expected.data["time"].to_numpy()[:-1]
    t_sol = p_solved.data["time"].to_numpy()[:-1]
    v_exp = p_expected.data["voltage_0"].to_numpy()[:-1]
    v_sol = p_solved.data["voltage_0"].to_numpy()[:-1]
    v_sol_interp = np.interp(t_exp, t_sol, v_sol)
    assert np.corrcoef(v_exp, v_sol_interp)[0, 1] > 0.999


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.spice
@pytest.mark.wires
@pytest.mark.probes
def test_spice_zener(tmp_path):
    fn = CASES_FOLDER + "zener/zener.fdtd.json"

    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p_expected = probe_from_fixture(tmp_path, "zener.fdtd_end_voltage_wire_V_10_10_12.dat")
    t_exp = p_expected.data["time"].to_numpy()[:-1]
    v_exp = p_expected.data["voltage_0"].to_numpy()[:-1]

    p_solved = Probe(solver.getSolvedProbeFolders("end_voltage_")[0])
    t_sol = p_solved.data["time"].to_numpy()[:-1]
    v_sol = p_solved.data["voltage_0"].to_numpy()[:-1]

    v_exp_interp = np.interp(t_sol, t_exp, v_exp)
    assert np.corrcoef(v_sol, v_exp_interp)[0, 1] > 0.999
