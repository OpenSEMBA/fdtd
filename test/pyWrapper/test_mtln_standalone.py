from utils import *


@no_mtln_skip
@pytest.mark.mtln
def test_paul_8_6_square(tmp_path):
    fn = CASES_FOLDER + 'paul/paul_8_6_square.fdtd.json'

    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver.run()
    
    p_expected = Probe(
        OUTPUTS_FOLDER+'paul_8_6_square.fdtd_start_voltage_wire_V_5_5_1.dat')

    probe_voltage = solver.getSolvedProbeFilenames("start_voltage")[0]
    probe_current = solver.getSolvedProbeFilenames("end_current")[0]
    probe_files = [probe_voltage, probe_current]
    p_solved = Probe(probe_files[0])

    assert np.allclose(p_expected.data.to_numpy()[:, 0:2], p_solved.data.to_numpy()[
                       :, 0:2], rtol=0.01, atol=0.2)



@no_mtln_skip
@pytest.mark.mtln
def test_paul_8_6_triangle(tmp_path):
    fn = CASES_FOLDER + 'paul/paul_8_6_triangle.fdtd.json'

    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver.run()
    
    p_expected = Probe(
        OUTPUTS_FOLDER+'paul_8_6_triangle.fdtd_start_voltage_wire_V_5_5_1.dat')

    probe_voltage = solver.getSolvedProbeFilenames("start_voltage")[0]
    probe_current = solver.getSolvedProbeFilenames("end_current")[0]
    probe_files = [probe_voltage, probe_current]
    p_solved = Probe(probe_files[0])

    assert np.allclose(p_expected.data.to_numpy()[:, 0:2], p_solved.data.to_numpy()[
                       :, 0:2], rtol=0.01, atol=0.5)


@no_mtln_skip
@pytest.mark.mtln
def test_paul_9_6(tmp_path):
    fn = CASES_FOLDER + 'paul/paul_9_6.fdtd.json'
    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver.run()
    
    p_expected = [Probe(OUTPUTS_FOLDER+'paul_9_6.fdtd_start_voltage_two_wires_V_5_5_1.dat'),
                  Probe(OUTPUTS_FOLDER+'paul_9_6.fdtd_end_voltage_two_wires_V_5_5_795.dat')]

    probe_voltage_left = solver.getSolvedProbeFilenames(
        "start_voltage")[0]
    probe_voltage_right = solver.getSolvedProbeFilenames("end_voltage")[
        0]
    probe_files = [probe_voltage_left, probe_voltage_right]

    p_solved = [Probe(probe_files[0]), Probe(probe_files[1])]

    for i in range(2):
        assert np.allclose(p_expected[i].data.to_numpy()[
                           :, :], p_solved[i].data.to_numpy()[:, :], rtol=0.01, atol=0.5)


@no_mtln_skip
@pytest.mark.mtln
def test_spice_multilines_opamp(tmp_path):
    fn = CASES_FOLDER + 'multilines_opamp/multilines_opamp.fdtd.json'

    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver.run()
    
    p_expected = [
        Probe(OUTPUTS_FOLDER+'multilines_opamp.fdtd_line_end_s2_V_5_5_102.dat')]

    probe_files = [solver.getSolvedProbeFilenames("line_end")[0]]

    p_solved = [Probe(probe_files[0]), Probe(probe_files[0])]

    assert np.allclose(p_expected[0].data.to_numpy()[
                       :-1, :], p_solved[0].data.to_numpy()[:-1, :], rtol=0.01, atol=0.05e-3)


@no_mtln_skip
@pytest.mark.mtln
def test_spice_connectors_diode(tmp_path):
    fn = CASES_FOLDER + 'spice_connectors/spice_connectors.fdtd.json'

    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver.run()
    
    p_expected = [Probe(OUTPUTS_FOLDER+'spice_connectors.fdtd_start_voltage_wire_V_10_10_8.dat'),
                  Probe(OUTPUTS_FOLDER+'spice_connectors.fdtd_end_voltage_wire_V_10_10_12.dat')]

    probe_voltage_left = solver.getSolvedProbeFilenames(
        "start_voltage_wire")[0]
    probe_voltage_right = solver.getSolvedProbeFilenames(
        "end_voltage_wire")[0]
    probe_files = [probe_voltage_left, probe_voltage_right]

    p_solved = [Probe(probe_files[0]), Probe(probe_files[1])]

    for i in range(2):
        assert np.allclose(p_expected[i].data.to_numpy()[
                           :-20, :], p_solved[i].data.to_numpy()[:-20, :], rtol=0.01, atol=0.05e-3)


@no_mtln_skip
@pytest.mark.mtln
def test_line_multiline_junction(tmp_path):
    fn = CASES_FOLDER + 'line_multiline_junction/line_multiline_junction.fdtd.json'
    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    
    solver.run()
    
    p_expected = [Probe(OUTPUTS_FOLDER+'line_multiline_junction.fdtd_s4_end_s4_V_5_5_159.dat'),
                  Probe(
                      OUTPUTS_FOLDER+'line_multiline_junction.fdtd_s5_end_s5_V_5_5_159.dat'),
                  Probe(OUTPUTS_FOLDER+'line_multiline_junction.fdtd_s2_start_s2_V_5_5_2.dat')]

    probe_s2 = solver.getSolvedProbeFilenames("s2_start")[0]
    probe_s4 = solver.getSolvedProbeFilenames("s4_end")[0]
    probe_s5 = solver.getSolvedProbeFilenames("s5_end")[0]
    probe_files = [probe_s4, probe_s5, probe_s2]

    for i in range(3):
        assert np.allclose(p_expected[i].data.to_numpy()[
                           :-20, :], Probe(probe_files[i]).data.to_numpy()[:-20, :], rtol=0.01, atol=5e-3)


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.codemodel
def test_spice_opamp_saturation(tmp_path):
    fn = CASES_FOLDER + 'opamp_saturation/opamp_saturation.fdtd.json'
    setNgspice(tmp_path)

    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver.run()
    
    p_expected = Probe(
        OUTPUTS_FOLDER+'opamp_saturation.fdtd_opamp_voltage_wire1_V_10_10_7.dat')
    p_solved = Probe(solver.getSolvedProbeFilenames(
        "opamp_voltage_wire1")[0])

    assert np.allclose(p_expected.data.to_numpy()[
                       :-5, :], p_solved.data.to_numpy()[:-5, :], rtol=0.01, atol=0.05e-3)


@no_mtln_skip
@pytest.mark.mtln
def test_spice_zener(tmp_path):
    fn = CASES_FOLDER + 'zener/zener.fdtd.json'

    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver.run()
    
    p_expected = Probe(
        OUTPUTS_FOLDER+'zener.fdtd_end_voltage_wire_V_10_10_12.dat')
    p_solved = Probe(solver.getSolvedProbeFilenames(
        "end_voltage_")[0])

    assert np.allclose(p_expected.data.to_numpy()[
                       :-5, :], p_solved.data.to_numpy()[:-5, :], rtol=0.01, atol=0.05e-3)
    
@no_mtln_skip
@pytest.mark.mtln
def test_mtln_sources(tmp_path):
    fn = CASES_FOLDER + 'sources/sources.fdtd.json'
    # terminal voltage source for wires & unshieldedmws : OK
    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path,
                  flags = ["-n", "1"])
    
    solver["materials"][0] = createWire(id = 1, r = 0.02)
    solver["mesh"]["elements"][3]["coordinateIds"] = [1]
    source = { "type" : "generator", 
              "field" : "voltage", 
              "elementIds" : [4], 
              "name" : "terminal_source_v",
              "magnitudeFile" : "sources.exc"}
    
    solver["sources"][0] = source
    solver.cleanUp()
    solver.run()
    assert solver.hasFinishedSuccessfully()

    # solver["materials"][0] = createUnshieldedWire(id = 1, lpul = 6.52188703e-08, cpul = 1.7060247700000001e-10)        
    # solver.cleanUp()
    # solver.run()
    # assert solver.hasFinishedSuccessfully()

    # # interior voltage source for wires & unshieldedmws : FAIL

    # solver["mesh"]["elements"][3]["coordinateIds"] = [2]
    # solver["materials"][0] = createWire(id = 1, r = 0.02)
    # solver.cleanUp()
    # solver.run()
    # assert (solver.hasFinishedSuccessfully() == False)

    # solver["materials"][0] = createUnshieldedWire(id = 1, lpul = 6.52188703e-08, cpul = 1.7060247700000001e-10)        
    # solver.cleanUp()
    # solver.run()
    # assert (solver.hasFinishedSuccessfully() == False)
    
    # # terminal voltage source for shieldedmws : OK
    # fn = CASES_FOLDER + 'sources/sources_shielded.fdtd.json'
    # solver = FDTD(input_filename=fn,
    #               path_to_exe=SEMBA_EXE,
    #               run_in_folder=tmp_path,
    #               flags = ["-n", "1"])
    
    # solver["materials"][0] = createWire(id = 1, r = 0.02)
    # solver["mesh"]["elements"][3]["coordinateIds"] = [1]
    # solver["sources"][0]["field"] = "voltage"
    # solver["sources"][0]["elemendIds"] = [4]
    # solver["sources"][0]["magnitudeFile"] = "source.exc"
    # solver.cleanUp()
    # solver.run()
    # assert solver.hasFinishedSuccessfully()
    
    # interior voltage source for shieldedmws : OK


    
    # inner voltage source: only for inner conductor of shieldedmws
    # solver = FDTD(input_filename=fn,
    #               path_to_exe=SEMBA_EXE,
    #               run_in_folder=tmp_path)
    # solver["mesh"]["elements"][3]["coordinateIds"] = [2]
    
    # solver["sources"][0]["field"] = "voltage"
    # solver["sources"][0]["elemendIds"] = [1]
    # solver["sources"][0]["magnitudeFile"] = "current_source_1A.exc"
    # solver.cleanUp()
    # solver.run()
    # assert solver.hasFinishedSuccessfully()
    
    # # interior current source
    # solver = FDTD(input_filename=fn,
    #               path_to_exe=SEMBA_EXE,
    #               run_in_folder=tmp_path)
    # solver["sources"][0]["field"] = "current"
    # solver["sources"][0]["elemendIds"] = [1]
    # solver["sources"][0]["magnitudeFile"] = "voltage_source_1V.exc"
    # solver.cleanUp()
    # solver.run()
    # assert solver.hasFinishedSuccessfully()
    
