from utils import *

# def test_paul_8_6_square(tmp_path):
#     case = 'paul_8_6_square'
#     makeCopy(tmp_path, EXCITATIONS_FOLDER+'coaxial_line_paul_8_6_0.25_square.exc')
#     makeCopy(tmp_path, CASE_FOLDER + case + '.fdtd.json')
#     fn = tmp_path._str + '/' + case + '.fdtd.json'

#     solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
#     solver.run()
#     assert solver.hasFinishedSuccessfully() == True

#     p_expected = Probe(OUTPUT_FOLDER+'paul_8_6_square.fdtd_start_voltage_bundle_wire_V_5_5_1.dat')

#     probe_voltage = solver.getSolvedProbeFilenames("start_voltage")[0]
#     probe_current = solver.getSolvedProbeFilenames("end_current")[0]
#     probe_files = [probe_voltage, probe_current]
#     p_solved = Probe(probe_files[0])

#     assert np.allclose(p_expected.df.to_numpy()[:,0:2], p_solved.df.to_numpy()[:,0:2], rtol = 0.01, atol=0.05)


# def test_paul_8_6_triangle(tmp_path):
#     case = 'paul_8_6_triangle'
#     makeCopy(tmp_path, EXCITATIONS_FOLDER+'coaxial_line_paul_8_6_0.05_triangle.exc')
#     makeCopy(tmp_path, CASE_FOLDER + case + '.fdtd.json')
#     fn = tmp_path._str + '/' + case + '.fdtd.json'

#     solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
#     solver.run()
#     assert solver.hasFinishedSuccessfully() == True

#     # p_expected = Probe(OUTPUT_FOLDER+'paul_8_6_triangle.fdtd_start_voltage_bundle_wire_V_5_5_1.dat')

#     probe_voltage = solver.getSolvedProbeFilenames("start_voltage")[0]
#     probe_current = solver.getSolvedProbeFilenames("end_current")[0]
#     probe_files = [probe_voltage, probe_current]
#     p_solved = Probe(probe_files[0])

#     times = [4.0, 5.9, 6.1, 8.0, 10.1, 12.0]
#     expected_voltages = [16.67, 12.5, -12.5, -25.0, 6.25, 12.5]
    
    
#     for t, v in zip(times, expected_voltages):
#         idx = np.argmin(np.abs(p_solved.df.to_numpy()[:, 0]-t*1e-6))
#         assert np.isclose(v, p_solved.df.to_numpy()[idx, 1], atol=0.5)


#     # assert np.allclose(p_expected.df.to_numpy()[:,0:2], p_solved.df.to_numpy()[:,0:2], rtol = 0.01, atol=0.5)

# def test_paul_9_6(tmp_path):
#     case = 'paul_9_6'
#     makeCopy(tmp_path, EXCITATIONS_FOLDER+'2_conductor_line_paul_9_6_pulse.exc')
#     makeCopy(tmp_path, CASE_FOLDER + case + '.fdtd.json')
#     fn = tmp_path._str + '/' + case + '.fdtd.json'

#     solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
#     solver.run()
#     assert solver.hasFinishedSuccessfully() == True

#     # p_expected = Probe(OUTPUT_FOLDER+'paul_8_6_triangle.fdtd_start_voltage_bundle_wire_V_5_5_1.dat')

#     probe_voltage = solver.getSolvedProbeFilenames("start_voltage")[0]
#     probe_current = solver.getSolvedProbeFilenames("end_current")[0]
#     probe_files = [probe_voltage, probe_current]
#     p_solved = Probe(probe_files[0])

#     assert np.allclose(p_expected.df.to_numpy()[:,0:2], p_solved.df.to_numpy()[:,0:2], rtol = 0.01, atol=0.05)


    
    
    
    

    