from utils import *

def test_paul_8_6_square(tmp_path):
    case = 'paul_8_6_square'
    makeCopy(tmp_path, EXCITATIONS_FOLDER+'coaxial_line_paul_8_6_0.25_square.exc')
    makeCopy(tmp_path, CASE_FOLDER + case + '.fdtd.json')
    fn = tmp_path._str + '/' + case + '.fdtd.json'

    solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
    solver.run()
    probe_voltage = solver.getSolvedProbeFilenames("start_voltage")[0]
    probe_current = solver.getSolvedProbeFilenames("end_current")[0]
    probe_files = [probe_voltage, probe_current]
    
    assert solver.hasFinishedSuccessfully() == True

    # p_expected = [Probe(OUTPUT_FOLDER+'paul_8_6_square.fdtd_start_voltage_bundle_wire_V_5_5_1.dat'),
    #               Probe(OUTPUT_FOLDER+'paul_8_6_square.fdtd_end_current_bundle_wire_I_5_5_100.dat')]
    
    p_solved = Probe(probe_files[0])
    start_times = [0.1, 4.1, 6.1, 8.1, 10.1, 12.1, 14.1, 16.1]
    end_times = [3.9, 5.9, 7.9, 9.9, 11.9, 13.9, 15.9, 18.9]
    expected_voltages = [25.0, -12.5, -37.5, -18.75, 18.75, 9.375, -9.375, -4.6875]

    for (start, end, v) in zip(start_times, end_times, expected_voltages):
        aux_time = 0.5*(start*1e-6 + end*1e-6)
        index = np.argmin(np.abs(p_solved.df.to_numpy()[:,0] - aux_time))
        assert np.all(np.isclose(p_solved.df.to_numpy()[index,1], v, atol=0.2))

    
    
    
    

    