from utils import *
import pytest

@pytest.mark.mtln
def test_shieldedPair(tmp_path):
    fn = CASE_FOLDER + 'shieldedPair/shieldedPair.fdtd.json'
    solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, \
        flags = ['-mtlnwires'], run_in_folder=tmp_path)
    solver.run()
    assert solver.hasFinishedSuccessfully() == True
    
    probe_files = ['shieldedPair.fdtd_wire_start_bundle_line_0_V_75_74_74.dat',
                   'shieldedPair.fdtd_wire_start_bundle_line_0_I_75_74_74.dat',
                   'shieldedPair.fdtd_wire_start_Wz_75_74_74_s4.dat',
                   'shieldedPair.fdtd_wire_end_Wz_75_71_74_s1.dat',
                   'shieldedPair.fdtd_wire_end_bundle_line_0_I_75_71_74.dat',
                   'shieldedPair.fdtd_wire_end_bundle_line_0_V_75_71_74.dat']

    p_expected = [Probe(OUTPUT_FOLDER+'shieldedPair.fdtd_wire_start_bundle_line_0_V_75_74_74.dat'),
                  Probe(OUTPUT_FOLDER+'shieldedPair.fdtd_wire_start_bundle_line_0_I_75_74_74.dat'),
                  Probe(OUTPUT_FOLDER+'shieldedPair.fdtd_wire_start_Wz_75_74_74_s4.dat'),
                  Probe(OUTPUT_FOLDER+'shieldedPair.fdtd_wire_end_Wz_75_71_74_s1.dat'),
                  Probe(OUTPUT_FOLDER+'shieldedPair.fdtd_wire_end_bundle_line_0_I_75_71_74.dat'),
                  Probe(OUTPUT_FOLDER+'shieldedPair.fdtd_wire_end_bundle_line_0_V_75_71_74.dat')]
    
    
    for i in [2,3]:
        p_solved = Probe(probe_files[i])
        assert np.allclose(p_expected[i].df.to_numpy()[:,0:3], p_solved.df.to_numpy()[:,0:3], rtol = 5e-2, atol=0.2)
    for i in [0,1,4,5]:
        p_solved = Probe(probe_files[i])
        assert np.allclose(p_expected[i].df.to_numpy()[:,0:4], p_solved.df.to_numpy()[:,0:4], rtol = 5e-2, atol=0.2)

@pytest.mark.mtln
def test_coated_antenna(tmp_path):
    fn = CASE_FOLDER + 'coated_antenna/coated_antenna.fdtd.json'

    solver = FDTD(\
        input_filename = fn, path_to_exe=SEMBA_EXE, \
        flags = ['-mtlnwires'], run_in_folder=tmp_path)
    solver.run()
    assert solver.hasFinishedSuccessfully() == True
    
    probe_current = solver.getSolvedProbeFilenames("mid_point_Wz")[0]
    probe_files = [probe_current]
    
    p_expected = Probe(OUTPUT_FOLDER+'coated_antenna.fdtd_mid_point_Wz_11_11_11_s2.dat')
    
    p_solved = Probe(probe_files[0])
    assert np.allclose(p_expected.df.to_numpy()[:,0], p_solved.df.to_numpy()[:,0], rtol = 0.0 , atol=10e-8)
    assert np.allclose(p_expected.df.to_numpy()[:,1], p_solved.df.to_numpy()[:,1], rtol = 0.0 , atol=10e-8)
    assert np.allclose(p_expected.df.to_numpy()[:,2], p_solved.df.to_numpy()[:,2], rtol = 0.0 , atol=10e-6)
    

def test_holland(tmp_path):
    fn = CASE_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    
    solver.run()
    assert solver.hasFinishedSuccessfully() == True
    
    probe_current = solver.getSolvedProbeFilenames("mid_point_Wz")[0]
    probe_files = [probe_current]
    
    p_expected = Probe(OUTPUT_FOLDER+'holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat')
    
    p_solved = Probe(probe_files[0])
    assert np.allclose(p_expected.df.to_numpy()[:,0:3], p_solved.df.to_numpy()[:,0:3], rtol = 1e-5, atol=1e-6)

@pytest.mark.mtln
def test_holland_mtln(tmp_path):
    fn = CASE_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags = ['-mtlnwires'], run_in_folder=tmp_path)
    
    solver.run()
    assert solver.hasFinishedSuccessfully() == True
    
    probe_current = solver.getSolvedProbeFilenames("mid_point_Wz")[0]
    probe_files = [probe_current]
    
    p_expected = Probe(OUTPUT_FOLDER+'holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat')
    
    p_solved = Probe(probe_files[0])
    assert np.allclose(p_expected.df.to_numpy()[:,0:3], p_solved.df.to_numpy()[:,0:3], rtol = 1e-5, atol=1e-6)


def test_towelHanger(tmp_path):
    fn = CASE_FOLDER + 'towelHanger/towelHanger.fdtd.json'
    solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE,\
        run_in_folder=tmp_path)
    solver.run()
    assert solver.hasFinishedSuccessfully() == True
    
    probe_files = [solver.getSolvedProbeFilenames("wire_start")[0], 
                   solver.getSolvedProbeFilenames("wire_mid")[0], 
                   solver.getSolvedProbeFilenames("wire_end")[0]]
    
    p_expected = [Probe(OUTPUT_FOLDER+'towelHanger.fdtd_wire_start_Wz_27_25_30_s1.dat'),
                  Probe(OUTPUT_FOLDER+'towelHanger.fdtd_wire_mid_Wx_35_25_32_s5.dat'),
                  Probe(OUTPUT_FOLDER+'towelHanger.fdtd_wire_end_Wz_43_25_30_s4.dat')]
    
    for i in range(3):
        p_solved = Probe(probe_files[i])
        assert np.allclose(p_expected[i].df.to_numpy()[:,0:3], p_solved.df.to_numpy()[:,0:3], rtol = 5e-2, atol=5e-2)
   
    
@pytest.mark.hdf5   
def test_sphere(tmp_path):    
    fn = CASE_FOLDER + 'sphere/sphere.fdtd.json'
    solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.input['general']['numberOfSteps'] = 20
    solver.input['probes'][0]['domain']['initialFrequency'] = 1e8
    solver.input['probes'][0]['domain']['finalFrequency'] = 1e9

    solver.run()
    assert solver.hasFinishedSuccessfully()
    
    far_field_probe_files = solver.getSolvedProbeFilenames("Far") # semba-fdtd seems to always use the name Far for "far field" probes.
    assert solver.hasFinishedSuccessfully() == True
    assert len(far_field_probe_files) == 1
    p = Probe(far_field_probe_files[0])
    assert p.type == 'farField'

    electric_field_movie_files = solver.getSolvedProbeFilenames("electric_field_movie")
    assert solver.hasFinishedSuccessfully() == True
    assert len(electric_field_movie_files) == 3
    p = Probe(electric_field_movie_files[0])
    assert p.type == 'movie'


def test_planewave_in_box(tmp_path):    
    fn = CASE_FOLDER + 'planewave/pw-in-box.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    
    solver.run()
    assert solver.hasFinishedSuccessfully()
    
    before = Probe(solver.getSolvedProbeFilenames("before")[0])
    inbox = Probe(solver.getSolvedProbeFilenames("inbox")[0])
    after = Probe(solver.getSolvedProbeFilenames("after")[0])
    
    np.allclose(inbox.df['field'], inbox.df['incident'])
    zeros = np.zeros_like(before.df['field'])
    np.allclose(before.df['field'], zeros)
    np.allclose(after.df['field'], zeros)

def test_planewave_with_periodic_boundaries(tmp_path):    
    fn = CASE_FOLDER + 'planewave/pw-with-periodic.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    
    solver.run()
    assert solver.hasFinishedSuccessfully()
    
    before = Probe(solver.getSolvedProbeFilenames("before")[0])
    inbox = Probe(solver.getSolvedProbeFilenames("inbox")[0])
    after = Probe(solver.getSolvedProbeFilenames("after")[0])
    
    np.allclose(inbox.df['field'], inbox.df['incident'])
    zeros = np.zeros_like(before.df['field'])
    np.allclose(before.df['field'], zeros)
    np.allclose(after.df['field'], zeros)