from utils import *
import pytest

def test_read_wire_probe():
    p = Probe(OUTPUT_FOLDER + 'holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat')
        
    assert p.case_name == 'holland1981'
    assert p.name == 'mid_point'
    assert p.type == 'wire'
    assert p.domainType == 'time'
    assert np.all(p.cell == np.array([11, 11, 12]))
    assert p.segment == 2
    
    assert len(p['time']) == 1001
    assert p['time'][0] == 0.0
    assert p['time'].iat[-1] == 0.2999999901276417E-007
    
    assert len(p['current']) == 1001
    assert p['current'][0] == 0.0
    assert p['current'].iat[-1] == -0.513576742E-004
    
def test_read_probe_from_NFDE():
    p = Probe(OUTPUT_FOLDER + 'fakeCurrentProbe_mid_point_Wz_11_11_11_s2.dat')
    
    assert p.type == 'wire'
    
def test_read_frequency_probe_from_NFDE():
    p = Probe(OUTPUT_FOLDER + 'edelcadfixZ_COR2_log__Wz_21_21_28_s10_df.dat')
    
    assert p.type == 'wire'
    assert p.domainType == 'frequency'
    assert np.all(p.cell == np.array([21, 21, 28]))
    assert p.segment == 10
    
        
def test_read_point_probe():
    p = Probe(OUTPUT_FOLDER + 'shieldingEffectiveness.fdtd_front_Ex_1_1_1.dat')
        
    assert p.case_name == 'shieldingEffectiveness'
    assert p.name == 'front'
    assert p.type == 'point'
    assert p.domainType == 'time'
    assert p.direction == 'x'
    assert p.field == 'E'
    assert np.all(p.cell == np.array([1, 1, 1]))
        
    assert len(p['time']) == 5193
    assert p['time'][0] == 0.0
    assert np.isclose(p['time'].iat[-1], 0.19997851853637005E-007)
    
    assert len(p['field']) == 5193
    assert p['field'][0] == 0.0
    assert np.isclose(p['field'].iat[-1], 0.120145380E+000)
    
    assert len(p['incident']) == 5193
    assert np.isclose(p['incident'][0], 0.134010895E-005)
    assert     p['incident'].iat[-1] == 0.0
  
  
def test_fdtd_set_new_folder_to_run(tmp_path):
    input = os.path.join(CASE_FOLDER, 'planewave', 'pw-in-box.fdtd.json')
    solver = FDTD(input, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    
    solver.input['general']['numberOfSteps'] = 1
    
    solver.run()
    assert solver.hasFinishedSuccessfully()
    
    
def test_fdtd_clean_up_after_run(tmp_path):
    input = CASE_FOLDER + 'planewave/pw-in-box.fdtd.json'
    solver = FDTD(input, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    
    solver.input['general']['numberOfSteps'] = 1
    
    solver.run()
    assert solver.hasFinishedSuccessfully()
    
    pn = solver.getSolvedProbeFilenames("inbox")
    assert os.path.isfile(pn[0])
    
    solver.cleanUp()
    assert not os.path.isfile(pn[0])
    
    
    