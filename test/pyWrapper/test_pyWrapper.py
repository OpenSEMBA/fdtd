from utils import *
import pytest

def test_read_wire_probe():
    p = Probe(OUTPUT_FOLDER + 'holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat')
        
    assert p.case_name == 'holland1981'
    assert p.name == 'mid_point'
    assert p.type == 'wire'
    assert np.all(p.cell == np.array([11, 11, 12]))
    assert p.segment_tag == 2
    
    assert len(p['time']) == 1001
    assert p['time'][0] == 0.0
    assert p['time'].iat[-1] == 0.2999999901276417E-007
    
    assert len(p['current']) == 1001
    assert p['current'][0] == 0.0
    assert p['current'].iat[-1] == -0.513576742E-004
    
    
def test_read_point_probe():
    p = Probe(OUTPUT_FOLDER + 'shieldingEffectiveness.fdtd_front_Ex_1_1_1.dat')
        
    assert p.case_name == 'shieldingEffectiveness'
    assert p.name == 'front'
    assert p.type == 'point'
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
  