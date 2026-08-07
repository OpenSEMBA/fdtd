# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'
OUTPUTS_FOLDER = '../../outputs/'
CASES_FOLDER = '../../cases/'
from pyWrapper import *
#%%
fn = CASES_FOLDER + 'paul/paul_8_6_triangle.fdtd.json'

solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE)
solver.run()

#%%
p_expected = Probe(
    OUTPUTS_FOLDER+'paul_8_6_triangle.fdtd_start_voltage_wire_V_5_5_1.dat')

probe_voltage = solver.getSolvedProbeFilenames("start_voltage")[0]
probe_current = solver.getSolvedProbeFilenames("end_current")[0]
probe_files = [probe_voltage, probe_current]
p_solved = Probe(probe_files[0])

#%%
plt.plot(p_expected['time'], p_expected['voltage_0'], label = 'reference')
plt.plot(p_solved['time'], p_solved['voltage_0'], '--',label = 'solved')
plt.xlabel('time [s]')
plt.ylabel('V [V]')
plt.legend()
#%%
solved = np.interp(p_expected['time'].to_numpy(), 
                    p_solved['time'].to_numpy(), 
                    p_solved['voltage_0'].to_numpy())
assert np.corrcoef(solved, p_expected['voltage_0'])[0,1] > 0.999

# %%
