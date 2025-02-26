# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants
from skrf.media import Freespace
from skrf.frequency import Frequency


import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

#####################################################
# %% Run solver
fn = 'zener.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()
#####################################################
# %% Postprocess
import csv
def readCSV(f):
    t, v = np.array([]), np.array([])
    csvFile = csv.reader(f)
    next(csvFile)
    i = 0
    for line in csvFile:
        if (i%5 == 0):
            t = np.append(t, float(line[0]))
            v = np.append(v, float(line[1]))
        i = i + 1
    return t, v

t_meas, v_meas = readCSV(open('./zener_measurement.csv'))
t_meas += 1.9700E-05

bulk = Probe(solver.getSolvedProbeFilenames("end_voltage")[0])
t = bulk['time'][:-1]
v = bulk['voltage_0'][:-1]+0.1

plt.figure()
plt.plot(t*1e6, v,'-', label='simulation')
plt.plot(t_meas*1e6, v_meas, 'k.', label='measurement')
plt.xlabel('Time (us)')
plt.ylabel('Voltage (V)')
plt.grid(which='both')
plt.xlim(0, t.to_numpy()[-1]*1e6)
plt.legend()


for n, time in enumerate(t_meas):
    if (time < t.to_numpy()[-1]):
        idx = (np.absolute(t.to_numpy()-time)).argmin()
        assert np.isclose(t_meas[n], t.to_numpy()[idx], rtol=1e-3)

# %%
