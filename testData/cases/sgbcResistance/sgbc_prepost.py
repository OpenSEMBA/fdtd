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
# %% Generate excitation and visualize
dt = 1e-11
tf = 20e-9
tramp = 5e-9
t = np.arange(0, tf, dt)
v = (t/tramp)*(t < tramp) + 1*(t >= tramp)

plt.figure()
plt.plot(t*1e9,v)
plt.xlabel('Time (ns)')

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = v
np.savetxt('ramp.exc', data)
 
#####################################################
# %% Run solver
fn = 'sgbcResistance.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()
#####################################################
# %% Postprocess
bulk = Probe(solver.getSolvedProbeFilenames("Bulk probe")[0])


t = bulk['time']
i = bulk['current']

plt.figure()
plt.plot(t*1e9, -i,'-', label='bulk probe')
plt.xlabel('Time (ns)')
plt.ylabel('Current (A)')
plt.grid(which='both')
plt.legend()

assert np.allclose(i.array[-101:-1], np.ones(100)*i.array[-100], rtol=1e-3)
assert np.allclose(-1/i.array[-101:-1], np.ones(100)*(50+45), rtol=0.05)


# %%
