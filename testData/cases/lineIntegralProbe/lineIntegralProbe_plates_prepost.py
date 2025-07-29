# %%
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import *

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *


#####################################################
# %% Generate excitation and visualize
dt = 0.25e-11
tf = 90e-9

t = np.arange(0, tf, dt)
v = 1/(1 + np.exp(0.3*(-t*1e9+25)))
plt.figure()
plt.plot(t*1e9,v)
plt.legend()
plt.xlabel('Time (ns)')

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = v
np.savetxt('lineIntegralProbe_plates.exc', data)

#####################################################
# %%
fn = 'lineIntegralProbe_plates.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

#####################################################
# %% Plot results
lineInt  = Probe(solver.getSolvedProbeFilenames("vprobe_LI_20_20_10")[0])
lineC = np.interp(t, lineInt['time'], lineInt['lineIntegral'])

plt.figure()
plt.plot(t*1e9, -lineC,  label = 'line integral')
plt.plot(t*1e9,v, '--', label='source')
plt.grid(which='both')
plt.legend()


# %%
