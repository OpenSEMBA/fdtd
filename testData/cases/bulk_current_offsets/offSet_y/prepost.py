# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import signal


import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../../build/bin/semba-fdtd'

from pyWrapper import *

# %% Run solver
f = 'offSet_y.fdtd.json'

solver = FDTD(input_filename = f, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()


# %% Visualizing initial values 

I_in = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=1)
time = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=0)

plt.figure()
plt.plot(time, I_in, label='Initial excitation, source 1')     
plt.grid(which='both')
plt.legend()
plt.show()

# %% Visualizing the current on bulk

probe1 = Probe(solver.getSolvedProbeFilenames("BulkCurrent1")[0])
probe2 = Probe(solver.getSolvedProbeFilenames("BulkCurrent2")[0])
probe3 = Probe(solver.getSolvedProbeFilenames("BulkCurrent3")[0])



plt.figure()
plt.plot(probe1["time"].to_numpy(), probe1["current"].to_numpy(), label='Current on bulk1', color='blue')
plt.plot(probe2["time"].to_numpy(), probe2["current"].to_numpy(), '--', label='Current on bulk2', color='red')
plt.plot(probe3["time"].to_numpy(), probe3["current"].to_numpy(), '--', label='Current on bulk3', color='green')
plt.grid(which='both')
plt.legend()
plt.show()

# %%
