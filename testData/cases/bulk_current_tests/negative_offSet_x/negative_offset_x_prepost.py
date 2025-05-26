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
f = 'offSet_negative_x.fdtd.json'

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

probeL = Probe(solver.getSolvedProbeFilenames("Bulk_left")[0])
probeR = Probe(solver.getSolvedProbeFilenames("Bulk_right")[0])
probeTotal = Probe(solver.getSolvedProbeFilenames("BulkTotal")[0])



plt.figure()
plt.plot(probeTotal["time"].to_numpy(), probeTotal["current"].to_numpy(), label='Total current on bulk', color='red')
plt.plot(probeL["time"].to_numpy(), probeL["current"].to_numpy(), '--', label='Current on left edge of bulk', color='blue')
plt.plot(probeR["time"].to_numpy(), probeR["current"].to_numpy(), '--', label='Current on right edge of bulk', color='green')
plt.grid(which='both')
plt.legend()
plt.show()

# %%
