# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build-rls/bin/semba-fdtd'

from pyWrapper import *

#####################################################
# %% Run solver
fn = 'tvs_diode.fdtd.json'
# setNgspice('/')
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
# solver.cleanUp()
# solver.run()

#%%
probe = Probe(solver.getSolvedProbeFilenames("crosstalk")[0])
fig, ax = plt.subplots(1,2, figsize=(10,5))
ax[0].plot(probe['time']*1e6, probe['voltage_0'], '.',label = 'conductor 0')
ax[1].plot(probe['time']*1e6, probe['voltage_1'], '.',label = 'conductor 1')

ax[0].grid(which='both')
ax[0].legend()
ax[0].set_xlabel('Time (us)')
ax[0].set_ylabel('V (V)')

ax[1].grid(which='both')
ax[1].legend()
ax[1].set_xlabel('Time (us)')
ax[1].set_ylabel('V (V)')
plt.show()


# %%
