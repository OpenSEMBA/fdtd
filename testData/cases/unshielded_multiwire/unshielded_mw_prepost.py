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
fn = 'unshielded_multiwire.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()


# %% Plot results

# plt.figure()
# plt.plot(t, f/(-1/iSGBC40_first-50),'--', label='SGBC cond = 40S/m, top')
# plt.plot(t, f/(-1/iSGBC40_last-50),'-.g', label='SGBC cond = 40S/m, bottom')
# plt.grid(which='both')
# plt.legend()
# plt.xlabel('Time (ns)')
# plt.ylabel('Conductivity (S/m)')
# plt.show()


