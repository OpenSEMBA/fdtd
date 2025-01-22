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
fn = 'sgbcOverlapping.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

#sgbc sigma= 20 first in materialAssociations
p = Probe(solver.getSolvedProbeFilenames("Bulk probe")[0])
t = p['time'].to_numpy()
iSGBC40_last = p['current'].to_numpy()

#sgbc sigma= 40 first in materialAssociations
solver['materialAssociations'][1]["materialId"] = 6
solver['materialAssociations'][2]["materialId"] = 2
solver.cleanUp()
solver.run()
iSGBC40_first = Probe(solver.getSolvedProbeFilenames("Bulk probe")[0])['current'].to_numpy()



#restore default values
solver['materialAssociations'][1]["materialId"] = 2
solver['materialAssociations'][2]["materialId"] = 6
solver.cleanUp()
solver.run()

# %% Plot results

l = 216e-3
h = 2e-3
w = 120e-3
f = l/(w*h)

plt.figure()
plt.plot(t, f/(-1/iSGBC40_first-50),'--', label='SGBC cond = 40S/m, top')
plt.plot(t, f/(-1/iSGBC40_last-50),'-.g', label='SGBC cond = 40S/m, bottom')
plt.grid(which='both')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('Conductivity (S/m)')
plt.show()


#####################################################
                                                                                                                                                                                                                                                                                                                                                                                                                                          # %% Run overlap mapvtk case
fn = 'sgbcOverlapping_vtk.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mapvtk'])
solver['general']["numberOfSteps"] = 1

# %%
