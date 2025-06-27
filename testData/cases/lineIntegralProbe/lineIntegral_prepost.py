# %%
import numpy as np
import matplotlib.pyplot as plt
import json

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

#####################################################
# %% Generate excitation and visualize
dt = 1e-10
tf = 10e-9
t = np.arange(0, tf, dt)
v = 1/(1 + np.exp(2.5*(-t*1e9+5)))

plt.figure()
plt.plot(t*1e9,v)
plt.xlabel('Time (ns)')

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = v
np.savetxt('lineIntegral.exc', data)

#####################################################
# %%
fn = 'lineIntegral.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

#####################################################
# %% Plot results

