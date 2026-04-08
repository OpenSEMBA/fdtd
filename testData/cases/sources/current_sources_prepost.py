# %%
import numpy as np
import matplotlib.pyplot as plt
import json

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build-rls/bin/semba-fdtd'
OUTPUTS_FOLDER = '../../outputs/'

from pyWrapper import *


#####################################################
# %% Generate excitation and visualize
dt = 3e-11
t0 = 1e-9
t1 = 30e-9
tf = 75e-9
A = 1.0
t = np.arange(0, tf, dt)
e = np.zeros(len(t))

for i, tv in enumerate(t):
    if (tv<t0):
        e[i] = 0.0
    elif (tv>=t0 and tv<=t1):
        e[i] = (tv-t0)*(A/(t1-t0))
    elif (tv>t1):
        e[i] = A

plt.figure()
plt.plot(t*1e9,e)
plt.xlabel('Time (ns)')

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = e
np.savetxt('current_source_1A.exc', data)


#####################################################
# %% Run solver
fn = 'sources_current.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver['materials'][0] = {"id":1,
                         "type": "unshieldedMultiwire",
                         "inductancePerMeter": [[6.5138e-7]],
                         "capacitancePerMeter": [[1.7069e-11]]
                         }
solver["general"]["timeStep"] = 1e-11
solver["general"]["numberOfSteps"] = 6000


solver.cleanUp()
solver.run()
#%%
probe_names = solver.getSolvedProbeFilenames("probe_end")
pendI = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

probe_names = solver.getSolvedProbeFilenames("probe_start")
pstI = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

probe_names = solver.getSolvedProbeFilenames("probe_end")
pendV = Probe(list(filter(lambda x: '_V_' in x, probe_names))[0])

probe_names = solver.getSolvedProbeFilenames("probe_start")
pstV = Probe(list(filter(lambda x: '_V_' in x, probe_names))[0])



#####################################################
# %% Plot results
plt.figure()
plt.plot(pstI['time']*1e9, pstI['current_0'], '*',label = 'probe start')
plt.plot(pendI['time']*1e9, pendI['current_0'], '--',label = 'probe end')
plt.grid(which='both')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('I (A)')
plt.show()
#%%
plt.figure()
plt.plot(pstV['time']*1e9, pstV['voltage_0'], '*',label = 'probe start')
plt.plot(pendV['time']*1e9, pendV['voltage_0'], '--',label = 'probe end')
plt.grid(which='both')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('V (V)')
plt.show()

# %%
