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
t0 = 2.5e-9
t1 = 30e-9
tf = 75e-9
A = 50.0
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
np.savetxt('voltage_source_50V.exc', data)


#####################################################
# %% Run solver
fn = 'sources_current.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver['materials'][0] = {"id":1,
                         "type": "wire",
                         "radius": 0.1e-3, 
                         "resistancePerMeter": 0.0, 
                         "inductancePerMeter": 0.0
                         }
solver["general"]["timeStep"] = 1e-11
solver["general"]["numberOfSteps"] = 7500

solver["materialAssociations"][1] = {"materialId": 4, "elementIds": [3]}
# solver["materialAssociations"][1] = {"name": "wire2","materialId": 1,"initialTerminalId": 3,"endTerminalId": 3,"elementIds": [7]}



solver["sources"][0]["elementIds"] = [1]
solver["sources"][0]["field"] = "voltage"
solver["sources"][0]["magnitudeFile"] = "voltage_source_50V.exc"
solver["mesh"]["elements"][0]["coordinateIds"] = [5]

solver["probes"][0]["elementIds"] = [4]
solver["probes"][0]["name"] = "probe_c8"
solver["probes"][0]["field"] = "voltage"
solver["probes"][0]["type"] = "wire"
solver["probes"][0]["domain"] = {"type" : "time"}

solver["probes"][1]["elementIds"] = [6]
solver["probes"][1]["name"] = "probe_end"
solver["probes"][1]["field"] = "voltage"
solver["probes"][1]["type"] = "wire"
solver["probes"][1]["domain"] = {"type" : "time"}

solver["probes"][2]["elementIds"] = [8]
solver["probes"][2]["name"] = "probe_c7"
solver["probes"][2]["field"] = "current"
solver["probes"][2]["type"] = "wire"
solver["probes"][2]["domain"] = {"type" : "time"}

solver["mesh"]["elements"][3]["coordinateIds"] = [8]
solver["mesh"]["elements"][5]["coordinateIds"] = [4]

# solver["probes"][2]["elementIds"] = [6]
# solver["probes"][2]["name"] = "probe_c8"
# solver["probes"][3]["elementIds"] = [6]
# solver["probes"][3]["name"] = "probe_c8"
# solver["mesh"]["elements"][5]["coordinateIds"] = [8]

solver.cleanUp()
solver.run()
#%%
probe_names = solver.getSolvedProbeFilenames("probe_end")
# p7I = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])
pendV = Probe(list(filter(lambda x: '_V_' in x, probe_names))[0])
probe_names = solver.getSolvedProbeFilenames("probe_c8")
# p8I = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])
p8V = Probe(list(filter(lambda x: '_V_' in x, probe_names))[0])
probe_names = solver.getSolvedProbeFilenames("probe_c7")
p7I = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])


#####################################################
# %% Plot results
# pf = "shieldedPair.fdtd_wire_start_bundle_line_0_I_75_74_74.dat"
# expected = Probe(OUTPUTS_FOLDER+pf)

plt.figure()
plt.plot(p7I['time']*1e9, p7I['current_0'], '.',label = 'probe7 (before source)')
# plt.plot(p8I['time']*1e9, p8I['current_0'], '--',label = 'probe8 (after source)')
plt.grid(which='both')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('I (A)')
plt.show()

plt.figure()
plt.plot(pendV['time']*1e9, pendV['voltage_0'], '.',label = 'line end')
plt.plot(p8V['time']*1e9, p8V['voltage_0'], '--',label = 'probe8 (after source)')
plt.plot(data[:,0]*1e9, data[:,1], '--',label = 'source')
# plt.plot(pendV['time']*1e9, 50*p8I['current_0'], '--',label = 'R*source I')
plt.grid(which='both')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('V (V)')
plt.show()
# %% Plot results


# %%
