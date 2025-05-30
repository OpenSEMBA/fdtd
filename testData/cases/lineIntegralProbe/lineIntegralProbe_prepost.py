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
tf = 200e-9
t = np.arange(0, tf, dt)
v = 1/(1 + np.exp(0.5*(-t*1e9+20)))

plt.figure()
plt.plot(t*1e9,v)
plt.xlabel('Time (ns)')

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = v
np.savetxt('lineIntegralProbe.exc', data)

#####################################################
# %%
fn = 'lineIntegralProbe.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

#####################################################
# %% Plot results

fprev1  = Probe(solver.getSolvedProbeFilenames("field_prev_1")[0])
fprev2  = Probe(solver.getSolvedProbeFilenames("field_prev_2")[0])
fprev3x  = Probe(solver.getSolvedProbeFilenames("field_prev_3_Ex")[0])
fprev3z  = Probe(solver.getSolvedProbeFilenames("field_prev_3_Ez")[0])
fpost1  = Probe(solver.getSolvedProbeFilenames("field_post_1")[0])
fpost2  = Probe(solver.getSolvedProbeFilenames("field_post_2")[0])
fpost3  = Probe(solver.getSolvedProbeFilenames("field_post_3")[0])
fu2  = Probe(solver.getSolvedProbeFilenames("field_u_2")[0])
fu3  = Probe(solver.getSolvedProbeFilenames("field_u_3")[0])
iprobe  = Probe(solver.getSolvedProbeFilenames("current")[0])

plt.figure()
# plt.plot(fprev1['time']*1e9, -(fprev1['field']+fprev2['field']+fprev3z['field'])*0.1, '--', label = 'edl from ground to wire, before R')
# plt.plot(fpost1['time']*1e9, -(fpost1['field']+fpost2['field']+fpost3['field'])*0.1, '-', label = 'edl from ground to wire, after R')
plt.plot(fpost1['time']*1e9, (-fprev3z['field']+fprev3x['field']+fu2['field']+fu3['field']+fpost3['field'])*0.1, '-', label = 'edl around R')

plt.xlabel('t [ns]')
plt.ylabel('V [V]')
plt.grid(which='both')
plt.legend()

plt.figure()
plt.plot(iprobe['time']*1e9, iprobe['current'], '--', label = 'current on wire')
plt.xlabel('t [ns]')
plt.ylabel('I [A]')
plt.grid(which='both')
plt.legend()


# %%
