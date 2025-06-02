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

# fprev1  = Probe(solver.getSolvedProbeFilenames("field_prev_1")[0])
# fprev2  = Probe(solver.getSolvedProbeFilenames("field_prev_2")[0])
# fprev3x  = Probe(solver.getSolvedProbeFilenames("field_prev_3_Ex")[0])
# fprev3z  = Probe(solver.getSolvedProbeFilenames("field_prev_3_Ez")[0])
# fpost1  = Probe(solver.getSolvedProbeFilenames("field_post_1")[0])
# fpost2  = Probe(solver.getSolvedProbeFilenames("field_post_2")[0])
# fpost3  = Probe(solver.getSolvedProbeFilenames("field_post_3")[0])
# fu2  = Probe(solver.getSolvedProbeFilenames("field_u_2")[0])
# fu3  = Probe(solver.getSolvedProbeFilenames("field_u_3")[0])
iprobe_12_10_13  = Probe(solver.getSolvedProbeFilenames("current_12_10_13_x")[0])
iprobe_7_10_12  = Probe(solver.getSolvedProbeFilenames("current_7_10_12_z")[0])


f7_10_10_z   =  Probe(solver.getSolvedProbeFilenames("field_7_10_10_z")[0])
f7_10_11_z   =  Probe(solver.getSolvedProbeFilenames("field_7_10_11_z")[0])
f7_10_12_z   =  Probe(solver.getSolvedProbeFilenames("field_7_10_12_z")[0])
f13_10_10_z  =  Probe(solver.getSolvedProbeFilenames("field_13_10_10_z")[0])
f13_10_11_z  =  Probe(solver.getSolvedProbeFilenames("field_13_10_11_z")[0])
f13_10_12_z  =  Probe(solver.getSolvedProbeFilenames("field_13_10_12_z")[0])

f7_10_13_x  =  Probe(solver.getSolvedProbeFilenames("field_7_10_13_x")[0])
f8_10_13_x  =  Probe(solver.getSolvedProbeFilenames("field_8_10_13_x")[0])
f9_10_13_x  =  Probe(solver.getSolvedProbeFilenames("field_9_10_13_x")[0])
f10_10_13_x  = Probe(solver.getSolvedProbeFilenames("field_10_10_13_x")[0])
f11_10_13_x  = Probe(solver.getSolvedProbeFilenames("field_11_10_13_x")[0])
f12_10_13_x  = Probe(solver.getSolvedProbeFilenames("field_12_10_13_x")[0])

plt.figure()
plt.plot(f7_10_10_z['time']*1e9, f7_10_10_z['field']*0.1, '-', label = '(7,10,10)z')
plt.plot(f7_10_11_z['time']*1e9, f7_10_11_z['field']*0.1, '-', label = '(7,10,11)z')
plt.plot(f7_10_12_z['time']*1e9, f7_10_12_z['field']*0.1, '-', label = '(7,10,12)z')
# plt.plot(f13_10_10_z['time']*1e9, f13_10_10_z['field']*0.1, '-', label = '(13,10,10)z')
# plt.plot(f13_10_11_z['time']*1e9, f13_10_11_z['field']*0.1, '-', label = '(13,10,11)z')
# plt.plot(f13_10_12_z['time']*1e9, f13_10_12_z['field']*0.1, '-', label = '(13,10,12)z')

plt.xlabel('t [ns]')
plt.ylabel('E*dl [V]')
plt.grid(which='both')
plt.legend()

plt.figure()
plt.plot(f7_10_13_x['time']*1e9, f7_10_13_x['field']*0.1, '--', label = '(7,10,13)x')
plt.plot(f8_10_13_x['time']*1e9, f8_10_13_x['field']*0.1, '--', label = '(8,10,13)x')
plt.plot(f9_10_13_x['time']*1e9, f9_10_13_x['field']*0.1, '--', label = '(9,10,13)x')
plt.plot(f10_10_13_x['time']*1e9, f10_10_13_x['field']*0.1, '--', label = '(10,10,13)x')
plt.plot(f11_10_13_x['time']*1e9, f11_10_13_x['field']*0.1, '--', label = '(11,10,13)x')
plt.plot(f12_10_13_x['time']*1e9, f12_10_13_x['field']*0.1, '--', label = '(12,10,13)x')

# plt.plot(f8_10_13_x['time']*1e9, (f8_10_13_x['field']+f9_10_13_x['field']+f10_10_13_x['field']+f11_10_13_x['field'])*0.1, '--', label = '(8-11,10,13)x')

plt.xlabel('t [ns]')
plt.ylabel('E*dl [V]')
plt.grid(which='both')
plt.legend()

plt.figure()
plt.plot(iprobe_12_10_13['time']*1e9, iprobe_12_10_13['current'], '--', label = 'I(12,10,13)x')
plt.plot(iprobe_7_10_12['time']*1e9, iprobe_7_10_12['current'], '--', label = 'I(7,10,12)z')
plt.xlabel('t [ns]')
plt.ylabel('I [A]')
plt.grid(which='both')
plt.legend()


# %%
