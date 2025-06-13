# %%
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import *

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

#####################################################
# %% Generate excitation and visualize
dt = 1e-10
tf = 600e-9
t = np.arange(0, tf, dt)
v = 1/(1 + np.exp(0.03*(-t*1e9+250)))
v2 = 0*(t<=100e-9) + ((t*1e9-100)/300) * (t>=100e-9)*(t<=400e-9) + 1*(t>=400e-9)
# vder = -t*1e9*0.5*np.exp(0.5*(-t*1e9+20))/(1 + np.exp(0.5*(-t*1e9+20)))**2
plt.figure()
plt.plot(t*1e9,v)
# plt.plot(t*1e9,v2)
# plt.plot(t*1e9,vder, label='der')
# plt.plot(t*1e9,v-0.15*vder, label='sum')
plt.legend()
plt.xlabel('Time (ns)')

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = v
np.savetxt('lineIntegralProbe.exc', data)

#####################################################

#%%
dt = t[1] - t[0]
fq  = fftfreq(len(t))/dt
vf1 = fft(v)
vf2 = fft(v2)
fmin = 1e6
fmax = 1e9
idx_min = (np.abs(fq - fmin)).argmin()
idx_max = (np.abs(fq - fmax)).argmin()
fv = fq[idx_min:idx_max]
vff1 = vf1[idx_min:idx_max]
# vff2 = vf2[idx_min:idx_max]
plt.figure()
plt.loglog((fv), np.abs(vff1),'-*', label='1')
# plt.loglog((fv), np.abs(vff2)/np.max(np.abs(vff2)),'-*', label='2')
plt.legend()
plt.show()

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
iprobe_7_10_11  = Probe(solver.getSolvedProbeFilenames("current_7_10_11_z")[0])


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

f7_10_10_x  =  Probe(solver.getSolvedProbeFilenames("field_7_10_10_x")[0])
f8_10_10_x  =  Probe(solver.getSolvedProbeFilenames("field_8_10_10_x")[0])
f9_10_10_x  =  Probe(solver.getSolvedProbeFilenames("field_9_10_10_x")[0])
f10_10_10_x  = Probe(solver.getSolvedProbeFilenames("field_10_10_10_x")[0])
f11_10_10_x  = Probe(solver.getSolvedProbeFilenames("field_11_10_10_x")[0])
f12_10_10_x  = Probe(solver.getSolvedProbeFilenames("field_12_10_10_x")[0])

fieldSum =  f7_10_10_z['field'] + f7_10_11_z['field'] + f7_10_12_z['field']\
          + f7_10_13_x['field'] + f8_10_13_x['field'] + f9_10_13_x['field']\
          + f10_10_13_x['field']+ f11_10_13_x['field']+ f12_10_13_x['field']\
          + f13_10_10_z['field']+ f13_10_11_z['field'] + f13_10_12_z['field']\
        #   + f7_10_10_x['field']+  f8_10_10_x['field']+ f9_10_10_x['field']\
        #   + f10_10_10_x['field']+ f11_10_10_x['field']+ f12_10_10_x['field']\

groundField = f7_10_10_x['field']+  f8_10_10_x['field']+ f9_10_10_x['field']\
            + f10_10_10_x['field']+ f11_10_10_x['field']+ f12_10_10_x['field']\

plt.figure()
plt.plot(f7_10_10_z['time']*1e9, fieldSum*0.1, '-', label = 'field Sum')
# plt.plot(f7_10_10_z['time']*1e9, groundField, '-', label = 'ground field')
plt.grid(which='both')
plt.legend()


plt.figure()
plt.plot(f7_10_10_z['time']*1e9, f7_10_10_z['field']*0.1, '-', label = '(7,10,10)z')
plt.plot(f7_10_11_z['time']*1e9, f7_10_11_z['field']*0.1, '-', label = '(7,10,11)z')
plt.plot(f7_10_12_z['time']*1e9, f7_10_12_z['field']*0.1, '-', label = '(7,10,12)z')
plt.plot(f13_10_10_z['time']*1e9, f13_10_10_z['field']*0.1, '-', label = '(13,10,10)z')
plt.plot(f13_10_11_z['time']*1e9, f13_10_11_z['field']*0.1, '-', label = '(13,10,11)z')
plt.plot(f13_10_12_z['time']*1e9, f13_10_12_z['field']*0.1, '-', label = '(13,10,12)z')

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
plt.plot(iprobe_7_10_12['time']*1e9, -iprobe_7_10_12['current'], '--', label = 'I(7,10,12)z')
plt.plot(iprobe_7_10_11['time']*1e9, iprobe_7_10_11['current'], '-.', label = 'I(7,10,11)z')
plt.axhline(y=1/(1*150), color='r', linestyle='-')
plt.xlabel('t [ns]')
plt.ylabel('I [A]')
plt.grid(which='both')
plt.legend()

# %% 
dt = iprobe_7_10_11['time'][1] - iprobe_7_10_11['time'][0]
fq  = fftfreq(len(iprobe_7_10_11['time']))/dt
If = fft(iprobe_7_10_11['current'])
fmin = 1e6
fmax = 1e9
idx_min = (np.abs(fq - fmin)).argmin()
idx_max = (np.abs(fq - fmax)).argmin()
f = fq[idx_min:idx_max]
Iff = If[idx_min:idx_max]
plt.figure()
plt.semilogx((f), np.abs(Iff),'-*')
plt.show()

# %%
Iff_v = np.interp(fv, f, Iff)

w = np.logspace(6,10,1000)
C = 0.65e-12
R = 200
Z = 1/(1/R + 1j*w*C)
plt.figure()
# plt.loglog((fv), np.abs(vff),'-.', label= 'V')
# plt.loglog((fv),  np.abs(Iff_v),'-.', label= 'I')
plt.loglog((fv),  np.abs(vff1/Iff_v),'-.', label= 'Z')
plt.loglog((w),  np.abs(Z),'-.', label= 'Z test')
plt.grid()
plt.legend()
plt.show()


# %%
