# %%
import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import *

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %% Generate excitation and visualize
fmin = 2e6; fmax = 6e9

#####################################################
# %% Generate excitation and visualize
dt = 0.25e-10
tf = 60e-9

t = np.arange(0, tf, dt)
v = 1/(1 + np.exp(0.3*(-t*1e9+25)))
# v2 = 0*(t<=100e-9) + ((t*1e9-100)/300) * (t>=100e-9)*(t<=400e-9) + 1*(t>=400e-9)
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
# vf2 = fft(v2)
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
#####################################################
# %%
fn = 'lineIntegralProbe.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

#####################################################
# %% Plot results
lineInt  = Probe(solver.getSolvedProbeFilenames("vprobe_LI_7_10_10")[0])
plt.figure()
plt.plot(lineInt['time']*1e9, lineInt['lineIntegral'], '-', label = 'line integral')
plt.grid(which='both')
plt.legend()

# closed  =   Probe(solver.getSolvedProbeFilenames("vprobe_LI_7_10_10")[0])
# left  =     Probe(solver.getSolvedProbeFilenames("vprobe_left_side_LI_7_10_10")[0])
# right  =    Probe(solver.getSolvedProbeFilenames("vprobe_right_side_LI_13_10_12")[0])
# above  =    Probe(solver.getSolvedProbeFilenames("vprobe_Rabove_LI_7_10_13")[0])
# ground  =   Probe(solver.getSolvedProbeFilenames("vprobe_ground_LI_7_10_10")[0])
# i      =    Probe(solver.getSolvedProbeFilenames("current_Wz_13_10_12")[0])
# plt.figure()
# plt.plot(above['time']*1e9, above['lineIntegral'], '-', label = 'line integral above R')
# plt.plot(left['time']*1e9, -left['lineIntegral'], '--', label = 'line integral left arm')
# plt.plot(right['time']*1e9, right['lineIntegral'], '-.', label = 'line integral right arm ')
# plt.plot(ground['time']*1e9, ground['lineIntegral'], '-.', label = 'line integral ground ')
# plt.grid(which='both')
# plt.legend()

# plt.figure()
# plt.plot(i['time']*1e9, i['current'], '-', label = 'I(7,10,13)x')
# plt.axhline(y=1/100, color='r', linestyle='-')
# plt.xlabel('t [ns]')
# plt.ylabel('I [A]')
# plt.grid(which='both')
# plt.legend()

#####################################################

# %% 
dt = iprobe_7_10_11['time'][1] - iprobe_7_10_11['time'][0]
fq  = fftfreq(len(iprobe_7_10_11['time']))/dt
If = fft(iprobe_7_10_11['current'])
idx_min = (np.abs(fq - fmin)).argmin()
idx_max = (np.abs(fq - fmax)).argmin()
f = fq[idx_min:idx_max]
Iff = If[idx_min:idx_max]
Iff_v = np.interp(fv, f, Iff)
plt.figure()
plt.semilogx((f), np.abs(Iff),'-*')
plt.semilogx((f), np.abs(Iff_v),'--')
plt.show()

# %%

w = np.logspace(6,10,1000)
C = 0.65e-12
R = 150
Z = 1/(1/R + 1j*w*C)
plt.figure()
# plt.loglog((fv), np.abs(vff),'-.', label= 'V')
# plt.loglog((fv),  np.abs(Iff_v),'-.', label= 'I')
plt.semilogx((fv),  np.abs(vff1/Iff_v),'-.', label= 'Z')
# plt.loglog((w),  np.abs(Z),'-.', label= 'Z test')
plt.ylim(100,200)
plt.grid()
plt.legend()
plt.show()


# %%
