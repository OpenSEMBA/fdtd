# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *


# %% Generate excitation and visualize
dt = 5e-11
w0 = 0.1e-7    # ~ 2 GHz bandwidth
# t0 = 10*w0
t0 = 3e-8
t = np.arange(0, t0+20*w0, dt)
f = np.exp( -np.power(t-t0,2)/ w0**2 )

plt.figure()
plt.plot(t,f)

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = f
np.savetxt('bidirectionalShield.exc', data)

fq = fftfreq(len(t))/dt
F = fft(f)
F = F[fq>=0]
fq = fq[fq>=0]
Fmax = np.max(np.abs(F))
Fcut = Fmax / np.sqrt(2.0)
idx = (np.abs(np.abs(F) - Fcut)).argmin()

plt.figure()
plt.semilogx(fq, np.abs(F),'.-')
plt.axvline(x = fq[idx], color = 'r')
plt.xlim(1e6,4e9)
plt.ylim(1e-9, 1e6)
plt.yscale('log')

plt.grid()

# %% Prepare solver
solver = FDTD(input_filename = 'bidirectionalShield.fdtd.json', path_to_exe=SEMBA_EXE, flags=['-mtlnwires'])
solver.cleanUp()
solver.run()
# %% Load probes

E_mid_x = Probe(solver.getSolvedProbeFilenames("electric_field_point_probe_mid_Ex_70_75_76")[0])
E_mid_y = Probe(solver.getSolvedProbeFilenames("electric_field_point_probe_mid_Ey_70_75_76")[0])
E_mid_z = Probe(solver.getSolvedProbeFilenames("electric_field_point_probe_mid_Ez_70_75_76")[0])

I1 = Probe(solver.getSolvedProbeFilenames("wire_1_end_bundle_cable_1_shield_I_70_70_74")[0])
I2 = Probe(solver.getSolvedProbeFilenames("wire_2_end_bundle_cable_2_shield_I_80_70_74")[0])

plt.plot(I1['time']*1e6, I1['current_1'], label='agressor cable, inner wire') 
plt.plot(I2['time']*1e6, 1e5*I2['current_0'], label=r'victim cable, shield (x$10^{5}$)') 
plt.plot(I2['time']*1e6, 1e7*I2['current_1'], '--',label=r'victim cable, inner wire (x$10^{7}$)') 
plt.legend()
plt.xlabel('t [us]')
plt.ylabel('I [A]')
plt.xlim(0,0.1)
plt.show()

plt.plot(E_mid_x['time']*1e6, E_mid_x['field'], label=r'E$_x$') 
plt.plot(E_mid_y['time']*1e6, E_mid_y['field'], label=r'E$_y$') 
plt.plot(E_mid_z['time']*1e6, E_mid_z['field'], label=r'E$_z$') 
plt.legend()
plt.title('Field at midpoint between wires')
plt.xlabel('t [us]')
plt.ylabel('|E| [V/m]')
plt.xlim(0,0.1)

# %%
