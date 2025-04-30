# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../../build/bin/semba-fdtd'

from pyWrapper import *

# # %% Generate excitation and visualize
# def generateGaussExcitation():
#     dt = 1e-12
#     w0 = 0.5e-8 
#     t0 = 10*w0
#     t = np.arange(0, t0+20*w0, dt)
#     f = np.exp( -np.power(t-t0,2)/ w0**2 )

#     data = np.zeros((len(t), 2))
#     data[:,0] = t
#     data[:,1] = f
#     np.savetxt('gauss.exc', data)

# %% Run solver
fn = 'current_bifurcation_lumped.fdtd.json'
# generateGaussExcitation()
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

# %% Visualizing initial values of currents and voltages

V_in = np.loadtxt("predefinedExcitation.1.exc", usecols=1)
time = np.loadtxt("predefinedExcitation.1.exc", usecols=0)
plt.figure()
plt.plot(time, V_in, label='Initial excitation voltage')    
plt.grid(which='both')
plt.legend()
plt.show()

InitialBulk_probe = Probe(solver.getSolvedProbeFilenames("Bulk Initial probe")[0])
TopBulk_probe = Probe(solver.getSolvedProbeFilenames("Bulk Top probe")[0])
BottomBulk_probe = Probe(solver.getSolvedProbeFilenames("Bulk Bottom probe")[0])

fig, axes = plt.subplots(1, 2, figsize=(15, 4), sharey=True)
axes[0].plot(InitialBulk_probe['time'], InitialBulk_probe['current'], label='Initial current', color='blue')
axes[0].set_title('Initial current')
axes[0].set_xlabel('Time')
axes[0].set_ylabel('Current')
axes[0].legend()
axes[0].grid(which='both')

axes[1].plot(TopBulk_probe['time'], TopBulk_probe['current'], label='Current on the top part of the circuit', color='green')
axes[1].plot(BottomBulk_probe['time'], BottomBulk_probe['current'], '--', label='Current on the bottom part of the circuit', color='red')
axes[1].set_title('Current on the circuit')
axes[1].set_xlabel('Time')
axes[1].legend()
axes[1].grid(which='both')

plt.tight_layout()
plt.show()

# %% Postprocess

Top_voltage = -2*np.loadtxt("current_bifurcation_alternative.fdtd_Top current_Wx_34_20_50_s3.dat", skiprows=1, usecols=2)
Bottom_voltage = -2*np.loadtxt("current_bifurcation_alternative.fdtd_Bottom current_Wx_34_20_30_s2.dat", skiprows=1, usecols=2)

plt.figure()
plt.plot(TopBulk_probe['time'], Top_voltage, label='Voltage on the top part of the circuit', color='green')
plt.plot(BottomBulk_probe['time'], Bottom_voltage, '--', label='Voltage on the bottom part of the circuit', color='red')
plt.title('Voltage on the circuit')
plt.xlabel('Time')
plt.ylabel('Voltage')
plt.legend()
plt.grid(which='both')
plt.show()


t = TopBulk_probe['time']
dt = t[1] - t[0]

fq = fftfreq(len(t))/dt

fmin = 1e4
fmax = 2e8
idx_min = (np.abs(fq - fmin)).argmin()
idx_max = (np.abs(fq - fmax)).argmin()

V_top_f = fft(Top_voltage)[idx_min:idx_max]
I_top_f = fft(TopBulk_probe['current'])[idx_min:idx_max]
Z_top_f = V_top_f / I_top_f

Z_bottom_f = fft(Bottom_voltage)[idx_min:idx_max]
I_bottom_f = fft(BottomBulk_probe['current'])[idx_min:idx_max]
Z_bottom_f = Z_bottom_f / I_bottom_f


plt.figure()
plt.plot(fq[idx_min:idx_max], np.abs(Z_top_f), '.-', label='Magnitude of impedance on the top part of the circuit')
plt.plot(fq[idx_min:idx_max], np.abs(Z_bottom_f), '.--', label='Magnitude of impedance on the bottom part of the circuit')
plt.xlabel('Frequency')
plt.ylabel('Impedance')
plt.grid(which='both')
plt.xlim(fmin, fmax)
plt.xscale('log')
plt.legend()
plt.show()

# %%
