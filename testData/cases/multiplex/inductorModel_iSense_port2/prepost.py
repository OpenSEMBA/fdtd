# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy import signal


import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../../build/bin/semba-fdtd'

from pyWrapper import *

def smooth_ramp():
    dt = 0.8e-12
    t_final = 20e-9
    t = np.arange(0, t_final, dt)

    f = np.zeros_like(t)

    t1 = 0.5e-9
    t2 = 0.75e-9
    k = 1 / t2  

    epsilon = 1e-15 # control of the smoothness of the ramp
    x = (t - t1)
    ramp_region = (t >= t1) & (t < t2)

    f[ramp_region] = k * x[ramp_region] * ( 1 + x[ramp_region] / np.sqrt(x[ramp_region]**2 + epsilon**2) ) / 2
    f[t >= t2] = k*x[ t >= t2]

    data = np.column_stack((t, f))
    np.savetxt('smoothRampExcitation.exc', data)


# %% Run solver
f = 'FullMontage.fdtd.json'
smooth_ramp()

solver = FDTD(input_filename = f, path_to_exe=SEMBA_EXE, flags=["-stableradholland"])
solver.cleanUp()
solver.run()
# %% Visualizing initial value of voltage adn the current on bulk

I_2_in = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=1)
time = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=0)

bulkProbe = Probe(solver.getSolvedProbeFilenames("Bulk probe_port2")[0])

plt.figure()
plt.plot(time, I_2_in, label='Initial current on nodal line. Port 2')
plt.plot(bulkProbe["time"].to_numpy(), bulkProbe["current"].to_numpy(), '--', label='Current on bulk probe. Port 2', color='orange')
plt.xlabel('Time [s]')
plt.ylabel('Ampere [A]')
plt.grid()
plt.legend()
plt.show()

# %% Visualizing the current on probe

probeWire = Probe(solver.getSolvedProbeFilenames("Wire probe_port1")[0])
I_1 = probeWire["current"].to_numpy()

plt.figure()
plt.plot(probeWire["time"].to_numpy(),I_1, label='Current on wire. Port 1', color='blue')
# plt.xlim((0, 2e-8))
# plt.ylim((-0.004, 0.0001))
plt.grid(which='both')
plt.legend()

plt.show()
# %% Calculation of mutual inductance without the DTFT

R_1 = 50
L_1 = 243e-9  # Inductance of the inductor in port 1.
i_1 = probeWire["current"].to_numpy()
i_2 = bulkProbe["current"].to_numpy()

mask_dI_2_neq_0 = np.abs(np.gradient(i_2, bulkProbe['time'])) != 0

M_theo = (R_1 * i_1[mask_dI_2_neq_0] + L_1 * np.gradient(i_1, probeWire['time'])[mask_dI_2_neq_0]) / np.gradient(i_2, bulkProbe['time'])[mask_dI_2_neq_0]

print(f"Estimated mutual inductance: {np.mean(M_theo) * 1e9:.2f} nH")

# %% Theoretical value of current on port 1 using the estimated mutual inductance

# num = [7.3e-9, 0]
# den = [L_1, R_1]
# system = signal.TransferFunction(num, den)
# tout, I_theo_port1, _ = signal.lsim(system, U=I_2_in, T=time)

# plt.figure()
# plt.plot(tout, I_theo_port1, label='Theoretical current on port 1', color='green')
# plt.plot(probeWire["time"].to_numpy(), np.abs(I_1), label='Current on wire. Port 1', color='blue')
# plt.xlabel('Time [s]')
# plt.ylabel('Ampere [A]')
# plt.xlim((0, 3e-9))
# plt.grid(which='both')
# plt.legend()
# plt.show()


# %% DTFT for mutual inductance estimation

R_1 = 0
L_1 = 243e-9  # Inductance of the inductor in port 1. Check the prepost on inductorModel_iSense_port1 for the value.

new_freqs = np.geomspace(1e4, 1e9, num=500)
t = probeWire["time"].to_numpy()
I_2 = np.interp(probeWire["time"].to_numpy(), time, I_2_in)

I_1_f = np.array([np.sum(I_1 * np.exp(-1j * 2 * np.pi * f * t)) for f in new_freqs])
I_2_f = np.array([np.sum(I_2 * np.exp(-1j * 2 * np.pi * f * t)) for f in new_freqs])

mask_I_2_f_neq_0 = I_2_f != 0

M_f = (R_1 + 1j * 2 * np.pi * new_freqs[mask_I_2_f_neq_0] * L_1) * I_1_f[mask_I_2_f_neq_0] / ( -1j * 2 * np.pi * new_freqs[mask_I_2_f_neq_0] * I_2_f[mask_I_2_f_neq_0] )

plt.figure()
plt.plot(new_freqs[mask_I_2_f_neq_0], np.abs(M_f), label='Mutual Inductance')
plt.xscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Mutual Inductance [H]')
plt.grid(which='both')
plt.legend()
plt.show()

M_mean = np.mean(np.abs(M_f[new_freqs < 1e6]))
print(f"Estimated mutual inductance: {M_mean * 1e9:.2f} nH")


# %%
