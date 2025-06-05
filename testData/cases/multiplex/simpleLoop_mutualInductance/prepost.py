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

# %% Generate smooth initial current on nodal line

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
f = 'nodalLine_nextTo_Loop_Port2.fdtd.json'
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
plt.ylim((-0.004, 3.1))
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

# %% Calculation without the DTFT

# R_1 = 50
# L_1 = 80e-9  # Inductance of the inductor in port 1.
# i_1 = probeWire["current"].to_numpy()
# i_2 = bulkProbe["current"].to_numpy()

# mask_dI_2_neq_0 = np.abs(np.gradient(i_2, bulkProbe['time'])) != 0

# M_theo = (R_1 * i_1[mask_dI_2_neq_0] + L_1 * np.gradient(i_1, probeWire['time'])[mask_dI_2_neq_0]) / np.gradient(i_2, bulkProbe['time'])[mask_dI_2_neq_0]

# print(f"Estimated mutual inductance: {np.mean(np.abs(M_theo)) * 1e9:.2f} nH")

# %% DTFT for mutual inductance estimation

R_1 = 50
L_1 = 80e-9  # Inductance of the inductor in port 1.

new_freqs = np.geomspace(1e4, 1e9, num=500)
t = probeWire["time"].to_numpy()
I_2 = np.interp(probeWire["time"].to_numpy(), time, I_2_in)

I_1_f = np.array([np.sum(I_1 * np.exp(-1j * 2 * np.pi * f * t)) for f in new_freqs])
I_2_f = np.array([np.sum(I_2 * np.exp(-1j * 2 * np.pi * f * t)) for f in new_freqs])

mask_I_2_f_neq_0 = I_2_f != 0

M_f = np.abs((R_1 + 1j * 2 * np.pi * new_freqs[mask_I_2_f_neq_0] * L_1) * I_1_f[mask_I_2_f_neq_0] / ( -1j * 2 * np.pi * new_freqs[mask_I_2_f_neq_0] * I_2_f[mask_I_2_f_neq_0] ))

low_freqs_mask = new_freqs[mask_I_2_f_neq_0] < 1e6

M_mean = np.mean(M_f[low_freqs_mask])
print(f"Estimated mutual inductance: {M_mean * 1e9:.2f} nH")

# %% Theoretical value of current on port 1 using the estimated mutual inductance

num = [2.03e-9, 0]
den = [80e-9, 50]
system = signal.TransferFunction(num, den)
tout, I_theo_port1, _ = signal.lsim(system, U=I_2_in, T=time)

plt.figure()
plt.plot(tout, -I_theo_port1, label='Theoretical current on port 1', color='green')
plt.plot(probeWire["time"].to_numpy(), I_1, label='Current on wire. Port 1', color='blue')
plt.xlabel('Time [s]')
plt.ylabel('Ampere [A]')
plt.xlim((0, 3e-9))
plt.grid(which='both')
plt.legend()
plt.show()

# %% Magnetic current comparison

magneticCurrent = Probe(solver.getSolvedProbeFilenames("Magnetic_Bulk")[0])

mu_0 = 4 * np.pi * 1e-7
l = 20e-3

theoretical_flux = mu_0 * l * I_2_in * np.log(2) / (2 * np.pi)

plt.figure()
plt.plot(magneticCurrent["time"].to_numpy(), magneticCurrent["current"].to_numpy(), label='Magnetic current on bulk probe', color='orange')
plt.plot(time, np.gradient(theoretical_flux, time), '--', label='Theoretical magnetic current', color='red')
plt.xlabel('Time [s]')
plt.xlim((0, 3e-9))
plt.grid(which='both')
plt.legend()
plt.show()



# %% Analysis for Port 1

f_p1 = 'nodalLine_nextTo_Loop_Port1.fdtd.json'

solver_p1 = FDTD(input_filename = f_p1, path_to_exe=SEMBA_EXE, flags=["-stableradholland"])
solver_p1.cleanUp()
solver_p1.run()

# %% Visualizing the current on probe for Port 1

V_in = np.loadtxt("predefinedExcitation.1.exc", usecols=1)
time_p1 = np.loadtxt("predefinedExcitation.1.exc", usecols=0)

plt.figure()
plt.plot(time_p1, V_in, label='Initial voltage on loop. Port 1')
plt.grid()
plt.xlabel('Time [s]')
plt.ylabel('Volt [V]')
plt.legend()
plt.show()

# %% Visualizing the current on probe for Port 1 and port 2	

probeWire_p1 = Probe(solver_p1.getSolvedProbeFilenames("Wire probe_port1")[0])
I_1_p1 = probeWire_p1["current"].to_numpy()

num = [1]
den = [80e-9, 50]
system = signal.TransferFunction(num, den)
tout, I_theo, _ = signal.lsim(system, U=V_in, T=time_p1)

plt.figure()
plt.plot(probeWire_p1["time"].to_numpy(), I_1_p1, '--', label='Current on wire. Port 1', color='blue')
plt.plot(tout, I_theo, label='Theoretical current on wire. Port 1', color='green')
plt.grid(which='both')
plt.xlabel('Time [s]')
plt.ylabel('Ampere [A]')
plt.legend()
plt.show()

probeBulk_p1 = Probe(solver_p1.getSolvedProbeFilenames("Bulk probe_port2")[0])
I_2_p1 = probeBulk_p1["current"].to_numpy()
plt.figure()
plt.plot(probeBulk_p1["time"].to_numpy(), I_2_p1, label='Current on bulk probe. Port 2', color='orange')
plt.grid(which='both')
plt.xlabel('Time [s]')
plt.ylabel('Ampere [A]')
plt.legend()
plt.show()


# %% Proper inductance estimation for Port 1

R_1 = 50

new_freqs = np.geomspace(1e4, 1e9, num=500)

V = np.interp(probeWire_p1["time"].to_numpy(), time_p1, V_in)
I = probeWire_p1["current"].to_numpy()
t = probeWire_p1["time"].to_numpy()

V_f = np.array([np.sum(V * np.exp(-1j * 2 * np.pi * f * t)) for f in new_freqs])
I_f = np.array([np.sum(I * np.exp(-1j * 2 * np.pi * f * t)) for f in new_freqs])

mask_I_f_neq_0 = I_f != 0

Z_f = V_f[mask_I_f_neq_0] / I_f[mask_I_f_neq_0]

L = np.abs(np.imag(Z_f))/(2 * np.pi * new_freqs[mask_I_f_neq_0])
low_freqs_mask = new_freqs[mask_I_f_neq_0] < 1e6
L_mean = np.mean(L[low_freqs_mask])



# %%
