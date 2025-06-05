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

# %% Run solver
f = 'FullMontage.fdtd.json'

solver = FDTD(input_filename = f, path_to_exe=SEMBA_EXE, flags=["-stableradholland"])
solver.cleanUp()
solver.run()
# %% Visualizing initial value of voltage

V_in = np.loadtxt("predefinedExcitation.1.exc", usecols=1)
time = np.loadtxt("predefinedExcitation.1.exc", usecols=0)

plt.figure()
plt.plot(time, V_in, label='Initial Voltage on source')
plt.xlabel('Time [s]')
plt.ylabel('Voltage [V]')
plt.grid()
plt.legend()
plt.show()

# %% Visualizing the current on probe

probeWire = Probe(solver.getSolvedProbeFilenames("Wire probe")[0])

plt.figure()
plt.plot(probeWire["time"].to_numpy(), probeWire["current"].to_numpy(), label='Current on wire', color='blue')
# plt.xlim((0, 2e-8))
# plt.ylim((-0.004, 0.0001))
plt.grid(which='both')
plt.legend()
plt.show()


# %% DTFT for impedance and cutoff frequency estimation

new_freqs = np.geomspace(1e4, 1e9, num=500)

V = np.interp(probeWire["time"].to_numpy(), time, V_in)
I = probeWire["current"].to_numpy()
t = probeWire["time"].to_numpy()

V_f = np.array([np.sum(V * np.exp(-1j * 2 * np.pi * f * t)) for f in new_freqs])
I_f = np.array([np.sum(I * np.exp(-1j * 2 * np.pi * f * t)) for f in new_freqs])

mask_I_f_neq_0 = I_f != 0

Z_f = V_f[mask_I_f_neq_0] / I_f[mask_I_f_neq_0]

H_Z = 20 * np.log10(np.abs(Z_f))
dec = np.log10(new_freqs)
m = (H_Z[1:] - H_Z[:-1]) / (dec[1:] - dec[:-1]) 

mask_20dB_decade = (18 <= m) & (m <= 22)

mag_ref = 20 * np.log10(np.abs(Z_f[np.where(m == m[mask_20dB_decade][0])]))
f_ref = new_freqs[np.where(m == m[mask_20dB_decade][0])]
mag_line_in_dB = mag_ref + 20 * np.log10(new_freqs / f_ref)
mag_line = 10 ** (mag_line_in_dB / 20)

plt.figure()
plt.plot(new_freqs[mask_I_f_neq_0], np.abs(Z_f), '.', label='Impedance in frequency domain using DTFT', color='red')
plt.semilogx(new_freqs, mag_line, '--', label='Reference line with 20dB/dec slope', color='black')
plt.axhline(y=50, color='green', linestyle='--', label='Reference line at 50 Ohm')
plt.xscale('log')
plt.yscale('log')
plt.ylim((30, 10e5))
plt.xlabel('Frequency [Hz]')
plt.ylabel('Impedance [Ohm]')
plt.grid(which='both')
plt.legend()
plt.show()

plt.figure()
plt.plot(new_freqs[mask_I_f_neq_0], np.abs(np.imag(Z_f)), '.', label='Complex impedance in frequency domain', color='red')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Frequency [Hz]')
plt.ylabel('Impedance [Ohm]')
plt.grid(which='both')
plt.legend()
plt.show()

# Inductance estimation
L = np.abs(np.imag(Z_f))/(2 * np.pi * new_freqs[mask_I_f_neq_0])

low_freqs_mask = new_freqs[mask_I_f_neq_0] < 1e6

L_mean = np.mean(L[low_freqs_mask])
print(f"Estimated inductance: {L_mean * 1e9:.2f} nH")

# %%
