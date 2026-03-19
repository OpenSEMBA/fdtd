# %%
import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# # %% Generate gaussian excitation
# # Gaussian peaks at t0=10 ns with spread (sigma) of 4 ns
# dt = 1e-11
# tf = 35e-9
# t = np.arange(0, tf, dt)

# t0    = 10e-9
# sigma = 4e-9
# gauss = np.exp(-((t - t0) / sigma) ** 2)

# exc_filename = 'gauss.exc'
# np.savetxt(exc_filename, np.column_stack((t, gauss)))

# plt.figure()
# plt.plot(t * 1e9, gauss)
# plt.xlabel('Time (ns)')
# plt.ylabel('Amplitude')
# plt.title('Gaussian excitation')
# plt.grid(which='both')
# plt.show()

# %% Run solver
# f = 'bulk_currents_X_oriented.fdtd.json'
f = 'bulk_currents_Y_oriented.fdtd.json'
# f = 'bulk_currents_Z_oriented.fdtd.json'

solver = FDTD(input_filename=f, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

# %% Show outputs
probe_LL = Probe(solver.getSolvedProbeFilenames("BC_LL")[0])
probe_LU = Probe(solver.getSolvedProbeFilenames("BC_LU")[0])
probe_UU = Probe(solver.getSolvedProbeFilenames("BC_UU")[0])
probe_UL = Probe(solver.getSolvedProbeFilenames("BC_UL")[0])

probe_time = probe_LL["time"].to_numpy()
exc_interp = np.interp(probe_time, t, gauss)

plt.figure()
plt.plot(probe_time, exc_interp, 'k-', label='Excitation (interp)')
plt.plot(probe_LL["time"].to_numpy(), probe_LL["current"].to_numpy(), ':', label='BC_LL')
plt.plot(probe_LU["time"].to_numpy(), probe_LU["current"].to_numpy(), ':', label='BC_LU')
plt.plot(probe_UU["time"].to_numpy(), probe_UU["current"].to_numpy(), ':', label='BC_UU')
plt.plot(probe_UL["time"].to_numpy(), probe_UL["current"].to_numpy(), ':', label='BC_UL')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')
plt.title('Bulk current – four probes')
plt.grid(which='both')
plt.legend()
plt.show()

# %% Checks
assert np.max(np.abs(probe_LL["current"].to_numpy())) < 2e-3
assert np.max(np.abs(probe_LU["current"].to_numpy())) < 2e-3
assert np.corrcoef(exc_interp, probe_UU["current"].to_numpy())[0, 1] > 0.9999
assert np.max(np.abs(probe_UL["current"].to_numpy())) < 2e-3

# %%