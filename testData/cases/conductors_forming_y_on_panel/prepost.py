# %%
import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %% Run solver
f = 'Conductors_50ohm_terminals.fdtd.json'

solver = FDTD(input_filename=f, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

# %% Read excitation waveform
exc_file = solver["sources"][0]["magnitudeFile"]
exc_data = np.loadtxt(exc_file)
exc_time = exc_data[:, 0]
exc_signal = exc_data[:, 1]

plt.figure()
plt.plot(exc_time * 1e9, exc_signal)
plt.xlabel('Time (ns)')
plt.ylabel('Amplitude')
plt.title('Excitation waveform')
plt.grid(which='both')
plt.show()

# %% Read wire current probes
probe_yplus  = Probe(solver.getSolvedProbeFilenames("curr_yplus_unc")[0])
probe_yminus = Probe(solver.getSolvedProbeFilenames("curr_yminus_unc")[0])
probe_joined = Probe(solver.getSolvedProbeFilenames("curr_joined")[0])

plt.figure()
plt.plot(probe_yplus["time"].to_numpy(),  probe_yplus["current"].to_numpy(),  label='curr_yplus_unc')
plt.plot(probe_yminus["time"].to_numpy(), probe_yminus["current"].to_numpy(), label='curr_yminus_unc')
plt.plot(probe_joined["time"].to_numpy(), probe_joined["current"].to_numpy(), label='curr_joined')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')
plt.title('Wire currents')
plt.grid(which='both')
plt.legend()
plt.show()

# %% Read bulk current probe
probe_bulk = Probe(solver.getSolvedProbeFilenames("BC")[0])

plt.figure()
plt.plot(probe_bulk["time"].to_numpy(), probe_bulk["current"].to_numpy(), label='BC')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')
plt.title('Bulk current')
plt.grid(which='both')
plt.legend()
plt.show()

# %%
