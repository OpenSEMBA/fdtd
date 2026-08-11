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


# %% Read bulk current probe
probe_bulk = Probe(solver.getSolvedProbeFilenames("BC")[0])
assert not probe_bulk["current"].isnull().any(), "probe_bulk current contains NaN values"

plt.figure()
plt.plot(exc_time, exc_signal, label='Source')
plt.plot(probe_bulk["time"].to_numpy(), probe_bulk["current"].to_numpy(), label='BC')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')
plt.title('Bulk current')
plt.grid(which='both')
plt.legend()
plt.show()



# %%
