# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import json
import os
from pathlib import Path
import scipy.constants
from skrf.media import Freespace
from skrf.frequency import Frequency

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %%
dt = 1e-12  # Time step
tf = 150e-9  # Final time
t = np.arange(0, tf, dt)

# Create a Gaussian pulse
t0 = 20e-9  # Center of the pulse
sigma = 5e-9  # Width of the pulse
v = 10.0 * np.exp(-(t - t0)**2 / (2 * sigma**2))

plt.figure(figsize=(10, 6))
plt.plot(t*1e9, v)
plt.xlabel('Time (ns)')
plt.ylabel('Voltage (V)')
plt.title('Source Voltage')
plt.grid(True)
plt.show()

# Save excitation data
data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = v
np.savetxt('step.exc', data)

# %%
fn = 'closedCircuit.fdtd.json'
solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

# %%
# resistance probe data
resistance_current = Probe(solver.getSolvedProbeFilenames("resistance_current")[0])

plt.plot(resistance_current['time']*1e9, resistance_current['current'], label='resistance')
plt.xlabel('Time (ns)')
plt.ylabel('Current (A)')
plt.title('Current Signals')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()