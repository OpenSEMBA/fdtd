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
# Load probe data
source_voltage = Probe(solver.getSolvedProbeFilenames("source_voltage")[0])
source_current = Probe(solver.getSolvedProbeFilenames("source_current")[0])
load_voltage = Probe(solver.getSolvedProbeFilenames("load_voltage")[0])
load_current = Probe(solver.getSolvedProbeFilenames("load_current")[0])

# Plot voltage signals
plt.figure(figsize=(12, 8))
plt.subplot(2, 1, 1)
plt.plot(source_voltage['time']*1e9, source_voltage['voltage'], label='Source')
plt.plot(load_voltage['time']*1e9, load_voltage['voltage'], label='Load')
plt.xlabel('Time (ns)')
plt.ylabel('Voltage (V)')
plt.title('Voltage Signals')
plt.legend()
plt.grid(True)

# Plot current signals
plt.subplot(2, 1, 2)
plt.plot(source_current['time']*1e9, source_current['current'], label='Source')
plt.plot(load_current['time']*1e9, load_current['current'], label='Load')
plt.xlabel('Time (ns)')
plt.ylabel('Current (A)')
plt.title('Current Signals')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

# %%
# Load electric field data
field_xz = Probe(solver.getSolvedProbeFilenames("electric_field_xz")[0])
field_xy = Probe(solver.getSolvedProbeFilenames("electric_field_xy")[0])

# Plot electric field magnitude over time
plt.figure(figsize=(12, 8))
plt.subplot(2, 1, 1)
plt.plot(field_xz['time']*1e9, np.abs(field_xz['field']))
plt.xlabel('Time (ns)')
plt.ylabel('|E| (V/m)')
plt.title('Electric Field Magnitude (XZ plane)')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(field_xy['time']*1e9, np.abs(field_xy['field']))
plt.xlabel('Time (ns)')
plt.ylabel('|E| (V/m)')
plt.title('Electric Field Magnitude (XY plane)')
plt.grid(True)

plt.tight_layout()
plt.show()