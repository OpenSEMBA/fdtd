# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import os

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %%
dt = 6e-9  # Time step
tf = 6e-7  # Final time
t = np.arange(0, tf, dt)

# Create a Gaussian pulse
t0 = 3e-7  # Center of the pulse
sigma = 3e-7  # Width of the pulse
v = 10.0 * np.exp(-(t - t0)**2 / (2 * sigma**2))

plt.figure(figsize=(10, 6))
plt.plot(t, v)
plt.xlabel('Time (s)')
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
solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, flags=['-mapvtk'])
solver.cleanUp()
solver.run()

# resistance probe data
wire_current = Probe(solver.getSolvedProbeFilenames("wire_current")[0])

plt.plot(wire_current['time'], wire_current['current'], label='wire')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')
plt.title('Current Signals')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
# %%

"""
# %%
# resistance probe data
resistance_current = Probe(solver.getSolvedProbeFilenames("resistance_current")[0])

plt.plot(resistance_current['time'], resistance_current['current'], label='resistance')
plt.xlabel('Time (s)')
plt.ylabel('Current (A)')
plt.title('Current Signals')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
# %%
"""