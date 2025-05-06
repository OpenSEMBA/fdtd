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

# %%
data = np.loadtxt('closedCircuit.fdtd_wire_current_Jx_105_99_99__105_100_100.dat', skiprows=1)
time = data[:, 0]
current = data[:, 1]

# Graficar los datos
plt.plot(time, current)
plt.xlabel('Tiempo (s)')
plt.ylabel('Corriente (A)')
plt.title('Respuesta de la sonda del cable')
plt.grid(True)
plt.show()
# %%
data = np.loadtxt('closedCircuit.fdtd_source_current_Jy_109_150_99__110_150_100.dat', skiprows=1)
time = data[:, 0]
current = data[:, 1]

# Graficar los datos
plt.plot(time, current)
plt.xlabel('Tiempo (s)')
plt.ylabel('Corriente (A)')
plt.title('Respuesta de la sonda del source')
plt.grid(True)
plt.show()

# %%
data = np.loadtxt('closedCircuit.fdtd_cube_current_Jy_99_150_99__100_150_100.dat', skiprows=1)
time = data[:, 0]
current = data[:, 1]

# Graficar los datos
plt.plot(time, current)
plt.xlabel('Tiempo (s)')
plt.ylabel('Corriente (A)')
plt.title('Respuesta de la sonda del source')
plt.grid(True)
plt.show()