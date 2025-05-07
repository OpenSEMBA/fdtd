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
# GAUSSIAN PULSE

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
np.savetxt('./step.exc', data)

# %%
# SIGMOID PULSE

t = np.linspace(0, 10e-8, 2000)

A = 10     
k = 30e7
t0 = 0.1e-8

sigmoid_raw = 1 / (1 + np.exp(-k * (t - t0)))
offset = 1 / (1 + np.exp(k * t0))
scaling = 1 / (1 - offset)
V = A * (sigmoid_raw - offset) * scaling

plt.figure(figsize=(8, 4))
plt.plot(t, V)
plt.xlabel('t (s)')
plt.ylabel('Voltage (V)')
plt.grid(True)
plt.tight_layout()
plt.show()

# Save excitation data
data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = V
np.savetxt('./step.exc', data)

# %%
fn = 'closedCircuit.fdtd.json'
solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, flags=['-mapvtk'])
solver.cleanUp()
solver.run()

# %%
data = np.loadtxt('closedCircuit.fdtd_cube_current_Jy_99_150_99__100_150_100.dat', skiprows=1)

x1 = data[:, 0]
y1 = data[:, 1]

data = np.loadtxt('closedCircuit.fdtd_wire_current_Jx_105_99_99__105_100_100.dat', skiprows=1)

x = data[:, 0]
y = data[:, 1]

data = np.loadtxt('closedCircuit.fdtd_source_current_Jy_109_150_99__110_150_100.dat', skiprows=1)

x2 = data[:, 0]
y2 = data[:, 1]

data = np.loadtxt('closedCircuit.fdtd_cube_current_Jy_98_150_98__101_150_101.dat', skiprows=1)

x3 = data[:, 0]
y3 = data[:, 1]

plt.figure(figsize=(8, 9))

plt.plot(x, y, color = 'red',label = 'left wire medition' )
plt.plot(x2, y2, color = 'green',label = 'source medition' )

plt.plot(x1, -y1, color = 'blue', label = 'resistance medition')
#plt.plot(x3, -y3, color = 'blue', label = 'perpendicular resistance medition')

plt.axvline(x=0.3e-8, color='orange', linestyle='-', label='1/c')
plt.axvline(x=3.6e-8, color='orange', linestyle='-', label='11/c')

plt.xlabel('time')
plt.ylabel('Current')
plt.grid(True)
plt.legend()

plt.xlim(right=1e-7)

plt.show()
# %%
