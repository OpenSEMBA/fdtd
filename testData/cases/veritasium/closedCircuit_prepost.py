# %% [markdown]
# # Veritasium Circuit Analysis
# 
# This script analyzes the simulation results of the circuit demonstrated in Veritasium's video about electricity speed.

# %% [markdown]
# ## 1. Import Libraries and Setup Environment

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

# Matplotlib configuration for better visualization
plt.style.use('default')  # Using default style instead of seaborn
plt.rcParams['figure.figsize'] = [12, 8]
plt.rcParams['axes.grid'] = True  # Enable grid by default
plt.rcParams['grid.alpha'] = 0.3  # Make grid slightly transparent

# %% [markdown]
# ## 2. Load and Visualize Geometry
# 
# First, we'll load the VTK file generated with the -mapvtk argument to visualize the circuit geometry.

# %%
def load_vtk_geometry(vtk_file):
    """Load and process the VTK geometry file"""
    # TODO: Implement VTK file loading
    pass

# Load geometry
geometry_file = "closedCircuit.vtk"
geometry = load_vtk_geometry(geometry_file)

def plot_geometry(geometry):
    """Visualize the circuit geometry"""
    # TODO: Implement visualization
    pass

# %% [markdown]
# ## 1. Generate and Visualize Excitation
# First, let's create a step function excitation for our voltage source

# %%
dt = 1e-12  # Time step
tf = 150e-9  # Final time
t = np.arange(0, tf, dt)
v = 10.0 * (t > 0)  # 10V step function

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

# %% [markdown]
# ## 2. Run Solver
# Now let's run the FDTD simulation

# %%
fn = 'closedCircuit.fdtd.json'
solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

# %% [markdown]
# ## 3. Load and Process Probe Data
# Let's analyze the voltage and current signals

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

# %% [markdown]
# ## 4. Calculate Propagation Speed
# Let's estimate the speed of propagation by finding the delay between source and load signals

# %%
def find_signal_delay(source_signal, load_signal, threshold=0.5):
    """Find the delay between source and load signals using a threshold crossing"""
    source_time = source_signal['time']
    source_data = source_signal['voltage']
    load_time = load_signal['time']
    load_data = load_signal['voltage']
    
    # Find first crossing of threshold
    source_crossing = source_time[source_data > threshold][0]
    load_crossing = load_time[load_data > threshold][0]
    
    return load_crossing - source_crossing

# Calculate delay and speed
delay = find_signal_delay(source_voltage, load_voltage)
distance = 300  # meters (cable length)
speed = distance / delay

print(f"Signal delay: {delay*1e9:.2f} ns")
print(f"Estimated propagation speed: {speed/1e8:.2f} × 10⁸ m/s")

# %% [markdown]
# ## 5. Analyze Electric Field
# Let's examine the electric field propagation

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

# %% [markdown]
# ## 6. Load and Visualize Field Movie
# Let's examine the field propagation animation

# %%
# Load movie data
movie_data = Probe(solver.getSolvedProbeFilenames("field_movie")[0])

def plot_field_frame(frame_data, frame_number):
    """Plot a specific frame of the electric field movie"""
    plt.figure(figsize=(10, 8))
    plt.imshow(np.abs(frame_data), cmap='viridis')
    plt.colorbar(label='|E| (V/m)')
    plt.title(f'Electric Field Magnitude - Frame {frame_number}')
    plt.show()

# Plot a few key frames
key_frames = [0, len(movie_data['field'])//4, len(movie_data['field'])//2, 
             3*len(movie_data['field'])//4, len(movie_data['field'])-1]

for frame in key_frames:
    plot_field_frame(movie_data['field'][frame], frame) 