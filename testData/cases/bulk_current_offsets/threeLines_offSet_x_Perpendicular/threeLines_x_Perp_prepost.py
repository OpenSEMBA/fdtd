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
f = 'threeLines.fdtd.json'
f_negative = 'threeLinesNegative.fdtd.json'
f_positive = 'threeLinesPositive.fdtd.json'

solver = FDTD(input_filename = f, path_to_exe=SEMBA_EXE)
solver_negative = FDTD(input_filename = f_negative, path_to_exe=SEMBA_EXE)
solver_positive = FDTD(input_filename = f_positive, path_to_exe=SEMBA_EXE)
solver.cleanUp()

solver.run()
solver_negative.run()
solver_positive.run()

# %% Visualizing initial values 

I_in_1 = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=1)
time_1 = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=0)

I_in_2 = np.loadtxt(solver["sources"][1]["magnitudeFile"], usecols=1)
time_2 = np.loadtxt(solver["sources"][1]["magnitudeFile"], usecols=0)

I_in_3 = np.loadtxt(solver["sources"][2]["magnitudeFile"], usecols=1)
time_3 = np.loadtxt(solver["sources"][2]["magnitudeFile"], usecols=0)

plt.figure()
plt.plot(time_1, I_in_1, label='Initial excitation, source 1')
plt.plot(time_2, I_in_2, '--', color='red', label='Initial excitation, source 2')   
plt.plot(time_3, I_in_3, '--', color='green', label='Initial excitation, source 3')     
plt.grid(which='both')
plt.legend()
plt.show()

# %% Visualizing the current on bulk

probe1 = Probe(solver.getSolvedProbeFilenames("Bulk probe1")[0])
probe2 = Probe(solver.getSolvedProbeFilenames("Bulk probe2")[0])
probe3 = Probe(solver.getSolvedProbeFilenames("Bulk probe3")[0])
probeTotal = Probe(solver.getSolvedProbeFilenames("Bulk probe total")[0])



plt.figure()
plt.plot(probe1["time"].to_numpy(), probe1["current"].to_numpy(), label='Current on bulk1', color='blue')
plt.plot(probe2["time"].to_numpy(), probe2["current"].to_numpy(), '--', label='Current on bulk2', color='red')
plt.plot(probe3["time"].to_numpy(), probe3["current"].to_numpy(), '--', label='Current on bulk3', color='green')
plt.grid(which='both')
plt.legend()
plt.show()

plt.figure()
plt.plot(probeTotal["time"].to_numpy(), probeTotal["current"].to_numpy(), label='Current on probe total', color='blue')
plt.plot(probe2["time"].to_numpy(), probe2["current"].to_numpy() + probe1["current"].to_numpy() + probe3["current"].to_numpy(), '--', label='Adition of current in bulk1, bulk2 and bulk3', color='red')
plt.grid(which='both')
plt.legend()
plt.show()

# %% Visualizing no currents in negative case

probe1_negative = Probe(solver_negative.getSolvedProbeFilenames("Bulk probe1")[0])
probe2_negative = Probe(solver_negative.getSolvedProbeFilenames("Bulk probe2")[0])
probe3_negative = Probe(solver_negative.getSolvedProbeFilenames("Bulk probe3")[0])

plt.figure()
plt.plot(probe1_negative["time"].to_numpy(), probe1_negative["current"].to_numpy(), label='Current on bulk1 negative', color='blue')
plt.plot(probe2_negative["time"].to_numpy(), probe2_negative["current"].to_numpy(), '--', label='Current on bulk2 negative', color='red')
plt.plot(probe3_negative["time"].to_numpy(), probe3_negative["current"].to_numpy(), '--', label='Current on bulk3 negative', color='green')
plt.grid(which='both')
plt.legend()
plt.show()

# %% Visualizing no currents in positive case

probe1_positive = Probe(solver_positive.getSolvedProbeFilenames("Bulk probe1")[0])
probe2_positive = Probe(solver_positive.getSolvedProbeFilenames("Bulk probe2")[0])
probe3_positive = Probe(solver_positive.getSolvedProbeFilenames("Bulk probe3")[0])

plt.figure()
plt.plot(probe1_positive["time"].to_numpy(), probe1_positive["current"].to_numpy(), label='Current on bulk1 positive', color='blue')
plt.plot(probe2_positive["time"].to_numpy(), probe2_positive["current"].to_numpy(), '--', label='Current on bulk2 positive', color='red')
plt.plot(probe3_positive["time"].to_numpy(), probe3_positive["current"].to_numpy(), '--', label='Current on bulk3 positive', color='green')
plt.grid(which='both')
plt.legend()
plt.show()
# %% Generating a new excitation file

# data = np.loadtxt('predefinedExcitation.1.exc')
# data[:, 1] *= 2
# np.savetxt('predefinedExcitation.2.exc', data, fmt='%.6e')
# %%
