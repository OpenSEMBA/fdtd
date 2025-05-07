# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../../build/bin/semba-fdtd'

from pyWrapper import *


# %% Run solver
fterminal = 'simple_loop_terminal.fdtd.json'
flumped = 'simple_loop_lumped.fdtd.json'
# generateGaussExcitation()
solver_terminal = FDTD(input_filename = fterminal, path_to_exe=SEMBA_EXE)
solver_lumped = FDTD(input_filename = flumped, path_to_exe=SEMBA_EXE)
solver_terminal.cleanUp()

solver_terminal.run()
solver_lumped.run()

# %% Visualizing initial values of currents and voltages

V_in = np.loadtxt("predefinedExcitation.1.exc", usecols=1)
time = np.loadtxt("predefinedExcitation.1.exc", usecols=0)
plt.figure()
plt.plot(time, V_in, label='Initial excitation voltage')    
plt.grid(which='both')
plt.legend()
plt.show()

InitialBulk_probe = Probe(solver_terminal.getSolvedProbeFilenames("Initial current")[0])
TerminalBulk_probe = Probe(solver_terminal.getSolvedProbeFilenames("Material current")[0])
LumpedBulk_probe = Probe(solver_lumped.getSolvedProbeFilenames("Material current")[0])

plt.figure()
plt.plot(InitialBulk_probe['time'].to_numpy(), InitialBulk_probe['current'].to_numpy(), label='Initial current on bulk', color='blue')
plt.plot(time, V_in/1000, '--', label='Theoretical current with a 1000 ohms resistor', color='black')
plt.title('Initial current')
plt.xlabel('Time')
plt.ylabel('Current')
plt.legend()
plt.grid(which='both')

plt.figure()
plt.plot(TerminalBulk_probe['time'].to_numpy(), TerminalBulk_probe['current'].to_numpy(), label='Current on the circuit with a terminal', color='green')
plt.plot(LumpedBulk_probe['time'].to_numpy(), LumpedBulk_probe['current'].to_numpy(), '--', label='Current on the with an equivalent lumped line', color='red')
plt.title('Current on the circuit')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')

plt.tight_layout()
plt.show()
# %%
