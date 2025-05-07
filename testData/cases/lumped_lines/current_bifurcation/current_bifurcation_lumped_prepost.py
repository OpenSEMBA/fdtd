# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../../build/bin/semba-fdtd'

from pyWrapper import *

# # %% Generate excitation and visualize
# def generateGaussExcitation():
#     dt = 1e-12
#     w0 = 0.5e-8 
#     t0 = 10*w0
#     t = np.arange(0, t0+20*w0, dt)
#     f = np.exp( -np.power(t-t0,2)/ w0**2 )

#     data = np.zeros((len(t), 2))
#     data[:,0] = t
#     data[:,1] = f
#     np.savetxt('gauss.exc', data)

# %% Run solver
fn = 'current_bifurcation_lumped.fdtd.json'
# generateGaussExcitation()
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

# %% Visualizing initial values of currents and voltages

V_in = np.loadtxt("predefinedExcitation.1.exc", usecols=1)
time = np.loadtxt("predefinedExcitation.1.exc", usecols=0)
plt.figure()
plt.plot(time, V_in, label='Initial excitation voltage')    
plt.grid(which='both')
plt.legend()
plt.show()

InitialBulk_probe = Probe(solver.getSolvedProbeFilenames("Bulk Initial probe")[0])
TopBulk_probe = Probe(solver.getSolvedProbeFilenames("Bulk Top probe")[0])
BottomBulk_probe = Probe(solver.getSolvedProbeFilenames("Bulk Bottom probe")[0])

plt.figure()
plt.plot(InitialBulk_probe['time'].to_numpy(), InitialBulk_probe['current'].to_numpy(), label='Initial current on bulk', color='blue')
plt.plot(time, V_in/500, '--', label='Theoretical current with two 1000 ohm resistors in parallel', color='black')
plt.title('Initial current')
plt.xlabel('Time')
plt.ylabel('Current')
plt.legend()
plt.grid(which='both')

plt.figure()
plt.plot(TopBulk_probe['time'].to_numpy(), TopBulk_probe['current'].to_numpy(), label='Current on resistor', color='green')
plt.plot(BottomBulk_probe['time'].to_numpy(), BottomBulk_probe['current'].to_numpy(), '--', label='Current on lumped line', color='red')
plt.title('Current on the circuit')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')

plt.tight_layout()
plt.show()

# %%
