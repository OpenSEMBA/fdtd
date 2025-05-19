# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
from scipy import signal

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../../build/bin/semba-fdtd'

from pyWrapper import *

# %% Generate excitation and visualize
def generateGaussExcitation():
    dt = 0.8e-12
    w0 = 1e-9 
    t0 = 10*w0
    t = np.arange(0, t0+20*w0, dt)
    f = np.exp( -np.power(t-t0,2)/ w0**2 )

    data = np.zeros((len(t), 2))
    data[:,0] = t
    data[:,1] = f
    np.savetxt('predefinedExcitation.1.exc', data)

# %% Generate excitation and visualize
def generateRampExcitation():
    dt = 0.8e-12
    t_final = 2e-8 
    t = np.arange(0, t_final, dt)

    f = np.zeros_like(t)

    t1 = 2e-9
    t2 = 16e-9

    ramp_region = (t >= t1) & (t < t2)
    f[ramp_region] = (t[ramp_region] - t1) / (t2 - t1)
    f[t >= t2] = 1

    # Guardar en archivo
    data = np.column_stack((t, f))
    np.savetxt('rampExcitation.exc', data)

# %% Run solver. 
# The file current_bifurcation_RC.fdtd.json contains an equivalent case to comparison but not effective.

# generateGaussExcitation()
generateRampExcitation()

# fterminal = 'current_bifurcation_RC.fdtd.json'
# solver_terminal = FDTD(input_filename = fterminal, path_to_exe=SEMBA_EXE)

flumped = 'simple_loop_lumped.fdtd.json'
solver_lumped = FDTD(input_filename = flumped, path_to_exe=SEMBA_EXE)

solver_lumped.cleanUp()

# solver_terminal.run()
solver_lumped.run()

# %% Visualizing initial values of currents and voltages

V_in = np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=1)
time = np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=0)
plt.figure()
plt.plot(time, V_in, label='Initial excitation voltage')    
plt.grid(which='both')
plt.legend()
plt.show()

#%%  Theoretical initial current
# InitialTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("Initial current")[0])
InitialLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("Initial current")[0])

R = solver_lumped.getMaterialProperties("lumped_RC")["resistance"]
C = solver_lumped.getMaterialProperties("lumped_RC")["capacitance"]
L = 1.65e-7     # Parasitic inductance of the system. See calculation at the lumped resistor prepost.

num = [R*C, 1]
den = [L*R*C, L, R]
system = signal.TransferFunction(num, den)
tout, I_theo, _ = signal.lsim(system, U=V_in, T=time)

plt.figure()
# plt.plot(InitialTerminal_probe['time'].to_numpy(), InitialTerminal_probe['current'].to_numpy(), label='On terminal case', color='green')
plt.plot(InitialLumped_probe['time'].to_numpy(), InitialLumped_probe['current'].to_numpy(), '--', label='On lumped case', color='red')
plt.plot(time, I_theo, '--', label='Theoretical current', color='black')
plt.title('Initial current')
plt.xlabel('Time')
plt.ylabel('Current')
plt.legend()
plt.grid(which='both')

#%% Comparison of currents between terminals and lumped
StartLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellStart")[0])
EndLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellEnd")[0])

AdjacentPostLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PostLumpedCell")[0])
# AdjacentPostTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("PostTerminalCell")[0])

AdjacentPreLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PreLumpedCell")[0])
# AdjacentPreTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("PreTerminalCell")[0])

plt.figure()
# plt.plot(AdjacentPreTerminalProbe['time'].to_numpy(), AdjacentPreTerminalProbe['current'].to_numpy(), label='PreTerminalCell', color='green')
plt.plot(AdjacentPreLumpedProbe['time'].to_numpy(), AdjacentPreLumpedProbe['current'].to_numpy(), '--', label='PreLumpedCell', color='red')
plt.plot(time, I_theo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit before the lumped/terminal cell')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(StartLumpedProbe['time'].to_numpy(), StartLumpedProbe['current'].to_numpy(), '--', label='LumpedCellStart', color='red')
plt.plot(EndLumpedProbe['time'].to_numpy(), EndLumpedProbe['current'].to_numpy(), '--', label='LumpedCellEnd', color='red')
plt.plot(time, I_theo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit at the start/end of the lumped cell')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
# plt.plot(AdjacentPostTerminalProbe['time'].to_numpy(), AdjacentPostTerminalProbe['current'].to_numpy(), label='PostTerminalCell', color='green')
plt.plot(AdjacentPostLumpedProbe['time'].to_numpy(), AdjacentPostLumpedProbe['current'].to_numpy(), '--', label='PostLumpedCell', color='red')
plt.plot(time, I_theo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit after the lumped/terminal cell')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()


# %%
