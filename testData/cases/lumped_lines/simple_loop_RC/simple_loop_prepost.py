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
    t2 = 6e-9

    ramp_region = (t >= t1) & (t < t2)
    f[ramp_region] = (t[ramp_region] - t1) / (t2 - t1)
    f[t >= t2] = 1

    # Guardar en archivo
    data = np.column_stack((t, f))
    np.savetxt('rampExcitation.exc', data)

# %% Run solver
fterminal = 'current_bifurcation_RC.fdtd.json'
flumped = 'simple_loop_lumped.fdtd.json'
# generateGaussExcitation()
generateRampExcitation()
solver_terminal = FDTD(input_filename = fterminal, path_to_exe=SEMBA_EXE)
solver_lumped = FDTD(input_filename = flumped, path_to_exe=SEMBA_EXE)
solver_lumped.cleanUp()

solver_terminal.run()
solver_lumped.run()

# %% Visualizing initial values of currents and voltages

V_in = np.loadtxt("rampExcitation.exc", usecols=1)
time = np.loadtxt("rampExcitation.exc", usecols=0)
plt.figure()
plt.plot(time, V_in, label='Initial excitation voltage')    
plt.grid(which='both')
plt.legend()
plt.show()

#%%  Theoretical initial current
InitialTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("Initial current")[0])
InitialLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("Initial current")[0])

R = 1000
C = 1e-12
L = 1.65e-7

num = [R*C, 1]
den = [L*R*C, L, R]
system = signal.TransferFunction(num, den)
tout, I_teo, _ = signal.lsim(system, U=V_in, T=time)

plt.figure()
# plt.plot(InitialTerminal_probe['time'].to_numpy(), InitialTerminal_probe['current'].to_numpy(), label='On terminal case', color='green')
plt.plot(InitialLumped_probe['time'].to_numpy(), InitialLumped_probe['current'].to_numpy(), '--', label='On lumped case', color='red')
plt.plot(time, I_teo, '--', label='Theoretical current', color='black')
plt.title('Initial current')
plt.xlabel('Time')
plt.ylabel('Current')
plt.legend()
plt.grid(which='both')

#%% Comparison of currents between terminals and lumped
StartLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("Start Material current")[0])
EndLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("End Material current")[0])

AfterLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("After Material current")[0])
AfterTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("After Material current")[0])

BeforeLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("Before Material current")[0])
BeforeTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("Before Material current")[0])

plt.figure()
# plt.plot(BeforeTerminal_probe['time'].to_numpy(), BeforeTerminal_probe['current'].to_numpy(), label='Before the terminal', color='green')
plt.plot(BeforeLumped_probe['time'].to_numpy(), BeforeLumped_probe['current'].to_numpy(), '--', label='Before the lumped line', color='red')
plt.plot(time, I_teo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(StartLumped_probe['time'].to_numpy(), StartLumped_probe['current'].to_numpy(), label='Start of the lumped line cell', color='blue')
plt.plot(EndLumped_probe['time'].to_numpy(), EndLumped_probe['current'].to_numpy(), '--', label='End of the lumped line cell', color='red')
plt.plot(time, I_teo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
# plt.plot(AfterTerminal_probe['time'].to_numpy(), AfterTerminal_probe['current'].to_numpy(), label='After the terminal', color='green')
plt.plot(AfterLumped_probe['time'].to_numpy(), AfterLumped_probe['current'].to_numpy(), '--', label='After the lumped line', color='red')
plt.plot(time, I_teo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()


# %%
