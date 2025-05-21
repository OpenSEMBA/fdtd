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

# %% Generate excitation and visualize
def generateRampExcitation():
    dt = 0.8e-12
    t_final = 20e-9 
    t = np.arange(0, t_final, dt)

    f = np.zeros_like(t)

    t1 = 1e-9
    t2 = 15e-9

    ramp_region = (t >= t1) & (t < t2)
    f[ramp_region] = (t[ramp_region] - t1) / (t2 - t1)
    f[t >= t2] = 1

    # Guardar en archivo
    data = np.column_stack((t, f))
    np.savetxt('rampExcitation.exc', data)


# %% Run solver
fterminal = 'simple_loop_terminal.fdtd.json'
flumped = 'simple_loop_lumped.fdtd.json'
generateRampExcitation()
solver_terminal = FDTD(input_filename = fterminal, path_to_exe=SEMBA_EXE)
solver_lumped = FDTD(input_filename = flumped, path_to_exe=SEMBA_EXE)
solver_terminal.cleanUp()

solver_terminal.run()
solver_lumped.run()

# %% Visualizing initial values of currents and voltages

V_in = np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=1)
time = np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=0)
plt.figure()
plt.plot(time, V_in, label='Initial excitation voltage')    
plt.grid(which='both')
plt.legend()
plt.show()

# %% Comparing initial current with theoretical

R = solver_lumped.getMaterialProperties("lumped_line")["resistance"]
L = 1.65e-7     # Parasitic inductance of the system. See calculation at the end. 

# # Theoretical current. Analytic solution.
# I_theo = np.zeros_like(time)
# t1 = 1e-9
# t2 = 15e-9
# ramp_region = (time >= t1) & (time < t2)
# I_theo[ramp_region] = (time[ramp_region] - t1) / (t2 - t1)/R  + (L/R**2)*(np.exp(-R*(time[ramp_region] - t1)/L) - 1)/(t2 - t1)
# I_theo[time >= t2] = 1/R + (L/R**2)*(np.exp(-R*(time[time >= t2] - t1)/L) - np.exp(-R*(time[time >= t2] - t2)/L))/(t2 - t1)

# Theoretical current using the transfer function and the initial voltage by laplace transform.
num = [1]
den = [L, R]
system = signal.TransferFunction(num, den)
tout, I_theo, _ = signal.lsim(system, U=V_in, T=time)

InitialTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("Initial current")[0])
InitialLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("Initial current")[0])

plt.figure()
plt.plot(InitialTerminal_probe['time'].to_numpy(), InitialTerminal_probe['current'].to_numpy(), label='Terminal case', color='green')
plt.plot(InitialLumped_probe['time'].to_numpy(), InitialLumped_probe['current'].to_numpy(), '--', label='Lumped case', color='red')
plt.plot(time, I_theo, '--', label='Theoretical current', color='black')
plt.title('Initial current at the source')
plt.xlabel('Time')
plt.ylabel('Current')
plt.legend()
plt.grid(which='both')

#%% Comparison of currents between terminals and lumped
StartTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("TerminalCellStart")[0])
StartLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellStart")[0])

EndTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("TerminalCellEnd")[0])
EndLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellEnd")[0])

AdjacentPostLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PostLumpedCell")[0])
AdjacentPostTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("PostTerminalCell")[0])

AdjacentPreLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PreLumpedCell")[0])
AdjacentPreTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("PreTerminalCell")[0])

plt.figure()
plt.plot(AdjacentPreTerminalProbe['time'].to_numpy(), AdjacentPreTerminalProbe['current'].to_numpy(), label='PreTerminalCell', color='green')
plt.plot(AdjacentPreLumpedProbe['time'].to_numpy(), AdjacentPreLumpedProbe['current'].to_numpy(), '--', label='PreLumpedCell', color='red')
plt.title('Current on the circuit before the lumped/terminal cell')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(StartTerminalProbe['time'].to_numpy(), StartTerminalProbe['current'].to_numpy(), label='TerminalCellStart', color='green')
plt.plot(StartLumpedProbe['time'].to_numpy(), StartLumpedProbe['current'].to_numpy(), '--', label='LumpedCellStart', color='red')
plt.plot(time, I_theo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit at the start of the lumped/terminal cell')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(EndTerminalProbe['time'].to_numpy(), EndTerminalProbe['current'].to_numpy(), label='TerminalCellEnd', color='green')
plt.plot(EndLumpedProbe['time'].to_numpy(), EndLumpedProbe['current'].to_numpy(), '--', label='LumpedCellEnd', color='red')
plt.plot(time, I_theo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit at the end of the lumped/terminal cell')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(AdjacentPostTerminalProbe['time'].to_numpy(), AdjacentPostTerminalProbe['current'].to_numpy(), label='PostTerminalCell', color='green')
plt.plot(AdjacentPostLumpedProbe['time'].to_numpy(), AdjacentPostLumpedProbe['current'].to_numpy(), '--', label='PostLumpedCell', color='red')
plt.title('Current on the circuit after the lumped/terminal cell')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()
# %% Calculation of the proper/parasitic inductance of the system. Only valid for the terminal case with R=0.

t1 = 1e-9       # Start of the ramp
t2 = 15e-9      # End of the ramp

time_probe = InitialTerminal_probe['time'].to_numpy()
V_in_interp_func = interp1d(time, V_in, bounds_error=False, fill_value="extrapolate")
V_in_interp = V_in_interp_func(time_probe)  
mask = (InitialTerminal_probe['time'] >= t1) & (InitialTerminal_probe['time'] <= t2)
current_slice = InitialTerminal_probe['current'][mask]
time_slice = InitialTerminal_probe['time'][mask]
vin_slice = V_in_interp[mask]
dI_dt = np.gradient(current_slice, time_slice)

with np.errstate(divide='ignore', invalid='ignore'):
    L_estimated = np.where(dI_dt != 0, vin_slice / dI_dt, np.nan)

L_mean = np.nanmean(L_estimated)
# %%
