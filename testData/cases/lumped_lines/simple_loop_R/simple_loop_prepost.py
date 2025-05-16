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

I_teo = np.zeros_like(time)
t1 = 1e-9
t2 = 15e-9
R = solver_lumped["materials"][1]["resistance"]
L = 1.65e-7

ramp_region = (time >= t1) & (time < t2)
I_teo[ramp_region] = (time[ramp_region] - t1) / (t2 - t1)/R  + (L/R**2)*(np.exp(-R*(time[ramp_region] - t1)/L) - 1)/(t2 - t1)
I_teo[time >= t2] = 1/R + (L/R**2)*(np.exp(-R*(time[time >= t2] - t1)/L) - np.exp(-R*(time[time >= t2] - t2)/L))/(t2 - t1)

# # Equivalent result for theoretical current but using scipy
# num = [1]
# den = [L, R]
# system = signal.TransferFunction(num, den)
# tout, I_out, _ = signal.lsim(system, U=V_in, T=time)

InitialTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("Initial current")[0])
InitialLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("Initial current")[0])

# # Calculation of the proper inductance of the system in case of R=0.

# time_probe = InitialTerminal_probe['time'].to_numpy()
# V_in_interp_func = interp1d(time, V_in, bounds_error=False, fill_value="extrapolate")
# V_in_interp = V_in_interp_func(time_probe)  
# mask = (InitialTerminal_probe['time'] >= t1) & (InitialTerminal_probe['time'] <= t2)
# current_slice = InitialTerminal_probe['current'][mask]
# time_slice = InitialTerminal_probe['time'][mask]
# vin_slice = V_in_interp[mask]
# dI_dt = np.gradient(current_slice, time_slice)

# with np.errstate(divide='ignore', invalid='ignore'):
#     L_estimated = np.where(dI_dt != 0, vin_slice / dI_dt, np.nan)

# L_mean = np.nanmean(L_estimated)

plt.figure()
plt.plot(InitialTerminal_probe['time'].to_numpy(), InitialTerminal_probe['current'].to_numpy(), label='On terminal case', color='green')
plt.plot(InitialLumped_probe['time'].to_numpy(), InitialLumped_probe['current'].to_numpy(), '--', label='On lumped case', color='red')
plt.plot(time, I_teo, '--', label='Theoretical current', color='black')
plt.title('Initial current')
plt.xlabel('Time')
plt.ylabel('Current')
plt.legend()
plt.grid(which='both')

#%% Comparison of currents between terminals and lumped
StartTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("Start Material current")[0])
StartLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("Start Material current")[0])

EndTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("End Material current")[0])
EndLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("End Material current")[0])

AfterLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("After Material current")[0])
AfterTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("After Material current")[0])

BeforeLumped_probe = Probe(solver_lumped.getSolvedProbeFilenames("Before Material current")[0])
BeforeTerminal_probe = Probe(solver_terminal.getSolvedProbeFilenames("Before Material current")[0])

plt.figure()
plt.plot(BeforeTerminal_probe['time'].to_numpy(), BeforeTerminal_probe['current'].to_numpy(), label='Before the terminal', color='green')
plt.plot(BeforeLumped_probe['time'].to_numpy(), BeforeLumped_probe['current'].to_numpy(), '--', label='Before the lumped line', color='red')
plt.title('Current on the circuit')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(StartTerminal_probe['time'].to_numpy(), StartTerminal_probe['current'].to_numpy(), label='Start of the terminal cell', color='green')
plt.plot(StartLumped_probe['time'].to_numpy(), StartLumped_probe['current'].to_numpy(), '--', label='Start of the lumped line cell', color='red')
# plt.plot(time, I_teo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(EndTerminal_probe['time'].to_numpy(), EndTerminal_probe['current'].to_numpy(), label='End of the terminal cell', color='green')
plt.plot(EndLumped_probe['time'].to_numpy(), EndLumped_probe['current'].to_numpy(), '--', label='End of the lumped line cell', color='red')
# plt.plot(time, I_teo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(AfterTerminal_probe['time'].to_numpy(), AfterTerminal_probe['current'].to_numpy(), label='After the terminal', color='green')
plt.plot(AfterLumped_probe['time'].to_numpy(), AfterLumped_probe['current'].to_numpy(), '--', label='After the lumped line', color='red')
plt.title('Current on the circuit')
plt.xlabel('Time')
plt.legend()
plt.grid(which='both')
plt.tight_layout()
plt.show()
# %%
