# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt

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


# %% Run solver
fterminal = 'simple_loop_terminal.fdtd.json'
flumped = 'simple_loop_lumped.fdtd.json'
generateGaussExcitation()
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


#%%  Theoretical initial current
InitialBulk_probe = Probe(solver_terminal.getSolvedProbeFilenames("Initial current")[0])

dt = time[1] - time[0]
fq = fftfreq(len(time))/dt

V_in_f = fft(V_in)
I_teo_f = V_in_f/(1000 + 1j*2*np.pi*fq*5.3e-6)
# I_teo_f = V_in_f/(1000 + 1j*2*np.pi*fq*5.3e-20)

I_teo = ifft(I_teo_f)

plt.figure()
plt.plot(InitialBulk_probe['time'].to_numpy(), InitialBulk_probe['current'].to_numpy(), label='Initial current on bulk', color='blue')
plt.plot(time, I_teo, '--', label='Theoretical current', color='black')
plt.title('Initial current')
plt.xlabel('Time')
plt.ylabel('Current')
plt.legend()
plt.grid(which='both')

#%% Comparison of currents between terminals and lumped
TerminalBulk_probe = Probe(solver_terminal.getSolvedProbeFilenames("Material current")[0])
LumpedBulk_probe = Probe(solver_lumped.getSolvedProbeFilenames("Material current")[0])

plt.figure()
plt.plot(TerminalBulk_probe['time'].to_numpy(), TerminalBulk_probe['current'].to_numpy(), label='Current on the circuit with a terminal', color='green')
plt.plot(LumpedBulk_probe['time'].to_numpy(), LumpedBulk_probe['current'].to_numpy(), '--', label='Current on the circuit with an equivalent lumped line', color='red')
# plt.plot(time, I_teo, '--', label='Theoretical current', color='black')
plt.title('Current on the circuit')
plt.xlabel('Time')
# plt.xlim((0, 2e-9))
plt.legend()
plt.grid(which='both')

plt.tight_layout()
plt.show()
# %% Frequency Domain

t = InitialBulk_probe['time']
dt = t[1] - t[0]

fq = fftfreq(len(t))/dt

fmin = 1e4
fmax = 6e8
idx_min = (np.abs(fq - fmin)).argmin()
idx_max = (np.abs(fq - fmax)).argmin()

I_initial_f = fft(InitialBulk_probe['current'])[idx_min:idx_max]
I_terminal_f = fft(TerminalBulk_probe['current'])[idx_min:idx_max]
# I_lumped_f = fft(LumpedBulk_probe['current'])[idx_min:idx_max]




plt.figure()
# plt.plot(fq[idx_min:idx_max], I_lumped_f, '.-', label='Current on the circuit with a lumped line')
plt.plot(fq[idx_min:idx_max], I_terminal_f, '.--', label='Current on the circuit with a terminal')
plt.xlabel('Frequency')
plt.grid(which='both')
plt.xlim(fmin, fmax)
# plt.ylim((-1, 5))
plt.xscale('log')
plt.legend()
plt.show()

# %%