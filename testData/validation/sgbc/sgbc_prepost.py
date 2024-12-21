# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %% Generate excitation and visualize
dx = 2.5e-3
cfl = 0.25
dt = cfl * np.sqrt(3*dx**3)/scipy.constants.speed_of_light

w0 = 0.187e-9
t0 = 10*w0
t = np.arange(0, t0+20*w0, dt)
f = np.exp( -np.power(t-t0,2)/ w0**2 )

plt.figure()
plt.plot(t,f)

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = f
np.savetxt('gauss_1GHz.exc', data)
 
fq = fftfreq(len(t))/dt
F = fft(f)
Fmax = np.max(np.abs(F))
Fcut = Fmax / np.sqrt(2.0)
idx = (np.abs(np.abs(F) - Fcut)).argmin()

plt.figure()
plt.plot(fq, np.abs(F),'.-')
plt.axvline(x = fq[idx], color = 'r')
plt.xlim(0,5e9)
plt.ylim(1e-9, 1e6)
plt.yscale('log')

plt.grid()

# %% Run solver
fn = 'shieldingEffectiveness.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.run()
assert solver.hasFinishedSuccessfully()

# %% Postprocess
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)

front = None
for pf in solver.getSolvedProbeFilenames("front"):
    if Probe(pf).direction == 'x':
        front = Probe(pf)
        break        

back = None 
for pf in solver.getSolvedProbeFilenames("back"):
    if Probe(pf).direction == 'x':
        back = Probe(pf)
        break        

plt.plot(front.df['time'], front.df['incident'])    
plt.plot(front.df['time'], front.df['field'])
plt.plot(back.df['time'], back.df['field'])


# %%
