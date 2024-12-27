# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants
from skrf.media import Freespace
from skrf.frequency import Frequency


import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %% Generate excitation and visualize
dx = 20e-3
cfl = 0.25
dt = cfl * np.sqrt(3*dx**3)/scipy.constants.speed_of_light

w0 = 0.1e-9    # ~ 2 GHz bandwidth
t0 = 10*w0
t = np.arange(0, t0+20*w0, dt)
f = np.exp( -np.power(t-t0,2)/ w0**2 )

plt.figure()
plt.plot(t,f)

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = f
np.savetxt('gauss.exc', data)
 
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
solver.cleanUp()
solver.run()
assert solver.hasFinishedSuccessfully()

# %% Postprocess
back = Probe(solver.getSolvedProbeFilenames("back")[0])
plt.plot(back.df['time'], back.df['incident'], label='incident')    
plt.plot(back.df['time'],  back.df['field'], label='back')
plt.legend()
plt.xlim(0,3e-9)


t = back.df['time']
dt = t[1] - t[0]
fq  = fftfreq(len(t))/dt
INC = fft(back.df['incident'])
BACK  = fft(back.df['field'])
S21 = BACK/INC

fmin = 8e6
fmax = 4e9
idx_min = (np.abs(fq - fmin)).argmin()
idx_max = (np.abs(fq - fmax)).argmin()
f = fq[idx_min:idx_max]
fdtd_s21 = S21[idx_min:idx_max]

freq = Frequency.from_f(f, unit='Hz')
air =  Freespace(freq)

sigma = 100
width = 10e-3

mat_ep_r = (1+sigma/(1j*freq.w*scipy.constants.epsilon_0))
conductive_material = Freespace(freq, ep_r=mat_ep_r)

slab = air.thru() ** conductive_material.line(width, unit='m') ** air.thru()

plt.figure()
plt.plot(f, 20*np.log10(np.abs(fdtd_s21)),'.-', label='FDTD')
plt.plot(f, 20*np.log10(np.abs(slab.s[:,0,1])),'.-', label='Analytical')
plt.xlim(fmin,fmax)
plt.grid(which='both')
plt.xscale('log')
plt.legend()

# %%
