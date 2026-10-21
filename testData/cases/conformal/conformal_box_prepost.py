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
dt = 1e-12
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
F = F[fq>=0]
fq = fq[fq>=0]
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
fn = 'conformal_box.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver["mesh"]["elements"][3]["subtype"] = "volume"
solver.cleanUp()
solver.run()

back_volume = Probe(solver.getSolvedProbeFolders("back")[0])
front_volume = Probe(solver.getSolvedProbeFolders("front")[0])

# %%
solver["mesh"]["elements"][3]["subtype"] = "surface"
solver.cleanUp()
solver.run()

back_surface = Probe(solver.getSolvedProbeFolders("back")[0])
front_surface = Probe(solver.getSolvedProbeFolders("front")[0])

t = back_surface['time']

# %% Postprocess

plt.figure()
plt.plot(t, back_surface['incident'], label="incident")
plt.plot(t, back_surface['field'], label="transmitted S")
plt.plot(t, back_volume['field'], '--', label="transmitted V")
plt.plot(t, front_surface['field'], label="reflected S")
plt.plot(t, front_volume['field'], '--', label="reflected V")
plt.grid(which='both')
plt.legend()

# %%
