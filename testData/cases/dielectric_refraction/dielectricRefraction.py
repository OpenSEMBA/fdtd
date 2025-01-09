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
# %%
dt = 1e-11
w0 = 0.5e-9
t0 = 10*w0
t = np.arange(0, t0+2*t0, dt)
f = np.exp( -np.power(t-t0,2)/ w0**2 )

#plt.figure()
#plt.plot(t,f)

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = f
np.savetxt('gauss.exc', data)

#%%
fq = fftfreq(len(t))/dt
F = fft(f)
#%%
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
# %%
fn = 'dielectricRefraction.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()
assert solver.hasFinishedSuccessfully()

#%%
outside = Probe(solver.getSolvedProbeFilenames("outside")[0]) 
plt.plot(outside.df['time'],  outside.df['field'], label='outside')
plt.legend()
plt.xlim(0,1e-7)
#%%
inside = Probe(solver.getSolvedProbeFilenames("inside")[0])   
plt.plot(inside.df['time'], inside.df['field'], label='inside')
plt.legend()
plt.xlim(0,1e-7)


# %%
getMaxTimeOutside = outside.df['field'].argmin()
getMaxTimeOutside
# %%
outside.df['field'][getMaxTimeOutside]
# %%
outside.df['time'][getMaxTimeOutside]
# %%
getMaxTimeInside = inside.df['field'].argmin()
inside.df['field'][getMaxTimeInside]
inside.df['time'][getMaxTimeInside]

# %%
