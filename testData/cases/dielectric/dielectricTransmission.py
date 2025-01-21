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
dt = 1e-12
w0 = 0.1e-8
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
fwave = np.sin(w0*t)

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = fwave
np.savetxt('wave.exc', data)

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
fn = 'dielectricTransmission.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

#%%
outside = Probe(solver.getSolvedProbeFilenames("outside")[0]) 
plt.plot(outside["time"].to_numpy(),  outside["field"].to_numpy(), label='outside')
plt.xlim(0, 5e-8)
plt.legend()
#%%
inside = Probe(solver.getSolvedProbeFilenames("inside")[0])   
plt.plot(inside["time"].to_numpy(), inside["field"].to_numpy(), label='inside')
plt.xlim(0, 5e-8)
plt.legend()

# %%
incidentPulseIndex = outside['field'].argmin()
reflectedPulseIndex = outside['field'].argmax()
refractedPulseIndex = inside['field'].argmin()

# %%
incidentPulse = outside['field'][incidentPulseIndex]
reflectedPulse = outside['field'][reflectedPulseIndex]
refractedPulse = inside['field'][refractedPulseIndex]
pulses = [
    incidentPulse,
    reflectedPulse,
    refractedPulse
]

incidentTime = outside["time"][incidentPulseIndex]
reflectedTime = outside["time"][reflectedPulseIndex]
refractedTime = inside["time"][refractedPulseIndex]

times = [
    incidentTime,
    reflectedTime,
    refractedTime
]

for i in range(len(pulses)):
    print(pulses[i])
    print(times[i])
    print()

# %%

print(pulses[0] + pulses[1] - pulses[2])
# %%

timeToDielectricSurface = (times[1] - times[0])/2 + times[0]

reflectedDelay = times[1] - timeToDielectricSurface
refractedDelay = times[2] - timeToDielectricSurface
print(reflectedDelay)
print(refractedDelay)
print(reflectedDelay/refractedDelay)
# %%
