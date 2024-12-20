# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants

# %%
dx = 2.5e-3
cfl = 0.25
dt = cfl * np.sqrt(3*dx**3)/scipy.constants.speed_of_light

t0 = 0.696e-9
w0 = 0.187e-9
t = np.arange(0, t0+20*w0, dt)
f = np.exp( -np.power(t-t0,2)/ w0**2 )
plt.plot(t,f)
data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = f
np.savetxt('gauss_1GHz.exc', data)
 
# %% 
fq = fftfreq(len(t))/dt
F = fft(f)
Fmax = np.max(np.abs(F))
Fcut = Fmax / np.sqrt(2.0)
idx = (np.abs(np.abs(F) - Fcut)).argmin()


plt.plot(fq, np.abs(F),'.-')
plt.axvline(x = fq[idx], color = 'r')
plt.xlim(0,5e9)
plt.yscale('log')

plt.grid()
