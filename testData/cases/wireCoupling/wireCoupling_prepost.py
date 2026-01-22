# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import json
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %% Generate excitation and visualize
dt = 1e-13
w0 = 0.1e-8 
t0 = 20*w0
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
plt.xlim(0,5e8)
plt.ylim(1e-9, 1e6)
plt.yscale('log')

plt.grid()

# %% Prepare solver
solver = FDTD(input_filename = 'wireCoupling.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run() 
i = Probe(solver.getSolvedProbeFilenames("wire_start_Wx")[0])
ic = Probe(solver.getSolvedProbeFilenames("wire_connection_Wx")[0])
e = Probe(solver.getSolvedProbeFilenames("wire_mid_Ex")[0])

solver = FDTD(input_filename = 'wireCoupling_H.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run() 
iH = Probe(solver.getSolvedProbeFilenames("wire_start_Wx")[0])
iHc = Probe(solver.getSolvedProbeFilenames("wire_connection_Wx")[0])
eH = Probe(solver.getSolvedProbeFilenames("wire_mid_Ex")[0])

solver = FDTD(input_filename = 'wireCoupling_V.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run() 
iV = Probe(solver.getSolvedProbeFilenames("wire_start_Wx")[0])
iVc = Probe(solver.getSolvedProbeFilenames("wire_connection_Wx")[0])
eV = Probe(solver.getSolvedProbeFilenames("wire_mid_Ex")[0])

# %% Postprocess

t = i['time']
dt = t[1] - t[0]
fq  = fftfreq(len(t))/dt
If = fft(i['current'])
IfH = fft(iH['current'])
IfV = fft(iV['current'])
Ifc = fft(ic['current'])
IfHc = fft(iHc['current'])
IfVc = fft(iVc['current'])

fmin = 1e6
fmax = 10e9
idx_min = (np.abs(fq - fmin)).argmin()
idx_max = (np.abs(fq - fmax)).argmin()
f = fq[idx_min:idx_max]
Iff = If[idx_min:idx_max]
IffH = IfH[idx_min:idx_max]
IffV = IfV[idx_min:idx_max]
Iffc = Ifc[idx_min:idx_max]
IffHc = IfHc[idx_min:idx_max]
IffVc = IfVc[idx_min:idx_max]

plt.figure()
plt.plot(f, 20*np.log10(np.abs(Iff)),'.-', label='no plane')
plt.plot(f, 20*np.log10(np.abs(Iffc)),'.-', label='no plane, on conn')
plt.grid(which='both')
plt.xlabel('f [Hz]')
plt.ylabel('I [dB]')
plt.xlim(0.5e8,5e8)
plt.ylim(-100,10)
plt.legend()

plt.figure()
plt.plot(f, 20*np.log10(np.abs(IffH)),'.--', label='H plane')
plt.plot(f, 20*np.log10(np.abs(IffHc)),'.--', label='H plane, on conn')
plt.grid(which='both')
plt.xlabel('f [Hz]')
plt.ylabel('I [dB]')
plt.xlim(0.5e8,5e8)
plt.ylim(-100,10)
plt.legend()

plt.figure()
plt.plot(f, 20*np.log10(np.abs(IffV)),'.--', label='V plane')
plt.plot(f, 20*np.log10(np.abs(IffVc)),'.--', label='V plane, on conn')
plt.grid(which='both')
plt.xlabel('f [Hz]')
plt.ylabel('I [dB]')
plt.xlim(0.5e8,5e8)
plt.ylim(-100,10)
plt.legend()

t = e['time']
dt = t[1] - t[0]
fq  = fftfreq(len(t))/dt
ef = fft(e['field'])
efH = fft(eH['field'])
efV = fft(eV['field'])

fmin = 1e6
fmax = 10e9
idx_min = (np.abs(fq - fmin)).argmin()
idx_max = (np.abs(fq - fmax)).argmin()
f = fq[idx_min:idx_max]
eff = ef[idx_min:idx_max]
effH = efH[idx_min:idx_max]
effV = efV[idx_min:idx_max]

plt.figure()
plt.plot(f, np.abs(eff),'.-', label='no plane')
plt.grid(which='both')
plt.xlabel('f [Hz]')
plt.ylabel('V [V/m]')
plt.xlim(0.5e8,5e8)
plt.legend()

plt.figure()
plt.plot(f, np.abs(effH),'.--', label='H plane')
plt.grid(which='both')
plt.xlabel('f [Hz]')
plt.ylabel('V [V/m]')
plt.xlim(0.5e8,5e8)
plt.legend()

plt.figure()
plt.plot(f, np.abs(effV),'.--', label='V plane')
plt.grid(which='both')
plt.xlabel('f [Hz]')
plt.ylabel('V [V/m]')
plt.xlim(0.5e8,5e8)
plt.legend()


# %%
