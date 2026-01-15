# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants

import pandas as pd

import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %% Read measurement files
solver = FDTD(input_filename = 'conformal_inductance_cylinder_conformal.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()
bulk_conf = Probe(solver.getSolvedProbeFilenames("BulkProbe")[0])

solver = FDTD(input_filename = 'conformal_inductance_cylinder_staircase.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()
bulk = Probe(solver.getSolvedProbeFilenames("BulkProbe")[0])

# %%
exc_file = "predefinedExcitation.1.exc"
exc = pd.read_csv(exc_file, sep='\\s+')
exc = exc.rename(columns={
    exc.columns[0]: 'time',
    exc.columns[1]: 'V'
})
fig, axs = plt.subplots(2, 1,figsize=(10, 10))
axs[0].plot(exc['time']*1e6, exc['V'], 'r', label = 'excitation')
axs[0].set_xlabel('time [us]')
axs[0].set_ylabel('V [V]')

axs[1].plot(bulk['time']*1e6, -bulk['current'], 'r', label = 'tessellator, staircase')
axs[1].plot(bulk_conf['time']*1e6, -bulk_conf['current'], 'b--', label = 'tessellator, conformal')
axs[1].legend()
axs[1].set_xlabel('time [us]')
axs[1].set_ylabel('I [A]')


# %% Plot discrete fourier transform of measured currents
new_freqs = np.geomspace(1e3, 1e7, num=100)
Ibulk = bulk["current"].to_numpy()
tbulk = bulk["time"].to_numpy()
dt_bulk = tbulk[1]-tbulk[0]
Ifbulk = dt_bulk*np.array([np.sum(Ibulk * np.exp(-1j * 2 * np.pi * f * tbulk)) for f in new_freqs])

Ibulk_conf = bulk_conf["current"].to_numpy()
tbulk_conf = bulk_conf["time"].to_numpy()
dt_bulk_conf = tbulk_conf[1]-tbulk_conf[0]
Ifbulk_conf = dt_bulk_conf*np.array([np.sum(Ibulk_conf * np.exp(-1j * 2 * np.pi * f * tbulk_conf)) for f in new_freqs])

Vexc = exc["V"].to_numpy()
texc = exc["time"].to_numpy()
dt_exc = texc[1]-texc[0]
Vfexc = dt_exc*np.array([np.sum(Vexc * np.exp(-1j * 2 * np.pi * f * texc)) for f in new_freqs])

plt.figure()
plt.loglog(new_freqs, np.abs(Vfexc), "b-", label = 'excitation')
plt.loglog(new_freqs, np.abs(Ifbulk), "r", label = 'current, tessellator staircase')
plt.loglog(new_freqs, np.abs(Ifbulk_conf), "k--", label = 'current, tessellator eP ')
plt.xlabel('f [Hz]')
plt.ylabel('Current')
plt.grid('both')
plt.legend()
plt.legend()

# %%
plt.figure()
plt.loglog(new_freqs, np.abs(Vfexc/Ifbulk), "b", label='tessellator, staircase')
plt.loglog(new_freqs, np.abs(Vfexc/Ifbulk_conf), "k.", label='tessellator, conformal')
plt.xlabel('f [Hz]')
plt.ylabel('|Z|')
plt.grid('both')

# %% 
assert np.allclose(np.abs(Vfexc/Ifbulk), np.abs(Vfexc/Ifbulk_conf), rtol=0.01, atol=0.001)


# %%
