# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import pandas as pd
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build-rls/bin/semba-fdtd'

from pyWrapper import *

#####################################################
# %% Run solver
fn = 'cable_panel.fdtd.json'
# setNgspice('/')
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
# solver.cleanUp()
# solver.run()
#%%

df = pd.read_csv('tvs_simulation_coupled_2.txt', sep = ' ')
#%%
probeV = Probe(solver.getSolvedProbeFilenames("volt")[0])
probeI = Probe(solver.getSolvedProbeFilenames("curr")[0])
fig, ax = plt.subplots(1,2, figsize=(10,5))
fig.suptitle("Culprit wire - V/I at load port")
ax[0].plot(probeV['time'][0:290000]*1e6, probeV['voltage_0'][0:290000], '.',label = 'full-wave')
ax[0].plot(df['time']*1e6, df['Vp1'], label = 'LTSspice')
ax[1].plot(probeI['time'][0:290000]*1e6, -probeI['current_0'][0:290000]*1e3, '.')

ax[0].grid(which='both')
ax[0].legend()
ax[0].set_xlabel('Time (us)')
ax[0].set_ylabel('V (V)')

ax[1].grid(which='both')
ax[1].legend()
ax[1].set_xlabel('Time (us)')
ax[1].set_ylabel('I (mA)')
plt.show()

fig, ax = plt.subplots(1,2, figsize=(10,5))
fig.suptitle("Victim wire - V/I at load port")
ax[0].plot(probeV['time'][0:290000]*1e6, probeV['voltage_1'][0:290000], '.',label = 'full-wave')
ax[0].plot(df['time']*1e6, df['Vp6'], label = 'LTSpice')
ax[1].plot(probeI['time'][0:290000]*1e6, -probeI['current_1'][0:290000]*1e6, '.')
ax[1].plot(df['time']*1e6, df['IL2']*1e6)

ax[0].grid(which='both')
ax[0].legend()
ax[0].set_xlabel('Time (us)')
ax[0].set_ylabel('V (V)')

ax[1].grid(which='both')
ax[1].legend()
ax[1].set_xlabel('Time (us)')
ax[1].set_ylabel('I (uA)')
plt.show()


# %%
x = np.linspace(0,10,10)
plt.semilogy(x, 10**x)
plt.show()

# %%
