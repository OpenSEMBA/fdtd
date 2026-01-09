# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import pandas as pd

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %%
solver = FDTD(input_filename = 'conformal_thin_strip.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()

# %%
exc = pd.read_csv("predefinedExcitation.1.exc", sep='\\s+')
exc = exc.rename(columns={
    exc.columns[0]: 'time',
    exc.columns[1]: 'V'
})
new_freqs = np.geomspace(1e3, 1e7, num=100)
Vexc = exc["V"].to_numpy()
texc = exc["time"].to_numpy()
dt_exc = texc[1]-texc[0]
Vfexc = dt_exc*np.array([np.sum(Vexc * np.exp(-1j * 2 * np.pi * f * texc)) for f in new_freqs])


bulk = Probe(solver.getSolvedProbeFilenames("BulkProbe")[0])
Ibulk = bulk["current"].to_numpy()
tbulk = bulk["time"].to_numpy()
dt_bulk = tbulk[1]-tbulk[0]
Ifbulk = dt_bulk*np.array([np.sum(Ibulk * np.exp(-1j * 2 * np.pi * f * tbulk)) for f in new_freqs])

#impedance comparison 
plt.loglog(new_freqs, np.abs(Vfexc/Ifbulk), label='conformal')
plt.ylabel(r'Z($\Omega$)')
plt.xlabel('f [Hz]')
