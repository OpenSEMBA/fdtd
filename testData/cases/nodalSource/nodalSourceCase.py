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

import pyWrapper as pW

# %%
fn = "NodalSourceTest.fdtd.json"
solver = pW.FDTD(input_filename=fn, path_to_exe=SEMBA_EXE)
solver.run()

# %%
resistanceBulkProbe = pW.Probe(solver.getSolvedProbeFilenames("Bulk probe Resistance")[0])
nodalBulkProbe = pW.Probe(solver.getSolvedProbeFilenames("Bulk probe Nodal Source")[0])
excitation = pW.ExcitationFile(excitation_filename=solver.getExcitationFile("predefinedExcitation")[0])

plt.figure()
plt.plot(resistanceBulkProbe['time'].to_numpy(), 
         resistanceBulkProbe['current'].to_numpy(), label='BP Current@resistance')
plt.plot(excitation.data['time'].to_numpy(), 
         excitation.data['value'].to_numpy(), label='excited current')
plt.plot(nodalBulkProbe['time'].to_numpy(),
         nodalBulkProbe['current'].to_numpy(), label='BP Current@nodal source')
plt.legend()


# %%
