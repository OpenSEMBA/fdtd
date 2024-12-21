# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

# %% Prepare solver
solver = FDTD(input_filename = 'pw-in-box.fdtd.json', path_to_exe=SEMBA_EXE)

# %% Run solver
solver.run()
assert solver.hasFinishedSuccessfully()

# %% Postprocess

before = Probe(solver.getSolvedProbeFilenames("before")[0])
inbox = Probe(solver.getSolvedProbeFilenames("inbox")[0])
after = Probe(solver.getSolvedProbeFilenames("after")[0])

plt.plot(before.df['time'], before.df['incident'], 'b.-', label='incident') 
plt.plot(before.df['time'], before.df['field'], 'r-', label='before')
plt.plot(inbox.df['time'], inbox.df['field'], 'g-', label='inbox')
plt.plot(after.df['time'], after.df['field'], 'k-', label='after')
plt.legend()


# %%