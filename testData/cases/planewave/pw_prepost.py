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
solver.cleanUp()
solver.run() 

# %% Postprocess
before = Probe(solver.getSolvedProbeFilenames("before")[0])
inbox = Probe(solver.getSolvedProbeFilenames("inbox")[0])
after = Probe(solver.getSolvedProbeFilenames("after")[0])

plt.plot(inbox['time'], inbox['incident'], 'b.', label='incident') 
plt.plot(before['time'], before['field'], 'r-', label='before')
plt.plot(inbox['time'], inbox['field'], 'g-', label='inbox')
plt.plot(after['time'], after['field'], 'k-', label='after')
plt.legend()
plt.xlim(0,4e-9)



# %%
# %% Prepare solver
solver = FDTD(input_filename = 'pw-in-box-with-movie.fdtd.json', path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run() 

# %% Postprocess
import h5py
fn = solver.getSolvedProbeFilenames("electric_field_movie")[2]
with h5py.File(fn, "r") as f:
    time_key = list(f.keys())[0]
    field_key = list(f.keys())[1]
    time_ds = f[time_key][()] 
    field_ds = f[field_key][()]
# %%
