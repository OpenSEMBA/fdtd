# %%
import numpy as np
import matplotlib.pyplot as plt
import json

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'
OUTPUTS_FOLDER = '../../outputs/'

from pyWrapper import *


#####################################################
# %% Run solver
fn = 'coated_antenna.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'])
solver.cleanUp()
solver.run()
probe_names = solver.getSolvedProbeFilenames("mid_point")
p = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])


#####################################################
# %% Plot results
pf = "coated_antenna.fdtd_mid_point_bundle_half_1_I_11_11_12.dat"
expected = Probe(OUTPUTS_FOLDER+pf)

plt.figure()
plt.plot(p['time'], -p['current_0'], '-',label = 'solved')
plt.plot(expected['time'], -expected['current_0'], '--',label = 'expected')

plt.grid(which='both')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('I (A)')
plt.show()


# %%
