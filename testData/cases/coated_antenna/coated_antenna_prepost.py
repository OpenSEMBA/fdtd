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

jcoated = json.load(open('./coated.json'))
jbare = json.load(open('./bare.json'))

coated_t, coated_i = np.array([]), np.array([])
for data in jcoated['datasetColl'][0]['data']:
    coated_t = np.append(coated_t, float(data['value'][0]))    
    coated_i = np.append(coated_i, float(data['value'][1]))    
bare_t, bare_i = np.array([]), np.array([])
for data in jbare['datasetColl'][0]['data']:
    bare_t = np.append(bare_t, float(data['value'][0]))    
    bare_i = np.append(bare_i, float(data['value'][1]))    


plt.figure()
plt.plot(bare_t*1e9, bare_i, 'k', label = r'$\varepsilon = 1.0$, paper')
plt.plot(coated_t*1e9, coated_i, 'b.', label = r'$\varepsilon = 3.2$, paper')
plt.plot(p['time']*1e9, -p['current_0'], '-',label = r'$\varepsilon = 3.2$, solved')
plt.plot(expected['time']*1e9, -expected['current_0'], '--',label = 'expected')

plt.grid(which='both')
plt.legend()
plt.xlabel('Time [ns]')
plt.ylabel('I [A]')
plt.xlim(0,12)
plt.show()


# %%
