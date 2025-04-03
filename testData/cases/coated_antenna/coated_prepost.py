# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import shutil 

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

 

#####################################################
# %% Run solver
fn = 'coated_antenna.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'])
# solver.cleanUp()
solver.run()

#####################################################
# %% Plot results

probe_names = solver.getSolvedProbeFilenames("mid_point")
mid = Probe(list(filter(lambda x: '_Wz_' in x, probe_names))[0])
probe_names = solver.getSolvedProbeFilenames("out_mid_point")
out = Probe(list(filter(lambda x: '_Wz_' in x, probe_names))[0])

jcoated = json.load(open('../../../tmp_cases/tmp_coated/coated.json'))
jbare = json.load(open('../../../tmp_cases/tmp_coated/bare.json'))

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
plt.plot(mid['time']*1e9, -mid['current'], label='solved')
plt.plot(out['time']*1e9, -out['current'], '--',label='expected')
plt.grid(which='both')
plt.xlim(0,12)
# plt.ylim(-1e-3,1e-3)
plt.legend()
plt.title('mid')
plt.xlabel('Time (ns)')
plt.ylabel('I (A)')
plt.show()
plt.figure()


# %%
