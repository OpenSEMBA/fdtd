# %%
import numpy as np
import matplotlib.pyplot as plt
import json

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build/bin/semba-fdtd'

from pyWrapper import *

 
#####################################################
# %% Run solver
fn = 'holland1981.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'])
solver.cleanUp()
solver.run()


#####################################################
# %% Plot results

jpaper12 = json.load(open('../../outputs/holland1981.fdtd_mid_point_Wz_11_11_12_s2_12.json'))
tp12, ip12 = np.array([]), np.array([])
for data in jpaper12['datasetColl'][0]['data']:
    tp12 = np.append(tp12, float(data['value'][0]))    
    ip12 = np.append(ip12, float(data['value'][1]))    

# probe_names = solver.getSolvedProbeFilenames("start_point")
# start = Probe(list(filter(lambda x: '_V_' in x, probe_names))[0])
# probe_names = solver.getSolvedProbeFilenames("end_point")
# end = Probe(list(filter(lambda x: '_V_' in x, probe_names))[0])

probe_names = solver.getSolvedProbeFilenames("mid_point")
mid = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

plt.figure()
plt.plot(mid['time']*1e9 - 3.05, mid['current_0'], label = 'mltnwires, mid')
# plt.plot(start['time']*1e9 - 3.05, start['voltage_0'], label = 'mltnwires, start')
# plt.plot(end['time']*1e9 - 3.05, end['voltage_0'], label = 'mltnwires, end')
plt.plot(tp12*1e9, ip12, 'k--', label = 'holland fig.12')
plt.grid(which='both')
plt.legend()
plt.title('mid')
plt.xlabel('Time (ns)')
plt.xlim(0,20)
# plt.ylim(-5e-3,5e-3)
plt.ylabel('I (A)')
plt.show()




#00 %%

# %%
