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

from pyWrapper import *

 
 
#####################################################
# %% Run solver
fn = 'unshielded_multiwire.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'])
solver.cleanUp()
solver.run()


# %% Plot results

j1 = json.load(open('berenger_multiwire_wire1.json'))
j2 = json.load(open('berenger_multiwire_wire2.json'))

jt1, ji1 = np.array([]), np.array([])
for data in j1['datasetColl'][0]['data']:
    jt1 = np.append(jt1, float(data['value'][0]))    
    ji1 = np.append(ji1, float(data['value'][1]))    
jt2, ji2 = np.array([]), np.array([])
for data in j2['datasetColl'][0]['data']:
    jt2 = np.append(jt2, float(data['value'][0]))    
    ji2 = np.append(ji2, float(data['value'][1]))    


probe_names = solver.getSolvedProbeFilenames("mid_point")
probe = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

plt.figure()
plt.plot(probe['time']*1e9-7, probe['current_1'], label='conductor 1')
plt.plot(probe['time']*1e9-7, probe['current_2'], label='conductor 2')
plt.grid(which='both')
plt.xlim(0,35)
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('I (A)')
plt.show()

# %% Plot output

