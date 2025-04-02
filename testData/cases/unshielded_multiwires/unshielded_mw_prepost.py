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
# %% Generate excitation and visualize
dt = 3e-11
tf = 50e-9
t = np.arange(0, tf, dt)
e = 1e5*(1-np.exp(-t*1e9/5))

plt.figure()
plt.plot(t*1e9,e)
plt.xlabel('Time (ns)')

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = e
np.savetxt('unshielded_50ns.exc', data)

# %%

p11 = 1/12.39e-12
p12 = 1/43.85e-12
p21 = 1/36.06e-12
p22 = 1/28.67e-12
mu0 = 4*np.pi*1e-7
eps0 = 8.85e-12
C = np.linalg.inv(np.array([[p11,p12],[p21, p22]]))
L = mu0*eps0*np.linalg.inv(C)
print(C)
print(L)
 
#####################################################
# %% Run solver
fn = 'unshielded_multiwire.fdtd.json'
shutil.copy('../../excitations/unshielded.exc','.')
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'])

solver["materialAssociations"][0]["materialId"] = 12
solver["materialAssociations"][0]["initialTerminalId"] = 4
solver["materialAssociations"][0]["endTerminalId"] = 4
solver["materialAssociations"][0]["elementIds"] = [5,6]
solver.cleanUp()
solver.run()

#####################################################
# %% Plot results

j1 = json.load(open('../../outputs/berenger_multiwire_wire1.json'))
j2 = json.load(open('../../outputs/berenger_multiwire_wire2.json'))

jt1, ji1 = np.array([]), np.array([])
for data in j1['datasetColl'][0]['data']:
    jt1 = np.append(jt1, float(data['value'][0]))    
    ji1 = np.append(ji1, float(data['value'][1]))    
jt2, ji2 = np.array([]), np.array([])
for data in j2['datasetColl'][0]['data']:
    jt2 = np.append(jt2, float(data['value'][0]))    
    ji2 = np.append(ji2, float(data['value'][1]))    

probe_names = solver.getSolvedProbeFilenames("mid_point")
mid_I = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])
probe_names = solver.getSolvedProbeFilenames("mid_point")
mid_W = Probe(list(filter(lambda x: '_Wz_' in x, probe_names))[0])

# probe_names = solver.getSolvedProbeFilenames("start_point")
# start = Probe(list(filter(lambda x: '_V_' in x, probe_names))[0])
# probe_names = solver.getSolvedProbeFilenames("end_point")
# end = Probe(list(filter(lambda x: '_V_' in x, probe_names))[0])


plt.figure()
plt.plot(mid_W['time']*1e9, -mid_W['current'], label='material 12 - wire 0 W')
plt.plot(mid_I['time']*1e9, mid_I['current_0'], label='material 12 - wire 0')
plt.plot(mid_I['time']*1e9, mid_I['current_1'], label='material 12 - wire 1')
plt.plot(mid_I['time']*1e9, mid_I['current_2'], label='material 12 - wire 2')
plt.plot(jt1*1e9+1, ji1, 'g.', label = 'wire1 - integral eq')
plt.plot(jt2*1e9+1, ji2, 'k.', label = 'wire2 - integral eq')
plt.grid(which='both')
plt.xlim(0,35)
plt.ylim(-160,300)
plt.legend()
plt.title('mid')
plt.xlabel('Time (ns)')
plt.ylabel('I (A)')
plt.show()
plt.figure()


# %%
