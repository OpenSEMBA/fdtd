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
fn = 'holland1981_unshielded.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'])
solver.cleanUp()
solver.run()
probe_names = solver.getSolvedProbeFilenames("mid_point")
mid_mtln = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

fn = 'holland1981.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()
probe_names = solver.getSolvedProbeFilenames("mid_point")
mid = Probe(list(filter(lambda x: '_Wz_' in x, probe_names))[0])
# %% Run solver

solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'],mpi_command='mpirun -np 1')
solver.cleanUp()
solver.run()
probe_names = solver.getSolvedProbeFilenames("mid_point")
mid_mpi1 = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])
# %% Run solver

solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'],mpi_command='mpirun -np 2')
solver.cleanUp()
solver.run()
probe_names = solver.getSolvedProbeFilenames("mid_point")
mid_mpi2 = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])
# %% Run solver

solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'],mpi_command='mpirun -np 3')
solver.cleanUp()
solver.run()
probe_names = solver.getSolvedProbeFilenames("mid_point")
mid_mpi3 = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])


#####################################################
# %% Plot results

jpaper12 = json.load(open('../../outputs/holland1981.fdtd_mid_point_Wz_11_11_12_s2_12.json'))
tp12, ip12 = np.array([]), np.array([])
for data in jpaper12['datasetColl'][0]['data']:
    tp12 = np.append(tp12, float(data['value'][0]))    
    ip12 = np.append(ip12, float(data['value'][1]))    

p_expected = Probe('../../outputs/holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat')

# p_solved['current_0'][118]
# 9.65666368e-06
# -p_expected['current'][118]
# 8.57409577e-06
i_inter = np.interp(mid['time']-3.05*1e-9, tp12, ip12)
assert np.allclose(
    i_inter, 
    -mid['current'], 
    rtol=1e-4, atol=5e-5)
assert np.allclose(
    i_inter, 
    mid_mtln['current_0'], 
    rtol=1e-4, atol=5e-5)


plt.figure()
plt.plot(mid['time']*1e9, -mid['current'], label = 'Wires')
plt.plot(mid_mtln['time']*1e9, mid_mtln['current_0'], label = 'MTLN')
plt.plot(mid['time']*1e9, i_inter, 'r--', label = 'holland fig.12 interp')
plt.grid(which='both')
plt.legend()
plt.title('mid')
plt.xlabel('Time (ns)')
plt.xlim(0,20)
plt.ylabel('I (A)')
plt.show()


# %% Plot results

jpaper12 = json.load(open('../../outputs/holland1981.fdtd_mid_point_Wz_11_11_12_s2_12.json'))
tp12, ip12 = np.array([]), np.array([])
for data in jpaper12['datasetColl'][0]['data']:
    tp12 = np.append(tp12, float(data['value'][0]))    
    ip12 = np.append(ip12, float(data['value'][1]))    


plt.figure()
plt.plot(mid['time']*1e9 - 3.05, mid['current_0'], label = 'No MPI')
plt.plot(mid_mpi1['time']*1e9 - 3.05, mid_mpi1['current_0'], '*', label = 'MPI np=1')
plt.plot(mid_mpi2['time']*1e9 - 3.05, mid_mpi2['current_0'], '--', label = 'MPI np=2')
plt.plot(mid_mpi3['time']*1e9 - 3.05, mid_mpi3['current_0'], '-.', label = 'MPI np=3')
plt.plot(tp12*1e9, ip12, 'k--', label = 'holland fig.12')
plt.grid(which='both')
plt.legend()
plt.title('mid')
plt.xlabel('Time (ns)')
plt.xlim(0,20)
plt.ylabel('I (A)')
plt.show()


# %%
