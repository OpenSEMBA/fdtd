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
fn = 'shieldedPair.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'])
solver.cleanUp()
solver.run()
probe_names = solver.getSolvedProbeFilenames("wire_start")
mid_nompi = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

fn = 'shieldedPair.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'],mpi_command='mpirun -np 2')
solver.cleanUp()
solver.run()
probe_names = solver.getSolvedProbeFilenames("wire_start")
mid_mpi2 = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

# %% Run solver
fn = 'shieldedPair.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE, flags=['-mtlnwires'],mpi_command='mpirun -np 3')
solver.cleanUp()
solver.run()
probe_names = solver.getSolvedProbeFilenames("wire_start")
mid_mpi3 = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])


#####################################################
# %% Plot results
pf = "shieldedPair.fdtd_wire_start_Wz_75_74_74_s4.dat"
expected = Probe(OUTPUTS_FOLDER+pf)
plt.figure()
plt.plot(mid_nompi['time'], -mid_nompi['current_0'], '.',label = 'no mpi')
plt.plot(mid_mpi2['time'], -mid_mpi2['current_0'], '-.',label = 'mpi 2')
plt.plot(expected['time'], -expected['current'], '--',label = 'expected')

plt.grid(which='both')
plt.legend()
plt.xlabel('Time (ns)')
plt.ylabel('I (A)')
plt.show()


# %%
