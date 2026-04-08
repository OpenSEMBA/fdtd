# %%
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))

from pyWrapper import *

case_dir = os.path.dirname(os.path.abspath(__file__))
SEMBA_MTLN_EXE    = os.path.normpath(os.path.join(case_dir, '../../../build-intel-rls/bin/semba-fdtd'))
SEMBA_NOMTLN_EXE  = os.path.normpath(os.path.join(case_dir, '../../../build-intel-rls-nomtln/bin/semba-fdtd'))
fn = os.path.join(case_dir, 'holland1981.fdtd.json')

# %% Run solvers
mtln_folder   = os.path.join(case_dir, 'mtln')
nomtln_folder = os.path.join(case_dir, 'nomtln')
os.makedirs(mtln_folder,   exist_ok=True)
os.makedirs(nomtln_folder, exist_ok=True)

solver_mtln = FDTD(input_filename=fn, path_to_exe=SEMBA_MTLN_EXE, run_in_folder=mtln_folder)
solver_mtln.cleanUp()
solver_mtln.run()
probe_files_mtln = [os.path.join(os.getcwd(), f)
                    for f in solver_mtln.getSolvedProbeFilenames("mid_point")]

solver_nomtln = FDTD(input_filename=fn, path_to_exe=SEMBA_NOMTLN_EXE, run_in_folder=nomtln_folder)
solver_nomtln.cleanUp()
solver_nomtln.run()
probe_files_nomtln = [os.path.join(os.getcwd(), f)
                      for f in solver_nomtln.getSolvedProbeFilenames("mid_point")]

# %% Plot results
plt.figure()
for pf in probe_files_mtln:
    p = Probe(pf)
    plt.plot(p['time'], p['current_0'], label=f'MTLN – {p.name}')
for pf in probe_files_nomtln:
    p = Probe(pf)
    plt.plot(p['time'], p['current'], linestyle='--', label=f'No MTLN – {p.name}')
plt.legend()
plt.tight_layout()
plt.show()
# %%
