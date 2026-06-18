# %%
import numpy as np
from numpy.fft import *
import matplotlib.pyplot as plt
import shutil 

import sys, os
from sys import platform
from os import environ as env

from resource import getrusage as resource_usage, RUSAGE_SELF
from time import time as timestamp

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '../../../build-rls-nomtln/bin/semba-fdtd'
OUTPUTS_FOLDER = '../../outputs/'
SPINIT_FOLDER = '../../spinit/'
from pyWrapper import *

def makeCopy(dest_dir, src_path):
    for file in glob.glob(src_path):
        src_file = file.split('/')[-1]
        dest_path = os.path.join(dest_dir, src_file)
        shutil.copy2(file, dest_path)

def setNgspice(tmp_path):
    if platform == "linux" or platform == "linux2":
        sys_name = "linux/"
        env["SPICE_SCRIPTS"] = "./"
    elif platform == "win32":
        sys_name = "windows/"
        env["SPICE_SCRIPTS"] = "./"

    makeCopy(tmp_path, SPINIT_FOLDER + sys_name + 'spinit')
    copyXSpiceModels(tmp_path, sys_name)
    # ngspice needs to read file 'spinit' to load code models needed by xspice
    # setSpiceScriptsFolder()


def copyXSpiceModels(temp_dir, sys_name):
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'analog.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'digital.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'spice2poly.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'table.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'xtradev.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'xtraevt.cm')

cwd = os.getcwd()
setNgspice(cwd)

#####################################################
# %% Run solver

fn = 'towelHanger.fdtd.json'
solver = FDTD(input_filename = fn, path_to_exe=SEMBA_EXE)
solver.cleanUp()
solver.run()
#####################################################
# %% Plot results

p_solved = [Probe(solver.getSolvedProbeFilenames("wire_start")[0]),
                Probe(solver.getSolvedProbeFilenames("wire_mid")[0]),
                Probe(solver.getSolvedProbeFilenames("wire_end")[0])]

p_expected = [Probe(OUTPUTS_FOLDER+'towelHanger.fdtd_wire_start_Wz_27_25_30_s1.dat'),
                Probe(OUTPUTS_FOLDER+'towelHanger.fdtd_wire_mid_Wx_35_25_32_s5.dat'),
                Probe(OUTPUTS_FOLDER+'towelHanger.fdtd_wire_end_Wz_43_25_30_s4.dat')]


plt.figure()
plt.plot(p_solved[0]['time']*1e9, p_solved[0]['current_0'], label='solved')
plt.plot(p_expected[0]['time']*1e9, p_expected[0]['current_0'], '--',label='expected')
plt.grid(which='both')
plt.xlabel('Time [ns]')
plt.ylabel('I [A]')
# plt.xlim(0,0.5)
plt.legend()


# %%
