from src_pyWrapper.pyWrapper import *
import shutil, glob, re
import json
import numpy as np
import matplotlib.pyplot as plt
<<<<<<< HEAD
import pyvista as pv

=======
from os import environ as env
from sys import platform
>>>>>>> origin/dev
# Use of absolute path to avoid conflicts when changing directory.
EXE_FOLDER = os.getcwd() + '/build/bin/'
SEMBA_EXE = os.getcwd() + '/build/bin/semba-fdtd'
SEMBA_EXE_INTEL_LLVM_RELEASE = os.getcwd() + '/build-ubuntu-intelLLVM-release/bin/semba-fdtd'
SEMBA_EXE_INTEL_LLVM_DEBUG   = os.getcwd() + '/build-ubuntu-intelLLVM-debug/bin/semba-fdtd'
TEST_DATA_FOLDER = os.getcwd() + '/testData/' 

CASE_FOLDER = TEST_DATA_FOLDER + 'cases/'
MODELS_FOLDER = TEST_DATA_FOLDER + 'models/'
EXCITATIONS_FOLDER = TEST_DATA_FOLDER + 'excitations/'
OUTPUT_FOLDER = TEST_DATA_FOLDER + 'outputs/'
SPINIT_FOLDER = TEST_DATA_FOLDER + 'spinit/'

def getCase(case):
    return json.load(open(CASE_FOLDER + case + '.fdtd.json'))

def makeCopy(temp_dir, file_path):
    file_name = file_path.split('/')[-1]
    temp_path = os.path.join(temp_dir, file_name)
    shutil.copy2(file_path, temp_path)



def copyInputFiles(temp_dir, input, excitation, executable):
    makeCopy(temp_dir, input)
    makeCopy(temp_dir, excitation)
    makeCopy(temp_dir, executable)

def getProbeFile(prefix, probe_name):
    extensions = ["dat", "xdmf", "h5", "bin"]
    
    probeFiles = []
    for ext in extensions:
        newFiles = ([x for x in glob.glob('*'+ext)  if re.match(prefix + '_' + probe_name + '.*'+ext,  x)])[0]
        probeFiles.append(newFiles)
    
    return  probeFiles

def countLinesInFile(probe_fn):
    with open(probe_fn, 'r') as f:
        return len(f.readlines())
    
def readWithoutHeader(file_name):
    with open(file_name) as f:
        _ = f.readline()
        rest = f.read()
    return rest
    
def compareFiles(expected_name, result_name):
    f_expected = readWithoutHeader(expected_name)
    f_result = readWithoutHeader(result_name)
    return f_expected == f_result            


def readSpiceFile(spice_file):
    t, val = np.array([]), np.array([])
    with open(spice_file) as f:
        for l in f:
            t = np.append(t, float(l.split()[0]))
            val = np.append(val, float(l.split()[1]))
    return t, val

def createFaceTagDictionary(vtkfile):
    ugrid = pv.UnstructuredGrid(vtkfile)
    quads = np.argwhere(ugrid.celltypes == 9) #[i][0]
    tags = np.array([])
    tagnumber_array = ugrid.cell_data['tagnumber']
    for i in range(quads[0][0], quads[-1][0]+1):
        tags = np.append(tags, tagnumber_array[i])

    unique, counts = np.unique(tags, return_counts=True)
    return dict(zip(unique, counts))

def setNgspice(tmp_path):
    if platform == "linux" or platform == "linux2":
        sys_name = "linux/"
        env["SPICE_SCRIPTS"] = "./"
    elif platform == "win32":
        sys_name = "windows/"    
        env["SPICE_SCRIPTS"] = "./"
    
    makeCopy(tmp_path, SPINIT_FOLDER + sys_name + 'spinit' )
    copyXSpiceModels(tmp_path, sys_name)
    #ngspice needs to read file 'spinit' to load code models needed by xspice
    # setSpiceScriptsFolder()

def copyXSpiceModels(temp_dir, sys_name):
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'analog.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'digital.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'spice2poly.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'table.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'xtradev.cm')
    makeCopy(temp_dir, SPINIT_FOLDER + sys_name + 'xtraevt.cm')
        