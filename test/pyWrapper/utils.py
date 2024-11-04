from src_pyWrapper.pyWrapper import *
import shutil, glob, re
import json
import numpy as np
import matplotlib.pyplot as plt

# Use of absolute path to avoid conflicts when changing directory.
SEMBA_EXE = os.getcwd() + '/build/bin/semba-fdtd'
SEMBA_EXE_INTEL_LLVM_RELEASE = os.getcwd() + '/build-ubuntu-intelLLVM-release/bin/semba-fdtd'
SEMBA_EXE_INTEL_LLVM_DEBUG   = os.getcwd() + '/build-ubuntu-intelLLVM-debug/bin/semba-fdtd'
TEST_DATA_FOLDER = os.getcwd() + '/testData/' 

CASE_FOLDER = TEST_DATA_FOLDER + 'cases/'
MODELS_FOLDER = TEST_DATA_FOLDER + 'models/'
EXCITATIONS_FOLDER = TEST_DATA_FOLDER + 'excitations/'
OUTPUT_FOLDER = TEST_DATA_FOLDER + 'outputs/'


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
    return  ([x for x in glob.glob('*dat') if re.match(prefix + '_' + probe_name + '.*dat',x)])[0]

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

