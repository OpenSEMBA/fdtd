from src_pyWrapper.pyWrapper import *
import shutil
import glob
import re
import json
import numpy as np
import pyvista as pv
from numpy.fft import *
import pytest
from os import environ as env
from sys import platform
from scipy.special import hankel2 as h
from scipy.special import h2vp as hp

no_mtln_skip = pytest.mark.skipif(
    os.getenv("SEMBA_FDTD_ENABLE_MTLN") == "OFF",
    reason="MTLN is not available",
)

no_hdf_skip = pytest.mark.skipif(
    os.getenv("SEMBA_FDTD_ENABLE_HDF") == "OFF",
    reason="HDF5 is not available",
)

no_mpi_skip = pytest.mark.skipif(
    os.getenv("SEMBA_FDTD_ENABLE_MPI") == "OFF",
    reason="MPI is not available",
)

# Use of absolute path to avoid conflicts when changing directory.
if platform == "linux":
    SEMBA_EXE = os.path.join(os.getcwd(), 'build', 'bin', 'semba-fdtd')
elif platform == "win32":
    SEMBA_EXE = os.path.join(os.getcwd(), 'build', 'bin', 'semba-fdtd.exe')

NGSPICE_DLL = os.path.join(os.getcwd(), 'precompiled_libraries', 'windows-intel', 'ngspice', 'ngspice.dll')
TEST_DATA_FOLDER = os.path.join(os.getcwd(), 'testData/')
CASES_FOLDER = os.path.join(TEST_DATA_FOLDER, 'cases/')
MODELS_FOLDER = os.path.join(TEST_DATA_FOLDER, 'models/')
EXCITATIONS_FOLDER = os.path.join(TEST_DATA_FOLDER, 'excitations/')
OUTPUTS_FOLDER = os.path.join(TEST_DATA_FOLDER, 'outputs/')
SPINIT_FOLDER = os.path.join(TEST_DATA_FOLDER, 'spinit/')
GEOMETRIES_FOLDER = os.path.join(TEST_DATA_FOLDER, 'geometries/')
PROBES_INPUT_EXAMPLE = os.path.join(TEST_DATA_FOLDER, 'input_examples/probes/')

def getCase(case):
    return json.load(open(CASES_FOLDER + case + '.fdtd.json'))


def makeCopy(dest_dir, src_path):
    for file in glob.glob(src_path):
        src_file = file.split('/')[-1]
        dest_path = os.path.join(dest_dir, src_file)
        shutil.copy2(file, dest_path)


def copyInputFiles(temp_dir, input, excitation, executable):
    makeCopy(temp_dir, input)
    makeCopy(temp_dir, excitation)
    makeCopy(temp_dir, executable)


def getProbeFile(prefix, probe_name):
    extensions = ["dat", "xdmf", "h5", "bin"]

    probeFiles = []
    for ext in extensions:
        newFiles = ([x for x in glob.glob('*'+ext)
                    if re.match(prefix + '_' + probe_name + '.*'+ext,  x)])[0]
        probeFiles.append(newFiles)

    return probeFiles


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


def createPropertyDictionary(vtkfile, celltype:int, property:str):
    ugrid = pv.UnstructuredGrid(vtkfile)
    objs = np.argwhere(ugrid.celltypes == celltype)  # [i][0]
    if len(objs) == 0:
        return dict()
    
    props = np.array([])
    prop_array = ugrid.cell_data[property]
    for i in range(objs[0][0], objs[-1][0]+1):
        props = np.append(props, prop_array[i])

    unique, counts = np.unique(props, return_counts=True)
    return dict(zip(unique, counts))


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

def RCS(fspace, radius):

    def h_sph(n,z):
        return np.sqrt(0.5*np.pi*z)*h(n+0.5,z)

    def h_sph_der(n,z):
        return 0.25*np.pi*(0.5*np.pi*z)**(-0.5)*h(n+0.5,z) + np.sqrt(0.5*np.pi*z)*hp(n+0.5,z)

    a = radius #radius
    f = fspace
    l = 3.0e8/f
    beta = 2*np.pi*f/3.0e8

    res = 0.0
    for n in range(1,50):
        res = res + (-1)**n*(2*n+1)/(h_sph_der(n,beta*a)*h_sph(n,beta*a))
    res = np.abs(res)**2*(l**2/(4*np.pi))
    return res