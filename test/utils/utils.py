from src_pyWrapper.pyWrapper import *
import os
import shutil
import glob
import json
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
import h5py
import numpy as np
import pyvista as pv
from numpy.fft import *
import pytest
from os import environ as env
from sys import platform
from scipy.special import hankel2 as h
from scipy.special import h2vp as hp
from test.utils.build_resolver import build_feature_enabled, solver_executable


mtln_skip = pytest.mark.skipif(
    build_feature_enabled("SEMBA_FDTD_ENABLE_MTLN"),
    reason="Binary compiled with MTLN. No tests wire wires",
)

no_mtln_skip = pytest.mark.skipif(
    not build_feature_enabled("SEMBA_FDTD_ENABLE_MTLN"),
    reason="MTLN is not available",
)

no_mpi_skip = pytest.mark.skipif(
    not build_feature_enabled("SEMBA_FDTD_ENABLE_MPI"),
    reason="MPI is not available",
)

def _default_semba_exe(project_root=None):
    return str(solver_executable(project_root))


# Use an absolute path to avoid conflicts when changing directory.
SEMBA_EXE = _default_semba_exe()

NGSPICE_DLL = os.path.join(
    os.getcwd(), 'precompiled_libraries', 'windows-intel', 'ngspice', 'ngspice.dll')
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


def createWire(id, r, rpul=0.0, lpul=0.0):
    return {"id": id,
            "type": "wire",
            "radius": r,
            "resistancePerMeter": rpul,
            "inductancePerMeter": lpul
            }


def createUnshieldedWire(id, lpul, cpul, rpul=0.0, gpul=0.0):
    return {
        "id": id,
        "type": "unshieldedMultiwire",
        "inductancePerMeter":  [[lpul]],
        "capacitancePerMeter": [[cpul]],
        "resistancePerMeter": [rpul],
        "conductancePerMeter": [gpul]
    }


def createPropertyDictionary(vtkfile, celltype: int, property: str):
    ugrid = pv.UnstructuredGrid(vtkfile)
    cell_mask = ugrid.celltypes == celltype
    if not np.any(cell_mask):
        return dict()

    props = np.asarray(ugrid.cell_data[property])[cell_mask]

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

    def h_sph(n, z):
        return np.sqrt(0.5*np.pi*z)*h(n+0.5, z)

    def h_sph_der(n, z):
        return 0.25*np.pi*(0.5*np.pi*z)**(-0.5)*h(n+0.5, z) + np.sqrt(0.5*np.pi*z)*hp(n+0.5, z)

    a = radius  # radius
    f = fspace
    l = 3.0e8/f
    beta = 2*np.pi*f/3.0e8

    res = 0.0
    for n in range(1, 50):
        res = res + (-1)**n*(2*n+1)/(h_sph_der(n, beta*a)*h_sph(n, beta*a))
    res = np.abs(res)**2*(l**2/(4*np.pi))
    return res


def dtft(signal, time, freqs):
    """Continuous-time DTFT approximation using trapezoidal integration."""
    signal = np.asarray(signal)
    time = np.asarray(time)
    res = np.empty_like(freqs, dtype=complex)
    trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
    for i, f in enumerate(freqs):
        res[i] = trapz(signal * np.exp(-1j * 2 * np.pi * f * time), time)
    return res


def compute_impedance(probe_path, time_exc, voltage_exc, freqs):
    probe = Probe(probe_path)
    time_I = probe["time"].to_numpy()
    current = probe["current_0"].to_numpy()
    V_interp = np.interp(time_I, time_exc, voltage_exc)
    I_f = dtft(current, time_I, freqs)
    V_f = dtft(V_interp, time_I, freqs)
    Z = V_f / I_f
    return time_I, current, Z


def probe_from_fixture(tmp_path, filename):
    """Expose a legacy flat fixture through the folder-based Probe API."""
    source = os.path.join(OUTPUTS_FOLDER, filename)
    folder = os.path.join(tmp_path, "expected-probes", os.path.splitext(filename)[0])
    os.makedirs(folder, exist_ok=True)
    shutil.copy2(source, os.path.join(folder, filename))
    return Probe(folder)

def corrcoef_on_common_time(t_ref, y_ref, t_cmp, y_cmp):
    t_ini = max(t_ref[0], t_cmp[0])
    t_end = min(t_ref[-1], t_cmp[-1])
    mask = (t_ref >= t_ini) & (t_ref <= t_end)
    t_common = t_ref[mask]
    y_ref_common = y_ref[mask]
    y_cmp_interp = np.interp(t_common, t_cmp, y_cmp)
    return np.corrcoef(y_ref_common, y_cmp_interp)[0, 1]

def check_values_are_comparable(arr):
    if (np.any(np.isnan(arr)) == True or
        np.any(np.isinf(arr)) == True or
        np.all(arr == 0.0) == True):
        assert False
    else:
        assert True


def assert_allclose_to_peak(actual, expected, peak_tolerance=1e-5, err_msg=""):
    # atol scaled to the peak of either array absorbs float32 round-off on near-zero values (an all-zero reference would give atol=0)
    scale = max(np.max(np.abs(expected)), np.max(np.abs(actual)))
    np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=peak_tolerance * scale, err_msg=err_msg)


def h5_datasets(h5_path):
    datasets = {}
    with h5py.File(h5_path, "r") as h5_file:
        h5_file.visititems(
            lambda name, item: datasets.__setitem__(name, item[()])
            if isinstance(item, h5py.Dataset)
            else None
        )
    return datasets


def assert_same_h5_datasets(h5_path, reference_h5_path, peak_tolerance=1e-5):
    datasets = h5_datasets(h5_path)
    reference = h5_datasets(reference_h5_path)
    assert datasets.keys() == reference.keys()
    for name in reference:
        assert_allclose_to_peak(datasets[name], reference[name], peak_tolerance, err_msg=name)


def assert_static_point_attributes(xdmf_path, attribute_names):
    """Require node attributes backed by one non-series HDF5 dataset."""
    xdmf_path = Path(xdmf_path)
    root = ET.parse(xdmf_path).getroot()
    with h5py.File(xdmf_path.with_suffix(".h5"), "r") as hdf5:
        for name in attribute_names:
            dataset_paths = []
            for attribute in root.findall(f'.//Attribute[@Name="{name}"]'):
                data_item = attribute.find("DataItem")
                assert data_item is not None
                assert data_item.get("ItemType") != "HyperSlab", (
                    f"Node attribute {name} must be static"
                )
                assert data_item.text is not None
                dataset_paths.append(data_item.text.strip().split(":", maxsplit=1)[1])

            assert dataset_paths, f"Missing static node attribute {name}"
            assert len(set(dataset_paths)) == 1
            assert hdf5[dataset_paths[0]].ndim in (1, 3)


def solver_command(case, *flags, mpi_procs=None, exe=None):
    # argv to run the solver on an input file; mpi_procs wraps it in mpirun -np N
    command = [str(exe or SEMBA_EXE), "-i", str(case), *map(str, flags)]
    return ["mpirun", "-np", str(mpi_procs), *command] if mpi_procs else command


def run_solver(case, *flags, mpi_procs=None, cwd=None, env=None, exe=None):
    return subprocess.run(solver_command(case, *flags, mpi_procs=mpi_procs, exe=exe),
                          cwd=cwd or os.path.dirname(str(case)), env=env,
                          capture_output=True, text=True, check=False)


def read_probe_samples(solver, name):
    dat = Probe.getDatFile(Probe(solver.getSolvedProbeFolders(name)[0]))
    samples = np.loadtxt(dat, skiprows=1)
    return samples[:, 0], samples[:, 1]
