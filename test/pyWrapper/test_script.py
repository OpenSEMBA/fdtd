from utils import *
from typing import Dict
import os
from sys import platform
from scipy import signal
import sys
from io import StringIO
import time
import matplotlib.pyplot as plt
import numpy as np

LOG_FILE = os.path.join(os.path.dirname(__file__), 'test_log.txt')

def printList(listObject: list, listName: str):
    print(f'List from {listName}')
    for element in listObject:
        print(str(element))
    print(f'End of list {listName}')

def read_xy_file(file_path):
    """
    Reads a file with two columns (x and y) and returns as numpy arrays.
    Skips the first row assuming it's a header.
    """
    data = np.loadtxt(file_path, delimiter=None, skiprows=1)
    if data.shape[1] != 2:
        raise ValueError(f"File {file_path} does not have two columns after skipping header.")
    x, y = data[:, 0], data[:, 1]
    return x, y

def plot_comparison(rls_folder, dbg_folder, output_folder):
    """
    Plot comparison between files in Release and Debug folders.
    Matches files by name.
    """
    rls_files = sorted(os.listdir(rls_folder))
    dbg_files = sorted(os.listdir(dbg_folder))
    skip_labels = ["Outputrequests", "Warnings", "Energy", "Report", "json", "paraviewfilters", ".exc"]
    for rls_file, dbg_file in zip(rls_files, dbg_files):
        if rls_file != dbg_file:
            print(f"Warning: file names do not match: {rls_file} vs {dbg_file}")
            continue

        if any(label in rls_file for label in skip_labels):
            print(f"Warning: file name to be skipped detected: {rls_file}")
            continue
        if any(label in dbg_file for label in skip_labels):
            print(f"Warning: file name to be skipped detected: {dbg_file}")
            continue        
        rls_path = os.path.join(rls_folder, rls_file)
        dbg_path = os.path.join(dbg_folder, dbg_file)
        
        x_rls, y_rls = read_xy_file(rls_path)
        x_dbg, y_dbg = read_xy_file(dbg_path)
        
        plt.figure(figsize=(8, 5))
        plt.plot(x_rls, y_rls, label='Release', linestyle='-', color='blue', alpha=0.7)
        plt.plot(x_dbg, y_dbg, label='Debug', linestyle='--', color='red', alpha=0.7)
        plt.title(f"Comparison: {rls_file}")
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

        # Save plot
        os.makedirs(output_folder, exist_ok=True)
        plt.savefig(os.path.join(output_folder, f"{rls_file}_comparison.png"))
        plt.close()

def test_compareOutput():
    # Carpeta raíz para la prueba
    base_test_folder = os.path.join(os.path.dirname(__file__), 'test')
    rlsFolder = os.path.join(base_test_folder, 'rls')
    dbgFolder = os.path.join(base_test_folder, 'dbg')
    
    # Crear carpetas si no existen
    os.makedirs(rlsFolder, exist_ok=True)
    os.makedirs(dbgFolder, exist_ok=True)
    
    # Archive old log if it exists
    if os.path.exists(LOG_FILE):
        old_log = os.path.join(os.path.dirname(LOG_FILE), 'test_log.old.txt')
        import shutil
        shutil.copy(LOG_FILE, old_log)
        print(f"Archived previous log to {old_log}")
    
    # Delete current log
    if os.path.exists(LOG_FILE):
        os.remove(LOG_FILE)
    
    # Redirect stdout to both console and log
    class Tee:
        def __init__(self, *files):
            self.files = files
        def write(self, obj):
            for f in self.files:
                f.write(obj)
                f.flush()
        def flush(self):
            for f in self.files:
                f.flush()
    
    with open(LOG_FILE, 'w') as log:
        sys.stdout = Tee(sys.stdout, log)
        
        fn = CASES_FOLDER + 'dielectric/dielectricTransmission.fdtd.json'
        
        # Release version
        solverrls = FDTD(
            input_filename=fn,
            path_to_exe=DEBUG_SEMBA_EXE,
            run_in_folder=rlsFolder
        )
        startTime = time.time()
        solverrls.run()
        endTime = time.time()
        print(f'Release version time spent {endTime-startTime}')
        
        # Debug version
        solverdbg = FDTD(
            input_filename=fn,
            path_to_exe=SEMBA_EXE,
            run_in_folder=dbgFolder
        )
        startTime = time.time()
        solverdbg.run()
        endTime = time.time()
        print(f'Debug version time spent {endTime-startTime}')
        
        rlsFiles = os.listdir(rlsFolder)
        printList(rlsFiles, 'Release')
        dbgFiles = os.listdir(dbgFolder)
        printList(dbgFiles, 'Debug')
        assert len(rlsFiles) == len(dbgFiles), 'Missing files'
        
        # Plot comparison
        plot_output_folder = os.path.join(base_test_folder, 'plots')
        plot_comparison(rlsFolder, dbgFolder, plot_output_folder)
        print(f"Plots saved in: {plot_output_folder}")
        
        sys.stdout = sys.__stdout__  # Restore stdout
