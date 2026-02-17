from utils import *
from typing import Dict
import os
from sys import platform
from scipy import signal
import sys
from io import StringIO
import time

LOG_FILE = os.path.join(os.path.dirname(__file__), 'test_log.txt')

def printList(listObject:list, listName:str):
        print(f'List from {listName}')
        for element in listObject:
                print(str(element))
        print(f'End of list {listName}')

def test_compareOutput(tmp_path):
        # Delete log file if it exists
        if os.path.exists(LOG_FILE):
                os.remove(LOG_FILE)
        
        # Redirect stdout to both console and log file
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
                print(os.path.join(tmp_path, 'rls'))
                fn = CASES_FOLDER + 'coated_antenna/coated_antenna.fdtd.json'
                rlsFolder = os.path.join(tmp_path, 'rls')
                solverrls = FDTD(
                        input_filename=fn,
                        path_to_exe=DEBUG_SEMBA_EXE,
                        run_in_folder=rlsFolder
                )
                startTime = time.time()
                solverrls.run()
                endTime = time.time()
                print(f'Release version time spent {endTime-startTime}')

                fn = CASES_FOLDER + 'coated_antenna/coated_antenna.fdtd.json'
                dbgFolder = os.path.join(tmp_path, 'dbg')
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


                sys.stdout = sys.__stdout__  # Restore stdout