import subprocess
import json
import os
import glob, re
import pandas as pd
import numpy as np

def positionStrToCell(pos_str):
    pos = pos_str.split('_')
    return np.array([int(pos[0]), int(pos[1]), int(pos[2])])

class Probe():
    
    def __init__(self, probe_filename):
        self.filename = probe_filename

        mtln_probe_tags = ['_V_','_I_']
        current_probe_tags = ['_Wx_', '_Wy_', '_Wz_']
        far_field_tag = ['_FF_']
        movie_tags = ['_ExC_', '_EyC_', '_EzC_', '_HxC_', '_HyC_', '_HzC_', '_ME_', '_MH_']
        all_tags = current_probe_tags + far_field_tag + movie_tags + mtln_probe_tags
      
        
        basename = os.path.basename(self.filename)
        self.case_name, basename_with_no_case_name = basename.split('.fdtd_')
        basename_with_no_case_name = os.path.splitext(basename_with_no_case_name)[0]
        
        
        for tag in all_tags:
            ids = [m.start() for m in re.finditer(tag, basename_with_no_case_name)]
            if len(ids) == 0:
                continue
            elif len(ids) == 1:
                if tag in current_probe_tags:
                    self.type = 'wire'
                    self.name, position_str = basename_with_no_case_name.split(tag)
                    self.cell = positionStrToCell(position_str)
                    self.segment_tag = int(position_str.split('_s')[1])
                    self.df = pd.read_csv(self.filename, sep='\s+')
                    self.df = self.df.rename(columns={'t': 'time', self.df.columns[1]: 'current'})
                    # self.df = self.df.rename(columns={'t': 'time', basename: 'current'})
                elif tag in far_field_tag:
                    self.type = 'farField'
                    self.name, positions_str = basename_with_no_case_name.split(tag)
                    init_str, end_str = positions_str.split('__')
                    self.cell_init = positionStrToCell(init_str)
                    self.cell_end = positionStrToCell(end_str)
                    self.df = pd.read_csv(self.filename, sep='\s+')
                elif tag in movie_tags:
                    self.type = 'movie'
                    self.name, positions_str = basename_with_no_case_name.split(tag)
                    init_str, end_str = positions_str.split('__')
                    self.cell_init = positionStrToCell(init_str)
                    self.cell_end = positionStrToCell(end_str)
                elif tag in mtln_probe_tags:
                    self.type ='mtln'
                    self.name, position_str = basename_with_no_case_name.split(tag)
                    self.cell = positionStrToCell(position_str)
                    self.df = pd.read_csv(self.filename, sep='\s+')
            else:
                raise ValueError("Unable to determine probe name or type for a probe with name:" + basename)
        try:
            self.type
            self.name
        except:
            raise ValueError('Unable to determine type for probe' + basename)
           
    
    def __getitem__(self, key):
        return self.df[key]

class FDTD():
    def __init__(self, input_filename, path_to_exe, flags = []):
        self.filename = input_filename
        self.path_to_exe = path_to_exe

        self.flags = flags

        self.folder = os.path.dirname(self.filename)
        self.case = os.path.basename(self.filename).split('.json')[0]
        self.hasRun = False
    
    def run(self):
        os.chdir(self.folder)
        case_name = self.case + ".json"
        self.output = subprocess.run([self.path_to_exe, "-i",case_name]+self.flags)
        self.hasRun = True
    
    def readJsonDict(self):
        with open(self.filename) as input_file:
            return json.load(input_file)
        
    def getSolvedProbeFilenames(self, probe_name):
        input_json = self.readJsonDict()
        if not "probes" in input_json:
            raise ValueError('Solver does not contain probes.')
        
        file_extensions = ('*.dat', '*.xdmf', '*.bin', '*.h5')
        probeFiles = []
        for ext in file_extensions:
            newProbes = [x for x in glob.glob(ext) if re.match(self.case + '_' + probe_name, x)]
            probeFiles.extend(newProbes)
            
        return probeFiles
    
    def getVTKMap(self):
        current_path = os.getcwd()
        folders = [item for item in os.listdir(current_path) if os.path.isdir(os.path.join(current_path, item))]
        if len(folders) != 1:
            return None
        mapFile = os.path.join(current_path,folders[0],folders[0]+"_1.vtk")
        assert os.path.isfile(mapFile)
        return mapFile
        
    def hasFinishedSuccessfully(self):
        if self.hasRun:
            if (self.output.returncode == 0):
                return True
            else:
                return False
        
        
