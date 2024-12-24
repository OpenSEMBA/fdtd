import subprocess
import json
import os
import shutil
import glob
import re
import pandas as pd
import numpy as np

DEFAULT_SEMBA_FDTD_PATH = '/build/bin/semba-fdtd'

class Probe():
    
    @staticmethod
    def _positionStrToCell(pos_str):
        pos = pos_str.split('_')
        return np.array([int(pos[0]), int(pos[1]), int(pos[2])])
    
    def __init__(self, probe_filename):
        self.filename = probe_filename

        mtln_probe_tags = ['_V_','_I_']
        current_probe_tags = ['_Wx_', '_Wy_', '_Wz_']
        point_probe_tags = ['_Ex_', '_Ey_', '_Ez_', '_Hx_', '_Hy_', '_Hz_']
        far_field_tag = ['_FF_']
        movie_tags = ['_ExC_', '_EyC_', '_EzC_', '_HxC_', '_HyC_', '_HzC_', '_ME_', '_MH_']
        
        all_tags = mtln_probe_tags \
            + current_probe_tags \
            + point_probe_tags \
            + far_field_tag \
            + movie_tags 
      
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
                    self.cell = self._positionStrToCell(position_str)
                    self.segment_tag = int(position_str.split('_s')[1])
                    self.df = pd.read_csv(self.filename, sep='\s+')
                    self.df = self.df.rename(columns={'t': 'time', self.df.columns[1]: 'current'})
                elif tag in point_probe_tags:
                    self.type = 'point'
                    self.name, position_str = basename_with_no_case_name.split(tag)
                    self.cell = self._positionStrToCell(position_str)
                    self.field = tag[1]
                    self.direction = tag[2]                       
                    self.df = pd.read_csv(self.filename, sep='\s+')
                    self.df = self.df.rename(columns = {
                        't': 'time', 
                        self.df.columns[1]: 'field',
                        self.df.columns[2]: 'incident'
                    })
                elif tag in far_field_tag:
                    self.type = 'farField'
                    self.name, positions_str = basename_with_no_case_name.split(tag)
                    init_str, end_str = positions_str.split('__')
                    self.cell_init = self._positionStrToCell(init_str)
                    self.cell_end  = self._positionStrToCell(end_str)
                    self.df = pd.read_csv(self.filename, sep='\s+')
                elif tag in movie_tags:
                    self.type = 'movie'
                    self.name, positions_str = basename_with_no_case_name.split(tag)
                    init_str, end_str = positions_str.split('__')
                    self.cell_init = self._positionStrToCell(init_str)
                    self.cell_end  = self._positionStrToCell(end_str)
                elif tag in mtln_probe_tags:
                    self.type ='mtln'
                    self.name, position_str = basename_with_no_case_name.split(tag)
                    self.cell = self._positionStrToCell(position_str)
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
    def __init__(self, input_filename, path_to_exe=None, flags = [], run_in_folder = None):
        self._setFilename(input_filename)
        
        if path_to_exe is None:
            self.path_to_exe = os.path.join(os.getcwd(), DEFAULT_SEMBA_FDTD_PATH)
        else:
            self.path_to_exe = path_to_exe
        assert os.path.isfile(self.path_to_exe)
        
        self.flags = flags
        self._hasRun = False
        
        if run_in_folder != None:
            self._setNewFolder(run_in_folder)
            
    
    def getFolder(self):
        res = os.path.dirname(self._filename)
        if len(res) == 0:
            return './'
        return res
        
    def getCaseName(self):
        return os.path.basename(self._filename).split('.json')[0]
    
    def _setFilename(self, newFilename):
        assert os.path.isfile(newFilename)
        self._filename = newFilename       
        self.input = json.load(open(self._filename))  
        
        
    def _getUsedFiles(self):
        res = []
        
        if 'sources' in self.input:
            for src in self.input['sources']:
                if 'magnitudeFile' in src:
                    res.append(src['magnitudeFile'])
                    
        if 'probes' in self.input:
            for src in self.input['probes']:
                if 'magnitudeFile' in src:
                    res.append(src['magnitudeFile'])
        return res
    
    def _setNewFolder(self, newFolder):
        assert os.path.isdir(newFolder)
        
        oldCaseFolder = self.getFolder()
        usedFiles = self._getUsedFiles()
        for usedFile in usedFiles:
            newFile = os.path.join(oldCaseFolder, usedFile)
            shutil.copy(newFile, newFolder)
        
        newFilename = shutil.copy(self._filename, newFolder)
        self._setFilename(newFilename)
    
    def run(self):
        if self.input != json.load(open(self._filename, 'r')):
            json.dump(self.input, open(self._filename,'w'))    
        
        os.chdir(self.getFolder())
        case_name = self.getCaseName() + ".json"
        self.output = subprocess.run([self.path_to_exe, "-i",case_name]+self.flags)
        self._hasRun = True
        
    def hasFinishedSuccessfully(self):
        if self._hasRun and (self.output.returncode == 0):
            return True
        else:
            return False
    
    def readJsonDict(self):
        with open(self._filename) as input_file:
            return json.load(input_file)
        
    def cleanUp(self):
        folder = self.getFolder()
        extensions = ('*.dat', '*.pl', '*.txt', '*.xdmf', '*.bin', '*.h5')
        for ext in extensions:
            files = glob.glob(folder + '/' + ext)
            for file in files:
                os.remove(file)
        
    def getSolvedProbeFilenames(self, probe_name):
        input_json = self.readJsonDict()
        if not "probes" in input_json:
            raise ValueError('Solver does not contain probes.')
        
        file_extensions = ('*.dat', '*.xdmf', '*.bin', '*.h5')
        probeFiles = []
        for ext in file_extensions:
            newProbes = [x for x in glob.glob(ext) if re.match(self.getCaseName() + '_' + probe_name, x)]
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
        

        
        
