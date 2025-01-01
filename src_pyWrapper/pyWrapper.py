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
    MTLN_PROBE_TAGS = ['_V_','_I_']
    CURRENT_PROBE_TAGS = ['_Wx_', '_Wy_', '_Wz_']
    POINT_PROBE_TAGS = ['_Ex_', '_Ey_', '_Ez_', '_Hx_', '_Hy_', '_Hz_']
    FAR_FIELD_TAG = ['_FF_']
    MOVIE_TAGS = ['_ExC_', '_EyC_', '_EzC_', '_HxC_', '_HyC_', '_HzC_', '_ME_', '_MH_']
    
    ALL_TAGS = MTLN_PROBE_TAGS \
        + CURRENT_PROBE_TAGS \
        + POINT_PROBE_TAGS \
        + FAR_FIELD_TAG \
        + MOVIE_TAGS 
       
    def __init__(self, probe_filename):
        assert os.path.isfile(probe_filename)
        
        # This initialization tries to infer all probe properties from the filename.
        self.filename = probe_filename
        self.case_name = self._getCaseNameFromFilename(probe_filename)
        self.name = self._getProbeNameFromFilename(probe_filename)
        self.domainType = self._getDomainTypeFromFilename(probe_filename)
        
        tag = self._getTagFromFilename(probe_filename)
        if tag not in Probe.MOVIE_TAGS:
            self.df = pd.read_csv(self.filename, sep='\\s+')
        
        position_str = self._getPositionStrFromFilename(probe_filename)
        if tag in Probe.CURRENT_PROBE_TAGS:
            self.type = 'wire'
            self.cell = self._positionStrToCell(position_str)
            self.segment = int(position_str.split('_s')[1])
            if self.domainType == 'time':
                self.df = self.df.rename(columns={\
                    't': 'time',
                    self.df.columns[1]: 'current'\
                })
            elif self.domainType == 'frequency':
                self.df = self.df.rename( columns={\
                    self.df.columns[0]: 'frequency', \
                    self.df.columns[1]: 'magnitude', \
                    self.df.columns[2]: 'phase'\
                })   
        elif tag in Probe.POINT_PROBE_TAGS:
            self.type = 'point'
            self.cell = self._positionStrToCell(position_str)
            self.field = tag[1]
            self.direction = tag[2]                       
            if self.domainType == 'time':
                self.df = self.df.rename(columns = {
                    't': 'time',
                    self.df.columns[1]: 'field',
                    self.df.columns[2]: 'incident'
                })
        elif tag in Probe.FAR_FIELD_TAG:
            self.type = 'farField'
            init_str, end_str = position_str.split('__')
            self.cell_init = self._positionStrToCell(init_str)
            self.cell_end  = self._positionStrToCell(end_str)
        elif tag in Probe.MOVIE_TAGS:
            self.type = 'movie'
            init_str, end_str = position_str.split('__')
            self.cell_init = self._positionStrToCell(init_str)
            self.cell_end  = self._positionStrToCell(end_str)
        elif tag in Probe.MTLN_PROBE_TAGS:
            self.type ='mtln'
            self.cell = self._positionStrToCell(position_str)
        else:
            raise ValueError("Unable to determine probe name or type for a probe with name:" + basename)
    
    def __getitem__(self, key):
        return self.df[key]
    
    @staticmethod
    def _getCaseNameFromFilename(fn):
        bn = os.path.basename(fn)
        if '.fdtd_' in bn:
           return bn.split('.fdtd_')[0]
        else:
            for tag in Probe.ALL_TAGS:
                if tag in bn:
                    return bn.split(tag)[0]
    
    @staticmethod
    def _getProbeNameFromFilename(fn):
        if '.fdtd_' in fn:
            tag = Probe._getTagFromFilename(fn)
            bn_without_case_name = fn.split('.fdtd_')[1]
            probe_name = bn_without_case_name.split(tag)[0]
            return probe_name
        else:
            return Probe._getCaseNameFromFilename(fn)
    
    @staticmethod
    def _getTagFromFilename(fn):
        for tag in Probe.ALL_TAGS:
            if tag in fn:
                return tag
        raise ValueError("Unable to determine probe tag")
    
    @staticmethod
    def _getPositionStrFromFilename(fn):
        bn_no_ext = os.path.splitext(os.path.basename(fn))[0]
        tag = Probe._getTagFromFilename(fn)
        if '_df' in bn_no_ext:
            bn_no_ext = bn_no_ext.replace('_df', '')
        position_str = bn_no_ext.split(tag)[1]
        return position_str
    
    @staticmethod
    def _getDomainTypeFromFilename(fn):
        if '_df.' in os.path.basename(fn):
            return 'frequency'    
        else:
            return 'time'
    
    @staticmethod
    def _positionStrToCell(pos_str):
        pos = pos_str.split('_')
        return np.array([int(pos[0]), int(pos[1]), int(pos[2])])
    

    
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
        
    def __getitem__(self, key):
        return self.input[key]
        
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
        

        
        
