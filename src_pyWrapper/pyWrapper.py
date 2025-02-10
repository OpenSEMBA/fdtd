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
    MTLN_PROBE_TAGS = ['_V_', '_I_']
    CURRENT_PROBE_TAGS = ['_Wx_', '_Wy_', '_Wz_']
    BULK_CURRENT_PROBE_TAGS = ['_Jx_', '_Jy_', '_Jz_']
    POINT_PROBE_TAGS = ['_Ex_', '_Ey_', '_Ez_', '_Hx_', '_Hy_', '_Hz_']
    FAR_FIELD_TAG = ['_FF_']
    MOVIE_TAGS = ['_ExC_', '_EyC_', '_EzC_',
                  '_HxC_', '_HyC_', '_HzC_', '_ME_', '_MH_']

    ALL_TAGS = MTLN_PROBE_TAGS \
        + CURRENT_PROBE_TAGS \
        + BULK_CURRENT_PROBE_TAGS \
        + POINT_PROBE_TAGS \
        + FAR_FIELD_TAG \
        + MOVIE_TAGS

    def __init__(self, probe_filename):
        if isinstance(probe_filename, os.PathLike):
            self.filename = probe_filename.as_posix()
        else:
            self.filename = probe_filename
        assert os.path.isfile(self.filename)

        # This initialization tries to infer all probe properties from the filename.
        self.case_name = self._getCaseNameFromFilename(self.filename)
        self.name = self._getProbeNameFromFilename(self.filename)
        self.domainType = self._getDomainTypeFromFilename(self.filename)

        tag = self._getTagFromFilename(self.filename)
        if tag not in Probe.MOVIE_TAGS:
            self.data = pd.read_csv(self.filename, sep='\\s+')

        position_str = self._getPositionStrFromFilename(self.filename)
        if tag in Probe.CURRENT_PROBE_TAGS:
            self.type = 'wire'
            self.field, self.direction = Probe._getFieldAndDirection(tag)
            self.cell = self._positionStrToCell(position_str)
            self.segment = int(position_str.split('_s')[1])
            if self.domainType == 'time':
                self.data = self.data.rename(columns={
                    't': 'time',
                    self.data.columns[1]: 'current'
                })
            elif self.domainType == 'frequency':
                self.data = self.data.rename(columns={
                    self.data.columns[0]: 'frequency',
                    self.data.columns[1]: 'magnitude',
                    self.data.columns[2]: 'phase'
                })
        elif tag in Probe.BULK_CURRENT_PROBE_TAGS:
            self.type = 'bulkCurrent'
            self.field, self.direction = Probe._getFieldAndDirection(tag)
            self.cell = self._positionStrToCell(position_str)
            if self.domainType == 'time':
                self.data = self.data.rename(columns={
                    't': 'time',
                    self.data.columns[1]: 'current'
                })
            elif self.domainType == 'frequency':
                self.data = self.data.rename(columns={
                    self.data.columns[0]: 'frequency',
                    self.data.columns[1]: 'magnitude',
                    self.data.columns[2]: 'phase'
                })
        elif tag in Probe.POINT_PROBE_TAGS:
            self.type = 'point'
            self.field, self.direction = Probe._getFieldAndDirection(tag)
            self.cell = self._positionStrToCell(position_str)

            if self.domainType == 'time':
                self.data = self.data.rename(columns={
                    't': 'time',
                    self.data.columns[1]: 'field'
                })
                if len(self.data.columns) == 3:
                    self.data = self.data.rename(columns={
                        self.data.columns[2]: 'incident'
                    })
        elif tag in Probe.BULK_CURRENT_PROBE_TAGS:
            self.type = 'bulkCurrent'
            self.field, self.direction = Probe._getFieldAndDirection(tag)
            self.cell_init, self.cell_end = \
                Probe._positionStrToTwoCells(position_str)

            if self.domainType == 'time':
                self.data = self.data.rename(columns={
                    't': 'time',
                    self.data.columns[1]: 'current'
                })
                if len(self.data.columns) == 3:
                    self.data = self.data.rename(columns={
                        self.data.columns[2]: 'incident'
                    })
        elif tag in Probe.FAR_FIELD_TAG:
            self.type = 'farField'
            self.cell_init, self.cell_end = \
                Probe._positionStrToTwoCells(position_str)
        elif tag in Probe.MOVIE_TAGS:
            self.type = 'movie'
            self.cell_init, self.cell_end = \
                Probe._positionStrToTwoCells(position_str)
        elif tag in Probe.MTLN_PROBE_TAGS:
            self.type = 'mtln'
            self.cell = self._positionStrToCell(position_str)
            if self.domainType == 'time':
                self.data = self.data.rename(columns={'t': 'time'})
            for n in range(self.data.shape[1]-1):
                if tag == '_V_':
                    self.data = self.data.rename(columns={
                        self.data.columns[n+1]: 'voltage_'+str(n)
                    })
                elif tag == '_I_':
                    self.data = self.data.rename(columns={
                        self.data.columns[n+1]: 'current_'+str(n)
                    })

        else:
            raise ValueError("Unable to determine probe type")

    def __getitem__(self, key):
        return self.data[key]

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
    def _positionStrToCell(pos_str: str):
        pos = pos_str.split('_')
        return np.array([int(pos[0]), int(pos[1]), int(pos[2])])

    @staticmethod
    def _positionStrToTwoCells(pos_str: str):
        init_str, end_str = pos_str.split('__')
        return Probe._positionStrToCell(init_str), \
            Probe._positionStrToCell(end_str)

    @staticmethod
    def _getFieldAndDirection(tag: str):
        return tag[1], tag[2]


class FDTD():
    def __init__(self, input_filename, path_to_exe=None,
                 flags=None, run_in_folder=None, mpi_command=None):

        self._setFilename(input_filename)

        if path_to_exe is None:
            semba_exe = \
                os.path.join(os.getcwd(), DEFAULT_SEMBA_FDTD_PATH)
        else:
            semba_exe = path_to_exe
        assert os.path.isfile(semba_exe)

        if mpi_command is None:
            mpi_command_parts = []
        else:
            mpi_command_parts = mpi_command.split()

        if flags is None:
            flags = []
        elif isinstance(flags, str):
            flags = flags.split()

        case_name = self.getCaseName() + ".json"
        self.run_command = \
            mpi_command_parts + [semba_exe] + ["-i", case_name] + flags

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

    def __getitem__(self, key):
        return self._input[key]

    def _setFilename(self, newFilename):
        if isinstance(newFilename, os.PathLike):
            self._filename = str(newFilename)
        else:
            self._filename = newFilename
        self._filename = newFilename
        self._input = json.load(open(self._filename))

    def getUsedFiles(self):
        res = []

        # Files used to define magnitudes.
        if 'sources' in self._input:
            for p in self._input['sources']:
                if 'magnitudeFile' in p:
                    res.append(p['magnitudeFile'])

        # Files used in transfer functions.
        if 'probes' in self._input:
            for p in self._input['probes']:
                if 'magnitudeFile' in p:
                    res.append(p['magnitudeFile'])

        # .model files in circuit type materials.
        if 'materials' in self._input:
            for p in self._input['materials']:
                if 'file' in p:
                    res.append(p['file'])

        # .model files in terminations.
        if 'materials' in self._input:
            for p in self._input['materials']:
                if 'terminations' in p:
                    for t in p['terminations']:
                        if 'file' in t:
                            res.append(t['file'])

        return res

    def _setNewFolder(self, newFolder):
        assert os.path.isdir(newFolder)

        oldCaseFolder = self.getFolder()
        usedFiles = self.getUsedFiles()
        for usedFile in usedFiles:
            newFile = os.path.join(oldCaseFolder, usedFile)
            shutil.copy(newFile, newFolder)

        newFilename = shutil.copy(self._filename, newFolder)
        self._setFilename(newFilename)

    def run(self):
        if self._input != json.load(open(self._filename, 'r')):
            json.dump(self._input, open(self._filename, 'w'))

        os.chdir(self.getFolder())
        self.output = subprocess.run(self.run_command)

        self._hasRun = True
        assert self.hasFinishedSuccessfully()

    def hasFinishedSuccessfully(self):
        if self._hasRun and (self.output.returncode == 0):
            return True
        else:
            return False

    def __getitem__(self, key):
        return self._input[key]

    def cleanUp(self):
        folder = self.getFolder()
        extensions = ('*.dat', '*.pl', '*.txt', '*.xdmf', '*.bin', '*.h5')
        for ext in extensions:
            files = glob.glob(folder + '/' + ext)
            for file in files:
                os.remove(file)

    def getSolvedProbeFilenames(self, probe_name):
        if not "probes" in self._input:
            raise ValueError('Solver does not contain probes.')

        file_extensions = ('*.dat', '*.xdmf', '*.bin', '*.h5')
        probeFiles = []
        for ext in file_extensions:
            newProbes = [x for x in glob.glob(ext) if re.match(
                self.getCaseName() + '_' + probe_name, x)]
            probeFiles.extend(newProbes)

        return sorted(probeFiles)

    def getVTKMap(self):
        current_path = os.getcwd()
        folders = [item for item in os.listdir(
            current_path) if os.path.isdir(os.path.join(current_path, item))]
        if len(folders) != 1:
            return None
        mapFile = os.path.join(current_path, folders[0], folders[0]+"_1.vtk")
        assert os.path.isfile(mapFile)
        return mapFile

    def getMaterialProperties(self, materialName):
        if 'materials' in self._input:
            for idx, element in enumerate(self._input['materials']):
                if element["name"] == materialName:
                    return self._input['materials'][idx]
