import subprocess
import json
import os
import shutil
import glob
import re
import pandas as pd
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from itertools import product
import copy
import matplotlib.pyplot as plt

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

class ExcitationFile():
    def __init__(self, excitation_filename):
        if isinstance(excitation_filename, os.PathLike):
            self.filename = excitation_filename.as_posix()
        else:
            self.filename = excitation_filename
        assert os.path.isfile(self.filename)

        self.data = pd.read_csv(self.filename, sep='\\s+', names=['time', 'value'])
    
    def plotExcitationFileResults(self, save_path=None) -> None:
        plt.figure(figsize=(8, 4))
        plt.plot(self.data['time'], self.data['value'], label='Excitation File Value')
        plt.xlabel('time')
        plt.ylabel('value')
        plt.title('Excitation File Results')
        plt.legend()
        plt.grid(True)
        plt.show()

        if save_path is None:
            base, _ = os.path.splitext(self.filename)
            save_path = base + '_plot.png'
        plt.savefig(save_path, bbox_inches='tight', dpi=150)
        
        plt.show()
        print(f"Plot saved to {save_path}")


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

        subfolders = [item for item in os.listdir(
            folder) if os.path.isdir(os.path.join(folder, item))]
        for f in subfolders:
            shutil.rmtree(f, ignore_errors=True)


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

    def getExcitationFile(self, excitation_file_name):
        file_extensions =('*.1.exc',)
        excitationFile = []
        for ext in file_extensions:
            newExcitationFile = [x for x in glob.glob(ext) if re.match(excitation_file_name, x)]
            excitationFile.extend(newExcitationFile)
        
        if ((len(excitationFile)) != 1):
            raise "Unexpected number of excitation Files found: {}".format(excitationFile)

        return excitationFile


    def getVTKMap(self):
        current_path = os.getcwd()
        folders = [item for item in os.listdir(
            current_path) if os.path.isdir(os.path.join(current_path, item))]
        if len(folders) == 0:
            return None
        for folder in folders:
            mapFile = os.path.join(current_path, folder, folder+"_1.vtk")
            if os.path.isfile(mapFile):
                return mapFile            
        raise ValueError("Unable to find mapvatk file")


    def getCurrentVTKMap(self):
        current_path = os.getcwd()
        folders = [item for item in os.listdir(
            current_path) if os.path.isdir(os.path.join(current_path, item))]
        if len(folders) != 1:
            return None
        for folder in folders:
            mapFile = os.path.join(current_path, folder, folder+"_1_current.vtk")
            if os.path.isfile(mapFile):
                return mapFile
        raise ValueError("Unable to find current vtk file")

    def getMaterialProperties(self, materialName):
        if 'materials' in self._input:
            for idx, element in enumerate(self._input['materials']):
                if element["name"] == materialName:
                    return self._input['materials'][idx]


class CaseMaker():
    # A collection of helper functions to create FDTD cases.
    def __init__(self):
        self.input = {}

    def __getitem__(self, key):
        return self.input[key]

    def setNumberOfTimeSteps(self, steps):
        if 'general' not in self.input:
            self.input['general'] = {}

        self.input['general']['numberOfSteps'] = steps

    def setGridFromVTK(self, path_to_grid):
        assert os.path.isfile(path_to_grid)

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(path_to_grid)
        reader.Update()

        points = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())

        if 'mesh' not in self.input:
            self.input['mesh'] = {}

        if 'grid' not in self.input['mesh']:
            self.input['mesh']['grid'] = {}

        grid = self.input['mesh']['grid']
        if 'steps' not in grid:
            grid['steps'] = {}

        steps = grid['steps']
        grid['numberOfCells'] = []
        grid['origin'] = []

        index = 0
        for x in ['x', 'y', 'z']:
            steps[x] = np.diff(np.unique(points[:, index])).tolist()
            grid['numberOfCells'].append(len(steps[x]))
            grid['origin'].append(float(points[0, index]))
            index += 1

    def setAllBoundaries(self, boundary_type):
        self.input['boundary'] = {
            "all": boundary_type
        }

    def addCellElementsFromVTK(self, path_to_vtk):
        assert os.path.isfile(path_to_vtk)

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(path_to_vtk)
        reader.Update()
        polyData = reader.GetOutput()

        vtkGroups = vtk_to_numpy(polyData.GetCellData().GetArray('group'))
        if len(np.unique(vtkGroups)) > 1:
            raise ValueError("Different groups are not supported.")

        # Stores in case input
        if 'mesh' not in self.input:
            self.input['mesh'] = {}

        if 'elements' not in self.input['mesh']:
            self.input['mesh']['elements'] = []
        elements = self.input['mesh']['elements']

        id = len(elements) + 1
        elements.append({
            "id": id,
            "type": "cell",
            "intervals": self._getCellsIntervals(polyData),
        })

        return id

    def addCellElementBox(self, boundingBox):
        ''' boundingBox is assumed to be in absolute coordinates.'''
        if 'mesh' not in self.input:
            self.input['mesh'] = {}

        if 'elements' not in self.input['mesh']:
            self.input['mesh']['elements'] = []
        elements = self.input['mesh']['elements']

        rel = self._absoluteToRelative(
            np.array([boundingBox[0], boundingBox[1]]))

        id = len(elements) + 1
        elements.append({
            "id": id,
            "type": "cell",
            "intervals": [
                [rel[0].tolist(), rel[1, :].tolist()]
            ]
        })

        return id

    def addNodeElement(self, position):
        ''' position is assumed to be in absolute coordinates.'''
        if isinstance(position, list):
            pos = np.array(position).reshape((1, 3))
        elif isinstance(position, np.ndarray):
            pos = position.reshape((1, 3))
        else:
            raise ValueError("Invalid type for position")
        relative = self._absoluteToRelative(pos)

        if 'mesh' not in self.input:
            self.input['mesh'] = {}

        if 'coordinates' not in self.input['mesh']:
            self.input['mesh']['coordinates'] = []

        coordId = len(self.input['mesh']['coordinates']) + 1
        self.input['mesh']['coordinates'].append({
            "id": coordId,
            "relativePosition": relative[0].tolist()
        })

        if 'elements' not in self.input['mesh']:
            self.input['mesh']['elements'] = []
        elementId = len(self.input['mesh']['elements']) + 1
        self.input['mesh']['elements'].append({
            "id": elementId,
            "type": "node",
            "coordinateIds": [coordId]
        })

        return elementId

    def addPECMaterial(self):
        if 'materials' not in self.input:
            self.input['materials'] = []

        id = len(self.input['materials']) + 1
        self.input['materials'].append({
            "id": id,
            "type": "pec",
        })

        return id

    def addMaterialAssociation(self, materialId: int, elementIds: list):
        if 'materials' not in self.input:
            raise ValueError("No materials defined.")
        if 'mesh' not in self.input:
            raise ValueError("No mesh defined.")
        if 'elements' not in self.input['mesh']:
            raise ValueError("No elements defined.")
        if len(elementIds) == 0:
            raise ValueError("No elements to associate.")

        if 'materialAssociations' not in self.input:
            self.input['materialAssociations'] = []

        self.input['materialAssociations'].append({
            'materialId': materialId,
            'elementIds': elementIds
        })

    def addPlanewaveSource(self, elementId, magnitudeFile, direction, polarization):
        assert os.path.isfile(magnitudeFile)

        if 'sources' not in self.input:
            self.input['sources'] = []

        self.input['sources'].append({
            'type': 'planewave',
            'magnitudeFile': magnitudeFile,
            'elementIds': [elementId],
            'direction': direction,
            'polarization': polarization
        })

    def addPointProbe(self, elementId, name):
        if 'probes' not in self.input:
            self.input['probes'] = []

        self.input['probes'].append({
            "name": name,
            "type": "point",
            "elementIds": [elementId]
        })

    def addFarFieldProbe(self, elementId, name, theta, phi, domain):
        if 'probes' not in self.input:
            self.input['probes'] = []

        self.input['probes'].append({
            "name": name,
            "type": "farField",
            "elementIds": [elementId],
            "theta": theta,
            "phi": phi,
            "domain": domain
        })

    def exportCase(self, case_name):
        with open(case_name + '.fdtd.json', 'w') as f:
            json.dump(self.input, f)

    def _buildGridLines(self):
        gridLines = []
        grid = self.input['mesh']['grid']
        index = 0
        for x in ['x', 'y', 'z']:
            newLines = np.cumsum(
                np.insert(grid['steps'][x], 0, grid['origin'][index]))
            index += 1
            gridLines.append(newLines)
        return gridLines

    def _relativeToAbsolute(self, relativeCoordinates):
        gridLines = self._buildGridLines()

        res = np.empty_like(relativeCoordinates, dtype=float)
        for x in range(3):
            res[:, x] = gridLines[x][relativeCoordinates[:, x]]

        return res

    def _absoluteToRelative(self, absoluteCoordinates):
        gridLines = self._buildGridLines()

        res = np.empty_like(absoluteCoordinates, dtype=int)
        for x in range(3):
            for i in range(absoluteCoordinates.shape[0]):
                array = gridLines[x]
                value = absoluteCoordinates[i, x]
                idx = np.searchsorted(array, value, side="left")
                if idx > 0 and (idx == len(array) or np.abs(value - array[idx-1]) < np.abs(value - array[idx])):
                    res[i, x] = idx-1
                else:
                    res[i, x] = idx

        if np.any(absoluteCoordinates - self._relativeToAbsolute(res) > 1e-10):
            raise ValueError(
                "Error in conversion from absolute to relative coordinates.")

        return res

    def _getCellsIntervals(self, polyData):
        relative = self._absoluteToRelative(
            vtk_to_numpy(polyData.GetPoints().GetData()))

        # Precounts
        numberOfIntervals = 0
        for i in range(polyData.GetNumberOfCells()):
            cellType = polyData.GetCellType(i)
            if cellType == vtk.VTK_QUAD or vtk.VTK_LINE:
                numberOfIntervals += 1
        res = [None] * numberOfIntervals

        # Writes as list
        c = 0
        for i in range(polyData.GetNumberOfCells()):
            cellType = polyData.GetCellType(i)
            if cellType == vtk.VTK_QUAD:
                p1 = polyData.GetCell(i).GetPointId(0)
                p2 = polyData.GetCell(i).GetPointId(2)
            elif cellType == vtk.VTK_LINE:
                p1 = polyData.GetCell(i).GetPointId(0)
                p2 = polyData.GetCell(i).GetPointId(1)
            else:
                continue
            res[c] = [relative[p1].tolist(), relative[p2].tolist()]
            c += 1

        return res
