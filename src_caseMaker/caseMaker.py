import json
import os
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np
from math import modf

class CaseMaker():
    # A collection of helper functions to create FDTD cases.
    def __init__(self, json_filename=""):
        self.input = {}
        if (json_filename != ""):
            self.input = json.load(open(json_filename))

    def __getitem__(self, key):
        return self.input[key]

    def setNumberOfTimeSteps(self, steps):
        if 'general' not in self.input:
            self.input['general'] = {}

        self.input['general']['numberOfSteps'] = steps

    def setTimeStep(self, time):
        if 'general' not in self.input:
            self.input['general'] = {}

        self.input['general']['timeStep'] = time

    def setGridFromVTK(self, path_to_grid):
        assert os.path.isfile(path_to_grid)

        with open(path_to_grid, 'r') as f:
            c = f.read(1)
            if (c == '<'):
                reader = vtk.vtkXMLPolyDataReader()
            else:
                reader = vtk.vtkUnstructuredGridReader()
                
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
            "all": {
                "type": boundary_type
            }
        }

    def addConformalVolumeFromVTK(self, path_to_vtk):
        assert os.path.isfile(path_to_vtk)

        with open(path_to_vtk, 'r') as f:
            c = f.read(1)
            if (c == '<'):
                reader = vtk.vtkXMLPolyDataReader()
            else:
                reader = vtk.vtkUnstructuredGridReader()

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
            "type": "conformalVolume",
            "intervals" : [],
            "triangles": self._getTriangles(polyData)
        })
        return id


    def addCellElementsFromVTK(self, path_to_vtk):
        assert os.path.isfile(path_to_vtk)

        with open(path_to_vtk, 'r') as f:
            c = f.read(1)
            if (c == '<'):
                reader = vtk.vtkXMLPolyDataReader()
            else:
                reader = vtk.vtkUnstructuredGridReader()

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
            "intervals": self._getCellsIntervals(polyData)
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
        for (idx, c) in enumerate(relativeCoordinates):
            for x in range(3):
                f, i = modf(c[x])[0], int(modf(c[x])[1])
                if (i == len(gridLines[x])-1) :
                    res[idx, x] = gridLines[x][i]
                else:
                    res[idx, x] = gridLines[x][i] + f*(gridLines[x][i+1]-gridLines[x][i])
                # res[:, x] = gridLines[x][relativeCoordinates[:, x]]

        return res

    def _filterIntegerCoordinates(self, absoluteCoordinates):
        nInts = 0
        for c in absoluteCoordinates:
            if all(v == True for v in [modf(c[x])[0] == 0 for x in range(3)]) == True:
                nInts = nInts + 1
        res = np.zeros([nInts,3])
        nInts = 0
        for c in absoluteCoordinates:
            if all(v == True for v in [modf(c[x])[0] == 0 for x in range(3)]) == True:
                res[nInts] = c
                nInts = nInts + 1
        return res
    
    def _absoluteToRelative(self, absoluteCoordinates, arg_type = int):
        gridLines = self._buildGridLines()

        if (arg_type == int):
            absoluteCoordinates = self._filterIntegerCoordinates(absoluteCoordinates)
        
        res = np.empty_like(absoluteCoordinates, dtype=arg_type)
        for x in range(3):
            for i in range(absoluteCoordinates.shape[0]):
                array = gridLines[x]
                value = absoluteCoordinates[i, x]
                idx = np.searchsorted(array, value, side="left")
                if (arg_type == int):
                    if idx > 0 and (idx == len(array) or np.abs(value - array[idx-1]) < np.abs(value - array[idx])):
                        res[i, x] = idx-1
                    else:
                        res[i,x] = idx
                else:
                    res[i, x] = idx-1 + (value-array[idx-1])/(array[idx]-array[idx-1])

        if np.any(absoluteCoordinates - self._relativeToAbsolute(res) > 1e-10):
            raise ValueError(
                "Error in conversion from absolute to relative coordinates.")

        return res

    def _getTriangles(self, polyData):
        relative = self._absoluteToRelative(
            vtk_to_numpy(polyData.GetPoints().GetData()), arg_type = float)

        numberOfTriangles = 0
        for i in range(polyData.GetNumberOfCells()):
            cellType = polyData.GetCellType(i)
            if cellType == vtk.VTK_TRIANGLE:
                numberOfTriangles += 1
        res = [None] * numberOfTriangles
        
        c = 0
        listOfCoords = []
        vtkToCoordId = {}
        if 'coordinates' not in self.input['mesh']:
            self.input['mesh']['coordinates'] = []

        for i in range(polyData.GetNumberOfCells()):
            cellType = polyData.GetCellType(i)
            if cellType == vtk.VTK_TRIANGLE:
                vtk_ids = [polyData.GetCell(i).GetPointId(0),
                          polyData.GetCell(i).GetPointId(1),
                          polyData.GetCell(i).GetPointId(2)]                 
                vertices = []
                for (id,p) in enumerate(vtk_ids):
                    if p not in listOfCoords:
                        listOfCoords.append(p)

                        coordId = len(self.input['mesh']['coordinates']) + 1
                        self.input['mesh']['coordinates'].append({
                            "id": coordId,
                            "relativePosition": relative[p].tolist()
                        })
                        
                        vtkToCoordId[p] = coordId
                    vertices.append(vtkToCoordId[p])
                res[c] = vertices
                    
                c = c+1

        return res            
      

    def _getCellsIntervals(self, polyData):
        relative = self._absoluteToRelative(
            vtk_to_numpy(polyData.GetPoints().GetData()))

        # Precounts
        numberOfIntervals = 0
        for i in range(polyData.GetNumberOfCells()):
            cellType = polyData.GetCellType(i)
            if cellType == (vtk.VTK_QUAD or vtk.VTK_LINE):
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
