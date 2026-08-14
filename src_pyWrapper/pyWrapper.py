import subprocess
import json
import os
import shutil
import glob
import re
from enum import Enum
import pandas as pd
import numpy as np

DEFAULT_SEMBA_FDTD_PATH = "build/bin/semba-fdtd"


class _ProbeType(str, Enum):
    WIRE = "wire"
    BULK_CURRENT = "bulkCurrent"
    LINE_INTEGRAL = "lineIntegral"
    POINT = "point"
    FAR_FIELD = "farField"
    MOVIE = "movie"
    MTLN = "mtln"


class _ProbeDomain(str, Enum):
    TIME = "time"
    FREQUENCY = "frequency"


_EXPECTED_DAT_COLUMNS = {
    (_ProbeType.WIRE, _ProbeDomain.TIME): (
        "t",
        "current",
        "delta_voltage",
        "plus_voltage",
        "minus_voltage",
        "voltage_difference",
    ),
    (_ProbeType.WIRE, _ProbeDomain.FREQUENCY): ("frequency", "magnitude", "phase"),
    (_ProbeType.BULK_CURRENT, _ProbeDomain.TIME): ("t", "current"),
    (_ProbeType.LINE_INTEGRAL, _ProbeDomain.TIME): ("t", "lineIntegral"),
    (_ProbeType.POINT, _ProbeDomain.TIME): ("t", "field"),
    (_ProbeType.POINT, _ProbeDomain.FREQUENCY): ("frequency", "real", "imaginary"),
    (_ProbeType.FAR_FIELD, _ProbeDomain.FREQUENCY): (
        "frequency",
        "Theta",
        "Phi",
        "Etheta_mod",
        "Etheta_phase",
        "Ephi_mod",
        "Ephi_phase",
        "RCS_arithmetic",
        "RCS_geometric",
    ),
    (_ProbeType.MTLN, _ProbeDomain.TIME): ("t", "{quantity}_0"),
    (_ProbeType.MTLN, _ProbeDomain.FREQUENCY): ("frequency", "{quantity}_0"),
}


class Probe:
    folder: str
    case_name: str
    name: str
    domainType: str
    type: str
    data: pd.DataFrame | None
    field: str
    direction: str
    cell: np.ndarray | None
    segment: int
    cell_init: np.ndarray
    cell_end: np.ndarray

    MTLN_PROBE_TAGS: list[str] = ["_V_", "_I_"]
    CURRENT_PROBE_TAGS: list[str] = ["_Wx_", "_Wy_", "_Wz_"]
    BULK_CURRENT_PROBE_TAGS: list[str] = ["_Jx_", "_Jy_", "_Jz_"]
    LINE_INTEGRAL_PROBE_TAG: list[str] = ["_LI_"]
    POINT_PROBE_TAGS: list[str] = ["_Ex_", "_Ey_", "_Ez_", "_Hx_", "_Hy_", "_Hz_"]
    FAR_FIELD_TAG: list[str] = ["_FF_"]
    MOVIE_TAGS: list[str] = ["_ExC_", "_EyC_", "_EzC_", "_HxC_", "_HyC_", "_HzC_", "_ME_", "_MH_"]

    ALL_TAGS: list[str] = (
        MTLN_PROBE_TAGS
        + CURRENT_PROBE_TAGS
        + BULK_CURRENT_PROBE_TAGS
        + LINE_INTEGRAL_PROBE_TAG
        + POINT_PROBE_TAGS
        + FAR_FIELD_TAG
        + MOVIE_TAGS
    )

    FILE_EXTENSIONS: list[str] = [".dat", ".xdmf", ".h5", ".bin"]
    DOMAIN_MARKERS: list[str] = ["_tm", "_fq", "_df"]

    def __init__(self, probe_folder):
        self.folder = os.path.abspath(os.fspath(probe_folder))
        assert os.path.isdir(self.folder), "Probe requires a probe output folder"

        self.case_name = self._getCaseNameFromFolder(self.folder)
        self.name = self._getProbeNameFromFolder(self.folder)
        self.domainType = self._getDomainTypeFromFolder(self.folder)
        self._domain = _ProbeDomain(self.domainType)
        self.data = None

        tag = self._getTagFromFolder(self.folder)

        position_str = self._getPositionStrFromFolder(self.folder)
        if tag in Probe.CURRENT_PROBE_TAGS:
            self.type = "wire"
            self.field, self.direction = Probe._getFieldAndDirection(tag)
            self.cell = self._positionStrToCell(position_str)
            self.segment = int(position_str.split("_s")[1])
        elif tag in Probe.BULK_CURRENT_PROBE_TAGS:
            self.type = "bulkCurrent"
            self.field, self.direction = Probe._getFieldAndDirection(tag)
            self.cell_init, self.cell_end = Probe._positionStrToTwoCells(position_str)
        elif tag in Probe.LINE_INTEGRAL_PROBE_TAG:
            self.type = "lineIntegral"
            self.field, self.direction = Probe._getFieldAndDirection(tag)
            self.cell = self._positionStrToCell(position_str) if position_str else None
        elif tag in Probe.POINT_PROBE_TAGS:
            self.type = "point"
            self.field, self.direction = Probe._getFieldAndDirection(tag)
            self.cell = self._positionStrToCell(position_str)
        elif tag in Probe.FAR_FIELD_TAG:
            self.type = "farField"
            self.cell_init, self.cell_end = Probe._positionStrToTwoCells(position_str)
        elif tag in Probe.MOVIE_TAGS:
            self.type = "movie"
            self.cell_init, self.cell_end = Probe._positionStrToTwoCells(position_str)
        elif tag in Probe.MTLN_PROBE_TAGS:
            self.type = "mtln"
            self.cell = self._positionStrToCell(position_str)
        else:
            raise ValueError("Unable to determine probe type")

        self._type = _ProbeType(self.type)
        text_file = self.getTextFile()
        if text_file is not None:
            self.data = pd.read_csv(text_file, sep="\\s+")
            self._normalise_data_columns(tag)

    def getExpectedOutputs(self):
        extensions = {
            _ProbeType.WIRE: [".dat", ".bin"],
            _ProbeType.BULK_CURRENT: [".dat", ".bin"],
            _ProbeType.LINE_INTEGRAL: [".dat", ".bin"],
            _ProbeType.POINT: [".dat", ".bin"],
            _ProbeType.FAR_FIELD: ["", ".bin"],
            _ProbeType.MOVIE: [".bin", ".xdmf", ".h5", "_geometry.xdmf", "_geometry.h5"],
            _ProbeType.MTLN: [".dat"],
        }[self._type]

        return sorted(
            os.path.join(self.folder, self._getFolderStem() + extension)
            for extension in extensions
            if os.path.isfile(os.path.join(self.folder, self._getFolderStem() + extension))
        )

    def getDatFile(self):
        if self._type in (_ProbeType.FAR_FIELD, _ProbeType.MOVIE):
            return None
        path = os.path.join(self.folder, self._getFolderStem() + ".dat")
        return path if os.path.isfile(path) else None

    def getTextFile(self):
        dat_file = self.getDatFile()
        if dat_file is not None:
            return dat_file
        if self._type == _ProbeType.FAR_FIELD:
            path = os.path.join(self.folder, self._getFolderStem())
            return path if os.path.isfile(path) else None
        return None

    def getXDMFFile(self):
        if self._type != _ProbeType.MOVIE:
            return None
        path = os.path.join(self.folder, self._getFolderStem() + ".xdmf")
        return path if os.path.isfile(path) else None

    def getBinFile(self):
        if self._type == _ProbeType.MTLN:
            return None
        path = os.path.join(self.folder, self._getFolderStem() + ".bin")
        return path if os.path.isfile(path) else None

    def getH5File(self):
        if self._type != _ProbeType.MOVIE:
            return None
        path = os.path.join(self.folder, self._getFolderStem() + ".h5")
        return path if os.path.isfile(path) else None

    def getExpectedColumns(self):
        try:
            columns = _EXPECTED_DAT_COLUMNS[(self._type, self._domain)]
        except KeyError as error:
            raise ValueError(
                f"Probe type {self._type.value} has no text output columns"
            ) from error

        if self._type == _ProbeType.MTLN:
            quantity = "voltage" if self._getTagFromFolder(self.folder) == "_V_" else "current"
            return [column.format(quantity=quantity) for column in columns]
        return list(columns)

    def _normalise_data_columns(self, tag):
        if self._type == _ProbeType.WIRE:
            if self._domain == _ProbeDomain.TIME:
                self.data = self.data.rename(columns={"t": "time", self.data.columns[1]: "current"})
            else:
                self.data = self.data.rename(
                    columns={
                        self.data.columns[0]: "frequency",
                        self.data.columns[1]: "magnitude",
                        self.data.columns[2]: "phase",
                    }
                )
        elif self._type == _ProbeType.BULK_CURRENT:
            self.data = self.data.rename(columns={"t": "time", self.data.columns[1]: "current"})
        elif self._type == _ProbeType.LINE_INTEGRAL:
            self.data = self.data.rename(columns={"t": "time", self.data.columns[1]: "lineIntegral"})
        elif self._type == _ProbeType.POINT and self._domain == _ProbeDomain.TIME:
            self.data = self.data.rename(columns={"t": "time", self.data.columns[1]: "field"})
            if len(self.data.columns) == 3:
                self.data = self.data.rename(columns={self.data.columns[2]: "incident"})
        elif self._type == _ProbeType.FAR_FIELD:
            self.data = self.data.rename(
                columns={
                    self.data.columns[0]: "freq",
                    self.data.columns[7]: "rcs_arit",
                    self.data.columns[8]: "rcs_geom",
                }
            )
        elif self._type == _ProbeType.MTLN:
            if self._domain == _ProbeDomain.TIME:
                self.data = self.data.rename(columns={"t": "time"})
            quantity = "voltage" if tag == "_V_" else "current"
            for index in range(self.data.shape[1] - 1):
                self.data = self.data.rename(
                    columns={self.data.columns[index + 1]: f"{quantity}_{index}"}
                )

    def __getitem__(self, key):
        if key not in self.data.columns:
            if key == "current" and "current_0" in self.data.columns:
                return self.data["current_0"]
            elif key == "current_0" and "current" in self.data.columns:
                return self.data["current"]

            if key == "voltage" and "voltage_0" in self.data.columns:
                return self.data["voltage_0"]
            elif key == "voltage_0" and "voltage" in self.data.columns:
                return self.data["voltage"]
        return self.data[key]

    @staticmethod
    def _getCaseNameFromFolder(folder):
        bn = os.path.basename(folder)
        if ".fdtd_" in bn:
            return bn.split(".fdtd_")[0]
        else:
            for tag in Probe.ALL_TAGS:
                if tag in bn:
                    return bn.split(tag)[0]

    @staticmethod
    def _getProbeNameFromFolder(folder):
        bn = os.path.basename(folder)
        if ".fdtd_" in bn:
            tag = Probe._getTagFromFolder(bn)
            bn_without_case_name = bn.split(".fdtd_", 1)[1]
            probe_name = bn_without_case_name.split(tag)[0]
            return probe_name
        else:
            return Probe._getCaseNameFromFolder(bn)

    @staticmethod
    def _getTagFromFolder(folder):
        bn = os.path.basename(folder)
        for tag in Probe.ALL_TAGS:
            if tag in bn:
                return tag
        raise ValueError("Unable to determine probe tag")

    def _getFolderStem(self):
        return os.path.basename(self.folder)

    @staticmethod
    def _getPositionStrFromFolder(folder):
        stem = os.path.basename(folder)
        tag = Probe._getTagFromFolder(folder)
        if tag in stem:
            position_str = stem.split(tag, 1)[1]
        elif stem.endswith(tag[:-1]):
            position_str = ""
        else:
            raise ValueError("Unable to determine probe position")
        for marker in Probe.DOMAIN_MARKERS:
            if position_str.endswith(marker):
                return position_str[: -len(marker)]
        return position_str

    @staticmethod
    def _getDomainTypeFromFolder(folder):
        if os.path.basename(folder).endswith(("_fq", "_df")):
            return "frequency"
        else:
            return "time"

    @staticmethod
    def _positionStrToCell(pos_str: str):
        pos = pos_str.split("_")
        return np.array([int(pos[0]), int(pos[1]), int(pos[2])])

    @staticmethod
    def _positionStrToTwoCells(pos_str: str):
        init_str, end_str = pos_str.split("__")
        return Probe._positionStrToCell(init_str), Probe._positionStrToCell(end_str)

    @staticmethod
    def _getFieldAndDirection(tag: str):
        return tag[1], tag[2]


class ExcitationFile:
    def __init__(self, excitation_filename):
        if isinstance(excitation_filename, os.PathLike):
            self.filename = excitation_filename.as_posix()
        else:
            self.filename = excitation_filename
        assert os.path.isfile(self.filename)

        self.data = pd.read_csv(self.filename, sep="\\s+", names=["time", "value"])

    def __getitem__(self, key):
        return self.data[key]


class FDTD:
    def __init__(
        self,
        input_filename,
        path_to_exe=None,
        flags=None,
        run_in_folder=None,
        mpi_command=None,
    ):

        self._setFilename(input_filename)

        if path_to_exe is None:
            semba_exe = os.path.join(os.getcwd(), DEFAULT_SEMBA_FDTD_PATH)
        else:
            semba_exe = path_to_exe
        assert os.path.isfile(semba_exe), f'Semba executable not found at: {semba_exe}'

        if mpi_command is None:
            mpi_command_parts = []
        else:
            mpi_command_parts = mpi_command.split()

        if flags is None:
            flags = []
        elif isinstance(flags, str):
            flags = flags.split()

        case_name = self.getCaseName() + ".json"
        self.run_command = mpi_command_parts + [semba_exe] + ["-i", case_name] + flags

        self._hasRun = False

        if run_in_folder != None:
            self._setNewFolder(run_in_folder)

    def getFolder(self):
        res = os.path.dirname(self._filename)
        if len(res) == 0:
            return "./"
        return res

    def getCaseName(self):
        return os.path.basename(self._filename).split(".json")[0]

    def __getitem__(self, key):
        return self._input[key]

    def _setFilename(self, newFilename):
        self._filename = os.path.abspath(os.fspath(newFilename))
        self._input = json.load(open(self._filename))

    def getUsedFiles(self):
        res = []

        # Files used to define magnitudes.
        if "sources" in self._input:
            for p in self._input["sources"]:
                if "magnitudeFile" in p:
                    res.append(p["magnitudeFile"])

        # Files used in transfer functions.
        if "probes" in self._input:
            for p in self._input["probes"]:
                if "magnitudeFile" in p:
                    res.append(p["magnitudeFile"])

        # .model files in terminations.
        if "materials" in self._input:
            for p in self._input["materials"]:
                if "terminations" in p:
                    for t in p["terminations"]:
                        if "file" in t and not t["file"] in res:
                            res.append(t["file"])

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
        if self._input != json.load(open(self._filename, "r")):
            json.dump(self._input, open(self._filename, "w"))

        os.chdir(self.getFolder())
        self.output = subprocess.run(self.run_command, capture_output=True)
        self._hasRun = True
        assert self.hasFinishedSuccessfully()

    def hasFinishedSuccessfully(self):
        if self._hasRun and (self.output.returncode == 0):
            return True
        else:
            if self.output.stdout:
                print(self.output.stdout.decode("utf-8"), end="")
            if self.output.stderr:
                print(self.output.stderr.decode("utf-8"), end="")
            return False

    def __getitem__(self, key):
        return self._input[key]

    def cleanUp(self):
        folder = self.getFolder()
        case_name = self.getCaseName()
        extensions = ("*.dat", "*.pl", "*.txt", "*.xdmf", "*.bin", "*.h5")
        for ext in extensions:
            files = glob.glob(os.path.join(folder, ext))
            for file in files:
                if os.path.basename(file).startswith(case_name):
                    os.remove(file)

        subfolders = [
            item
            for item in os.listdir(folder)
            if os.path.isdir(os.path.join(folder, item))
        ]
        for f in subfolders:
            if f.startswith(case_name):
                shutil.rmtree(os.path.join(folder, f), ignore_errors=True)

    def getSolvedProbeFilenames(self, probe_name, include_binary=False):
        if not "probes" in self._input:
            raise ValueError("Solver does not contain probes.")

        file_extensions = (".dat", ".xdmf", ".h5", "")
        if include_binary:
            file_extensions += (".bin",)
        search_root = os.path.abspath(self.getFolder())
        probe_prefix = self.getCaseName() + "_" + probe_name
        probeFiles = []
        for path in glob.glob(os.path.join(search_root, "**", "*"), recursive=True):
            if not os.path.isfile(path):
                continue
            basename = os.path.basename(path)
            if not basename.startswith(probe_prefix):
                continue
            if basename.endswith("_geometry.xdmf") or basename.endswith("_geometry.h5"):
                continue
            extension = next(
                (value for value in Probe.FILE_EXTENSIONS if basename.endswith(value)),
                "",
            )
            artifact_suffix = os.path.splitext(
                basename[len(probe_prefix) :].lstrip("_")
            )[1]
            if not extension and artifact_suffix:
                continue
            if not extension and not any(
                tag in basename[len(self.getCaseName()) + 1 :]
                for tag in Probe.FAR_FIELD_TAG
            ):
                continue
            if extension not in file_extensions:
                continue
            probeFiles.append(os.path.abspath(path))

        return sorted(probeFiles)

    def getExcitationFile(self, excitation_file_name):
        file_extensions = ("*.exc",)
        excitationFile = []
        for ext in file_extensions:
            newExcitationFile = [
                x for x in glob.glob(ext) if re.match(excitation_file_name, x)
            ]
            excitationFile.extend(newExcitationFile)

        if (len(excitationFile)) != 1:
            raise ValueError(
                "Unexpected number of excitation Files found: {}".format(excitationFile)
            )

        return excitationFile

    def getVTKMap(self):
        map_files = sorted(glob.glob(os.path.join("**", "*.vtu"), recursive=True))
        if not map_files:
            map_files = sorted(glob.glob(os.path.join("**", "*_1.vtk"), recursive=True))
        return map_files[0] if map_files else None

    def getCurrentVTKMap(self):
        current_maps = sorted(
            glob.glob(os.path.join("**", "*current*.vtu"), recursive=True)
        )
        if not current_maps:
            current_maps = sorted(
                glob.glob(os.path.join("**", "*_1_current.vtk"), recursive=True)
            )
        return current_maps[0] if current_maps else None

    def getCurrentMovie(self):
        current_movies = sorted(
            glob.glob(os.path.join("**", "*current_movie*", "*.xdmf"), recursive=True)
        )
        return current_movies[0] if current_movies else None

    def getCurrentMovieGeometry(self):
        geometries = sorted(
            glob.glob(os.path.join("**", "*current_movie*", "*_geometry.xdmf"), recursive=True)
        )
        return geometries[0] if geometries else None

    def getMaterialProperties(self, materialName):
        if "materials" in self._input:
            for idx, element in enumerate(self._input["materials"]):
                if element.get("name") == materialName:
                    return self._input["materials"][idx]
