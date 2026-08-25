from test.utils.utils import *
from pathlib import Path
import re


def assert_current_movie_has_tag_number(solver, expected_tags=None):
    import h5py

    xdmf_path = Path(solver.getCurrentMovie())
    hdf_path = xdmf_path.with_suffix(".h5")
    xdmf_contents = xdmf_path.read_text()

    assert xdmf_path.is_file()
    assert hdf_path.is_file()
    tag_match = re.search(
        r'<Attribute Name="tagnumber".*?/attributes/(a\d+)/values',
        xdmf_contents,
        re.DOTALL,
    )
    assert tag_match is not None
    with h5py.File(hdf_path, "r") as hdf_file:
        tag_values = hdf_file[f"attributes/{tag_match.group(1)}/values"][0]
        for coordinates, expected_tag in (expected_tags or {}).items():
            assert tag_values[coordinates] == expected_tag


def current_movie_geometry_tag_counts(solver):
    import h5py

    geometry_path = Path(solver.getCurrentMovieGeometry())
    assert geometry_path.is_file()
    with h5py.File(geometry_path.with_suffix(".h5"), "r") as hdf_file:
        tags = hdf_file["attributes/a0001/values"][()]
    values, counts = np.unique(tags, return_counts=True)
    return dict(zip(values, counts))


def current_movie_geometry_all_tag_counts(solver):
    import h5py

    geometry_path = Path(solver.getCurrentMovieGeometry())
    attribute_ids = re.findall(
        r'<Attribute Name="tagnumber".*?/attributes/(a\d+)/values',
        geometry_path.read_text(),
        re.DOTALL,
    )
    with h5py.File(geometry_path.with_suffix(".h5"), "r") as hdf_file:
        tags = np.concatenate(
            [hdf_file[f"attributes/{attribute_id}/values"][()] for attribute_id in attribute_ids]
        )
    values, counts = np.unique(tags, return_counts=True)
    return dict(zip(values, counts))


@pytest.mark.planewave
def test_fdtd_set_new_folder_to_run(tmp_path):
    input_filename = os.path.join(CASES_FOLDER, "planewave", "pw-in-box.fdtd.json")
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver["general"]["numberOfSteps"] = 1

    solver.run()


@pytest.mark.planewave
def test_fdtd_with_string_args(tmp_path):
    input_filename = os.path.join(CASES_FOLDER, "planewave", "pw-in-box.fdtd.json")
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path, flags="-h")
    solver["general"]["numberOfSteps"] = 1

    solver.run()


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.planewave
def test_fdtd_with_mpi_run(tmp_path):
    input_filename = os.path.join(CASES_FOLDER, "planewave", "pw-in-box.fdtd.json")
    solver = FDTD(
        input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-h"],
        mpi_command="mpirun -np 2",
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()


@pytest.mark.planewave
def test_fdtd_clean_up_after_run(tmp_path):
    input = CASES_FOLDER + "planewave/pw-in-box.fdtd.json"
    solver = FDTD(input, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver["general"]["numberOfSteps"] = 1

    solver.run()

    solved_probe_folders = solver.getSolvedProbeFolders("inbox")
    assert os.path.isdir(solved_probe_folders[0])

    solver.cleanUp()

    assert not os.path.exists(solved_probe_folders[0])


@pytest.mark.planewave
def test_fdtd_probe_folders_include_binary_artifacts(tmp_path):
    input = CASES_FOLDER + "planewave/pw-in-box.fdtd.json"
    solver = FDTD(input, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 1
    solver.run()

    probe_folders = solver.getSolvedProbeFolders("inbox")

    assert len(probe_folders) == 1
    probe = Probe(probe_folders[0])
    assert probe.getDatFile() is not None
    assert probe.getBinFile() is not None


@pytest.mark.planewave
def test_fdtd_clean_up_does_not_delete_other_cases_files(tmp_path):
    input = CASES_FOLDER + "planewave/pw-in-box.fdtd.json"
    solver = FDTD(input, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    case_name = solver.getCaseName()
    other_case_name = "other_case.fdtd"

    own_file = os.path.join(str(tmp_path), case_name + "_probe_Ex_1_2_3.dat")
    other_file = os.path.join(str(tmp_path), other_case_name + "_probe_Ex_1_2_3.dat")

    open(own_file, "w").close()
    open(other_file, "w").close()

    solver.cleanUp()

    assert not os.path.isfile(own_file)
    assert os.path.isfile(other_file)


@pytest.mark.wires
@pytest.mark.termination
def test_holland_case_checking_number_of_outputs_single_wire(tmp_path):
    input_filename = CASES_FOLDER + "holland/holland1981.fdtd.json"
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    number_of_steps = 10
    solver["general"]["numberOfSteps"] = number_of_steps

    solver["materials"][0] = createWire(id=1, r=0.02)
    solver.run()

    probe_files = solver.getSolvedProbeFolders("mid_point")

    assert len(probe_files) == 1
    output_probe = Probe(probe_files[0])
    assert len(output_probe["current"]) == 11


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.termination
def test_holland_case_checking_number_of_outputs_wire(tmp_path):
    input_filename = CASES_FOLDER + "holland/holland1981.fdtd.json"
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    number_of_steps = 10
    solver["general"]["numberOfSteps"] = number_of_steps

    solver["materials"][0] = createWire(id=1, r=0.02)
    solver.run()

    probe_files = solver.getSolvedProbeFolders("mid_point")

    assert len(probe_files) == 1
    output_probe = Probe(probe_files[0])
    assert len(output_probe["current"]) == 11


@no_mtln_skip
@pytest.mark.mtln
@pytest.mark.wires
@pytest.mark.multiwire
@pytest.mark.termination
def test_holland_case_checking_number_of_outputs_unshielded(tmp_path):
    input_filename = CASES_FOLDER + "holland/holland1981.fdtd.json"
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    number_of_steps = 10
    solver["general"]["numberOfSteps"] = number_of_steps
    solver["materials"][0] = createUnshieldedWire(
        id=1, lpul=6.52188703e-08, cpul=1.7060247700000001e-10
    )
    expected_output_filename = "holland1981.fdtd_mid_point_single_wire_I_11_11_12.dat"
    solver.run()

    probe_files = solver.getSolvedProbeFolders("mid_point")

    assert len(probe_files) == 1
    output_probe = Probe(probe_files[0])
    assert len(output_probe["current"]) == 11


@mtln_skip
@pytest.mark.wires
@pytest.mark.probes
@pytest.mark.termination
def test_towel_hanger_case_creates_output_probes(tmp_path):
    input_filename = CASES_FOLDER + "towelHanger/towelHanger.fdtd.json"
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    probe_start = solver.getSolvedProbeFolders("wire_start")
    probe_mid = solver.getSolvedProbeFolders("wire_mid")
    probe_end = solver.getSolvedProbeFolders("wire_end")

    assert len(probe_start) == 1
    assert len(probe_mid) == 1
    assert len(probe_end) == 1

    assert "towelHanger.fdtd_wire_start_Wz_27_25_30_s3" == Path(probe_start[0]).name
    assert "towelHanger.fdtd_wire_mid_Wx_35_25_32_s13" == Path(probe_mid[0]).name
    assert "towelHanger.fdtd_wire_end_Wz_43_25_30_s22" == Path(probe_end[0]).name
    assert countLinesInFile(Probe(probe_start[0]).getDatFile()) == 3
    assert countLinesInFile(Probe(probe_mid[0]).getDatFile()) == 3
    assert countLinesInFile(Probe(probe_end[0]).getDatFile()) == 3


@pytest.mark.wires
def test_simple_cabin_initialization(tmp_path):
    input_filename = CASES_FOLDER + "simple_cabin/simple_cabin.fdtd.json"
    solver = FDTD(input_filename=input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    assert solver.hasFinishedSuccessfully()


@pytest.mark.probes
@pytest.mark.hdf
@pytest.mark.farfield
def test_sphere_case_with_far_field_probe_launches(tmp_path):
    input_filename = CASES_FOLDER + "sphere/sphere.fdtd.json"
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 1
    solver["probes"][0]["domain"]["numberOfFrequencies"] = 1
    solver["mesh"]["grid"]["numberOfCells"] = [4, 4, 4]
    solver["mesh"]["elements"] = solver["mesh"]["elements"][:2]
    solver["mesh"]["elements"][0]["intervals"] = [[[0, 0, 0], [4, 4, 4]]]
    solver["mesh"]["elements"][1]["intervals"] = [[[1, 1, 1], [3, 3, 3]]]
    solver["materialAssociations"][:] = []

    solver.run()

    far_field_files = solver.getSolvedProbeFolders("farfield")
    assert len(far_field_files) == 1
    far_field_probe = Probe(far_field_files[0])
    assert far_field_probe.case_name == "sphere"
    assert far_field_probe.type == "farField"
    assert np.all(far_field_probe.cell_init == np.array([1, 1, 1]))

    movie_folders = solver.getSolvedProbeFolders("electric_field_movie")
    assert len(movie_folders) == 1
    movie_probe = Probe(movie_folders[0])
    assert movie_probe.getXDMFFile() is not None
    assert movie_probe.case_name == "sphere"
    assert movie_probe.type == "movie"
    assert np.all(movie_probe.cell_init == np.array([1, 1, 1]))


@pytest.mark.conformal
@pytest.mark.vtk
def test_fill_conformal_vtk_sphere(tmp_path):
    input_filename = CASES_FOLDER + "conformal/conformal_sphere_1mm_rcs_delta.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtk_map_filename = solver.getVTKMap()
    assert os.path.isfile(vtk_map_filename)

    line_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=3, property="mediatype"
    )

    face_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=9, property="mediatype"
    )

    assert line_media_dict[0.5] == 12  # PEC line
    assert line_media_dict[2004] == 24  # Conformal line

    assert face_media_dict[0] == 6  # PEC surface
    assert face_media_dict[1005] == 24  # Conformal PEC surface
    assert face_media_dict[1006] == 24  # Conformal PEC surface


@pytest.mark.conformal
@pytest.mark.vtk
def test_fill_conformal_fL_0_005_vtk_large_sphere(tmp_path):
    input_filename = CASES_FOLDER + "conformal/conformal_fL_sphere_rcs.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtk_map_filename = solver.getVTKMap()
    assert os.path.isfile(vtk_map_filename)

    line_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=3, property="mediatype"
    )

    assert -0.5 not in line_media_dict.keys()

    face_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=9, property="mediatype"
    )

    assert -1 not in face_media_dict.keys()


@pytest.mark.conformal
@pytest.mark.vtk
def test_fill_conformal_fL_0_15_vtk_large_sphere(tmp_path):
    input_filename = CASES_FOLDER + "conformal/conformal_fL_0.15_sphere_rcs.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtk_map_filename = solver.getVTKMap()
    assert os.path.isfile(vtk_map_filename)

    line_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=3, property="mediatype"
    )

    assert -0.5 not in line_media_dict.keys()

    face_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=9, property="mediatype"
    )

    assert -1 not in face_media_dict.keys()


@pytest.mark.conformal
@pytest.mark.vtk
def test_fill_slanted_vtk_large_sphere(tmp_path):
    input_filename = CASES_FOLDER + "conformal/slanted_sphere_rcs.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtk_map_filename = solver.getVTKMap()
    assert os.path.isfile(vtk_map_filename)

    line_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=3, property="mediatype"
    )

    assert -0.5 not in line_media_dict.keys()

    face_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=9, property="mediatype"
    )

    assert -1 not in face_media_dict.keys()


@pytest.mark.conformal
@pytest.mark.vtk
def test_fill_conformal_vtk_corner(tmp_path):
    #          /|
    #        5  |
    #      / |\ |
    #    3___|_4|_______
    #    |   | ||_______|______
    #    |   | |        |      /
    #    |    6|        |    /
    #    |  / \|        |  /
    #    1/____2________|/

    input_filename = CASES_FOLDER + "conformal/conformal_corner.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtk_map_filename = solver.getVTKMap()
    assert os.path.isfile(vtk_map_filename)

    face_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=9, property="mediatype"
    )

    assert 0 not in face_media_dict.keys()
    assert face_media_dict[1005] == 2  # Conformal PEC surface #1
    assert face_media_dict[1006] == 2  # Conformal PEC surface #2

    line_media_dict = createPropertyDictionary(
        vtk_map_filename, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 1  # PEC line
    assert line_media_dict[2004] == 4  # Conformal line #1


@pytest.mark.probes
@pytest.mark.movie
def test_movie_with_frequency_domain(tmp_path):
    input_filename = CASES_FOLDER + "observation/movieFrequency.fdtd.json"
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 100
    solver["general"]["timeStep"] = 2.0e-9
    solver["probes"][0]["domain"]["numberOfFrequencies"] = 100

    solver.run()

    movie_folders = solver.getSolvedProbeFolders("movie_electric")
    assert len(movie_folders) == 1
    movie_probe = Probe(movie_folders[0])
    assert movie_probe.getXDMFFile() is not None
    assert movie_probe.case_name == "movieFrequency"
    assert movie_probe.type == "movie"
    assert np.all(movie_probe.cell_init == np.array([1, 1, 1]))


@pytest.mark.probes
@pytest.mark.movie
def test_movie_with_time_domain(tmp_path):
    input_filename = CASES_FOLDER + "observation/movieTime.fdtd.json"
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 1
    solver["probes"][0]["domain"]["samplingPeriod"] = 1e-9

    solver.run()

    movie_folders = solver.getSolvedProbeFolders("movie_electric")
    assert len(movie_folders) == 1
    movie_probe = Probe(movie_folders[0])
    assert movie_probe.getXDMFFile() is not None
    assert movie_probe.case_name == "movieTime"
    assert movie_probe.type == "movie"
    assert np.all(movie_probe.cell_init == np.array([1, 1, 1]))


@pytest.mark.sgbc
@pytest.mark.vtk
def test_three_surfaces(tmp_path):
    input_filename = CASES_FOLDER + "observation/three_surfaces.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert face_tag_dict[64] == 4
    assert face_tag_dict[128] == 4
    assert face_tag_dict[192] == 4

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 8
    assert line_tag_dict[128] == 6
    assert line_tag_dict[192] == 4

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert face_media_dict[0] == 4  # PEC surface
    assert face_media_dict[304] == 4  # SGBC surface
    assert face_media_dict[305] == 4  # SGBC surface

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 8  # PEC line
    assert line_media_dict[3.5] == 10  # SGBC line


@pytest.mark.sgbc
@pytest.mark.probes
def test_three_surfaces_Jprobe(tmp_path):
    input_filename = CASES_FOLDER + "observation/three_surfaces_Jprobe.fdtd.json"
    solver = FDTD(input_filename=input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    assert_current_movie_has_tag_number(
        solver,
        {
            (3, 3, 3): 64,
            (3, 5, 4): 64,
            (0, 0, 0): 0,
        },
    )
    tag_counts = current_movie_geometry_all_tag_counts(solver)
    assert tag_counts[64] == 12
    assert tag_counts[128] == 10
    assert tag_counts[192] == 8


@pytest.mark.sgbc
@pytest.mark.probes
def test_three_one_cell_surfaces_Jprobe(tmp_path):
    input_filename = CASES_FOLDER + "observation/three_one_cell_surfaces_Jprobe.fdtd.json"
    solver = FDTD(input_filename=input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    assert_current_movie_has_tag_number(
        solver,
        {
            (3, 3, 3): 64,
            (3, 4, 3): 64,
            (0, 0, 0): 0,
        },
    )
    tag_counts = current_movie_geometry_tag_counts(solver)
    assert tag_counts[64] == 4
    assert tag_counts[128] == 3
    assert tag_counts[192] == 2


@pytest.mark.probes
def test_one_cell_PEC_surface_Jprobe(tmp_path):
    input_filename = CASES_FOLDER + "observation/one_cell_surface_Jprobe.fdtd.json"
    solver = FDTD(input_filename=input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver["general"]["numberOfSteps"] = 1
    solver["materialAssociations"][0]["materialId"] = 1

    for probe_component in ["x", "y", "z"]:
        solver["probes"][0]["component"] = probe_component
        solver.cleanUp()
        solver.run()
        expected_tag = 0 if probe_component == "z" else 64
        assert_current_movie_has_tag_number(solver, {(3, 3, 3): expected_tag, (0, 0, 0): 0})
        assert current_movie_geometry_all_tag_counts(solver)[64] == 5


def test_one_cell_SGBC_surface_Jprobe(tmp_path):
    input_filename = CASES_FOLDER + "observation/one_cell_surface_Jprobe.fdtd.json"
    solver = FDTD(input_filename=input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver["general"]["numberOfSteps"] = 1
    solver["materialAssociations"][0]["materialId"] = 2

    for probe_component in ["x", "y", "z"]:
        solver["probes"][0]["component"] = probe_component
        solver.cleanUp()
        solver.run()
        assert_current_movie_has_tag_number(solver, {(3, 3, 3): 0, (0, 0, 0): 0})
        assert current_movie_geometry_all_tag_counts(solver)[64] == 5


def test_1_volume(tmp_path):
    input_filename = CASES_FOLDER + "observation/pec_volume.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert face_tag_dict[64] == 36

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert len(line_tag_dict) == 0

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert face_media_dict[0] == 36  # PEC surface

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert len(line_media_dict) == 0


def test_2_volumes(tmp_path):
    input_filename = CASES_FOLDER + "observation/pec_volumes.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert face_tag_dict[64] == 36
    assert face_tag_dict[128] == 36

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert len(line_tag_dict) == 0

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert face_media_dict[0] == 72  # PEC surface

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert len(line_media_dict) == 0


def test_1_line(tmp_path):
    input_filename = CASES_FOLDER + "observation/pec_line.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line


def test_volume_and_surfaces(tmp_path):
    input_filename = CASES_FOLDER + "observation/volume_and_surfaces.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert face_tag_dict[64] == 6
    assert face_tag_dict[128] == 1
    assert face_tag_dict[192] == 1

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 1
    assert line_tag_dict[128] == 4
    assert line_tag_dict[192] == 3

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert face_media_dict[16] == 1  # PMC surface
    assert face_media_dict[0] == 6  # PEC surface
    assert face_media_dict[305] == 1  # SGBC surface

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[16.5] == 4  # PMC line
    assert line_media_dict[0.5] == 1  # PEC line
    assert line_media_dict[3.5] == 3  # SGBC line


def test_count_bug(tmp_path):
    input_filename = CASES_FOLDER + "observation/count_bug.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1
    solver.run()

    solver["materialAssociations"][0]["materialId"] = 3
    solver["materialAssociations"][1]["materialId"] = 1
    solver["materialAssociations"][2]["materialId"] = 3
    solver.cleanUp()
    solver.run()

    solver["materialAssociations"][0]["materialId"] = 3
    solver["materialAssociations"][1]["materialId"] = 3
    solver["materialAssociations"][2]["materialId"] = 1
    solver.cleanUp()
    solver.run()


@mtln_skip
def test_wires(tmp_path):
    input_filename = CASES_FOLDER + "observation/wires.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 4
    assert line_tag_dict[128] == 6
    assert line_tag_dict[192] == 4

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[7] == 7
    assert line_media_dict[10] == 6
    assert line_media_dict[21] == 1


@mtln_skip
def test_wires_collision_Jprobe(tmp_path):
    input_filename = CASES_FOLDER + "observation/wires_collision_Jprobe.fdtd.json"
    solver = FDTD(input_filename=input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    assert_current_movie_has_tag_number(solver, {(3, 3, 3): 0, (0, 0, 0): 0})
    tag_counts = current_movie_geometry_tag_counts(solver)
    assert tag_counts[64] == 2
    assert tag_counts[128] == 4
    assert tag_counts[192] == 6
    assert tag_counts[256] == 4


@mtln_skip
def test_wires_collision(tmp_path):
    input_filename = CASES_FOLDER + "observation/wires_collision.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2
    assert line_tag_dict[128] == 4
    assert line_tag_dict[192] == 6
    assert line_tag_dict[256] == 4

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[7] == 6
    assert line_media_dict[8] == 1
    assert line_media_dict[10] == 6
    assert line_media_dict[21] == 1
    assert line_media_dict[0.5] == 2  # PEC line


@mtln_skip
def test_wire_x_collision_y(tmp_path):
    input_filename = CASES_FOLDER + "observation/wire_x_collision_y.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  # Wire w/o collision
    assert line_media_dict[8] == 1  # Wire touching non wire
    assert line_media_dict[10] == 2  # Wire extreme


@mtln_skip
def test_wire_x_collision_y_Jprobe(tmp_path):
    input_filename = CASES_FOLDER + "observation/wire_x_collision_y_Jprobe.fdtd.json"
    solver = FDTD(input_filename=input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    assert_current_movie_has_tag_number(solver, {(3, 3, 3): 128, (0, 0, 0): 0})
    tag_counts = current_movie_geometry_tag_counts(solver)
    assert tag_counts[64] == 2
    assert tag_counts[128] == 4


@mtln_skip
def test_wire_x_collision_z(tmp_path):
    input_filename = CASES_FOLDER + "observation/wire_x_collision_z.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  # Wire w/o collision
    assert line_media_dict[8] == 1  # Wire touching non wire
    assert line_media_dict[10] == 2  # Wire extreme


@mtln_skip
def test_wire_x_long_collision_z(tmp_path):
    input_filename = CASES_FOLDER + "observation/wire_x_long_collision_z.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 8  # Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 3  # Wire w/o collision
    assert line_media_dict[8] == 1  # Wire touching non wire
    assert line_media_dict[10] == 2  # Wire extreme
    assert line_media_dict[21] == 2  # Wire extreme


@mtln_skip
def test_wire_y_collision_x(tmp_path):
    input_filename = CASES_FOLDER + "observation/wire_y_collision_x.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  # Wire w/o collision
    assert line_media_dict[8] == 1  # Wire touching non wire
    assert line_media_dict[10] == 2  # Wire extreme


@mtln_skip
def test_wire_y_collision_z(tmp_path):
    input_filename = CASES_FOLDER + "observation/wire_y_collision_z.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  # Wire w/o collision
    assert line_media_dict[8] == 1  # Wire touching non wire
    assert line_media_dict[10] == 2  # Wire extreme


@mtln_skip
def test_wire_z_collision_x(tmp_path):
    input_filename = CASES_FOLDER + "observation/wire_z_collision_x.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  # Wire w/o collision
    assert line_media_dict[8] == 1  # Wire touching non wire
    assert line_media_dict[10] == 2  # Wire extreme


@mtln_skip
def test_wire_z_collision_y(tmp_path):
    input_filename = CASES_FOLDER + "observation/wire_z_collision_y.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  # Wire w/o collision
    assert line_media_dict[8] == 1  # Wire touching non wire
    assert line_media_dict[10] == 2  # Wire extreme


@no_mtln_skip
def test_multiwire_z_collision_y(tmp_path):
    input_filename = CASES_FOLDER + "observation/multiwire_z_collision_y.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  # MultiWire w/o collision
    assert line_media_dict[13] == 1  # Multiwire touching non multiwire
    assert line_media_dict[14] == 2  # Multiwire extreme


@no_mtln_skip
def test_multiwire_z_collision_x(tmp_path):
    input_filename = CASES_FOLDER + "observation/multiwire_z_collision_x.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  # MultiWire w/o collision
    assert line_media_dict[13] == 1  # Multiwire touching non multiwire
    assert line_media_dict[14] == 2  # Multiwire extreme


@no_mtln_skip
def test_multiwire_y_collision_x(tmp_path):
    input_filename = CASES_FOLDER + "observation/multiwire_y_collision_x.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  # MultiWire w/o collision
    assert line_media_dict[13] == 1  # Multiwire touching non multiwire
    assert line_media_dict[14] == 2  # Multiwire extreme


@no_mtln_skip
def test_multiwire_y_collision_z(tmp_path):
    input_filename = CASES_FOLDER + "observation/multiwire_y_collision_z.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  # MultiWire w/o collision
    assert line_media_dict[13] == 1  # Multiwire touching non multiwire
    assert line_media_dict[14] == 2  # Multiwire extreme


@no_mtln_skip
def test_multiwire_x_collision_y(tmp_path):
    input_filename = CASES_FOLDER + "observation/multiwire_x_collision_y.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  # MultiWire w/o collision
    assert line_media_dict[13] == 1  # Multiwire touching non multiwire
    assert line_media_dict[14] == 2  # Multiwire extreme


@no_mtln_skip
def test_multiwire_x_collision_z(tmp_path):
    input_filename = CASES_FOLDER + "observation/multiwire_x_collision_z.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 4  # Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  # MultiWire w/o collision
    assert line_media_dict[13] == 1  # Multiwire touching non multiwire
    assert line_media_dict[14] == 2  # Multiwire extreme


@no_mtln_skip
def test_multiwire_x_long_collision_z(tmp_path):
    input_filename = CASES_FOLDER + "observation/multiwire_x_long_collision_z.fdtd.json"
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=tmp_path,
        flags=["-mapvtk"],
    )
    solver["general"]["numberOfSteps"] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="tagnumber"
    )
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="tagnumber"
    )
    assert line_tag_dict[64] == 2  # PEC
    assert line_tag_dict[128] == 8  # Multiwire:
    # 3 segments w/o collision, 2 of them are also intermediate
    # 1 segment adjacent something not multiwire, 2 extremes,

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property="mediatype"
    )
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property="mediatype"
    )
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 3  # MultiWire w/o collision
    assert line_media_dict[13] == 1  # Multiwire touching non multiwire
    assert line_media_dict[14] == 2  # Multiwire extreme
    assert line_media_dict[61] == 2  # Intermediate multiwire segment


def test_can_assign_same_surface_impedance_to_multiple_geometries(tmp_path):
    input_filename = CASES_FOLDER + "multipleAssigments/multipleSurfaceImpedance.fdtd.json"

    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    probe_files = solver.getSolvedProbeFolders("BulkProbeEntry")
    assert len(probe_files) == 1
    assert Probe(probe_files[0]) is not None


def test_can_assign_same_dielectric_material_to_multiple_geometries(tmp_path):
    input_filename = CASES_FOLDER + "multipleAssigments/multipleDielectricMaterial.fdtd.json"

    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    probe_files = solver.getSolvedProbeFolders("BulkProbeEntry")
    assert len(probe_files) == 1
    assert Probe(probe_files[0]) is not None


def test_can_execute_fdtd_from_folder_with_spaces_and_can_process_additional_arguments(
    tmp_path,
):
    folder_with_spaces: str = os.path.join(tmp_path, "spaced bin")
    os.mkdir(folder_with_spaces)
    if platform == "win32":
        shutil.copy2(NGSPICE_DLL, folder_with_spaces)

    semba_executable = SEMBA_EXE.split(os.path.sep)[-1]
    executable_path: str = os.path.join(folder_with_spaces, semba_executable)
    shutil.copy2(SEMBA_EXE, executable_path)
    print(executable_path)

    input_filename = CASES_FOLDER + "dielectric/dielectricTransmission.fdtd.json"
    solver = FDTD(input_filename, path_to_exe=executable_path, run_in_folder=tmp_path)
    solver.run()
    probe_files = solver.getSolvedProbeFolders("outside")
    assert len(probe_files) == 1
    assert Probe(probe_files[0]) is not None
    assert solver.getVTKMap()[0] is not None


def test_bulk_current_outputs(tmp_path):
    # This test uses bulk_probe_cases_over_nodal_source.fdtd from input_examples as input.
    # Verifies all kind of bulk probes are recognised and setted properly by checking outputFile format.
    input_filename = PROBES_INPUT_EXAMPLE + "bulk_probe_cases_over_nodal_source.fdtd.json"
    solver = FDTD(input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    bulk_x_plane_filenames = solver.getSolvedProbeFolders("BulkXPlane")
    bulk_y_plane_filenames = solver.getSolvedProbeFolders("BulkYPlane")
    bulk_z_plane_filenames = solver.getSolvedProbeFolders("BulkZPlane")
    bulk_y_point_filenames = solver.getSolvedProbeFolders("BulkYPoint")
    bulk_z_volume_filenames = solver.getSolvedProbeFolders("BulkZVolume")

    assert len(bulk_x_plane_filenames) == 1
    assert len(bulk_y_plane_filenames) == 1
    assert len(bulk_z_plane_filenames) == 1
    assert len(bulk_y_point_filenames) == 1
    assert len(bulk_z_volume_filenames) == 10

    bulk_x_plane_probe = Probe(bulk_x_plane_filenames[0])
    bulk_y_plane_probe = Probe(bulk_y_plane_filenames[0])
    bulk_z_plane_probe = Probe(bulk_z_plane_filenames[0])
    bulk_y_point_probe = Probe(bulk_y_point_filenames[0])
    bulk_z_volume_probe = Probe(bulk_z_volume_filenames[0])

    assert bulk_x_plane_probe.direction == "x"
    assert bulk_y_plane_probe.direction == "y"
    assert bulk_z_plane_probe.direction == "z"
    assert bulk_y_point_probe.direction == "y"


def test_wires_vtk(tmp_path):
    input_filename = CASES_FOLDER + "wire_vtk/wire_vtk.fdtd.json"

    solver = FDTD(input_filename=input_filename, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    vtkmapfile = solver.getVTKMap()
    reader = pv.get_reader(vtkmapfile)
    mesh = reader.read()
    assert mesh.n_cells != 0
