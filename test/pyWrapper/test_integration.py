from utils import *

def test_holland_case_checking_number_of_outputs_single_wire(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    number_of_steps = 10
    solver['general']['numberOfSteps'] = number_of_steps

    solver['materials'][0] = createWire(id = 1, r = 0.02)
    solver.run()

    probe_files = solver.getSolvedProbeFilenames("mid_point")

    assert len(probe_files) == 1
    p = Probe(probe_files[0])
    assert len(p['current']) == 10

@no_mtln_skip
def test_holland_case_checking_number_of_outputs_wire(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    number_of_steps = 10
    solver['general']['numberOfSteps'] = number_of_steps

    solver['materials'][0] = createWire(id = 1, r = 0.02)
    outfile = 'holland1981.fdtd_mid_point_single_wire_I_11_11_12.dat'
    solver.run()

    probe_files = solver.getSolvedProbeFilenames("mid_point")

    assert len(probe_files) == 1
    p = Probe(probe_files[0])
    assert len(p['current']) == 10

@no_mtln_skip
def test_holland_case_checking_number_of_outputs_unshielded(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    number_of_steps = 10
    solver['general']['numberOfSteps'] = number_of_steps
    solver['materials'][0] = createUnshieldedWire(id = 1, lpul = 6.52188703e-08, cpul = 1.7060247700000001e-10)        
    outfile = 'holland1981.fdtd_mid_point_single_wire_I_11_11_12.dat'
    solver.run()

    probe_files = solver.getSolvedProbeFilenames("mid_point")

    assert len(probe_files) == 1
    p = Probe(probe_files[0])
    assert len(p['current']) == 10


@mtln_skip
def test_towel_hanger_case_creates_output_probes(tmp_path):
    fn = CASES_FOLDER + 'towelHanger/towelHanger.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver['general']['numberOfSteps'] = 1

    solver.run()

    probe_start = solver.getSolvedProbeFilenames("wire_start")
    probe_mid = solver.getSolvedProbeFilenames("wire_mid")
    probe_end = solver.getSolvedProbeFilenames("wire_end")

    assert len(probe_start) == 1
    assert len(probe_mid) == 1
    assert len(probe_end) == 1

    assert 'towelHanger.fdtd_wire_start_Wz_27_25_30_s3.dat' == probe_start[0]
    assert 'towelHanger.fdtd_wire_mid_Wx_35_25_32_s13.dat' == probe_mid[0]
    assert 'towelHanger.fdtd_wire_end_Wz_43_25_30_s22.dat' == probe_end[0]
    assert countLinesInFile(probe_start[0]) == 3
    assert countLinesInFile(probe_mid[0]) == 3
    assert countLinesInFile(probe_end[0]) == 3

@no_mpi_skip
def test_airplane_case_with_mpi(tmp_path):
    fn = CASES_FOLDER + 'airplane/airplane.fdtd.json'
    solver = FDTD(fn, 
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, 
                  flags=['-mapvtk'],
                  mpi_command='mpirun -np 2')
    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)


def test_sphere_case_with_far_field_probe_launches(tmp_path):
    fn = CASES_FOLDER + 'sphere/sphere.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver['general']['numberOfSteps'] = 1
    solver['probes'][0]['domain']['numberOfFrequencies'] = 100

    solver.run()

    p = Probe(solver.getSolvedProbeFilenames("far")[0])
    assert p.case_name == 'sphere'
    assert p.type == 'farField'
    assert np.all(p.cell_init == np.array([2, 2, 2]))

    p = Probe(solver.getSolvedProbeFilenames("electric_field_movie")[0])
    assert p.case_name == 'sphere'
    assert p.type == 'movie'
    assert np.all(p.cell_init == np.array([2, 2, 2]))

def test_fill_conformal_vtk_sphere(tmp_path):
    fn = CASES_FOLDER + 'conformal/conformal_sphere_1mm_rcs_delta.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')

    assert line_media_dict[0.5] == 12  # PEC line
    assert line_media_dict[2004] == 24  # Conformal line

    assert face_media_dict[0] == 6  # PEC surface
    assert face_media_dict[1005] == 24  # Conformal PEC surface
    assert face_media_dict[1006] == 24  # Conformal PEC surface

def test_fill_conformal_fL_0_005_vtk_large_sphere(tmp_path):
    fn = CASES_FOLDER + 'conformal/conformal_fL_sphere_rcs.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    assert -0.5 not in line_media_dict.keys()

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')

    assert -1 not in face_media_dict.keys()

def test_fill_conformal_fL_0_15_vtk_large_sphere(tmp_path):
    fn = CASES_FOLDER + 'conformal/conformal_fL_0.15_sphere_rcs.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    assert -0.5 not in line_media_dict.keys()

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')

    assert -1 not in face_media_dict.keys()

def test_fill_slanted_vtk_large_sphere(tmp_path):
    fn = CASES_FOLDER + 'conformal/slanted_sphere_rcs.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    assert -0.5 not in line_media_dict.keys()

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')

    assert -1 not in face_media_dict.keys()

    
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
    
    
    
    fn = CASES_FOLDER + 'conformal/conformal_corner.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    
    assert(0 not in face_media_dict.keys())
    assert face_media_dict[1005] == 2  # Conformal PEC surface #1
    assert face_media_dict[1006] == 2  # Conformal PEC surface #2

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 1  # PEC line
    assert line_media_dict[2004] == 4  # Conformal line #1
    
def test_movie_with_frequency_domain(tmp_path):
    fn = CASES_FOLDER + 'observation/movieFrequency.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver['general']['numberOfSteps'] = 100
    solver['general']['timeStep'] = 2.0e-9
    solver['probes'][0]['domain']['numberOfFrequencies'] = 100

    solver.run()

    p = Probe(solver.getSolvedProbeFilenames("movie_electric")[0])
    assert p.case_name == 'movieFrequency'
    assert p.type == 'movie'
    assert np.all(p.cell_init == np.array([1, 1, 1]))

def test_movie_with_time_domain(tmp_path):
    fn = CASES_FOLDER + 'observation/movieTime.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver['general']['numberOfSteps'] = 1
    solver['probes'][0]['domain']['samplingPeriod'] = 1e-9

    solver.run()

    p = Probe(solver.getSolvedProbeFilenames("movie_electric")[0])
    assert p.case_name == 'movieTime'
    assert p.type == 'movie'
    assert np.all(p.cell_init == np.array([1, 1, 1]))




def test_three_surfaces(tmp_path):
    fn = CASES_FOLDER + 'observation/three_surfaces.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert face_tag_dict[64] == 4
    assert face_tag_dict[128] == 4
    assert face_tag_dict[192] == 4

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 8
    assert line_tag_dict[128] == 4
    assert line_tag_dict[192] == 4

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert face_media_dict[0] == 4  # PEC surface
    assert face_media_dict[304] == 4  # SGBC surface
    assert face_media_dict[305] == 4  # SGBC surface

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 8  # PEC line
    assert line_media_dict[3.5] == 8  # SGBC line

def test_three_surfaces_Jprobe(tmp_path):
    fn = CASES_FOLDER + 'observation/three_surfaces_Jprobe.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getCurrentVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert face_tag_dict[0] == 725

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 12
    assert line_tag_dict[128] == 10
    assert line_tag_dict[192] == 8

def test_three_one_cell_surfaces_Jprobe(tmp_path):
    fn = CASES_FOLDER + 'observation/three_one_cell_surfaces_Jprobe.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getCurrentVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert face_tag_dict[0] == 728

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 4
    assert line_tag_dict[128] == 3
    assert line_tag_dict[192] == 2

def test_one_cell_PEC_surface_Jprobe(tmp_path):
    fn = CASES_FOLDER + 'observation/one_cell_surface_Jprobe.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)

    solver['general']['numberOfSteps'] = 1
    solver['materialAssociations'][0]['materialId'] = 1

    expected_face_tags = {
        "x": 729,
        "y": 729,
        "z": 728
    }

    for x in ["x", "y", "z"]:
        solver["probes"][0]["component"] = x
        solver.cleanUp()
        solver.run()
        vtkmapfile = solver.getCurrentVTKMap()
        assert os.path.isfile(vtkmapfile)

        face_tag_dict = createPropertyDictionary(vtkmapfile, celltype=9, property='tagnumber')
        assert face_tag_dict[0] == expected_face_tags[x]

        line_tag_dict = createPropertyDictionary(vtkmapfile, celltype=3, property='tagnumber')
        assert line_tag_dict[64] == 4

def test_one_cell_SGBC_surface_Jprobe(tmp_path):
    fn = CASES_FOLDER + 'observation/one_cell_surface_Jprobe.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)

    solver['general']['numberOfSteps'] = 1
    solver['materialAssociations'][0]['materialId'] = 2
    

    solver["probes"][0]["component"] = "x"
    solver.cleanUp()
    solver.run()
    vtkmapfile = solver.getCurrentVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(vtkmapfile, celltype=9, property='tagnumber')
    assert face_tag_dict[0] == 729

    line_tag_dict = createPropertyDictionary(vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 4

    solver["probes"][0]["component"] = "y"
    solver.cleanUp()
    solver.run()
    vtkmapfile = solver.getCurrentVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(vtkmapfile, celltype=9, property='tagnumber')
    assert face_tag_dict[0] == 729

    line_tag_dict = createPropertyDictionary(vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 4

    solver["probes"][0]["component"] = "z"
    solver.cleanUp()
    solver.run()
    vtkmapfile = solver.getCurrentVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(vtkmapfile, celltype=9, property='tagnumber')
    assert face_tag_dict[0] == 728

    line_tag_dict = createPropertyDictionary(vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 4




def test_1_volume(tmp_path):
    fn = CASES_FOLDER + 'observation/pec_volume.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert face_tag_dict[64] == 36

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert len(line_tag_dict) == 0

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert face_media_dict[0] == 36  # PEC surface

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert len(line_media_dict) == 0


def test_2_volumes(tmp_path):
    fn = CASES_FOLDER + 'observation/pec_volumes.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert face_tag_dict[64] == 36
    assert face_tag_dict[128] == 36

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert len(line_tag_dict) == 0

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert face_media_dict[0] == 72  # PEC surface

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert len(line_media_dict) == 0


def test_1_line(tmp_path):
    fn = CASES_FOLDER + 'observation/pec_line.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line


def test_volume_and_surfaces(tmp_path):
    fn = CASES_FOLDER + 'observation/volume_and_surfaces.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert face_tag_dict[64] == 6
    assert face_tag_dict[128] == 1
    assert face_tag_dict[192] == 1

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 1
    assert line_tag_dict[128] == 4
    assert line_tag_dict[192] == 3

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert face_media_dict[-1] == 1  # PEC surface
    assert face_media_dict[0] == 6  # PEC surface
    assert face_media_dict[305] == 1  # SGBC surface

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[-0.5] == 4  # PMC line
    assert line_media_dict[0.5] == 1  # PEC line
    assert line_media_dict[3.5] == 3  # SGBC line


def test_count_bug(tmp_path):
    fn = CASES_FOLDER + 'observation/count_bug.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1
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
    fn = CASES_FOLDER + 'observation/wires.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 4
    assert line_tag_dict[128] == 6
    assert line_tag_dict[192] == 4

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[7] == 7  
    assert line_media_dict[10] == 6 
    assert line_media_dict[21] == 1 

@mtln_skip
def test_wires_collision_Jprobe(tmp_path):
    fn = CASES_FOLDER + 'observation/wires_collision_Jprobe.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getCurrentVTKMap()
    
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 1
    assert face_tag_dict[0] == 729

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert len(line_tag_dict) == 4
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Wire1
    assert line_tag_dict[192] == 6 #Wire2
    assert line_tag_dict[256] == 4 #Wire3


@mtln_skip
def test_wires_collision(tmp_path):
    fn = CASES_FOLDER + 'observation/wires_collision.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2
    assert line_tag_dict[128] == 4
    assert line_tag_dict[192] == 6
    assert line_tag_dict[256] == 4

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[7] == 6  
    assert line_media_dict[8] == 1  
    assert line_media_dict[10] == 6 
    assert line_media_dict[21] == 1 
    assert line_media_dict[0.5] == 2  # PEC line

@mtln_skip
def test_wire_x_collision_y(tmp_path):
    fn = CASES_FOLDER + 'observation/wire_x_collision_y.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  #Wire w/o collision
    assert line_media_dict[8] == 1  #Wire touching non wire
    assert line_media_dict[10] == 2 #Wire extreme

@mtln_skip
def test_wire_x_collision_y_Jprobe(tmp_path):
    fn = CASES_FOLDER + 'observation/wire_x_collision_y_Jprobe.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getCurrentVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert line_tag_dict[0] ==  729

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Wire



@mtln_skip
def test_wire_x_collision_z(tmp_path):
    fn = CASES_FOLDER + 'observation/wire_x_collision_z.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  #Wire w/o collision
    assert line_media_dict[8] == 1  #Wire touching non wire
    assert line_media_dict[10] == 2 #Wire extreme

@mtln_skip
def test_wire_x_long_collision_z(tmp_path):
    fn = CASES_FOLDER + 'observation/wire_x_long_collision_z.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 8 #Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 3  #Wire w/o collision
    assert line_media_dict[8] == 1  #Wire touching non wire
    assert line_media_dict[10] == 2 #Wire extreme
    assert line_media_dict[21] == 2 #Wire extreme

@mtln_skip
def test_wire_y_collision_x(tmp_path):
    fn = CASES_FOLDER + 'observation/wire_y_collision_x.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  #Wire w/o collision
    assert line_media_dict[8] == 1  #Wire touching non wire
    assert line_media_dict[10] == 2 #Wire extreme

@mtln_skip
def test_wire_y_collision_z(tmp_path):
    fn = CASES_FOLDER + 'observation/wire_y_collision_z.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  #Wire w/o collision
    assert line_media_dict[8] == 1  #Wire touching non wire
    assert line_media_dict[10] == 2 #Wire extreme

@mtln_skip
def test_wire_z_collision_x(tmp_path):
    fn = CASES_FOLDER + 'observation/wire_z_collision_x.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  #Wire w/o collision
    assert line_media_dict[8] == 1  #Wire touching non wire
    assert line_media_dict[10] == 2 #Wire extreme

@mtln_skip
def test_wire_z_collision_y(tmp_path):
    fn = CASES_FOLDER + 'observation/wire_z_collision_y.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Wire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[7] == 1  #Wire w/o collision
    assert line_media_dict[8] == 1  #Wire touching non wire
    assert line_media_dict[10] == 2 #Wire extreme

@no_mtln_skip
def test_multiwire_z_collision_y(tmp_path):
    fn = CASES_FOLDER + 'observation/multiwire_z_collision_y.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  #MultiWire w/o collision
    assert line_media_dict[13] == 1  #Multiwire touching non multiwire
    assert line_media_dict[14] == 2 #Multiwire extreme

@no_mtln_skip
def test_multiwire_z_collision_x(tmp_path):
    fn = CASES_FOLDER + 'observation/multiwire_z_collision_x.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  #MultiWire w/o collision
    assert line_media_dict[13] == 1  #Multiwire touching non multiwire
    assert line_media_dict[14] == 2 #Multiwire extreme

@no_mtln_skip
def test_multiwire_y_collision_x(tmp_path):
    fn = CASES_FOLDER + 'observation/multiwire_y_collision_x.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  #MultiWire w/o collision
    assert line_media_dict[13] == 1  #Multiwire touching non multiwire
    assert line_media_dict[14] == 2 #Multiwire extreme

@no_mtln_skip
def test_multiwire_y_collision_z(tmp_path):
    fn = CASES_FOLDER + 'observation/multiwire_y_collision_z.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  #MultiWire w/o collision
    assert line_media_dict[13] == 1  #Multiwire touching non multiwire
    assert line_media_dict[14] == 2 #Multiwire extreme

@no_mtln_skip
def test_multiwire_x_collision_y(tmp_path):
    fn = CASES_FOLDER + 'observation/multiwire_x_collision_y.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  #MultiWire w/o collision
    assert line_media_dict[13] == 1  #Multiwire touching non multiwire
    assert line_media_dict[14] == 2 #Multiwire extreme

@no_mtln_skip
def test_multiwire_x_collision_z(tmp_path):
    fn = CASES_FOLDER + 'observation/multiwire_x_collision_z.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 4 #Multiwire

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 1  #MultiWire w/o collision
    assert line_media_dict[13] == 1  #Multiwire touching non multiwire
    assert line_media_dict[14] == 2 #Multiwire extreme

@no_mtln_skip
def test_multiwire_x_long_collision_z(tmp_path):
    fn = CASES_FOLDER + 'observation/multiwire_x_long_collision_z.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')

    face_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='tagnumber')
    assert len(face_tag_dict) == 0

    line_tag_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='tagnumber')
    assert line_tag_dict[64] == 2 #PEC
    assert line_tag_dict[128] == 8 #Multiwire: 
    # 3 segments w/o collision, 2 of them are also intermediate
    # 1 segment adjacent something not multiwire, 2 extremes, 

    face_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=9, property='mediatype')
    assert len(face_media_dict) == 0

    line_media_dict = createPropertyDictionary(
        vtkmapfile, celltype=3, property='mediatype')
    assert line_media_dict[0.5] == 2  # PEC line
    assert line_media_dict[12] == 3  #MultiWire w/o collision
    assert line_media_dict[13] == 1  #Multiwire touching non multiwire
    assert line_media_dict[14] == 2  #Multiwire extreme
    assert line_media_dict[61] == 2  #Intermediate multiwire segment

def test_can_assign_same_surface_impedance_to_multiple_geometries(tmp_path):
    fn = CASES_FOLDER + 'multipleAssigments/multipleSurfaceImpedance.fdtd.json'

    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    assert (Probe(solver.getSolvedProbeFilenames("BulkProbeEntry")[0]) is not None)

def test_can_assign_same_dielectric_material_to_multiple_geometries(tmp_path):
    fn = CASES_FOLDER + 'multipleAssigments/multipleDielectricMaterial.fdtd.json'

    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    assert (Probe(solver.getSolvedProbeFilenames("BulkProbeEntry")[0]) is not None)

def test_can_execute_fdtd_from_folder_with_spaces_and_can_process_additional_arguments(tmp_path):
    projectRoot = os.getcwd()
    folderWitSpaces: str  = os.path.join(tmp_path, "spaced bin")
    os.mkdir(folderWitSpaces)
    if platform == 'win32':
        shutil.copy2(NGSPICE_DLL, folderWitSpaces)
 
    sembaExecutable = SEMBA_EXE.split(os.path.sep)[-1]
    pathToExe: str = os.path.join(folderWitSpaces, sembaExecutable)
    shutil.copy2(SEMBA_EXE, pathToExe)
    print(pathToExe)
    
    fn = CASES_FOLDER + "dielectric/dielectricTransmission.fdtd.json"
    solver = FDTD(fn, path_to_exe=pathToExe, run_in_folder=tmp_path)
    solver.run()
    assert (Probe(solver.getSolvedProbeFilenames("outside")[0]) is not None)
    assert (solver.getVTKMap()[0] is not None)

def test_bulk_current_outputs(tmp_path):
    # This test uses bulk_probe_cases_over_nodal_source.fdtd from input_examples as input.
    # Verifies all kind of bulk probes are recognised and setted properly by checking outputFile format.
    fn = PROBES_INPUT_EXAMPLE + 'bulk_probe_cases_over_nodal_source.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    bulkXPlaneFiles = solver.getSolvedProbeFilenames("BulkXPlane") 
    bulkYPlaneFiles = solver.getSolvedProbeFilenames("BulkYPlane") 
    bulkZPlaneFiles = solver.getSolvedProbeFilenames("BulkZPlane") 
    bulkYPointFiles = solver.getSolvedProbeFilenames("BulkYPoint") 
    bulkZVolumeFiles = solver.getSolvedProbeFilenames("BulkZVolume") 

    assert len(bulkXPlaneFiles) == 1
    assert len(bulkYPlaneFiles) == 1
    assert len(bulkZPlaneFiles) == 1
    assert len(bulkYPointFiles) == 1
    assert len(bulkZVolumeFiles) == 10

    probeBulkXPlane = Probe(bulkXPlaneFiles[0])
    probeBulkYPlane = Probe(bulkYPlaneFiles[0])
    probeBulkZPlane = Probe(bulkZPlaneFiles[0])
    probeBulkYPoint = Probe(bulkYPointFiles[0])
    probeBulkZVolume = Probe(bulkZVolumeFiles[0])

    assert probeBulkXPlane.direction == 'x'
    assert probeBulkYPlane.direction == 'y'
    assert probeBulkZPlane.direction == 'z'
    assert probeBulkYPoint.direction == 'y'

