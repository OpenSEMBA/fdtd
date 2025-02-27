from utils import *


def test_holland_case_checking_number_of_outputs(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    number_of_steps = 10
    solver['general']['numberOfSteps'] = number_of_steps

    solver.run()

    probe_files = solver.getSolvedProbeFilenames("mid_point")

    assert len(probe_files) == 1
    assert 'holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat' == probe_files[0]
    assert countLinesInFile(probe_files[0]) == number_of_steps + 2


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

    assert 'towelHanger.fdtd_wire_start_Wz_27_25_30_s1.dat' == probe_start[0]
    assert 'towelHanger.fdtd_wire_mid_Wx_35_25_32_s5.dat' == probe_mid[0]
    assert 'towelHanger.fdtd_wire_end_Wz_43_25_30_s4.dat' == probe_end[0]
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

    p = Probe(solver.getSolvedProbeFilenames("Far")[0])
    assert p.case_name == 'sphere'
    assert p.type == 'farField'
    assert np.all(p.cell_init == np.array([2, 2, 2]))

    p = Probe(solver.getSolvedProbeFilenames("electric_field_movie")[0])
    assert p.case_name == 'sphere'
    assert p.type == 'movie'
    assert np.all(p.cell_init == np.array([2, 2, 2]))


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

@pytest.mark.xfail(reason="feature not implemented")
def test_conformal_delay(tmp_path):
    fn = CASES_FOLDER + 'observation/conformal.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])

    solver.run()    
    source = solver.getSolvedProbeFilenames("source")
    t = source.data['time']
    inc = source.data['incident']
    back = source.data['field']
    #assert delay time corresponds
