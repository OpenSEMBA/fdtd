from utils import *
import pytest


def test_holland_case_checking_number_of_outputs(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    number_of_steps = 10
    solver['general']['numberOfSteps'] = number_of_steps

    solver.run()
    assert solver.hasFinishedSuccessfully()

    probe_files = solver.getSolvedProbeFilenames("mid_point")

    assert solver.hasFinishedSuccessfully() == True
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

    assert solver.hasFinishedSuccessfully() == True
    assert len(probe_start) == 1
    assert len(probe_mid) == 1
    assert len(probe_end) == 1

    assert 'towelHanger.fdtd_wire_start_Wz_27_25_30_s1.dat' == probe_start[0]
    assert 'towelHanger.fdtd_wire_mid_Wx_35_25_32_s5.dat' == probe_mid[0]
    assert 'towelHanger.fdtd_wire_end_Wz_43_25_30_s4.dat' == probe_end[0]
    assert countLinesInFile(probe_start[0]) == 3
    assert countLinesInFile(probe_mid[0]) == 3
    assert countLinesInFile(probe_end[0]) == 3


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


def test_sgbc_with_mapvtk_checking_tagnumbers(tmp_path):
    fn = CASES_FOLDER + 'sgbc/sgbc.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path, flags=['-mapvtk'])
    solver['general']['numberOfSteps'] = 1

    solver.run()
    assert solver.hasFinishedSuccessfully() == True

    vtkmapfile = solver.getVTKMap()
    assert os.path.isfile(vtkmapfile)

    d = createFaceTagDictionary(vtkmapfile)
    assert d[64] == 4
    assert d[128] == 4
    assert d[192] == 4
