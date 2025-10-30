from utils import *
from src_caseMaker.caseMaker import *


def test_read_wire_probe():
    p = Probe(OUTPUTS_FOLDER + 'holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat')

    assert p.case_name == 'holland1981'
    assert p.name == 'mid_point'
    assert p.type == 'wire'
    assert p.domainType == 'time'
    assert np.all(p.cell == np.array([11, 11, 12]))
    assert p.segment == 2

    assert len(p['time']) == 1001
    assert p['time'][0] == 0.0
    assert p['time'].iat[-1] == 0.2999999901276417E-007

    assert len(p['current']) == 1001
    assert p['current'][0] == 0.0
    assert p['current'].iat[-1] == -0.513576742E-004


def test_read_probe_from_NFDE():
    p = Probe(OUTPUTS_FOLDER + 'fakeCurrentProbe_mid_point_Wz_11_11_11_s2.dat')

    assert p.type == 'wire'


def test_read_frequency_probe_from_NFDE():
    p = Probe(OUTPUTS_FOLDER + 'edelcadfixZ_COR2_log__Wz_21_21_28_s10_df.dat')

    assert p.type == 'wire'
    assert p.domainType == 'frequency'
    assert np.all(p.cell == np.array([21, 21, 28]))
    assert p.segment == 10


def test_read_point_probe():
    p = Probe(OUTPUTS_FOLDER + 'shieldingEffectiveness.fdtd_front_Ex_1_1_1.dat')

    assert p.case_name == 'shieldingEffectiveness'
    assert p.name == 'front'
    assert p.type == 'point'
    assert p.domainType == 'time'
    assert p.direction == 'x'
    assert p.field == 'E'
    assert np.all(p.cell == np.array([1, 1, 1]))

    assert len(p['time']) == 5193
    assert p['time'][0] == 0.0
    assert np.isclose(p['time'].iat[-1], 0.19997851853637005E-007)

    assert len(p['field']) == 5193
    assert p['field'][0] == 0.0
    assert np.isclose(p['field'].iat[-1], 0.120145380E+000)

    assert len(p['incident']) == 5193
    assert np.isclose(p['incident'][0], 0.134010895E-005)
    assert p['incident'].iat[-1] == 0.0


def test_read_point_probe_without_planewave():
    p = Probe(OUTPUTS_FOLDER + 'twoWires.fdtd_ProbeEnd_Ey_25_13_5.dat')

    assert p.case_name == 'twoWires'
    assert p.name == 'ProbeEnd'
    assert p.type == 'point'
    assert p.domainType == 'time'
    assert p.direction == 'y'
    assert p.field == 'E'
    assert np.all(p.cell == np.array([25, 13, 5]))


def test_read_bulk_current_probe():
    p = Probe(OUTPUTS_FOLDER +
              'twoWires.fdtd_Bulk probe_Jx_15_11_13__15_13_17.dat')

    assert p.case_name == 'twoWires'
    assert p.name == 'Bulk probe'
    assert p.type == 'bulkCurrent'
    assert p.domainType == 'time'
    assert p.direction == 'x'


def test_fdtd_set_new_folder_to_run(tmp_path):
    input = os.path.join(CASES_FOLDER, 'planewave', 'pw-in-box.fdtd.json')
    solver = FDTD(input, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver['general']['numberOfSteps'] = 1

    solver.run()


def test_fdtd_with_string_args(tmp_path):
    input = os.path.join(CASES_FOLDER, 'planewave', 'pw-in-box.fdtd.json')
    solver = FDTD(input,
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path,
                  flags='-h')
    solver['general']['numberOfSteps'] = 1

    solver.run()


@no_mpi_skip
def test_fdtd_with_mpi_run(tmp_path):
    input = os.path.join(CASES_FOLDER, 'planewave', 'pw-in-box.fdtd.json')
    solver = FDTD(input,
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path,
                  flags=['-h'],
                  mpi_command='mpirun -np 2')
    solver['general']['numberOfSteps'] = 1

    solver.run()


def test_fdtd_clean_up_after_run(tmp_path):
    input = CASES_FOLDER + 'planewave/pw-in-box.fdtd.json'
    solver = FDTD(input, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver['general']['numberOfSteps'] = 1

    solver.run()

    pn = solver.getSolvedProbeFilenames("inbox")
    assert os.path.isfile(pn[0])

    solver.cleanUp()

    assert not os.path.isfile(pn[0])


def test_fdtd_get_used_files():
    fn = CASES_FOLDER + 'multilines_opamp/multilines_opamp.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE)

    used_files = solver.getUsedFiles()

    assert len(used_files) == 2
    assert used_files[0] == 'spice_4port_pulse_start_75.exc'
    assert used_files[1] == 'opamp.model'


def test_casemaker_setGridFromVTK():
    cm = CaseMaker()
    cm.setGridFromVTK(GEOMETRIES_FOLDER + 'sphere.grid.vtp')

    grid = cm['mesh']['grid']
    
    steps = grid['steps']
    assert np.allclose(np.array(steps['x']), 2.0)
    assert np.allclose(np.array(steps['y']), 2.0)
    assert np.allclose(np.array(steps['z']), 2.0)
    assert np.array(steps['x']).size == 100
    assert np.array(steps['y']).size == 100
    assert np.array(steps['z']).size == 100

    assert grid['numberOfCells'] == [100, 100, 100]
    assert grid['origin'] == [-100.0, -100.0, -100.0]


def test_casemaker_sphere_rcs_case(tmp_path):
    os.chdir(tmp_path)
    cm = CaseMaker()

    cm.setNumberOfTimeSteps(1)
    cm.setTimeStep(1e-12)

    cm.setAllBoundaries("mur")

    cm.setGridFromVTK(GEOMETRIES_FOLDER + "sphere.grid.vtp")
    sphereId = cm.addCellElementsFromVTK(
        GEOMETRIES_FOLDER + "buggy_sphere.str.vtp")

    pecId = cm.addPECMaterial()
    cm.addMaterialAssociation(pecId, [sphereId])

    planewaveBoxId = cm.addCellElementBox(
        [[-75.0, -75.0, -75.0], [75.0, 75.0, 75.0]])
    direction = {"theta": np.pi/2, "phi": 0.0}
    polarization = {"theta": np.pi/2, "phi": np.pi/2}
    dt = 1e-12
    w0 = 0.1e-9    # ~ 2 GHz bandwidth
    t0 = 10*w0
    t = np.arange(0, t0+20*w0, dt)
    data = np.empty((len(t), 2))
    data[:,0] = t
    data[:,1] = np.exp( -np.power(t-t0,2)/ w0**2 )
    np.savetxt('gauss.exc', data)
    cm.addPlanewaveSource(planewaveBoxId, 'gauss.exc', direction, polarization)

    pointProbeNodeId = cm.addNodeElement([-65.0, 0.0, 0.0])
    cm.addPointProbe(pointProbeNodeId, name="front")

    n2ffBoxId = cm.addCellElementBox(
        [[-85.0, -85.0, -85.0], [85.0, 85.0, 85.0]])
    theta = {"initial": np.pi/2, "final": np.pi/2, "step": 0.0}
    phi = {"initial": np.pi, "final": np.pi, "step": 0.0}
    domain = {
        "type": "frequency",
        "initialFrequency": 10e6,
        "finalFrequency": 1e9,
        "numberOfFrequencies": 10,
        "frequencySpacing": "logarithmic"
    }
    cm.addFarFieldProbe(n2ffBoxId, "n2ff", theta, phi, domain)

    case_name = 'sphere_rcs'
    cm.exportCase(case_name)

    solver = FDTD(
        case_name + ".fdtd.json",
        path_to_exe=SEMBA_EXE)
    solver.run()
    assert solver.hasFinishedSuccessfully()