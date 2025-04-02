from utils import *
from typing import Dict
import os
from sys import platform

@no_mtln_skip
@pytest.mark.mtln
def test_shieldedPair(tmp_path):
    fn = CASES_FOLDER + 'shieldedPair/shieldedPair.fdtd.json'
    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'],
                  run_in_folder=tmp_path)
    solver.run()

    probe_files = ['shieldedPair.fdtd_wire_start_bundle_line_0_V_75_74_74.dat',
                   'shieldedPair.fdtd_wire_start_bundle_line_0_I_75_74_74.dat',
                   'shieldedPair.fdtd_wire_start_Wz_75_74_74_s4.dat',
                   'shieldedPair.fdtd_wire_end_Wz_75_71_74_s1.dat',
                   'shieldedPair.fdtd_wire_end_bundle_line_0_I_75_71_74.dat',
                   'shieldedPair.fdtd_wire_end_bundle_line_0_V_75_71_74.dat']

    p_expected = []
    for pf in probe_files:
        p_expected.append(Probe(OUTPUTS_FOLDER+pf))

    for i in [2, 3]:
        p_solved = Probe(probe_files[i])
        assert np.allclose(p_expected[i].data.to_numpy()[:, 0:3], p_solved.data.to_numpy()[
                           :, 0:3], rtol=5e-2, atol=0.2)
    for i in [0, 1, 4, 5]:
        p_solved = Probe(probe_files[i])
        assert np.allclose(p_expected[i].data.to_numpy()[:, 0:4], p_solved.data.to_numpy()[
                           :, 0:4], rtol=5e-2, atol=0.2)


@no_mtln_skip
@pytest.mark.mtln
def test_coated_antenna(tmp_path):
    fn = CASES_FOLDER + 'coated_antenna/coated_antenna.fdtd.json'

    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        flags=['-mtlnwires'],
        run_in_folder=tmp_path)
    solver.run()

    probe_current = solver.getSolvedProbeFilenames("mid_point_Wz")[0]
    probe_files = [probe_current]

    p_expected = Probe(
        OUTPUTS_FOLDER+'coated_antenna.fdtd_mid_point_Wz_11_11_11_s2.dat')

    p_solved = Probe(probe_files[0])
    assert np.allclose(
        p_expected.data.to_numpy()[:, 0],
        p_solved.data.to_numpy()[:, 0],
        rtol=0.0, atol=10e-8)
    assert np.allclose(
        p_expected.data.to_numpy()[:, 1],
        p_solved.data.to_numpy()[:, 1],
        rtol=0.0, atol=10e-8)
    assert np.allclose(
        p_expected.data.to_numpy()[:, 2],
        p_solved.data.to_numpy()[:, 2],
        rtol=0.0, atol=10e-6)


def test_holland(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(input_filename=fn, 
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)

    solver.run()

    probe_current = solver.getSolvedProbeFilenames("mid_point_Wz")[0]
    probe_files = [probe_current]

    p_expected = Probe(
        OUTPUTS_FOLDER+'holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat')

    p_solved = Probe(probe_files[0])
    assert np.allclose(
        p_expected.data.to_numpy()[:, 0:3], 
        p_solved.data.to_numpy()[:, 0:3], 
        rtol=1e-5, atol=1e-6)


@no_mtln_skip
@pytest.mark.mtln
def test_holland_mtln(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'], run_in_folder=tmp_path)

    solver.run()

    probe_current = solver.getSolvedProbeFilenames("mid_point_Wz")[0]
    probe_files = [probe_current]

    p_expected = Probe(
        OUTPUTS_FOLDER+'holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat')

    p_solved = Probe(probe_files[0])
    assert np.allclose(
        p_expected.data.to_numpy()[:, 0:3], 
        p_solved.data.to_numpy()[:, 0:3], 
        rtol=1e-5, atol=1e-6)

@no_mtln_skip
@pytest.mark.mtln
def test_unshielded_multiwires(tmp_path):
    fn = CASES_FOLDER + 'unshielded_multiwires/unshielded_multiwires.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'], run_in_folder=tmp_path)

    solver.run()

    probe_names = solver.getSolvedProbeFilenames("mid_point")
    p_solved = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

    p_expected = Probe(
        OUTPUTS_FOLDER+'unshielded_multiwires.fdtd_mid_point_bundle_single_inner_wire_passthrough_I_2_11_14.dat')

    assert np.allclose(
        p_expected.data.to_numpy()[:, 0:3], 
        p_solved.data.to_numpy()[:, 0:3], 
        rtol=1e-5, atol=1e-6)


def test_towelHanger(tmp_path):
    fn = CASES_FOLDER + 'towelHanger/towelHanger.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver.run()

    probe_files = [solver.getSolvedProbeFilenames("wire_start")[0],
                   solver.getSolvedProbeFilenames("wire_mid")[0],
                   solver.getSolvedProbeFilenames("wire_end")[0]]

    p_expected = [Probe(OUTPUTS_FOLDER+'towelHanger.fdtd_wire_start_Wz_27_25_30_s1.dat'),
                  Probe(OUTPUTS_FOLDER +
                        'towelHanger.fdtd_wire_mid_Wx_35_25_32_s5.dat'),
                  Probe(OUTPUTS_FOLDER+'towelHanger.fdtd_wire_end_Wz_43_25_30_s4.dat')]

    for i in range(3):
        p_solved = Probe(probe_files[i])
        assert np.allclose(p_expected[i].data.to_numpy()[:, 0:3], p_solved.data.to_numpy()[
                           :, 0:3], rtol=5e-2, atol=5e-2)


@no_hdf_skip
@pytest.mark.hdf
def test_sphere(tmp_path):
    fn = CASES_FOLDER + 'sphere/sphere.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)
    solver['general']['numberOfSteps'] = 20
    solver['probes'][0]['domain']['initialFrequency'] = 1e8
    solver['probes'][0]['domain']['finalFrequency'] = 1e9

    solver.run()

    # semba-fdtd seems to always use the name Far for "far field" probes.
    far_field_probe_files = solver.getSolvedProbeFilenames("Far")
    assert len(far_field_probe_files) == 1
    p = Probe(far_field_probe_files[0])
    assert p.type == 'farField'

    electric_field_movie_files = solver.getSolvedProbeFilenames(
        "electric_field_movie")
    assert len(electric_field_movie_files) == 3
    p = Probe(electric_field_movie_files[0])
    assert p.type == 'movie'


@no_hdf_skip
@pytest.mark.hdf
def test_movie_in_planewave_in_box(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-in-box-with-movie.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    movie_files = solver.getSolvedProbeFilenames("electric_field_movie")
    h5file = [f for f in movie_files if f.endswith('.h5')][0]
    with h5py.File(h5file, "r") as f:
        time_key = list(f.keys())[0]
        field_key = list(f.keys())[1]
        time_ds = f[time_key][()]
        field_ds = f[field_key][()]

    assert np.isclose(np.max(field_ds), 1.0, rtol=1e-2)
    assert np.min(field_ds) == 0.0


def test_planewave_in_box(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-in-box.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    before = Probe(solver.getSolvedProbeFilenames("before")[0])
    inbox = Probe(solver.getSolvedProbeFilenames("inbox")[0])
    after = Probe(solver.getSolvedProbeFilenames("after")[0])

    np.allclose(inbox.data['field'], inbox.data['incident'])
    zeros = np.zeros_like(before.data['field'])
    np.allclose(before.data['field'], zeros)
    np.allclose(after.data['field'], zeros)


def test_planewave_with_periodic_boundaries(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-with-periodic.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()
    
    before = Probe(solver.getSolvedProbeFilenames("before")[0])
    inbox = Probe(solver.getSolvedProbeFilenames("inbox")[0])
    after = Probe(solver.getSolvedProbeFilenames("after")[0])

    np.allclose(inbox.data['field'], inbox.data['incident'])
    zeros = np.zeros_like(before.data['field'])
    np.allclose(before.data['field'], zeros)
    np.allclose(after.data['field'], zeros)


def test_sgbc_shielding_effectiveness(tmp_path):
    fn = CASES_FOLDER + 'sgbcShieldingEffectiveness/shieldingEffectiveness.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    # FDTD results
    back = Probe(solver.getSolvedProbeFilenames("back")[0])

    t = back.data['time']
    dt = t[1] - t[0]
    fq = fftfreq(len(t))/dt
    INC = fft(back.data['incident'])
    BACK = fft(back.data['field'])
    S21 = BACK/INC

    fmin = 8e6
    fmax = 1e9
    idx_min = (np.abs(fq - fmin)).argmin()
    idx_max = (np.abs(fq - fmax)).argmin()
    f = fq[idx_min:idx_max]
    fdtd_s21 = S21[idx_min:idx_max]

    # Analytical results
    from skrf.media import Freespace
    from skrf.frequency import Frequency
    import scipy.constants

    freq = Frequency.from_f(f, unit='Hz')
    air = Freespace(freq)

    sigma = 100
    width = 10e-3
    mat_ep_r = (1+sigma/(1j*freq.w*scipy.constants.epsilon_0))
    conductive_material = Freespace(freq, ep_r=mat_ep_r)

    slab = air.thru() ** conductive_material.line(width, unit='m') ** air.thru()

    # For debugging only.
    # plt.figure()
    # plt.plot(f, 20*np.log10(np.abs(fdtd_s21)),'.-', label='FDTD S21')
    # plt.plot(f, 20*np.log10(np.abs(slab.s[:,0,1])),'.-', label='Analytical S21')
    # plt.grid(which='both')
    # plt.xlim(f[0], f[-1])
    # plt.xscale('log')
    # plt.legend()
    # plt.show()

    fdtd_s21_db = 20*np.log10(np.abs(fdtd_s21))
    anal_s21_db = 20*np.log10(np.abs(slab.s[:, 0, 1]))

    assert np.allclose(fdtd_s21_db, anal_s21_db, rtol=0.05)

def test_sgbc_structured_resistance(tmp_path):
    fn = CASES_FOLDER + 'sgbcResistance/sgbcResistance.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    i = Probe(solver.getSolvedProbeFilenames("Bulk probe")[0]).data['current']
    assert np.allclose(i.array[-101:-1], np.ones(100)*i.array[-100], rtol=1e-3)
    assert np.allclose(-1/i.array[-101:-1], np.ones(100)*(50+45), rtol=0.05)


def test_pec_overlapping_sgbcs(tmp_path):
    """ Test that PEC surfaces overlapping SGBC surfaces prioritize PEC.
    """
    fn = CASES_FOLDER + 'sgbcOverlapping/sgbcOverlapping.fdtd.json'

    # Runs case without overlap.
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    p = Probe(solver.getSolvedProbeFilenames("Bulk probe")[0])
    t = p['time'].to_numpy()
    iSGBC = p['current'].to_numpy()

    # Adds current SGBC elements as PEC. Now both are defined over same surface.
    sgbcElementIds = solver["materialAssociations"][1]["elementIds"]
    solver['materialAssociations'][0]["elementIds"].extend(sgbcElementIds)
    solver.cleanUp()
    solver.run()
    iPEC = Probe(solver.getSolvedProbeFilenames("Bulk probe")[0])['current'].to_numpy()

    
    # For debugging only.
    # plt.figure()
    # plt.plot(t, iSGBC,'.-', label='SGBC case')
    # plt.plot(t, iPEC,'.-', label='PEC overlapping')
    # plt.grid(which='both')
    # plt.legend()
    # plt.show()

    
    # Checks values are different due to PEC prioritization.
    assert np.all(np.greater(np.abs(iPEC[1000:]), np.abs(iSGBC[1000:])))

def test_sgbc_overlapping_sgbc(tmp_path):
    """ Test that SGBC surfaces overlapping SGBC surfaces prioritize first in MatAss.
    """
    fn = CASES_FOLDER + 'sgbcOverlapping/sgbcOverlapping.fdtd.json'

    # Runs case without overlap.
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    # Changes materialId in first SGBC in MatAss to material with larger conductivity.
    solver['materialAssociations'][1]["materialId"] = 6
    solver.cleanUp()
    solver.run()
    p = Probe(solver.getSolvedProbeFilenames("Bulk probe")[0])

    t = p['time'].to_numpy()
    iSGBC_top = p['current'].to_numpy()

    # Changes materialId in second SGBC in MatAss to material with larger conductivity.
    solver['materialAssociations'][1]["materialId"] = 2
    solver['materialAssociations'][2]["materialId"] = 6
    solver.cleanUp()
    solver.run()
    iSGBC_bottom = Probe(solver.getSolvedProbeFilenames("Bulk probe")[0])['current'].to_numpy()

    
    # For debugging only.
    # plt.figure()
    # plt.plot(t, iSGBC_top,'.-', label='SGBC sigma = 40 S/m, top')
    # plt.plot(t, iSGBC_bottom,'.-', label='SGBC sigma = 20 S/m, bottom')
    # plt.grid(which='both')
    # plt.legend()
    # plt.show()

    
    # Checks values are different due to prioritization of first written.
    assert np.all(np.greater(np.abs(iSGBC_top[1000:]), np.abs(iSGBC_bottom[1000:])))

def test_dielectric_transmission(tmp_path):
    _FIELD_TOLERANCE = 0.05

    def getIncidentField(probe:Probe) -> Dict:
        idx = probe["field"].argmin()
        time = probe["time"][idx]
        value = probe["field"][idx]
        return {"time":time, "value": value}
    
    def getReflectedField(probe:Probe) -> Dict:
        idx = probe["field"].argmax()
        time = probe["time"][idx]
        value = probe["field"][idx]
        return {"time":time, "value": value}
    
    def getTransmittedField(probe:Probe) -> Dict:
        idx = probe["field"].argmin()
        time = probe["time"][idx]
        value = probe["field"][idx]
        return {"time":time, "value": value}
    
    def getReflectedDelay(incidentTime:float, reflectedTime:float):
        timeToSurface:float = ((reflectedTime-incidentTime)/2) + incidentTime
        reflectedDelay:float = reflectedTime - timeToSurface
        return reflectedDelay
        
    def getTransmittedDelay(incidentTime:float, reflectedTime:float, transmittedTime:float):
        timeToSurface:float = ((reflectedTime-incidentTime)/2) + incidentTime
        transmitedDelay = transmittedTime - timeToSurface
        return transmitedDelay

    fn = CASES_FOLDER + "dielectric/dielectricTransmission.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    relativePermittivity = solver.getMaterialProperties('DielectricMaterial')["relativePermittivity"]
    materialRelativeImpedance = np.sqrt(1/relativePermittivity)
    
    expectedReflectedCoeff = (materialRelativeImpedance - 1) / (materialRelativeImpedance + 1)
    expectedtransmittedCoeff = (1 + expectedReflectedCoeff)
    expectedDelayRatio = 1/np.sqrt(relativePermittivity)

    outsideProbe = Probe(solver.getSolvedProbeFilenames("outside")[0])
    insideProbe = Probe(solver.getSolvedProbeFilenames("inside")[0])


    incidentField = getIncidentField(outsideProbe)
    reflectedField = getReflectedField(outsideProbe)
    transmittedField = getTransmittedField(insideProbe)
 
    assert (incidentField['value'] - transmittedField['value'] + reflectedField['value']) < _FIELD_TOLERANCE
    assert np.allclose(reflectedField["value"]/incidentField["value"], expectedReflectedCoeff, rtol=_FIELD_TOLERANCE)
    assert np.allclose(transmittedField["value"]/incidentField["value"], expectedtransmittedCoeff, rtol=_FIELD_TOLERANCE)

    reflectedDelay:float = getReflectedDelay(incidentField['time'], reflectedField['time'])
    transmitedDelay:float = getTransmittedDelay(incidentField['time'], reflectedField['time'], transmittedField['time'])

    assert np.allclose(reflectedDelay/transmitedDelay, expectedDelayRatio, rtol=_FIELD_TOLERANCE)

    
def test_rectilinear_mode(tmp_path):
    _FIELD_TOLERANCE = 4
    _TIME_TOLERANCE = 4

    def getPeakPulse(probe:Probe) -> Dict:
        idx = probe["field"].argmax()
        time = probe["time"][idx]
        value = probe["field"][idx]
        return {"time":time, "value": value}
    
    rectilinearModeFile = CASES_FOLDER + "rectilinear_mode/rectilinearMode.fdtd.json"
    noRectilinearModeFile = CASES_FOLDER + "rectilinear_mode/noRectilinearMode.fdtd.json"

    rectilinearModeFolder = os.path.join(tmp_path, 'rectilinear')
    noRectilinearModeFolder = os.path.join(tmp_path, 'noRectilinear')

    os.mkdir(rectilinearModeFolder)
    os.mkdir(noRectilinearModeFolder)

    solverRectilinear = FDTD(rectilinearModeFile, path_to_exe=SEMBA_EXE, run_in_folder=rectilinearModeFolder)
    solverRectilinear.run()
    rectilinearFrontProbe = Probe(solverRectilinear.getSolvedProbeFilenames("Front probe")[0])
    rectilinearVertexProbe = Probe(solverRectilinear.getSolvedProbeFilenames("Vertex probe")[0])


    solverNoRectilinear = FDTD(noRectilinearModeFile, path_to_exe=SEMBA_EXE, run_in_folder=noRectilinearModeFolder)
    solverNoRectilinear.run()
    noRectilinearFrontProbe = Probe(solverNoRectilinear.getSolvedProbeFilenames("Front probe")[0])
    noRectilinearVertexProbe = Probe(solverNoRectilinear.getSolvedProbeFilenames("Vertex probe")[0])

    np.testing.assert_almost_equal(getPeakPulse(rectilinearFrontProbe)['value'], getPeakPulse(noRectilinearFrontProbe)['value'], decimal=_FIELD_TOLERANCE)
    np.testing.assert_almost_equal(getPeakPulse(rectilinearFrontProbe)['time'], getPeakPulse(noRectilinearFrontProbe)['time'], decimal=_TIME_TOLERANCE)
    np.testing.assert_almost_equal(getPeakPulse(rectilinearVertexProbe)['value'], getPeakPulse(noRectilinearVertexProbe)['value'], decimal=_FIELD_TOLERANCE)
    np.testing.assert_almost_equal(getPeakPulse(rectilinearVertexProbe)['time'], getPeakPulse(noRectilinearVertexProbe)['time'], decimal=_TIME_TOLERANCE)
    
def testCanExecuteFDTDFromFolderWithSpacesAndCanProcessAdditionalArguments(tmp_path):
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
    