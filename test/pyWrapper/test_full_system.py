from utils import *
from typing import Dict
import os
from sys import platform
from scipy import signal

def test_lineIntegralProbe(tmp_path):
    fn = CASES_FOLDER + 'lineIntegralProbe/lineIntegralProbe_plates.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    
    pf = 'lineIntegralProbe_plates.fdtd_vprobe_LI_20_20_10.dat'
    li_probe  = Probe(solver.getSolvedProbeFilenames("vprobe_LI_20_20_10")[0])
    expected  = Probe(OUTPUTS_FOLDER+pf)
    np.allclose(li_probe['lineIntegral'].to_numpy(), expected['lineIntegral'].to_numpy(), rtol =0.01 , atol=0.01)


@no_mtln_skip
@pytest.mark.mtln
def test_shieldedPair(tmp_path):
    fn = CASES_FOLDER + 'shieldedPair/shieldedPair.fdtd.json'
    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'],
                  run_in_folder=tmp_path)
    solver.run()

    probe_files = ['shieldedPair.fdtd_wire_start_bundle_line_out_V_75_74_74.dat',
                   'shieldedPair.fdtd_wire_start_bundle_line_out_I_75_74_74.dat',
                   'shieldedPair.fdtd_wire_end_bundle_line_out_I_75_71_74.dat',
                   'shieldedPair.fdtd_wire_end_bundle_line_out_V_75_71_74.dat']

    p_expected = []
    for pf in probe_files:
        p_expected.append(Probe(OUTPUTS_FOLDER+pf))

    for i in [0, 1, 2, 3]:
        p_solved = Probe(probe_files[i])
        assert np.allclose(p_expected[i].data.to_numpy()[:, 0:4], p_solved.data.to_numpy()[
                           :, 0:4], rtol=5e-2, atol=0.2)

@no_mtln_skip
@no_mpi_skip
@pytest.mark.mtln
@pytest.mark.mpi
def test_shieldedPair_mpi(tmp_path):
    fn = CASES_FOLDER + 'shieldedPair/shieldedPair.fdtd.json'
    solver = FDTD(input_filename=fn,
                  path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'],
                  mpi_command='mpirun -np 2',
                  run_in_folder=tmp_path)
    solver.run()

    probe_files = ['shieldedPair.fdtd_wire_start_bundle_line_out_V_75_74_74.dat',
                   'shieldedPair.fdtd_wire_start_bundle_line_out_I_75_74_74.dat',
                   'shieldedPair.fdtd_wire_end_bundle_line_out_I_75_71_74.dat',
                   'shieldedPair.fdtd_wire_end_bundle_line_out_V_75_71_74.dat']

    p_expected = []
    for pf in probe_files:
        p_expected.append(Probe(OUTPUTS_FOLDER+pf))

    for i in [0, 1, 2, 3]:
        p_solved = Probe(probe_files[i])
        assert np.allclose(p_expected[i].data.to_numpy()[:, 0:4], p_solved.data.to_numpy()[
                           :, 0:4], rtol=5e-2, atol=0.2)


@no_mtln_skip
@pytest.mark.mtln
def test_coated_antenna(tmp_path):
    """ Test for a coated antenna with MTLN wires reproducing Fig. 2 in:
        A. Rubio Bretones, R. Gomez Martin, A. Salinas and I. Sanchez, 
        "Time domain analysis of dielectric coated wire scatterers and antennas," 
        Proceedings of MELECON '94. Mediterranean Electrotechnical Conference,
        Antalya, Turkey, 1994, pp. 1174-1176 vol.3, doi: 10.1109/MELCON.1994.380859.
    """
    fn = CASES_FOLDER + 'coated_antenna/coated_antenna.fdtd.json'

    solver = FDTD(
        input_filename=fn,
        path_to_exe=SEMBA_EXE,
        flags=['-mtlnwires'],
        run_in_folder=tmp_path)
    solver.run()

    probe_current = solver.getSolvedProbeFilenames("mid_point")[0]
    probe_files = [probe_current]

    p_expected = Probe(
        OUTPUTS_FOLDER+'coated_antenna.fdtd_mid_point_bundle_half_1_I_11_11_12.dat')

    p_solved = Probe(probe_files[0])
    assert np.allclose(
        p_expected['time'].to_numpy(),
        p_solved['time'].to_numpy(),
        rtol=0.0, atol=10e-8)
    assert np.allclose(
        p_expected['current_0'].to_numpy(),
        p_solved['current_0'].to_numpy(),
        rtol=0.0, atol=10e-8)


def test_holland(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981.fdtd.json'
    solver = FDTD(input_filename=fn, 
                  path_to_exe=SEMBA_EXE,
                  run_in_folder=tmp_path)

    solver.run()

    probe_current = solver.getSolvedProbeFilenames("mid_point_Wz")[0]
    probe_files = [probe_current]
    p_solved = Probe(probe_files[0])

    expected_f = json.load(open(OUTPUTS_FOLDER+'holland1981.fdtd_mid_point_Wz_11_11_12_s2_12.json'))
    expected_t, expected_i = np.array([]), np.array([])
    for data in expected_f['datasetColl'][0]['data']:
        expected_t = np.append(expected_t, float(data['value'][0]))    
        expected_i = np.append(expected_i, float(data['value'][1]))    

    expected_i_interp = np.interp(p_solved['time']-3.05*1e-9, expected_t, expected_i)

    assert np.allclose(
        expected_i_interp, 
        -p_solved['current'], 
        rtol=1e-4, atol=5e-5)


@no_mtln_skip
@pytest.mark.mtln
def test_holland_mtln(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981_unshielded.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'], run_in_folder=tmp_path)

    solver.run()

    probe_current = solver.getSolvedProbeFilenames("mid_point_bundle_single_unshielded_multiwire_I")[0]
    probe_files = [probe_current]
    p_solved = Probe(probe_files[0])

    expected_f = json.load(open(OUTPUTS_FOLDER+'holland1981.fdtd_mid_point_Wz_11_11_12_s2_12.json'))
    expected_t, expected_i = np.array([]), np.array([])
    for data in expected_f['datasetColl'][0]['data']:
        expected_t = np.append(expected_t, float(data['value'][0]))    
        expected_i = np.append(expected_i, float(data['value'][1]))    

    expected_i_interp = np.interp(p_solved['time']-3.05*1e-9, expected_t, expected_i)

    assert np.allclose(
        expected_i_interp, 
        p_solved['current_0'], 
        rtol=1e-4, atol=5e-5)

@no_mtln_skip
@no_mpi_skip
@pytest.mark.mtln
@pytest.mark.mpi
def test_holland_mtln_mpi(tmp_path):
    fn = CASES_FOLDER + 'holland/holland1981_unshielded.fdtd.json'
    # no mpi
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'] ,run_in_folder=tmp_path)
    solver.run()
    probe_names = solver.getSolvedProbeFilenames("mid_point")
    probe_mid_no_mpi = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

    # no mpi -np 1
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'], mpi_command='mpirun -np 1',run_in_folder=tmp_path)
    solver.cleanUp()
    solver.run()
    probe_mid_mpi_1 = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

    # no mpi -np 2
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'], mpi_command='mpirun -np 2',run_in_folder=tmp_path)
    solver.cleanUp()
    solver.run()
    probe_mid_mpi_2 = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

    # # no mpi -np 3
    # solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
    #               flags=['-mtlnwires'], mpi_command='mpirun -np 3',run_in_folder=tmp_path)
    # solver.cleanUp()
    # solver.run()
    # probe_mid_mpi_3 = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

    expected_f = json.load(open(OUTPUTS_FOLDER+'holland1981.fdtd_mid_point_Wz_11_11_12_s2_12.json'))
    expected_t, expected_i = np.array([]), np.array([])
    for data in expected_f['datasetColl'][0]['data']:
        expected_t = np.append(expected_t, float(data['value'][0]))    
        expected_i = np.append(expected_i, float(data['value'][1]))    

    expected_i_interp = np.interp(probe_mid_no_mpi['time']-3.05*1e-9, expected_t, expected_i)

    p_expected = Probe(
        OUTPUTS_FOLDER+'holland1981.fdtd_mid_point_Wz_11_11_12_s2.dat')

    assert np.allclose(
        expected_i_interp, 
        probe_mid_no_mpi['current_0'], 
        rtol=1e-4, atol=5e-5)

    assert np.allclose(
        expected_i_interp, 
        probe_mid_mpi_1['current_0'], 
        rtol=1e-4, atol=5e-5)

    assert np.allclose(
        expected_i_interp, 
        probe_mid_mpi_2['current_0'], 
        rtol=1e-4, atol=5e-5)

    # assert np.allclose(
    #     expected_i_interp, 
    #     probe_mid_mpi_3['current_0'], 
    #     rtol=1e-4, atol=5e-5)


@no_mtln_skip
@pytest.mark.mtln
def test_unshielded_multiwires(tmp_path):
    fn = CASES_FOLDER + 'unshielded_multiwires/unshielded_multiwires_berenger.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE,
                  flags=['-mtlnwires'], run_in_folder=tmp_path)

    solver.run()

    probe_names = solver.getSolvedProbeFilenames("mid_point")
    p_solved = Probe(list(filter(lambda x: '_I_' in x, probe_names))[0])

    p_expected = Probe(
        OUTPUTS_FOLDER+'unshielded_multiwires_berenger.fdtd_mid_point_bundle_unshielded_two_wire_I_2_11_14.dat')

    assert np.allclose(
        p_expected['current_0'], 
        p_solved['current_0'], 
        rtol=1e-3, atol=0.01)
    assert np.allclose(
        p_expected['current_1'], 
        p_solved['current_1'], 
        rtol=1e-3, atol=0.01)


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

    far_field_probe_files = solver.getSolvedProbeFilenames("far")
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

    assert np.corrcoef(inbox.data['field'].to_numpy(), inbox.data['incident'].to_numpy())[0, 1] > 0.999
    zeros = np.zeros_like(before.data['field'])
    assert np.allclose(before.data['field'].to_numpy(), zeros, atol=5e-4)
    assert np.allclose(after.data['field'].to_numpy(), zeros, atol=5e-4)


def test_planewave_with_periodic_boundaries(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-with-periodic.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()
    
    before = Probe(solver.getSolvedProbeFilenames("before")[0])
    inbox = Probe(solver.getSolvedProbeFilenames("inbox")[0])
    after = Probe(solver.getSolvedProbeFilenames("after")[0])

    assert np.corrcoef(inbox.data['field'].to_numpy(), inbox.data['incident'].to_numpy())[0, 1] > 0.999
    zeros = np.zeros_like(before.data['field'])
    assert np.allclose(before.data['field'].to_numpy(), zeros, atol=1.5e-3)
    assert np.allclose(after.data['field'].to_numpy(), zeros, atol=1.5e-3)


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


def test_nodal_source(tmp_path):
    fn = CASES_FOLDER + "nodalSource/nodalSource.fdtd.json"
    assert (os.path.isfile(fn))
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    
    resistanceBulkProbe = Probe( \
    solver.getSolvedProbeFilenames("Bulk probe Resistance")[0])
    nodalBulkProbe = Probe( \
        solver.getSolvedProbeFilenames("Bulk probe Nodal Source")[0])
    excitation = ExcitationFile( \
        excitation_filename=solver.getExcitationFile("predefinedExcitation")[0])

    # For debugging.
    # plt.figure()
    # plt.plot(resistanceBulkProbe['time'].to_numpy(), 
    #         resistanceBulkProbe['current'].to_numpy(), label='BP Current@resistance')
    # plt.plot(excitation.data['time'].to_numpy(), 
    #         excitation.data['value'].to_numpy(), label='excited current')
    # plt.plot(nodalBulkProbe['time'].to_numpy(),
    #         -nodalBulkProbe['current'].to_numpy(), label='BP Current@nodal source')
    # plt.legend()

    exc = np.interp(nodalBulkProbe['time'].to_numpy(), 
                    excitation.data['time'].to_numpy(), 
                    excitation.data['value'].to_numpy())
    assert np.corrcoef(exc, -nodalBulkProbe['current'])[0,1] > 0.999
    assert np.corrcoef(-nodalBulkProbe['current'], resistanceBulkProbe['current'])[0,1] > 0.998


def testCanAssignSameSurfaceImpedanceToMultipleGeometries(tmp_path):
    fn = CASES_FOLDER + 'multipleAssigments/multipleSurfaceImpedance.fdtd.json'

    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    assert (Probe(solver.getSolvedProbeFilenames("BulkProbeEntry")[0]) is not None)

def testCanAssignSameDielectricMaterialToMultipleGeometries(tmp_path):
    fn = CASES_FOLDER + 'multipleAssigments/multipleDielectricMaterial.fdtd.json'

    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    assert (Probe(solver.getSolvedProbeFilenames("BulkProbeEntry")[0]) is not None)

def test_lumped_resistor(tmp_path):
    # This test validates the behavior of lumped resistor materials in a simplified circuit.
    # The circuit consists of a 40mm x 40mm simple loop with a lumped resistor line inserted along one edge.
    # Due to the geometry, the circuit naturally exhibits a parasitic inductance of approximately 1.65e-7 H.
    #
    # Current probes are placed at the start and end of the lumped line, as well as a few cells before and after it.
    # These measurements are used to evaluate the accuracy of the lumped material model.
    #
    # For validation, the results are compared against two reference cases:
    # 1. A simple loop circuit with the same dimensions where a terminal  is inserted in place of the lumped line, 
    #    using the same resistance of the lumped line.
    # 2. Theoretical current response calculated using Laplace transforms from the initial pulse excitation.
    #
    # For better interaction with the case, the user can go to the file: testData/cases/lumped_lines/simple_loop_R/simple_loop_prepost.py

    fn_lumped = CASES_FOLDER + 'lumped_lines/simple_loop_R/simple_loop_lumped.fdtd.json'
    fn_terminal = CASES_FOLDER + 'lumped_lines/simple_loop_R/simple_loop_terminal.fdtd.json'
    
    solver_lumped = FDTD(fn_lumped, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver_terminal = FDTD(fn_terminal, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver_lumped.run()
    solver_terminal.run()

    StartTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("TerminalCellStart")[0])
    StartLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellStart")[0])

    EndTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("TerminalCellEnd")[0])
    EndLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellEnd")[0])

    AdjacentPostLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PostLumpedCell")[0])
    AdjacentPostTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("PostTerminalCell")[0])

    AdjacentPreLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PreLumpedCell")[0])
    AdjacentPreTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("PreTerminalCell")[0])

    assert np.corrcoef(StartLumpedProbe['current'].to_numpy(), StartTerminalProbe['current'].to_numpy())[0, 1] > 0.999
    assert np.corrcoef(EndLumpedProbe['current'].to_numpy(), EndTerminalProbe['current'].to_numpy())[0, 1] > 0.999
    assert np.corrcoef(AdjacentPostLumpedProbe['current'].to_numpy(), AdjacentPostTerminalProbe['current'].to_numpy())[0, 1] > 0.999
    assert np.corrcoef(AdjacentPreLumpedProbe['current'].to_numpy(), AdjacentPreTerminalProbe['current'].to_numpy())[0, 1] > 0.999

    R = solver_lumped.getMaterialProperties("lumped_line")["resistance"]
    L = 1.65e-7     # parasitic inductance mentioned above

    num = [1]
    den = [L, R]
    system = signal.TransferFunction(num, den)
    tout, I_out, _ = signal.lsim(system, 
                                 U=np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=1), 
                                 T=np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=0))
    
    I_theo = np.interp(AdjacentPreLumpedProbe['time'], tout, I_out)
    
    assert np.corrcoef(AdjacentPostLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(AdjacentPreLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(StartLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(EndLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999

def test_lumped_capacitor(tmp_path):
    # This test validates the behavior of lumped capacitor materials in a simplified circuit. The lumped capacitor 
    # can be modeled as a capacitor in parallel with a resistor.
    # The circuit consists of a 40mm x 40mm simple loop with a lumped capacitor line inserted along one edge.
    # Due to the geometry, the circuit naturally exhibits a parasitic inductance of approximately 1.65e-7 H.
    #
    # Current probes are placed at the start and end of the lumped line, as well as a few cells before and after it.
    # These measurements are used to evaluate the accuracy of the lumped material model.
    #
    # For validation, the results are compared against only one reference case:
    # 1. Theoretical current response calculated using Laplace transforms from the initial pulse excitation.
    #
    # For better interaction with the case, the user can go to the file: testData/cases/lumped_lines/simple_loop_RC/simple_loop_prepost.py

    fn_lumped = CASES_FOLDER + 'lumped_lines/simple_loop_RC/simple_loop_lumped.fdtd.json'
    solver_lumped = FDTD(fn_lumped, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver_lumped.run()


    StartLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellStart")[0])
    EndLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellEnd")[0])
    AdjacentPostLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PostLumpedCell")[0])
    AdjacentPreLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PreLumpedCell")[0])

    R = solver_lumped.getMaterialProperties("lumped_RC")["resistance"]
    C = solver_lumped.getMaterialProperties("lumped_RC")["capacitance"]
    L = 1.65e-7 # parasitic inductance mentioned above

    num = [R*C, 1]
    den = [L*R*C, L, R]
    system = signal.TransferFunction(num, den)
    tout, I_out, _ = signal.lsim(system, 
                                 U=np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=1), 
                                 T=np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=0))
    
    I_theo = np.interp(AdjacentPreLumpedProbe['time'], tout, I_out)
    
    assert np.corrcoef(AdjacentPostLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(AdjacentPreLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(StartLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(EndLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999

def test_lumped_inductor(tmp_path):
    # This test validates the behavior of lumped inductor materials in a simplified circuit. The lumped inductor
    # can be modeled as an inductor in series with a resistor.
    # The circuit consists of a 40mm x 40mm simple loop with a lumped resistor line inserted along one edge.
    # Due to the geometry, the circuit naturally exhibits a parasitic inductance of approximately 1.65e-7 H.
    #
    # Current probes are placed at the start and end of the lumped line, as well as a few cells before and after it.
    # These measurements are used to evaluate the accuracy of the lumped material model.
    #
    # For validation, the results are compared against two reference cases:
    # 1. A simple loop circuit with the same dimensions where a terminal  is inserted in place of the lumped line, 
    #    using the same resistance and inductance of the lumped line.
    # 2. Theoretical current response calculated using Laplace transforms from the initial pulse excitation.
    #
    # For better interaction with the case, the user can go to the file: testData/cases/lumped_lines/simple_loop_RL/simple_loop_prepost.py

    fn_lumped = CASES_FOLDER + 'lumped_lines/simple_loop_RL/simple_loop_lumped.fdtd.json'
    fn_terminal = CASES_FOLDER + 'lumped_lines/simple_loop_RL/simple_loop_terminal.fdtd.json'
    
    solver_lumped = FDTD(fn_lumped, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver_terminal = FDTD(fn_terminal, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver_lumped.run()
    solver_terminal.run()

    StartTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("TerminalCellStart")[0])
    StartLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellStart")[0])

    EndTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("TerminalCellEnd")[0])
    EndLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("LumpedCellEnd")[0])

    AdjacentPostLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PostLumpedCell")[0])
    AdjacentPostTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("PostTerminalCell")[0])

    AdjacentPreLumpedProbe = Probe(solver_lumped.getSolvedProbeFilenames("PreLumpedCell")[0])
    AdjacentPreTerminalProbe = Probe(solver_terminal.getSolvedProbeFilenames("PreTerminalCell")[0])

    assert np.corrcoef(StartLumpedProbe['current'].to_numpy(), StartTerminalProbe['current'].to_numpy())[0, 1] > 0.999
    assert np.corrcoef(EndLumpedProbe['current'].to_numpy(), EndTerminalProbe['current'].to_numpy())[0, 1] > 0.999
    assert np.corrcoef(AdjacentPostLumpedProbe['current'].to_numpy(), AdjacentPostTerminalProbe['current'].to_numpy())[0, 1] > 0.999
    assert np.corrcoef(AdjacentPreLumpedProbe['current'].to_numpy(), AdjacentPreTerminalProbe['current'].to_numpy())[0, 1] > 0.999

    R = solver_lumped.getMaterialProperties("lumped_RL")["resistance"]
    L = solver_lumped.getMaterialProperties("lumped_RL")["inductance"] + 1.65e-7 # inductance of the lumped line + parasitic inductance mentioned above

    num = [1]
    den = [L, R]
    system = signal.TransferFunction(num, den)
    tout, I_out, _ = signal.lsim(system, 
                                 U=np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=1), 
                                 T=np.loadtxt(solver_lumped["sources"][0]["magnitudeFile"], usecols=0))
    
    I_theo = np.interp(AdjacentPreLumpedProbe['time'], tout, I_out)
    
    assert np.corrcoef(AdjacentPostLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(AdjacentPreLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(StartLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(EndLumpedProbe['current'].to_numpy(), I_theo)[0, 1] > 0.999

def test_lumped_resistor_parallel_terminal_resistor(tmp_path):
    # This test verifies current splitting behavior in a parallel resistive configuration.
    # The setup consists of a 40mm x 40mm circuit with two parallel elements:
    # - A lumped resistor line
    # - A terminal with the same resistance
    #
    # Current measurements are taken at three key locations:
    # 1. At the source (input) (Bulk Initial probe)
    # 2. At the terminal branch (Bulk Top probe)
    # 3. At the lumped resistor line branch (Bulk Bottom probe)
    #
    # The goal is to validate that the current divides between the two parallel resistive paths, i.e., the
    # current at the terminal branch plus the lumped resistor line branch should be equal to the current at the source.
    # The test also compares the current at the source with the theoretical current response calculated using Laplace transforms.
    #
    # For better interaction with the case, the user can go to the file: testData/cases/lumped_lines/current_bifurcation/current_bifurcation_lumped_prepost.py

    fn = CASES_FOLDER + 'lumped_lines/current_bifurcation/current_bifurcation_lumped.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()


    InitialBulk_probe = Probe(solver.getSolvedProbeFilenames("Bulk Initial probe")[0])
    TopBulk_probe = Probe(solver.getSolvedProbeFilenames("Bulk Top probe")[0])
    BottomBulk_probe = Probe(solver.getSolvedProbeFilenames("Bulk Bottom probe")[0])

    R_lumped = solver.getMaterialProperties("lumped_resistor")["resistance"]
    R_terminal = solver.getMaterialProperties("Terminal_R")["terminations"][0]["resistance"]
    R = 1/(1/R_lumped + 1/R_terminal)  
    L = 1.65e-7 # parasitic inductance mentioned above

    num = [1]
    den = [L, R]
    system = signal.TransferFunction(num, den)
    tout, I_out, _ = signal.lsim(system, 
                                 U=np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=1), 
                                 T=np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=0))
    
    I_theo = np.interp(InitialBulk_probe['time'], tout, I_out)
    
    assert np.corrcoef(TopBulk_probe['current'].to_numpy() + BottomBulk_probe['current'].to_numpy(), I_theo)[0, 1] > 0.999
    assert np.corrcoef(InitialBulk_probe['current'].to_numpy(), I_theo)[0, 1] > 0.999

def test_offset_normal_in_x(tmp_path):
    # This test verifies the positive offset along the normal vector (in x-direction) respect to the bulk plane
    # used to measure the current values of the system.
    # The setup consists in a polyline with points [(0 mm,0 mm,0 mm), (20 mm,0 mm,0 mm), (20 mm,-20 mm,0 mm)]
    # as nodal source and three bulk planes defined at x=18 mm, x=20 mm and x=22 mm.
    # The test checks that only the plane defined at x=18 mm has non-zero current values.

    fn = CASES_FOLDER + 'bulk_current_tests/offSet_x/offSet_x.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    I_in = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=1)
    time = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=0)
    
    probe_at_x_18 = Probe(solver.getSolvedProbeFilenames("BulkCurrent1")[0])
    probe_at_x_20 = Probe(solver.getSolvedProbeFilenames("BulkCurrent2")[0])
    probe_at_x_22 = Probe(solver.getSolvedProbeFilenames("BulkCurrent3")[0])
    
    I_interp = np.interp(
        probe_at_x_18['time'].to_numpy(),
        time,
        I_in
    )

    assert np.corrcoef(probe_at_x_18['current'].to_numpy(), I_interp)[0, 1] > 0.999
    assert np.allclose(probe_at_x_20['current'].to_numpy(), 0.0, atol=1.5e-3)
    assert np.allclose(probe_at_x_22['current'].to_numpy(), 0.0, atol=1.5e-3)

def test_offset_normal_in_y(tmp_path):
    # This test verifies the positive offset along the normal vector (in y-direction) respect to the bulk plane
    # used to measure the current values of the system.
    # The setup consists in a polyline with points [(0 mm,0 mm,0 mm), (0 mm,0 mm,20 mm), (0 mm,-20 mm,20 mm)]
    # as nodal source and three bulk planes defined at y=-2 mm, y=0 mm and y=2 mm.
    # The test checks that only the plane defined at y=-2 mm has non-zero current values.

    fn = CASES_FOLDER + 'bulk_current_tests/offSet_y/offSet_y.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    I_in = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=1)
    time = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=0)
    
    probe_at_y_m2 = Probe(solver.getSolvedProbeFilenames("BulkCurrent1")[0])
    probe_at_y_0 = Probe(solver.getSolvedProbeFilenames("BulkCurrent2")[0])
    probe_at_y_2 = Probe(solver.getSolvedProbeFilenames("BulkCurrent3")[0])
    
    I_interp = np.interp(
        probe_at_y_m2['time'].to_numpy(),
        time,
        I_in
    )

    assert np.corrcoef(np.abs(probe_at_y_m2['current'].to_numpy()), I_interp)[0, 1] > 0.999
    assert np.allclose(probe_at_y_0['current'].to_numpy(), 0.0, atol=1.5e-3)
    assert np.allclose(probe_at_y_2['current'].to_numpy(), 0.0, atol=1.5e-3)

def test_offset_normal_in_z(tmp_path):
    # This test verifies the positive offset along the normal vector (in z-direction) respect to the bulk plane
    # used to measure the current values of the system.
    # The setup consists in a polyline with points [(0 mm,0 mm,0 mm), (0 mm,0 mm,20 mm), (0 mm,-20 mm,20 mm)]
    # as nodal source and three bulk planes defined at z=18 mm, z=20 mm and z=22 mm.
    # The test checks that only the plane defined at z=18 mm has non-zero current values.

    fn = CASES_FOLDER + 'bulk_current_tests/offSet_z/offSet_z.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    I_in = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=1)
    time = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=0)
    
    probe_at_z_18 = Probe(solver.getSolvedProbeFilenames("BulkCurrent1")[0])
    probe_at_z_20 = Probe(solver.getSolvedProbeFilenames("BulkCurrent2")[0])
    probe_at_z_22 = Probe(solver.getSolvedProbeFilenames("BulkCurrent3")[0])
    
    I_interp = np.interp(
        probe_at_z_18['time'].to_numpy(),
        time,
        I_in
    )

    assert np.corrcoef(probe_at_z_18['current'].to_numpy(), I_interp)[0, 1] > 0.999
    assert np.allclose(probe_at_z_20['current'].to_numpy(), 0.0, atol=1.5e-3)
    assert np.allclose(probe_at_z_22['current'].to_numpy(), 0.0, atol=1.5e-3)

def test_offset_perpendicular_in_x(tmp_path):
    # This test verifies the negative offset presented in the y and z directions when the bulk plane is defined
    # with a normal vector in the x-direction.
    # The setup consists in three nodal sources defined as follow:
    #   - Line [(0 mm,18 mm,20 mm), (50 mm,18 mm,20 mm)]. Gaussian pulse with 1 A amplitude.
    #   - Line [(0 mm,20 mm,20 mm), (50 mm,20 mm,20 mm)]. Gaussian pulse with 2 A amplitude.
    #   - Line [(0 mm,18 mm,18 mm), (50 mm,18 mm,18 mm)]. Gaussian pulse with 3 A amplitude.
    #
    # The nodal lines are separeted by one cell in y and z directions. We define four bulk planes:
    #   - Plane at x=30 mm, (16 mm, 20 mm) and (20 mm, 24 mm) in y and z directions. For nodal source 1.
    #   - Plane at x=30 mm, (20 mm, 24 mm) and (20 mm, 24 mm) in y and z directions. For nodal source 2.
    #   - Plane at x=30 mm, (16 mm, 20 mm) and (16 mm, 20 mm) in y and z directions. For nodal source 3.
    #   - Plane at x=30 mm, (16 mm, 24 mm) and (16 mm, 24 mm) in y and z directions. For the total current.
    #
    # The test checks that the bulk planes only measure the current values of the respective nodal sources
    # and if we move one negative cell in the y and z directions, the current values are zero for the bulk planes
    # 1 and 3 respectively; similar behavior if we move one positive cell in the y and z directions for  
    # the bulk planes 1 and 2. This proves that the bulk have a negative offset in the y and z directions.

    fn = CASES_FOLDER + 'bulk_current_tests/threeLines_offSet_x_Perpendicular/threeLines.fdtd.json'
    fn_negative = CASES_FOLDER + 'bulk_current_tests/threeLines_offSet_x_Perpendicular/threeLinesNegative.fdtd.json'
    fn_positive = CASES_FOLDER + 'bulk_current_tests/threeLines_offSet_x_Perpendicular/threeLinesPositive.fdtd.json'
    
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver_negative = FDTD(input_filename = fn_negative, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver_positive = FDTD(input_filename = fn_positive, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()
    solver_negative.run()
    solver_positive.run()

    probe1 = Probe(solver.getSolvedProbeFilenames("Bulk probe1")[0])
    probe2 = Probe(solver.getSolvedProbeFilenames("Bulk probe2")[0])
    probe3 = Probe(solver.getSolvedProbeFilenames("Bulk probe3")[0])
    probeTotal = Probe(solver.getSolvedProbeFilenames("Bulk probe total")[0])

    assert np.corrcoef(probe1['current'].to_numpy() + probe2['current'].to_numpy() + probe3['current'].to_numpy(),
                          probeTotal['current'].to_numpy())[0, 1] > 0.999
    
    probe1_negative = Probe(solver_negative.getSolvedProbeFilenames("Bulk probe1")[0])
    probe2_negative = Probe(solver_negative.getSolvedProbeFilenames("Bulk probe2")[0])
    probe3_negative = Probe(solver_negative.getSolvedProbeFilenames("Bulk probe3")[0])

    I_in_2 = np.loadtxt(solver["sources"][1]["magnitudeFile"], usecols=1)
    time_2 = np.loadtxt(solver["sources"][1]["magnitudeFile"], usecols=0)

    I_2_interp = np.interp(
        probe2_negative['time'].to_numpy(),
        time_2,
        I_in_2
    )

    assert np.corrcoef(probe2_negative['current'].to_numpy(), I_2_interp)[0, 1] > 0.999
    assert np.allclose(probe1_negative['current'].to_numpy(), 0.0, atol=1.5e-2)
    assert np.allclose(probe3_negative['current'].to_numpy(), 0.0, atol=1.5e-2)

    probe1_positive = Probe(solver_positive.getSolvedProbeFilenames("Bulk probe1")[0])
    probe2_positive = Probe(solver_positive.getSolvedProbeFilenames("Bulk probe2")[0])
    probe3_positive = Probe(solver_positive.getSolvedProbeFilenames("Bulk probe3")[0])

    I_in_3 = np.loadtxt(solver["sources"][2]["magnitudeFile"], usecols=1)
    time_3 = np.loadtxt(solver["sources"][2]["magnitudeFile"], usecols=0)

    I_3_interp = np.interp(
        probe3_positive['time'].to_numpy(),
        time_3,
        I_in_3
    )

    assert np.corrcoef(probe3_positive['current'].to_numpy(), I_3_interp)[0, 1] > 0.999
    assert np.allclose(probe1_positive['current'].to_numpy(), 0.0, atol=1.5e-2)
    assert np.allclose(probe2_positive['current'].to_numpy(), 0.0, atol=1.5e-2)

def test_negative_offset_in_x(tmp_path):
    # Following the previous test, we have seen that the bulk surfaces has negative offsets in the directions
    # perpendicular to the normal vector of the bulk surface. The previous test checks the negative offset in the
    # y and z directions when the normal vector is in the x-direction. Now we check the negative offset in the
    # x-direction when the normal vector is in the y-direction.
    # The setup consists in a nodal source placed in a line defined by [(0 mm,0 mm,0 mm), (0 mm,50 mm,0 mm)]. 
    # The nodal source is a Gaussian pulse with 1 A amplitude. And three bulk planes defined at:
    #  - Plane at y=36 mm, (-4 mm, -2 mm) and (0 mm, 2 mm) in x and z directions.
    #  - Plane at y=38 mm, (0 mm, -2 mm) and (4 mm, 2 mm) in x and z directions.
    #  - Plane at y=40 mm, (-2 mm, -2 mm) and (2 mm, 2 mm) in x and z directions.
    #
    # All the bulk planes measure the current values of the nodal source, however, the first only captures the
    # current values at the right edge of the plane, similarly, the second captures the current values at the
    # left edge of the plane. This test checks that the second and third bulk planes measure correctly the current values
    # while the first bulk plane has zero current values. This proves that the bulk have a negative offset in the x-direction.

    fn = CASES_FOLDER + 'bulk_current_tests/negative_offSet_x/offSet_negative_x.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    I_in = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=1)
    time = np.loadtxt(solver["sources"][0]["magnitudeFile"], usecols=0)
    
    probeR = Probe(solver.getSolvedProbeFilenames("Bulk_right")[0])
    probeL = Probe(solver.getSolvedProbeFilenames("Bulk_left")[0])
    probeTotal = Probe(solver.getSolvedProbeFilenames("BulkTotal")[0])
    
    I_interp = np.interp(
        probeTotal['time'].to_numpy(),
        time,
        I_in
    )

    assert np.corrcoef(probeTotal['current'].to_numpy(), I_interp)[0, 1] > 0.999
    assert np.corrcoef(probeL['current'].to_numpy(), I_interp)[0, 1] > 0.999
    assert np.allclose(probeR['current'].to_numpy(), 0.0, atol=3e-3)