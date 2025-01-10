from utils import *


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
        assert np.allclose(p_expected[i].df.to_numpy()[:, 0:3], p_solved._recordedData.to_numpy()[
                           :, 0:3], rtol=5e-2, atol=0.2)
    for i in [0, 1, 4, 5]:
        p_solved = Probe(probe_files[i])
        assert np.allclose(p_expected[i].df.to_numpy()[:, 0:4], p_solved._recordedData.to_numpy()[
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
        p_expected._recordedData.to_numpy()[:, 0],
        p_solved._recordedData.to_numpy()[:, 0],
        rtol=0.0, atol=10e-8)
    assert np.allclose(
        p_expected._recordedData.to_numpy()[:, 1],
        p_solved._recordedData.to_numpy()[:, 1],
        rtol=0.0, atol=10e-8)
    assert np.allclose(
        p_expected._recordedData.to_numpy()[:, 2],
        p_solved._recordedData.to_numpy()[:, 2],
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
        p_expected._recordedData.to_numpy()[:, 0:3], 
        p_solved._recordedData.to_numpy()[:, 0:3], 
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
        p_expected._recordedData.to_numpy()[:, 0:3], 
        p_solved._recordedData.to_numpy()[:, 0:3], 
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
        assert np.allclose(p_expected[i]._recordedData.to_numpy()[:, 0:3], p_solved._recordedData.to_numpy()[
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

    h5file = solver.getSolvedProbeFilenames("electric_field_movie")[2]
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

    np.allclose(inbox._recordedData['field'], inbox._recordedData['incident'])
    zeros = np.zeros_like(before._recordedData['field'])
    np.allclose(before._recordedData['field'], zeros)
    np.allclose(after._recordedData['field'], zeros)


def test_planewave_with_periodic_boundaries(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-with-periodic.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()
    
    before = Probe(solver.getSolvedProbeFilenames("before")[0])
    inbox = Probe(solver.getSolvedProbeFilenames("inbox")[0])
    after = Probe(solver.getSolvedProbeFilenames("after")[0])

    np.allclose(inbox._recordedData['field'], inbox._recordedData['incident'])
    zeros = np.zeros_like(before._recordedData['field'])
    np.allclose(before._recordedData['field'], zeros)
    np.allclose(after._recordedData['field'], zeros)


def test_sgbc_shielding_effectiveness(tmp_path):
    fn = CASES_FOLDER + 'sgbcShieldingEffectiveness/shieldingEffectiveness.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    # FDTD results
    back = Probe(solver.getSolvedProbeFilenames("back")[0])

    t = back._recordedData['time']
    dt = t[1] - t[0]
    fq = fftfreq(len(t))/dt
    INC = fft(back._recordedData['incident'])
    BACK = fft(back._recordedData['field'])
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
