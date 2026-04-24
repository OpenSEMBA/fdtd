from utils import *


def test_point_probe_regression_planewave_in_box(tmp_path):
    fn = CASES_FOLDER + 'planewave/pw-in-box.fdtd.json'
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)

    solver.run()

    before = Probe(solver.getSolvedProbeFilenames("before")[0])
    inbox = Probe(solver.getSolvedProbeFilenames("inbox")[0])
    after = Probe(solver.getSolvedProbeFilenames("after")[0])

    assert np.corrcoef(
        inbox.data['field'].to_numpy(),
        inbox.data['incident'].to_numpy(),
    )[0, 1] > 0.999

    zeros = np.zeros_like(before.data['field'])
    assert np.allclose(before.data['field'].to_numpy(), zeros, atol=5e-4)
    assert np.allclose(after.data['field'].to_numpy(), zeros, atol=5e-4)


@mtln_skip
def test_wire_probe_regression_towel_hanger(tmp_path):
    fn = CASES_FOLDER + 'towelHanger/towelHanger.fdtd.json'
    solver = FDTD(input_filename=fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()

    p_solved = [
        Probe(solver.getSolvedProbeFilenames("wire_start")[0]),
        Probe(solver.getSolvedProbeFilenames("wire_mid")[0]),
        Probe(solver.getSolvedProbeFilenames("wire_end")[0]),
    ]

    p_expected = [
        Probe(OUTPUTS_FOLDER + 'towelHanger.fdtd_wire_start_Wz_27_25_30_s1.dat'),
        Probe(OUTPUTS_FOLDER + 'towelHanger.fdtd_wire_mid_Wx_35_25_32_s5.dat'),
        Probe(OUTPUTS_FOLDER + 'towelHanger.fdtd_wire_end_Wz_43_25_30_s4.dat'),
    ]

    for i in range(3):
        solved = np.interp(
            p_expected[i]['time'].to_numpy(),
            p_solved[i]['time'].to_numpy(),
            p_solved[i]['current_0'].to_numpy(),
        )
        assert np.corrcoef(solved, p_expected[i]['current_0'])[0, 1] > 0.999
