from test.utils.utils import *
from scipy.constants import c, epsilon_0
from scipy.signal import find_peaks

CAVITY_SIZE = (0.20, 0.16, 0.12)
FREQS = np.linspace(1e7, 7e8, 200)


def _read_probe(solver, name):
    dat = Probe.getDatFile(Probe(solver.getSolvedProbeFolders(name)[0]))
    samples = np.loadtxt(dat, skiprows=1)
    return samples[:, 0], samples[:, 1]


def _cavity_modes(n_modes):
    a, b, h = CAVITY_SIZE
    modes = [
        c / 2 * np.sqrt((i / a) ** 2 + (j / b) ** 2 + (k / h) ** 2)
        for i in range(3) for j in range(3) for k in range(3)
        if i > 0 and j > 0  # modes with Ez, as seen by the z probe and z source
    ]
    return np.sort(modes)[:n_modes]


@pytest.mark.slow
def test_pec_cavity_resonances_and_late_time_stability(tmp_path):
    solver = FDTD(CASES_FOLDER + "stability/pec_cavity.fdtd.json",
                  path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    t, ez = _read_probe(solver, "ez")

    spectrum = np.abs(np.fft.rfft(ez * np.hanning(len(ez))))
    freqs = np.fft.rfftfreq(len(ez), t[1] - t[0])
    peaks, _ = find_peaks(spectrum, height=0.05 * spectrum.max())
    found = freqs[peaks]
    for expected in _cavity_modes(4):  # first four modes, below the grid dispersion error
        assert np.min(np.abs(found - expected)) < 0.01 * expected

    # lossless and source-free after the pulse: no late-time growth
    quarter = len(ez) // 4
    assert np.abs(ez[-quarter:]).max() <= 1.05 * np.abs(ez[quarter:2 * quarter]).max()


@pytest.mark.slow
def test_lossy_cavity_decay_matches_conductivity(tmp_path):
    fn = CASES_FOLDER + "stability/lossy_cavity.fdtd.json"
    solver = FDTD(fn, path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    t, ez = _read_probe(solver, "ez")
    sigma = solver.getMaterialProperties("lossy")["electricConductivity"]

    # envelope decay of a uniformly filled cavity: exp(-sigma t / (2 eps0))
    win = 400
    n = len(ez) // win
    env = np.abs(ez[: n * win]).reshape(n, win).max(axis=1)
    tt = t[: n * win].reshape(n, win).mean(axis=1)
    sel = tt > 40e-9  # after the pulse has left the source
    rate = -np.polyfit(tt[sel], np.log(env[sel]), 1)[0]
    assert np.isclose(rate, sigma / (2 * epsilon_0), rtol=0.1)


def test_cpml_reflection_below_minus_50_db(tmp_path):
    traces = {}
    for size in ("small", "large"):
        folder = tmp_path / size
        folder.mkdir()
        solver = FDTD(CASES_FOLDER + f"stability/cpml_{size}.fdtd.json",
                      path_to_exe=SEMBA_EXE, run_in_folder=folder)
        solver.run()
        traces[size] = _read_probe(solver, "ez")[1]

    # the large domain is reflection-free within the simulated time, so it is the reference
    n = min(len(v) for v in traces.values())
    error = np.abs(traces["small"][:n] - traces["large"][:n]).max()
    assert 20 * np.log10(error / np.abs(traces["large"]).max()) < -50


@pytest.mark.slow
@pytest.mark.parametrize("cfl", [0.99, 1.1])  # the solver brings dt back to the CFL limit when asked for more
def test_courant_limit_stays_bounded(tmp_path, cfl):
    solver = FDTD(CASES_FOLDER + "stability/courant_cavity.fdtd.json",
                  path_to_exe=SEMBA_EXE, run_in_folder=tmp_path, flags=f"-cfl {cfl}")
    solver.run()
    ez = _read_probe(solver, "ez")[1]
    half = len(ez) // 2
    assert np.all(np.isfinite(ez))
    assert np.abs(ez[half:]).max() <= 1.05 * np.abs(ez[:half]).max()


@pytest.mark.slow
def test_dielectric_interface_cavity_resonances(tmp_path):
    from scipy.optimize import brentq

    solver = FDTD(CASES_FOLDER + "stability/dielectric_cavity.fdtd.json",
                  path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    t, ez = _read_probe(solver, "ez")

    # half-filled 1D cavity: tan(k1 a)/k1 + tan(k2 b)/k2 = 0, k1 = w/c, k2 = n w/c
    a = b = 0.20
    n = 2.0
    g = lambda f: (lambda k: np.sin(k * a) * n * np.cos(n * k * b) + np.cos(k * a) * np.sin(n * k * b))(2 * np.pi * f / c)
    grid = np.linspace(0.2e9, 1.6e9, 4000)
    roots = [brentq(g, x0, x1) for x0, x1 in zip(grid[:-1], grid[1:]) if g(x0) * g(x1) < 0][:3]

    spectrum = np.abs(np.fft.rfft(ez * np.hanning(len(ez))))
    freqs = np.fft.rfftfreq(len(ez), t[1] - t[0])
    peaks, _ = find_peaks(spectrum, height=0.02 * spectrum.max())
    for expected in roots:
        assert np.min(np.abs(freqs[peaks] - expected)) < 0.03 * expected
    quarter = len(ez) // 4
    assert np.abs(ez[-quarter:]).max() <= 1.05 * np.abs(ez[quarter:2 * quarter]).max()


@pytest.mark.slow
def test_cpml_grazing_long_run_decays(tmp_path):
    solver = FDTD(CASES_FOLDER + "stability/cpml_grazing.fdtd.json",
                  path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    ez = _read_probe(solver, "ez")[1]
    assert np.all(np.isfinite(ez))
    assert np.abs(ez[-len(ez) // 4:]).max() < 1e-3 * np.abs(ez).max()


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.slow
def test_cavity_mpi_matches_single_rank(tmp_path):
    traces = {}
    for ranks in (1, 2):
        folder = tmp_path / f"np{ranks}"
        folder.mkdir()
        solver = FDTD(CASES_FOLDER + "stability/courant_cavity.fdtd.json", path_to_exe=SEMBA_EXE,
                      mpi_command=f"mpirun -np {ranks}", run_in_folder=folder)
        solver.run()
        traces[ranks] = _read_probe(solver, "ez")[1]
    assert np.allclose(traces[1], traces[2], rtol=1e-4, atol=1e-4 * np.abs(traces[1]).max())


@pytest.mark.farfield
@pytest.mark.slow  # ~9 min
def test_staircase_sphere_rcs_matches_mie(tmp_path):
    solver = FDTD(CASES_FOLDER + "stability/mie_sphere.fdtd.json",
                  path_to_exe=SEMBA_EXE, run_in_folder=tmp_path)
    solver.run()
    assert solver.hasFinishedSuccessfully()

    far = Probe(solver.getSolvedProbeFolders("n2f")[0])
    sel = (far.data["Theta"] == 90.0) & (far.data["Phi"] == 0.0)
    rcs_sim = np.interp(FREQS, far.data.loc[sel, "freq"], far.data.loc[sel, "rcs_geom"])
    rcs_mie = RCS(fspace=FREQS, radius=0.5)

    band = slice(5, 150)  # 27-530 MHz; 20 cells per radius (10 cells gave corr 0.962 here)
    assert np.corrcoef(rcs_mie[band], rcs_sim[band])[0, 1] > 0.98
    assert np.median(np.abs(rcs_sim[band] / rcs_mie[band] - 1)) < 0.15
