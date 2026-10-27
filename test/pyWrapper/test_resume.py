from test.utils.utils import *
import filecmp
import re
import signal
import subprocess
import time

MOVIES = (
    "electric_field_movie_ME_0_0_0__29_29_9",
    "electric_field_movie_low_ME_0_0_0__29_29_1",
    "electric_field_movie_high_ME_0_0_8__29_29_9",
)


def _run(case, *flags):
    return subprocess.run([SEMBA_EXE, "-i", case, *flags], cwd=os.path.dirname(case),
                          capture_output=True, text=True)


def _set_steps(case, steps):
    data = json.load(open(case))
    data["general"]["numberOfSteps"] = steps
    json.dump(data, open(case, "w"))


@pytest.mark.slow
@pytest.mark.hdf
@pytest.mark.movie
def test_resume_after_sigkill_keeps_previous_movies(tmp_path):
    case_folder = os.path.join(CASES_FOLDER, "planewave")
    for name in ("pw-in-box-with-partial-movies.fdtd.json", "gauss_1GHz.exc"):
        shutil.copy(os.path.join(case_folder, name), tmp_path)
    case = str(tmp_path / "case.fdtd.json")
    os.rename(tmp_path / "pw-in-box-with-partial-movies.fdtd.json", case)

    # far more steps than can finish: the run is only ever stopped by the kill
    _set_steps(case, 100_000_000)
    solver = subprocess.Popen(
        [SEMBA_EXE, "-i", case, "-flush", "1"],  # shortest restart-fields flush interval, one minute
        cwd=tmp_path, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    try:
        deadline = time.time() + 180
        while not list(tmp_path.glob("*.fields")):
            assert solver.poll() is None, "solver ended before writing restart fields"
            assert time.time() < deadline, "restart fields were not written"
            time.sleep(1)
        time.sleep(2)  # let the writer finish the file
    finally:
        solver.send_signal(signal.SIGKILL)
        solver.wait()

    # a resume cannot end before the saved step; the failed attempt reports it
    _set_steps(case, 1)
    probe = _run(case, "-r")
    saved_step = int(re.search(r"Last executed time step\s+(\d+)", probe.stdout).group(1))

    _set_steps(case, saved_step + 200)
    resumed = _run(case, "-r")
    assert resumed.returncode == 0, resumed.stdout[-2000:]

    import h5py

    for movie in MOVIES:
        folder = tmp_path / f"case.fdtd_{movie}"
        previous = folder / "before_resume_1" / f"{folder.name}.h5"
        assert previous.is_file()
        assert (folder / f"{folder.name}.h5").is_file()
        for path in (previous, folder / f"{folder.name}.h5"):
            with h5py.File(path, "r"):
                pass


@pytest.mark.slow
@pytest.mark.hdf
def test_resume_after_sigkill_reproduces_uninterrupted_probes(tmp_path):
    reference = tmp_path / "reference"
    interrupted = tmp_path / "interrupted"
    case_folder = os.path.join(CASES_FOLDER, "planewave")
    for folder in (reference, interrupted):
        folder.mkdir()
        for name in ("pw-in-box.fdtd.json", "gauss_1GHz.exc"):
            shutil.copy(os.path.join(case_folder, name), folder)
    case_reference = str(reference / "pw-in-box.fdtd.json")
    case_interrupted = str(interrupted / "pw-in-box.fdtd.json")

    _set_steps(case_interrupted, 100_000_000)
    solver = subprocess.Popen(
        [SEMBA_EXE, "-i", case_interrupted, "-flush", "1"],
        cwd=interrupted, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    try:
        deadline = time.time() + 180
        while not list(interrupted.glob("*.fields")):
            assert solver.poll() is None, "solver ended before writing restart fields"
            assert time.time() < deadline, "restart fields were not written"
            time.sleep(1)
        time.sleep(2)  # let the writer finish the file
    finally:
        solver.send_signal(signal.SIGKILL)
        solver.wait()

    _set_steps(case_interrupted, 1)
    probe = _run(case_interrupted, "-r", "-n", "1")
    saved_step = int(re.search(r"Last executed time step\s+(\d+)", probe.stdout).group(1))
    steps = saved_step + 2000

    uninterrupted = _run(case_reference, "-n", str(steps))
    assert uninterrupted.returncode == 0, uninterrupted.stdout[-2000:]
    resumed = _run(case_interrupted, "-r", "-n", str(steps))
    assert resumed.returncode == 0, resumed.stdout[-2000:]

    # Energy is reported on wall-clock intervals, so only the deterministic outputs are compared
    for name in ("before_Ex_3_3_1_tm", "after_Ex_3_3_5_tm", "inbox_Ex_3_3_3_tm", "inbox_Ex_3_3_3_fq"):
        expected = reference / f"pw-in-box.fdtd_{name}.dat"
        assert filecmp.cmp(expected, interrupted / expected.name, shallow=False), name
