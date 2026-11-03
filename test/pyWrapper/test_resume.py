from test.utils.utils import *
import filecmp
import re
import signal
import subprocess
import sys
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
        solver.kill()  # SIGKILL on POSIX, TerminateProcess on Windows (signal.SIGKILL does not exist there)
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
        solver.kill()  # SIGKILL on POSIX, TerminateProcess on Windows (signal.SIGKILL does not exist there)
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


def _flush_case(tmp_path, *, steps=20, minutes=None, additional=None):
    case_folder = os.path.join(CASES_FOLDER, "planewave")
    for name in ("pw-in-box.fdtd.json", "gauss_1GHz.exc"):
        shutil.copy(os.path.join(case_folder, name), tmp_path)
    case = str(tmp_path / "pw-in-box.fdtd.json")
    data = json.load(open(case))
    data["general"]["numberOfSteps"] = steps
    if minutes is not None:
        data["general"]["restartFlushMinutes"] = minutes
    if additional is not None:
        data["general"]["additionalArguments"] = additional
    json.dump(data, open(case, "w"))
    return case


@pytest.mark.parametrize(
    "minutes, additional, flags, expected",
    [
        (7, None, [], 7),  # json general.restartFlushMinutes alone
        (7, "-flush 3", [], 3),  # additionalArguments overrides the json field
        (7, "-flush 3", ["-flush", "2"], 2),  # console flag overrides both
    ],
)
def test_restart_flush_priority(tmp_path, minutes, additional, flags, expected):
    case = _flush_case(tmp_path, minutes=minutes, additional=additional)
    run = _run(case, *flags)
    assert run.returncode == 0, run.stdout[-2000:]
    assert f"Restart FIELDS flush: every {expected} minutes" in run.stdout


def test_restart_flush_negative_is_rejected(tmp_path):
    case = _flush_case(tmp_path, minutes=-1)
    run = _run(case)
    assert run.returncode != 0
    assert "restartFlushMinutes cannot be negative" in run.stdout + run.stderr


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif(sys.platform == "win32", reason="kills the job's process group, POSIX only")
def test_resume_rejects_different_process_count_then_resumes_with_same(tmp_path):
    case = _flush_case(tmp_path, steps=100_000_000, minutes=1)
    solver = subprocess.Popen(
        ["mpirun", "-np", "2", SEMBA_EXE, "-i", case],
        cwd=tmp_path, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True,
    )
    try:
        deadline = time.time() + 300
        while len(list(tmp_path.glob("*.fields"))) < 2:
            assert solver.poll() is None, "solver ended before writing restart fields"
            assert time.time() < deadline, "restart fields were not written by both ranks"
            time.sleep(1)
        time.sleep(2)  # let the writers finish the files
    finally:
        os.killpg(solver.pid, signal.SIGKILL)  # only the job launched above
        solver.wait()

    _set_steps(case, 1)
    wrong = subprocess.run(["mpirun", "-np", "1", SEMBA_EXE, "-i", case, "-r"], cwd=tmp_path,
                           capture_output=True, text=True)
    assert wrong.returncode != 0
    assert "mpirun -n 2" in wrong.stdout + wrong.stderr  # the saved process count is named

    # same count as saved: the probe reports the checkpointed step, as in the tests above
    probe = subprocess.run(["mpirun", "-np", "2", SEMBA_EXE, "-i", case, "-r"], cwd=tmp_path,
                           capture_output=True, text=True)
    saved_step = int(re.search(r"Last executed time step\s+(\d+)", probe.stdout).group(1))

    _set_steps(case, saved_step + 200)
    resumed = subprocess.run(["mpirun", "-np", "2", SEMBA_EXE, "-i", case, "-r"], cwd=tmp_path,
                             capture_output=True, text=True)
    assert resumed.returncode == 0, resumed.stdout[-2000:]


def _mpi_run(case, *flags):
    return subprocess.run(["mpirun", "-np", "2", SEMBA_EXE, "-i", case, *flags], cwd=os.path.dirname(case),
                          capture_output=True, text=True)


def _assert_resumed_frames_match_reference(resumed_h5, reference_h5):
    # every frame the resumed run wrote must equal the uninterrupted frame recorded at the same time
    import h5py
    import numpy as np

    with h5py.File(reference_h5, "r") as ref, h5py.File(resumed_h5, "r") as res:
        index = {round(float(t), 15): i for i, t in enumerate(ref["series/values"][...])}
        ref_frames = ref["attributes/a0001/values"]
        res_frames = res["attributes/a0001/values"]
        for j, t in enumerate(res["series/values"][...]):
            i = index[round(float(t), 15)]
            assert np.array_equal(res_frames[j], ref_frames[i]), f"frame at t={t} differs after resume"


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.hdf
@pytest.mark.movie
@pytest.mark.farfield
@pytest.mark.slow
@pytest.mark.skipif(sys.platform == "win32", reason="kills the job's process group, POSIX only")
def test_resume_planewave_sphere_with_farfield_point_probe_and_movies_matches_uninterrupted(tmp_path):
    reference = tmp_path / "reference"
    interrupted = tmp_path / "interrupted"
    case_folder = os.path.join(CASES_FOLDER, "sphere")
    for folder in (reference, interrupted):
        folder.mkdir()
        for name in ("sphere_restart.fdtd.json", "gauss.exc"):
            shutil.copy(os.path.join(case_folder, name), folder)
    case_reference = str(reference / "sphere_restart.fdtd.json")
    case_interrupted = str(interrupted / "sphere_restart.fdtd.json")

    _set_steps(case_interrupted, 100_000_000)
    solver = subprocess.Popen(
        ["mpirun", "-np", "2", SEMBA_EXE, "-i", case_interrupted, "-flush", "1"],
        cwd=interrupted, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, start_new_session=True,
    )
    try:
        deadline = time.time() + 600
        while len(list(interrupted.glob("*.fields"))) < 2:
            assert solver.poll() is None, "solver ended before writing restart fields"
            assert time.time() < deadline, "restart fields were not written by both ranks"
            time.sleep(1)
        time.sleep(3)  # let both writers finish their files
    finally:
        os.killpg(solver.pid, signal.SIGKILL)  # only the job launched above
        solver.wait()

    _set_steps(case_interrupted, 1)
    probe = _mpi_run(case_interrupted, "-r")
    saved_step = int(re.search(r"Last executed time step\s+(\d+)", probe.stdout + probe.stderr).group(1))
    steps = saved_step + 2000

    _set_steps(case_reference, steps)
    uninterrupted = _mpi_run(case_reference)
    assert uninterrupted.returncode == 0, uninterrupted.stdout[-2000:]

    _set_steps(case_interrupted, steps)
    resumed = _mpi_run(case_interrupted, "-r")
    assert resumed.returncode == 0, resumed.stdout[-2000:]

    # the far-field log and the point probes are deterministic; the energy file is written on wall-clock intervals
    deterministic = list(reference.glob("*_point_probe_*_tm.dat")) + list(reference.glob("*farfield_log*.dat"))
    assert len(deterministic) == 4, [p.name for p in deterministic]
    for expected in deterministic:
        assert filecmp.cmp(expected, interrupted / expected.name, shallow=False), expected.name

    movies = [p for p in reference.iterdir() if p.is_dir() and "_movie_" in p.name]
    assert len(movies) == 2, [p.name for p in movies]
    for movie in movies:
        resumed_h5 = interrupted / movie.name / f"{movie.name}.h5"
        assert resumed_h5.is_file()
        _assert_resumed_frames_match_reference(resumed_h5, movie / f"{movie.name}.h5")
