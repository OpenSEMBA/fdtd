from test.utils.utils import *
import re
import signal
import sys
import time
from pathlib import Path

import numpy as np
from src_pyWrapper.pyWrapper import FDTD

MOVIES: tuple[str, ...] = (
    "electric_field_movie_ME_0_0_0__29_29_9",
    "electric_field_movie_low_ME_0_0_0__29_29_1",
    "electric_field_movie_high_ME_0_0_8__29_29_9",
)


def _fdtd(case: str | Path, flags: list[str] | None = None, mpi_procs: int | None = None) -> FDTD:
    mpi_command = f"mpirun -np {mpi_procs}" if mpi_procs else None
    return FDTD(case, path_to_exe=SEMBA_EXE, flags=flags, mpi_command=mpi_command)


def _output(solver: FDTD) -> str:
    return (solver.output.stdout + solver.output.stderr).decode()


def _saved_step(solver: FDTD) -> int:
    # a resume cannot end before the saved step, so the run fails and reports the checkpointed step
    with pytest.raises(AssertionError):
        solver.run()
    return int(re.search(r"Last executed time step\s+(\d+)", _output(solver)).group(1))


def _is_number(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def _assert_probe_file_close(expected: Path, actual: Path, peak_tolerance: float = 1e-5) -> None:
    # text probe files: header lines must match exactly, numeric rows within a tolerance scaled to the peak
    lines_expected = open(expected).read().splitlines()
    lines_actual = open(actual).read().splitlines()
    assert len(lines_expected) == len(lines_actual), expected.name
    rows_expected: list[list[float]] = []
    rows_actual: list[list[float]] = []
    for line_expected, line_actual in zip(lines_expected, lines_actual):
        tokens = line_expected.split()
        if tokens and all(_is_number(t) for t in tokens):
            rows_expected.append([float(t) for t in tokens])
            rows_actual.append([float(t) for t in line_actual.split()])
        else:
            assert line_expected == line_actual, expected.name
    if rows_expected:
        ref = np.array(rows_expected)
        atol = peak_tolerance * np.max(np.abs(ref))
        np.testing.assert_allclose(np.array(rows_actual), ref, rtol=1e-5, atol=atol, err_msg=expected.name)


@pytest.mark.slow
@pytest.mark.hdf
@pytest.mark.movie
def test_resume_after_sigkill_keeps_previous_movies(tmp_path: Path) -> None:
    case_folder = os.path.join(CASES_FOLDER, "planewave")
    for name in ("pw-in-box-with-partial-movies.fdtd.json", "gauss_1GHz.exc"):
        shutil.copy(os.path.join(case_folder, name), tmp_path)
    case = str(tmp_path / "case.fdtd.json")
    os.rename(tmp_path / "pw-in-box-with-partial-movies.fdtd.json", case)

    # far more steps than can finish: the run is only ever stopped by the kill
    set_steps(case, 100_000_000)
    solver = _fdtd(case, ["-flush", "1"]).start()  # shortest restart-fields flush interval, one minute
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

    set_steps(case, 1)
    saved_step = _saved_step(_fdtd(case, ["-r"]))

    set_steps(case, saved_step + 200)
    _fdtd(case, ["-r"]).run()

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
def test_resume_after_sigkill_reproduces_uninterrupted_probes(tmp_path: Path) -> None:
    reference = tmp_path / "reference"
    interrupted = tmp_path / "interrupted"
    case_folder = os.path.join(CASES_FOLDER, "planewave")
    for folder in (reference, interrupted):
        folder.mkdir()
        for name in ("pw-in-box.fdtd.json", "gauss_1GHz.exc"):
            shutil.copy(os.path.join(case_folder, name), folder)
    case_reference = str(reference / "pw-in-box.fdtd.json")
    case_interrupted = str(interrupted / "pw-in-box.fdtd.json")

    set_steps(case_interrupted, 100_000_000)
    solver = _fdtd(case_interrupted, ["-flush", "1"]).start()
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

    set_steps(case_interrupted, 1)
    saved_step = _saved_step(_fdtd(case_interrupted, ["-r", "-n", "1"]))
    steps = saved_step + 2000

    uninterrupted = _fdtd(case_reference, ["-n", str(steps)])
    uninterrupted.run()
    _fdtd(case_interrupted, ["-r", "-n", str(steps)]).run()

    # Energy is reported on wall-clock intervals, so only the deterministic outputs are compared
    for name in ("before_Ex_3_3_1_tm", "after_Ex_3_3_5_tm", "inbox_Ex_3_3_3_tm", "inbox_Ex_3_3_3_fq"):
        expected = reference / f"pw-in-box.fdtd_{name}.dat"
        _assert_probe_file_close(expected, interrupted / expected.name)


def _flush_case(
    tmp_path: Path,
    *,
    steps: int = 20,
    minutes: int | None = None,
    additional: str | None = None,
) -> str:
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
def test_restart_flush_priority(
    tmp_path: Path, minutes: int, additional: str | None, flags: list[str], expected: int
) -> None:
    case = _flush_case(tmp_path, minutes=minutes, additional=additional)
    solver = _fdtd(case, flags)
    solver.run()
    assert f"Restart FIELDS flush: every {expected} minutes" in _output(solver)


def test_restart_flush_negative_is_rejected(tmp_path: Path) -> None:
    case = _flush_case(tmp_path, minutes=-1)
    solver = _fdtd(case)
    with pytest.raises(AssertionError):
        solver.run()
    assert "restartFlushMinutes cannot be negative" in _output(solver)


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.slow
@pytest.mark.skipif(sys.platform == "win32", reason="kills the job's process group, POSIX only")
def test_resume_rejects_different_process_count_then_resumes_with_same(tmp_path: Path) -> None:
    case = _flush_case(tmp_path, steps=100_000_000, minutes=1)
    solver = _fdtd(case, mpi_procs=2).start(start_new_session=True)
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

    set_steps(case, 1)
    wrong = _fdtd(case, ["-r"], mpi_procs=1)
    with pytest.raises(AssertionError):
        wrong.run()
    assert "mpirun -n 2" in _output(wrong)  # the saved process count is named

    # same count as saved: the probe reports the checkpointed step, as in the tests above
    saved_step = _saved_step(_fdtd(case, ["-r"], mpi_procs=2))

    set_steps(case, saved_step + 200)
    _fdtd(case, ["-r"], mpi_procs=2).run()


def _assert_resumed_frames_match_reference(
    resumed_h5: Path, reference_h5: Path, peak_tolerance: float = 1e-5
) -> None:
    # every frame the resumed run wrote must match the uninterrupted frame recorded at the same time
    import h5py

    with h5py.File(reference_h5, "r") as ref, h5py.File(resumed_h5, "r") as res:
        index = {round(float(t), 15): i for i, t in enumerate(ref["series/values"][...])}
        ref_frames = ref["attributes/a0001/values"]
        res_frames = res["attributes/a0001/values"]
        for j, t in enumerate(res["series/values"][...]):
            i = index[round(float(t), 15)]
            assert_allclose_to_peak(res_frames[j], ref_frames[i], peak_tolerance,
                                    err_msg=f"frame at t={t} differs after resume")


@no_mpi_skip
@pytest.mark.mpi
@pytest.mark.hdf
@pytest.mark.movie
@pytest.mark.farfield
@pytest.mark.slow
@pytest.mark.skipif(sys.platform == "win32", reason="kills the job's process group, POSIX only")
def test_resume_planewave_sphere_with_farfield_point_probe_and_movies_matches_uninterrupted(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference"
    interrupted = tmp_path / "interrupted"
    case_folder = os.path.join(CASES_FOLDER, "sphere")
    for folder in (reference, interrupted):
        folder.mkdir()
        for name in ("sphere_restart.fdtd.json", "gauss.exc"):
            shutil.copy(os.path.join(case_folder, name), folder)
    case_reference = str(reference / "sphere_restart.fdtd.json")
    case_interrupted = str(interrupted / "sphere_restart.fdtd.json")

    set_steps(case_interrupted, 100_000_000)
    solver = _fdtd(case_interrupted, ["-flush", "1"], mpi_procs=2).start(start_new_session=True)
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

    set_steps(case_interrupted, 1)
    saved_step = _saved_step(_fdtd(case_interrupted, ["-r"], mpi_procs=2))
    steps = saved_step + 2000

    set_steps(case_reference, steps)
    _fdtd(case_reference, mpi_procs=2).run()

    set_steps(case_interrupted, steps)
    _fdtd(case_interrupted, ["-r"], mpi_procs=2).run()

    # the far-field log and the point probes are deterministic; the energy file is written on wall-clock intervals
    deterministic = list(reference.glob("*_point_probe_*_tm.dat")) + list(reference.glob("*farfield_log*.dat"))
    assert len(deterministic) == 4, [p.name for p in deterministic]
    for expected in deterministic:
        _assert_probe_file_close(expected, interrupted / expected.name)

    movies = [p for p in reference.iterdir() if p.is_dir() and "_movie_" in p.name]
    assert len(movies) == 2, [p.name for p in movies]
    for movie in movies:
        resumed_h5 = interrupted / movie.name / f"{movie.name}.h5"
        assert resumed_h5.is_file()
        _assert_resumed_frames_match_reference(resumed_h5, movie / f"{movie.name}.h5")
