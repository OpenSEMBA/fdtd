from utils import *


legacy_exe = os.getenv(
    "SEMBA_LEGACY_EXE",
    os.path.join(os.getcwd(), "build-legacy", "bin", "semba-fdtd"),
)

missing_legacy_exe_skip = pytest.mark.skipif(
    not os.path.isfile(legacy_exe),
    reason="Legacy executable not found. Set SEMBA_LEGACY_EXE or build build-legacy/bin/semba-fdtd",
)


def _run_solver(case_file, exe_path, run_dir):
    solver = FDTD(input_filename=case_file, path_to_exe=exe_path, run_in_folder=run_dir)
    solver.run()
    return solver


def _assert_probe_columns_equal(probe_legacy, probe_new, columns, atol=1e-8, rtol=1e-5):
    legacy_time = probe_legacy["time"].to_numpy()
    new_time = probe_new["time"].to_numpy()

    for col in columns:
        new_vals = np.interp(new_time, legacy_time, probe_legacy[col].to_numpy())
        assert np.allclose(new_vals, probe_new[col].to_numpy(), atol=atol, rtol=rtol)


def _compare_case_probes(case_file, probe_names, columns_by_probe, tmp_path):
    legacy_dir = tmp_path / "legacy"
    new_dir = tmp_path / "new"
    legacy_dir.mkdir()
    new_dir.mkdir()

    solver_legacy = _run_solver(case_file, legacy_exe, legacy_dir)
    solver_new = _run_solver(case_file, SEMBA_EXE, new_dir)

    for probe_name in probe_names:
        legacy_probe = Probe(solver_legacy.getSolvedProbeFilenames(probe_name)[0])
        new_probe = Probe(solver_new.getSolvedProbeFilenames(probe_name)[0])
        _assert_probe_columns_equal(
            legacy_probe,
            new_probe,
            columns_by_probe[probe_name],
        )


@missing_legacy_exe_skip
def test_point_probe_output_equivalence_planewave_in_box(tmp_path):
    case_file = CASES_FOLDER + "planewave/pw-in-box.fdtd.json"
    _compare_case_probes(
        case_file=case_file,
        probe_names=["before", "inbox", "after"],
        columns_by_probe={
            "before": ["field"],
            "inbox": ["field", "incident"],
            "after": ["field"],
        },
        tmp_path=tmp_path,
    )


@mtln_skip
@missing_legacy_exe_skip
def test_wire_probe_output_equivalence_towel_hanger(tmp_path):
    case_file = CASES_FOLDER + "towelHanger/towelHanger.fdtd.json"
    _compare_case_probes(
        case_file=case_file,
        probe_names=["wire_start", "wire_mid", "wire_end"],
        columns_by_probe={
            "wire_start": ["current_0"],
            "wire_mid": ["current_0"],
            "wire_end": ["current_0"],
        },
        tmp_path=tmp_path,
    )
