# %% [markdown]
# # semba-fdtd Performance Benchmark (Cybonera)
#
# Runs the three cybonera cases via `pyWrapper.FDTD` and writes a single
# Markdown report with command, timing and memory usage.

# %% Configuration and setup
import os
import sys
import glob
import datetime
import threading
import time
import subprocess

try:
    import psutil
    _HAVE_PSUTIL = True
except ImportError:
    _HAVE_PSUTIL = False


def _find_repo_root(start_dir):
    """Find repo root by walking up until src_pyWrapper/pyWrapper.py exists."""
    cur = os.path.abspath(start_dir)
    while True:
        candidate = os.path.join(cur, 'src_pyWrapper', 'pyWrapper.py')
        if os.path.isfile(candidate):
            return cur
        parent = os.path.dirname(cur)
        if parent == cur:
            return os.path.abspath(start_dir)
        cur = parent


try:
    _BASE_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _BASE_DIR = os.getcwd()

ROOT = _find_repo_root(_BASE_DIR)
PYWRAPPER_DIR = os.path.join(ROOT, 'src_pyWrapper')
if PYWRAPPER_DIR not in sys.path:
    sys.path.append(PYWRAPPER_DIR)

import pyWrapper

SEMBA_EXE = os.path.join(ROOT, 'build', 'bin', 'semba-fdtd')
FLAGS = []
MPI_COMMAND = None
BENCHMARK_NAME = 'cybonera'
CYBONERA_DIR = os.path.join(ROOT, 'tmp_cases', 'cybonera')


def _fmt_memory(bytes_val):
    if bytes_val is None:
        return 'not measured (psutil not available)'
    for unit, threshold in [('GiB', 1024**3), ('MiB', 1024**2), ('KiB', 1024)]:
        if bytes_val >= threshold:
            return f"{bytes_val / threshold:.2f} {unit} ({bytes_val:,} bytes)"
    return f"{bytes_val} bytes"


def _fmt_time(seconds):
    h = int(seconds // 3600)
    m = int((seconds % 3600) // 60)
    s = seconds % 60
    parts = []
    if h:
        parts.append(f"{h} h")
    if m:
        parts.append(f"{m} min")
    parts.append(f"{s:.3f} s")
    return ' '.join(parts)


def _poll_memory(proc, interval, results):
    peak = 0
    try:
        ps = psutil.Process(proc.pid)
        while proc.poll() is None:
            try:
                rss = ps.memory_info().rss
                for child in ps.children(recursive=True):
                    try:
                        rss += child.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                if rss > peak:
                    peak = rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
            time.sleep(interval)
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        pass
    results['peak_rss'] = peak


def _run_with_peak_memory(command, cwd=None):
    peak_result = {'peak_rss': None}
    proc = subprocess.Popen(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )

    mem_thread = None
    if _HAVE_PSUTIL:
        mem_thread = threading.Thread(
            target=_poll_memory, args=(proc, 0.05, peak_result), daemon=True
        )
        mem_thread.start()

    stdout_bytes, stderr_bytes = proc.communicate()
    if mem_thread is not None:
        mem_thread.join(timeout=2)

    cp = subprocess.CompletedProcess(
        args=command,
        returncode=proc.returncode,
        stdout=stdout_bytes,
        stderr=stderr_bytes,
    )
    return cp, peak_result['peak_rss']


CASE_FILES = sorted(glob.glob(os.path.join(CYBONERA_DIR, '*.fdtd.json')))

if not os.path.isfile(SEMBA_EXE):
    raise FileNotFoundError(f"semba-fdtd binary not found: {SEMBA_EXE}")

print(f"Root       : {ROOT}")
print(f"Binary     : {SEMBA_EXE}")
print(f"Cases found: {len(CASE_FILES)}")
for case in CASE_FILES:
    print(f"  - {case}")


# %% Run cases with pyWrapper
results = []

for case in CASE_FILES:
    case_abs = os.path.abspath(case)
    case_name = os.path.basename(case_abs)

    print(f"\nRunning {case_name} ...")

    solver = pyWrapper.FDTD(
        input_filename=case_abs,
        path_to_exe=SEMBA_EXE,
        flags=FLAGS,
        mpi_command=MPI_COMMAND,
    )

    solver.cleanUp()

    patched_state = {'peak_rss': None}
    original_run = pyWrapper.subprocess.run

    def _patched_run(cmd, capture_output=False):
        completed, peak = _run_with_peak_memory(cmd)
        patched_state['peak_rss'] = peak
        return completed

    pyWrapper.subprocess.run = _patched_run

    t_start = time.perf_counter()
    status = 'ok'
    error_message = ''
    return_code = 0
    try:
        solver.run()
        if hasattr(solver, 'output') and solver.output is not None:
            return_code = solver.output.returncode
    except Exception as exc:
        status = 'failed'
        error_message = str(exc)
        if hasattr(solver, 'output') and solver.output is not None:
            return_code = solver.output.returncode
        else:
            return_code = -1
    finally:
        pyWrapper.subprocess.run = original_run
    t_end = time.perf_counter()

    wall_time_s = t_end - t_start
    peak_rss = patched_state['peak_rss']
    full_command = ' '.join(solver.run_command)

    print(
        f"  status={status} time={wall_time_s:.2f}s "
        f"peak={_fmt_memory(peak_rss)}"
    )

    results.append({
        'case': case_abs,
        'binary': SEMBA_EXE,
        'flags': ' '.join(FLAGS) if FLAGS else '(none)',
        'mpi_command': MPI_COMMAND if MPI_COMMAND else '(none)',
        'command': full_command,
        'wall_time_s': wall_time_s,
        'peak_rss': peak_rss,
        'return_code': return_code,
        'status': status,
        'error': error_message,
    })


# %% Generate consolidated markdown report
total_time = sum(r['wall_time_s'] for r in results)
max_peak = max(
    (r['peak_rss'] for r in results if r['peak_rss'] is not None),
    default=None,
)

report_lines = [
    '# semba-fdtd Performance Benchmark Report (cybonera)',
    '',
    f"**Generated:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
    '',
    '## Summary',
    '',
    '| Metric | Value |',
    '|--------|-------|',
    f"| Number of cases | {len(results)} |",
    f"| Total wall-clock time | {_fmt_time(total_time)} |",
    f"| Max peak memory (RSS) across cases | {_fmt_memory(max_peak)} |",
    '',
    '## Per-case results',
    '',
    '| Case | Binary and options | Total execution time | Memory consumption (peak RSS) | Exit code | Status |',
    '|------|---------------------|----------------------|-------------------------------|-----------|--------|',
]

for r in results:
    report_lines.append(
        f"| `{r['case']}` | `{r['command']}` | {_fmt_time(r['wall_time_s'])} | "
        f"{_fmt_memory(r['peak_rss'])} | {r['return_code']} | {r['status']} |"
    )

for r in results:
    if r['status'] != 'ok':
        report_lines += [
            '',
            f"### Error details: `{os.path.basename(r['case'])}`",
            '',
            '```text',
            r['error'] if r['error'] else '(no message)',
            '```',
        ]

report_md = '\n'.join(report_lines)
print(report_md)


# %% Save report
report_filename = BENCHMARK_NAME + '_benchmark_report.md'
with open(report_filename, 'w') as f:
    f.write(report_md + '\n')

print(f"\nReport saved to: {os.path.abspath(report_filename)}")
