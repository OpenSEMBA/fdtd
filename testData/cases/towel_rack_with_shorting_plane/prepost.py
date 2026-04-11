# %%
import numpy as np
import sys, os
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants
from skrf.media import Freespace
from skrf.frequency import Frequency

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../', 'src_pyWrapper'))
SEMBA_EXE = '/home/luis/ugrfdtd/publico/build-rls/bin/semba-fdtd'
from pyWrapper import *

case_dir = os.path.dirname(os.path.abspath(__file__))


# Source Intel oneAPI compiler + MPI environment variables into the current
# process so that any subprocess (e.g. the FDTD solver) inherits them.
INTEL_SETVARS = '/opt/intel/oneapi/setvars.sh'
import subprocess as _sp
_result = _sp.run(
    ['bash', '-c',
     f'source {INTEL_SETVARS} --force compiler mpi > /dev/null 2>&1 && env -0'],
    capture_output=True, text=True
)
for _line in _result.stdout.split('\0'):
    if '=' in _line:
        _k, _v = _line.split('=', 1)
        os.environ[_k] = _v



def dtft(signal, time, freqs):
	"""Continuous-time DTFT approximation using trapezoidal integration."""
	signal = np.asarray(signal)
	time = np.asarray(time)
	res = np.empty_like(freqs, dtype=complex)
	trapz = getattr(np, "trapezoid", getattr(np, "trapz", None)) 
	for i, f in enumerate(freqs):
		res[i] = trapz(signal * np.exp(-1j * 2 * np.pi * f * time), time)
	return res

def compute_impedance(probe_path, time_exc, voltage_exc, freqs):
	probe = Probe(probe_path)
	time_I = probe["time"].to_numpy()
	current = probe["current_0"].to_numpy()
	V_interp = np.interp(time_I, time_exc, voltage_exc)
	I_f = dtft(current, time_I, freqs)
	V_f = dtft(V_interp, time_I, freqs)
	Z = V_f / I_f
	return time_I, current, Z

# %% Generate excitation and visualize
dt = 1e-12
w0 = 0.1e-8    # ~ 200 MHz bandwidth
t0 = 10*w0
t = np.arange(0, t0+20*w0, dt)
f = np.exp( -np.power(t-t0,2)/ w0**2 )

plt.figure()
plt.plot(t,f)

data = np.zeros((len(t), 2))
data[:,0] = t
data[:,1] = f
np.savetxt(os.path.join(case_dir, 'gauss.exc'), data)
 
fq = fftfreq(len(t))/dt
F = fft(f)
F = F[fq>=0]
fq = fq[fq>=0]
Fmax = np.max(np.abs(F))
Fcut = Fmax / np.sqrt(2.0)
idx = (np.abs(np.abs(F) - Fcut)).argmin()

plt.figure()
plt.plot(fq, np.abs(F),'.-')
plt.axvline(x = fq[idx], color = 'r')
plt.xlim(0,1e9)
plt.ylim(1e-9, 1e6)
plt.yscale('log')

plt.grid()

# %% Run simulation
CASE_NAME = "towel_rack_with_shorting_plane.fdtd.json"
FOLDER_WITH_SHORTING_PLANE = "with_shorting_plane"
FOLDER_WITHOUT_SHORTING_PLANE = "without_shorting_plane"

os.makedirs(os.path.join(case_dir, FOLDER_WITH_SHORTING_PLANE), exist_ok=True)
solver = FDTD(
	input_filename = os.path.join(case_dir, CASE_NAME),
	path_to_exe=SEMBA_EXE,
	run_in_folder=os.path.join(case_dir, FOLDER_WITH_SHORTING_PLANE)
)
solver.cleanUp()
solver.run()

os.makedirs(os.path.join(case_dir, FOLDER_WITHOUT_SHORTING_PLANE), exist_ok=True)
solver = FDTD(
	input_filename = os.path.join(case_dir, CASE_NAME),
	path_to_exe=SEMBA_EXE,
	run_in_folder=os.path.join(case_dir, FOLDER_WITHOUT_SHORTING_PLANE)
)
solver['materialAssociations'][0]['elementIds'] = [1]
solver.cleanUp()
solver.run()
PROBE_NAME = solver.getSolvedProbeFilenames("Wire probe")[0]

# %% Postprocessing
excitation_filename = os.path.join(case_dir, 'gauss.exc')
exc = ExcitationFile(excitation_filename)

time_exc = exc["time"].to_numpy()
voltage_exc = exc["value"].to_numpy()

freqs = np.geomspace(1e3, 1e9, 61)


probe_with_path    = os.path.join(case_dir, FOLDER_WITH_SHORTING_PLANE,    PROBE_NAME)
probe_without_path = os.path.join(case_dir, FOLDER_WITHOUT_SHORTING_PLANE, PROBE_NAME)

time_I_w,  current_w,  Z_in_w  = compute_impedance(probe_with_path,    time_exc, voltage_exc, freqs)
time_I_wo, current_wo, Z_in_wo = compute_impedance(probe_without_path, time_exc, voltage_exc, freqs)

# %% Time-domain comparison
plt.figure()
plt.plot(time_exc, voltage_exc, label="Excitation")
plt.xlabel("Time [s]")
plt.ylabel("Excitation [V]")
plt.grid(which="both")
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(time_I_w,  current_w,  label="With shorting plane")
plt.plot(time_I_wo, current_wo, label="Without shorting plane")
plt.xlabel("Time [s]")
plt.ylabel("Current [A]")
plt.grid(which="both")
plt.legend()
plt.tight_layout()

# %% Impedance comparison
plt.figure()
plt.semilogx(freqs, 20 * np.log10(np.abs(Z_in_w)),  label="With shorting plane")
plt.semilogx(freqs, 20 * np.log10(np.abs(Z_in_wo)), label="Without shorting plane")
plt.xlabel("Frequency [Hz]")
plt.ylabel("|Z(j2πf)| [dB]")
plt.grid(which="both")
plt.legend()
plt.tight_layout()

print("=== END ===")