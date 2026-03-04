# %%
import numpy as np
import sys, os
from numpy.fft import *
import matplotlib.pyplot as plt
import scipy.constants
from skrf.media import Freespace
from skrf.frequency import Frequency

sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../publico/', 'src_pyWrapper'))
SEMBA_EXE = '/home/luis/ugrfdtd/build-intel-rls/bin/ugrfdtd'
from pyWrapper import *

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
np.savetxt('gauss.exc', data)
 
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
solver = FDTD(
	input_filename = 'towel_rack_with_shorting_plane.fdtd.json', 
	path_to_exe=SEMBA_EXE
)
solver.cleanUp()
solver.run()

# %% Postprocessing
case_dir = os.path.dirname(__file__)

excitation_filename = os.path.join(case_dir, 'gauss.exc')
exc = ExcitationFile(excitation_filename)

time_exc = exc["time"].to_numpy()
excitation = exc["value"].to_numpy()


def compute_impedance(probe_path: str):
	"""Compute time signal and input impedance Z_in(f) for a given probe file."""
	probe = Probe(os.path.join(case_dir, probe_path))
	time_I = probe["time"].to_numpy()
	current = probe["current_0"].to_numpy()

	# Interpolate excitation voltage onto this current probe time grid
	V_interp = np.interp(time_I, time_exc, excitation)

	# DTFT of current and voltage on this time grid
	I_f = dtft(current, time_I, freqs)
	V_f = dtft(V_interp, time_I, freqs)

	Z_in = V_f / I_f
	return time_I, current, Z_in

# DTFT 
f_min = 1e3
f_max = 1e9
points_per_decade = 10
num_points = int((np.log10(f_max) - np.log10(f_min)) * points_per_decade) + 1
freqs = np.logspace(np.log10(f_min), np.log10(f_max), num=num_points)


def dtft(signal, time, freqs):
	"""Continuous-time DTFT approximation using trapezoidal integration."""
	signal = np.asarray(signal)
	time = np.asarray(time)
	res = np.empty_like(freqs, dtype=complex)
	for i, f in enumerate(freqs):
		res[i] = np.trapz(signal * np.exp(-1j * 2 * np.pi * f * time), time)
	return res

# Compute impedances for both simulations
# probe_without = 'towel_rack_without_shorting_plane.fdtd_Wire probe_Wx_72_44_44_s33.dat'
probe_with = 'towel_rack_with_shorting_plane.fdtd_Wire probe_Cable_I_65_30_30.dat'

# time_I_wo, current_wo, Z_in_wo = compute_impedance(probe_without)
time_I_w, current_w, Z_in_w = compute_impedance(probe_with)

#  Time-domain comparison
plt.figure()
plt.plot(time_exc, excitation, label="Excitation")
plt.xlabel("Time [s]")
plt.ylabel("Excitation [V]")
plt.grid(which="both")
plt.legend()
plt.tight_layout()

plt.figure()
plt.plot(time_I_w, current_w, label="Current with shorting plane")
plt.xlabel("Time [s]")
plt.ylabel("Current [A]")
plt.grid(which="both")
plt.legend()
plt.tight_layout()


# %% Impedance comparison (magnitude and phase)
plt.figure()
plt.semilogx(freqs, 20 * np.log10(np.abs(Z_in_w)), label="With shorting plane")
plt.ylabel("|Z(j2πf)| [dB]")
plt.grid(which="both")
plt.legend()

plt.figure()
plt.semilogx(freqs, np.angle(Z_in_w), label="With shorting plane")
plt.ylabel("Phase Z(j2πf)")
plt.grid(which="both")
plt.legend()
plt.show()

print("=== END ===")