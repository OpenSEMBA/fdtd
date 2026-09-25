# %% Setup
from pathlib import Path
import os
import sys

import matplotlib.pyplot as plt
import numpy as np


CASE_DIR = Path(__file__).resolve().parent
REPOSITORY_DIR = CASE_DIR.parents[2]
sys.path.insert(0, str(REPOSITORY_DIR / 'src_pyWrapper'))

from pyWrapper import ExcitationFile, FDTD, Probe


SEMBA_EXE = Path(
    os.environ.get('SEMBA_FDTD_EXECUTABLE', REPOSITORY_DIR / 'build/bin/semba-fdtd')
).resolve()
if not SEMBA_EXE.is_file():
    raise FileNotFoundError(
        'Set SEMBA_FDTD_EXECUTABLE or build the solver at ' f'{SEMBA_EXE}'
    )

RUN_DIR = CASE_DIR / 'run'
RUN_DIR.mkdir(exist_ok=True)

CASES = (
    ('Staircased', CASE_DIR / 'coaxial_staircased.fdtd.json'),
    ('Conformal', CASE_DIR / 'coaxial_conformal.fdtd.json'),
)
FREQUENCIES = np.geomspace(1e6, 100e6, 101)
RESISTANCE_BAND = (1e6, 10e6)
INDUCTANCE_BAND = (10e6, 100e6)


# %% Helpers
def dtft(signal, time, frequencies):
    """Approximate the continuous Fourier transform by trapezoidal integration."""
    signal = np.asarray(signal)
    time = np.asarray(time)
    transform = np.empty_like(frequencies, dtype=complex)
    trapezoid = getattr(np, 'trapezoid', np.trapz)
    for index, frequency in enumerate(frequencies):
        transform[index] = trapezoid(
            signal * np.exp(-2j * np.pi * frequency * time), time
        )
    return transform


def in_frequency_band(frequencies, band):
    return (frequencies >= band[0]) & (frequencies <= band[1])


# %% Run both cases and read the input current
measurements = {}
for label, input_filename in CASES:
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=RUN_DIR,
    )
    solver.cleanUp()
    solver.run()

    probe_filenames = solver.getSolvedProbeFilenames('Input current')
    if len(probe_filenames) != 1:
        raise RuntimeError(
            f'Expected one Input current output for {label}, got {probe_filenames}'
        )
    probe = Probe(probe_filenames[0])
    measurements[label] = {
        'time': probe['time'].to_numpy(),
        'current': probe['current'].to_numpy(),
    }


# %% Read the excitation voltage
excitation = ExcitationFile(CASE_DIR / 'gauss.exc')
voltage_time = excitation['time'].to_numpy()
voltage = excitation['value'].to_numpy()


# %% Time-domain comparison
figure, (voltage_axis, current_axis) = plt.subplots(2, 1, sharex=True, figsize=(9, 7))
voltage_axis.plot(voltage_time * 1e9, voltage, color='black', label='Generator voltage')
voltage_axis.set_ylabel('Voltage [V]')
voltage_axis.grid()
voltage_axis.legend()

for label, measurement in measurements.items():
    current_axis.plot(measurement['time'] * 1e9, measurement['current'], label=label)
current_axis.set_xlabel('Time [ns]')
current_axis.set_ylabel('Current [A]')
current_axis.grid()
current_axis.legend()
figure.tight_layout()


# %% Frequency-domain input impedance and R/L estimates
resistance_mask = in_frequency_band(FREQUENCIES, RESISTANCE_BAND)
inductance_mask = in_frequency_band(FREQUENCIES, INDUCTANCE_BAND)
angular_frequency = 2 * np.pi * FREQUENCIES[inductance_mask]

figure, components_axis = plt.subplots(figsize=(9, 4.5))
estimates = {}
for label, measurement in measurements.items():
    voltage_at_current_time = np.interp(
        measurement['time'], voltage_time, voltage
    )
    impedance = dtft(voltage_at_current_time, measurement['time'], FREQUENCIES) / dtft(
        measurement['current'], measurement['time'], FREQUENCIES
    )
    resistance = np.mean(np.real(impedance[resistance_mask]))
    inductance = np.dot(angular_frequency, np.imag(impedance[inductance_mask])) / np.dot(
        angular_frequency, angular_frequency
    )
    estimates[label] = (resistance, inductance)

    components_axis.semilogx(
        FREQUENCIES * 1e-6, np.real(impedance), label=f'{label}: Re(Z_in)'
    )
    components_axis.semilogx(
        FREQUENCIES * 1e-6, np.imag(impedance), '--', label=f'{label}: Im(Z_in)'
    )
    components_axis.semilogx(
        FREQUENCIES[inductance_mask] * 1e-6,
        angular_frequency * inductance,
        ':',
        label=f'{label}: 2πfL fit',
    )

components_axis.set_xlabel('Frequency [MHz]')
components_axis.set_ylabel('Z_in [ohm]')
components_axis.grid(which='both')
components_axis.legend(ncol=3)
figure.tight_layout()

for label, (resistance, inductance) in estimates.items():
    print(
        f'{label}: R = {resistance:.6g} ohm '
        f'({RESISTANCE_BAND[0] * 1e-6:g}-{RESISTANCE_BAND[1] * 1e-6:g} MHz), '
        f'L = {inductance * 1e9:.6g} nH '
        f'({INDUCTANCE_BAND[0] * 1e-6:g}-{INDUCTANCE_BAND[1] * 1e-6:g} MHz)'
    )

# %%
