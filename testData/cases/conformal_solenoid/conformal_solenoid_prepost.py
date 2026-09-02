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
    ('Solenoid', CASE_DIR / 'solenoid.fdtd.json', 'BC'),
    ('45 deg staircased', CASE_DIR / 'solenoid_45deg_with_staircased.fdtd.json', 'Bulk probe'),
    ('45 deg conformal', CASE_DIR / 'solenoid_45deg_with_conformal.fdtd.json', 'Bulk probe'),
)


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


def input_impedance(time, current, voltage_time, voltage):
    voltage_at_current_time = np.interp(time, voltage_time, voltage)
    return dtft(voltage_at_current_time, time, FREQUENCIES) / dtft(
        current, time, FREQUENCIES
    )


# %% Run all cases in one output folder
measurements = {}
for label, input_filename, probe_name in CASES:
    solver = FDTD(
        input_filename=input_filename,
        path_to_exe=SEMBA_EXE,
        run_in_folder=RUN_DIR,
    )
    solver.cleanUp()
    solver.run()

    probe = Probe(solver.getSolvedProbeFilenames(probe_name)[0])
    measurements[label] = {
        'time': probe['time'].to_numpy(),
        # The bulk-current probe normal is opposite to source-current direction.
        'current': -probe['current'].to_numpy(),
    }


# %% Time-domain comparison
excitation = ExcitationFile(CASE_DIR / 'predefinedExcitation.1.exc')
voltage_time = excitation['time'].to_numpy()
voltage = excitation['value'].to_numpy()

figure, (voltage_axis, current_axis) = plt.subplots(2, 1, sharex=True, figsize=(9, 7))
voltage_axis.plot(voltage_time * 1e9, voltage, color='black', label='Generator voltage')
voltage_axis.set_ylabel('Voltage [V]')
voltage_axis.grid()
voltage_axis.legend()
voltage_axis.set_xlim(0, 4)
voltage_axis.plot(voltage_time * 1e9, voltage, color='black', label='Generator voltage')

for label, measurement in measurements.items():
    current_axis.plot(measurement['time'] * 1e9, measurement['current'], label=label)
current_axis.set_xlabel('Time [ns]')
current_axis.set_ylabel('Current [A]')
current_axis.grid()
current_axis.legend()
figure.tight_layout()


# %% Frequency-domain input impedance and R/L estimates
FREQUENCIES = np.geomspace(1e6, 100e6, 101)
RESISTANCE_BAND = (1e6, 10e6)
INDUCTANCE_BAND = (10e6, 100e6)


def in_frequency_band(frequencies, band):
    return (frequencies >= band[0]) & (frequencies <= band[1])


figure, components_axis = plt.subplots(figsize=(9, 4.5))
estimates = {}
for label, measurement in measurements.items():
    impedance = input_impedance(
        measurement['time'], measurement['current'], voltage_time, voltage
    )
    resistance_mask = in_frequency_band(FREQUENCIES, RESISTANCE_BAND)
    inductance_mask = in_frequency_band(FREQUENCIES, INDUCTANCE_BAND)
    resistance = np.mean(np.real(impedance[resistance_mask]))
    angular_frequency = 2 * np.pi * FREQUENCIES[inductance_mask]
    inductance = np.dot(angular_frequency, np.imag(impedance[inductance_mask])) / np.dot(
        angular_frequency, angular_frequency
    )
    estimates[label] = (resistance, inductance)

    components_axis.semilogx(
        FREQUENCIES * 1e-6, np.real(impedance), label=f'{label}: Re(Z)'
    )
    components_axis.semilogx(
        FREQUENCIES * 1e-6, np.imag(impedance), '--', label=f'{label}: Im(Z)'
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
