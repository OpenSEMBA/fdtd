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
    ('Conformal capacitor', CASE_DIR / 'capacitor.fdtd.json', 'Wire probe'),
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

    probe = Probe(solver.getSolvedProbeFolders(probe_name)[0])
    measurements[label] = {
        'time': probe['time'].to_numpy(),
        # The bulk-current probe normal is opposite to source-current direction.
        'current': -probe['current'].to_numpy(),
    }


# %% Time-domain comparison
excitation = ExcitationFile(CASE_DIR / 'gauss_1ghz.exc')
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


# %% Frequency-domain input impedance and R/C estimates
FREQUENCIES = np.geomspace(1e6, 1e9, 201)
CUTOFF_FREQUENCY = 100e6
RESISTANCE_BAND = (300e6, 1e9)  # resistor-dominated plateau, Re(Z) ~ R
CAPACITANCE_BAND = (1e6, 10e6)  # capacitor-dominated region, Im(Z) ~ -1/(2*pi*f*C)


def in_frequency_band(frequencies, band):
    return (frequencies >= band[0]) & (frequencies <= band[1])


figure, (magnitude_axis, phase_axis) = plt.subplots(2, 1, sharex=True, figsize=(9, 7))
estimates = {}
for label, measurement in measurements.items():
    impedance = input_impedance(
        measurement['time'], measurement['current'], voltage_time, voltage
    )
    resistance_mask = in_frequency_band(FREQUENCIES, RESISTANCE_BAND)
    capacitance_mask = in_frequency_band(FREQUENCIES, CAPACITANCE_BAND)
    resistance = np.mean(np.real(impedance[resistance_mask]))
    angular_frequency = 2 * np.pi * FREQUENCIES[capacitance_mask]
    # Least-squares fit of Im(Z) = -1/(omega*C) against omega.
    capacitance = -np.dot(angular_frequency, angular_frequency) / np.dot(
        angular_frequency, np.imag(impedance[capacitance_mask])
    )
    cutoff = 1 / (2 * np.pi * resistance * capacitance)
    estimates[label] = (resistance, capacitance, cutoff)

    magnitude_axis.loglog(FREQUENCIES * 1e-6, np.abs(impedance), label=label)
    phase_axis.semilogx(
        FREQUENCIES * 1e-6, np.degrees(np.angle(impedance)), label=label
    )

magnitude_axis.axhline(
    estimates[list(estimates)[0]][0], color='gray', linestyle='--', label='Fitted R'
)
for axis in (magnitude_axis, phase_axis):
    axis.axvline(
        CUTOFF_FREQUENCY * 1e-6, color='gray', linestyle=':', label='100 MHz cutoff'
    )
    axis.grid(which='both')
    axis.legend()

magnitude_axis.set_ylabel('|Z_in| [ohm]')
phase_axis.set_xlabel('Frequency [MHz]')
phase_axis.set_ylabel('phase(Z_in) [deg]')
figure.tight_layout()

for label, (resistance, capacitance, cutoff) in estimates.items():
    print(
        f'{label}: R = {resistance:.6g} ohm '
        f'({RESISTANCE_BAND[0] * 1e-6:g}-{RESISTANCE_BAND[1] * 1e-6:g} MHz), '
        f'C = {capacitance * 1e12:.6g} pF '
        f'({CAPACITANCE_BAND[0] * 1e-6:g}-{CAPACITANCE_BAND[1] * 1e-6:g} MHz), '
        f'cutoff = {cutoff * 1e-6:.6g} MHz'
    )

# %%
