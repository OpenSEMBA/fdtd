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

EPSILON_0 = 8.8541878128e-12
INPUT_FILENAME = CASE_DIR / 'capacitor_charge.fdtd.json'
PROBE_NAME = 'Point probe'


# %% Helpers
def plate_geometry(solver):
    """Derive plate area and gap from the two constant-z 'cell' surfaces in the mesh."""
    steps = solver['mesh']['grid']['steps']
    dx, dy, dz = steps['x'][0], steps['y'][0], steps['z'][0]

    plate_z_indices = []
    area = None
    for element in solver['mesh']['elements']:
        if element.get('type') != 'cell':
            continue
        for (ax, ay, az), (bx, by, bz) in element['intervals']:
            if az != bz:
                continue
            plate_z_indices.append(az)
            area = abs(bx - ax) * dx * abs(by - ay) * dy

    gap = abs(plate_z_indices[0] - plate_z_indices[1]) * dz
    return area, gap


# %% Run the case
solver = FDTD(
    input_filename=CASE_DIR / 'capacitor_charge.fdtd.json',
    path_to_exe=SEMBA_EXE,
    flags='-mapvtk -ignoresamplingerrors',
    run_in_folder=RUN_DIR,
)
solver.cleanUp()
solver.run()

ez_probe = next(
    probe
    for probe in map(Probe, solver.getSolvedProbeFolders(PROBE_NAME))
    if probe.field == 'E' and probe.direction == 'z'
)
time = ez_probe['time'].to_numpy()
ez = ez_probe['field'].to_numpy()


# %% Injected charge from the source's magnitude file
excitation = ExcitationFile(CASE_DIR / 'predefinedExcitation.1.exc')
current_time = excitation['time'].to_numpy()
current_value = excitation['value'].to_numpy()
current_at_probe_time = np.interp(time, current_time, current_value)
trapezoid = getattr(np, 'trapezoid', np.trapz)
charge = trapezoid(current_at_probe_time, time)


# %% Compare measured and theoretical capacitance
area, gap = plate_geometry(solver)
ez_final = ez[-1]
voltage = ez_final * gap
measured_capacitance = charge / np.abs(voltage)
theoretical_capacitance = EPSILON_0 * area / gap
relative_error = abs(measured_capacitance - theoretical_capacitance) / theoretical_capacitance

print(
    f'Measured C    = {measured_capacitance * 1e12:.6g} pF\n'
    f'Theoretical C = {theoretical_capacitance * 1e12:.6g} pF '
    f'(eps0 * {area * 1e6:.6g} mm^2 / {gap * 1e3:.6g} mm)\n'
    f'Relative error = {relative_error * 100:.3g} %'
)


# %% Time-domain diagnostics
figure, (current_axis, field_axis) = plt.subplots(2, 1, sharex=True, figsize=(9, 7))
current_axis.plot(current_time * 1e9, current_value, color='black')
current_axis.set_ylabel('Injected current [A]')
current_axis.grid()

field_axis.plot(time * 1e9, ez)
field_axis.set_xlabel('Time [ns]')
field_axis.set_ylabel('Ez [V/m]')
field_axis.grid()
field_axis.legend()
figure.tight_layout()

# %%
