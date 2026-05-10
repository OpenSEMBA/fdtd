---
name: emc-engineer
description: Writing .fdtd.json inputs, launching simulations via pyWrapper FDTD class, analyzing probe outputs (.dat files) in time/frequency domains, generating excitation signals with NumPy, and common EMC analysis patterns like shielding effectiveness and radar cross section
---

## When to use

- Writing or modifying `.fdtd.json` input files for new simulation cases
- Launching simulations and monitoring progress
- Analyzing probe output data (`.dat` files) in time and frequency domains
- Creating pre/post-processing Python scripts with the pyWrapper interface
- Comparing FDTD results with analytical models or measurements
- Setting up excitation signals (`.exc` files) for time-domain sources

## Key Files to Read

### Core Solver
- `src_main_pub/launcher.F90` — Entry point (15 lines)
- `src_main_pub/semba_fdtd.F90` — Main type: init/launch/end lifecycle
- `src_main_pub/preprocess_geom.F90` — Geometry parsing, mesh building, material assignment
- `src_main_pub/postprocess.F90` — DTFT, transfer functions, resampling
- `src_main_pub/observation.F90` — Probe data collection and storage

### Python Interface
- `src_pyWrapper/pyWrapper.py` — FDTD, Probe, ExcitationFile classes

### Documentation
- `doc/fdtdjson.md` — Complete `.fdtd.json` input format specification (876 lines)
- `doc/tutorials/veritasium/veritasium.md` — End-to-end workflow tutorial
- `doc/mtln.md` — MTLN solver documentation

### Example Cases
- `testData/input_examples/` — 23 `.fdtd.json` files demonstrating various features
- `testData/cases/` — 44 full simulation case directories
- `testData/excitations/` — Source waveform definitions
- `testData/netlists/` — SPICE netlist files for MTLN
- `testData/cases/planewave/pw_prepost.py` — Simple pre/post example
- `testData/cases/sgbcShieldingEffectiveness/sgbc_prepost.py` — EMC analysis with analytical comparison
- `testData/cases/holland/holland_prepost.py` — Wire case comparison

## Simulation Workflow

### Step 1: Write the Input File

The primary input format is `.fdtd.json` with these top-level keys:

| Key | Required | Description |
|-----|----------|-------------|
| `general` | Yes | `timeStep`, `numberOfSteps`, optional `mtlnProblem`, `additionalArguments` |
| `boundary` | No (defaults to `mur` on all sides) | Per-face: `pec`, `pmc`, `periodic`, `mur`, `pml` |
| `mesh` | Yes | Contains `grid`, `coordinates`, `elements` |
| `materials` | No | Material definitions |
| `materialAssociations` | No | Links materials to mesh elements |
| `sources` | No | Electromagnetic sources |
| `probes` | No | Output sensors |

### Mesh Definition

```json
"mesh": {
  "grid": {
    "numberOfCells": [nx, ny, nz],
    "steps": {"x": [dx], "y": [dy], "z": [dz]},
    "origin": [0, 0, 0]
  },
  "coordinates": {
    "name": {"relativePosition": [cellIndex, fractionalX, fractionalY, fractionalZ]}
  },
  "elements": [
    {"type": "node", "position": [x, y, z], "id": materialID},
    {"type": "polyline", "points": [[x1,y1,z1], ...], "id": materialID},
    {"type": "cell", "intervals": [[[ax,ay,az],[bx,by,bz]], ...], "id": materialID}
  ]
}
```

**Interval convention for `cell` elements:** An interval `[[ax,ay,az],[bx,by,bz]]` defines a region `[ax,bx) x [ay,by) x [az,bz)`. Depending on how many dimensions differ:
- 0 dimensions differ = point
- 1 dimension differs = oriented line
- 2 dimensions differ = oriented surface
- 3 dimensions differ = volume

### Step 2: Generate Excitation Signals (if needed)

For time-domain sources, generate `.exc` files using NumPy:

```python
import numpy as np

# Gaussian pulse
dt = 1e-12
t0 = 20 * dt
w0 = 5 * dt
t = np.arange(0, t0 + 20 * w0, dt)
f = np.exp(-(t - t0)**2 / w0**2)
np.savetxt('gauss.exc', np.column_stack((t, f)))

# Sigmoid step
t = np.linspace(0, 10e-8, 2000)
V = A * (sigmoid_raw - offset) * scaling
np.savetxt('./step.exc', np.column_stack((t, V)))
```

### Step 3: Launch the Simulation

**Direct command line:**
```bash
semba-fdtd -i CASE_NAME.fdtd.json
mpirun -n NPROCS semba-fdtd -i CASE_NAME.fdtd.json
```

**Via Python pyWrapper:**
```python
from src_pyWrapper.pyWrapper import FDTD

solver = FDTD(
    input_filename='my_case.fdtd.json',
    path_to_exe='/path/to/semba-fdtd',
    flags=['-mapvtk'],          # Optional: generate VTK geometry map
    mpi_command='mpirun -n 4'   # Optional: MPI launch
)
solver.cleanUp()    # Remove old output files
solver.run()        # Execute the solver
print(f"Return code: {solver.returncode}")
```

**Key command-line flags:**
| Flag | Purpose |
|------|---------|
| `-mapvtk` | Generate VTK geometry map files for ParaView |
| `-sgbc` | Enable SGBC (multi-layer surface) solver |
| `-r` | Resume from a previous checkpoint |
| `-stoch` / `-nostoch` | Stochastic simulation toggle |
| `-cfl <value>` | Override CFL factor |

### Step 4: Monitor Progress

Check `SEMBA_FDTD_temp.log` for progress. The solver writes timing info periodically.

### Step 5: Analyze Results

**Output file types:**

| Extension | Description | Example |
|-----------|-------------|---------|
| `.dat` | Time-domain probe data (`time field_value`) | `case_name.fdtd_probeName_Ex_5_10_3.dat` |
| `_df.dat` | Frequency-domain (`freq magnitude phase`) | `case_name.fdtd_probeName_Ex_5_10_3_df.dat` |
| `_tr.dat` | Transfer function in dB (`freq 20*log10(mag) phase`) | `case_name.fdtd_probeName_Ex_5_10_3_tr.dat` |
| `.xdmf` + `.h5` | 3D field movies (HDF5 binary + XDMF metadata) | `case_name.xdmf` |
| `.vtk` | Geometry map for ParaView | `case_name_1.vtk` |
| `.old` | Restart field checkpoint | `case_name.fields.old` |

**Load and analyze with pyWrapper:**
```python
from src_pyWrapper.pyWrapper import Probe, FDTD

solver = FDTD(input_filename='pw-in-box.fdtd.json', path_to_exe=SEMBA_EXE)
solver.run()

# Load time-domain probe
before = Probe(solver.getSolvedProbeFilenames("before")[0])
print(before.dataframe)        # pandas DataFrame
print(before.name)             # Extracted probe name
print(before.field_column)     # Field column name
print(before.cell_positions)   # Cell positions

# Plot time domain
import matplotlib.pyplot as plt
plt.plot(before['time'], before['field'], label='Ex component')
plt.xlabel('Time (s)')
plt.ylabel('Field (V/m)')
plt.legend()
plt.show()

# Load frequency-domain probe
freq_probe = Probe(solver.getSolvedProbeFilenames("before")[0].replace('.dat', '_df.dat'))
print(freq_probe.dataframe)    # Columns: frequency, magnitude, phase
```

**Manual analysis:**
```python
import numpy as np
import matplotlib.pyplot as plt

# Time-domain data
data = np.loadtxt('case_name.fdtd_probe_Ex_5_10_3.dat', skiprows=1)
time = data[:, 0]
field = data[:, 1]
plt.plot(time, field)

# Frequency-domain data
freq_data = np.loadtxt('case_name.fdtd_probe_Ex_5_10_3_df.dat', skiprows=1)
freq = freq_data[:, 0]
magnitude = freq_data[:, 1]
phase = freq_data[:, 2]

# Transfer function
tr_data = np.loadtxt('case_name.fdtd_probe_Ex_5_10_3_tr.dat', skiprows=1)
plt.plot(tr_data[:, 0], tr_data[:, 1])  # dB vs frequency

# FFT-based postprocessing (e.g., shielding effectiveness)
INC = np.fft.fft(back['incident'])
S transmitted = np.fft.fft(back['field']) / INC
fq = np.fft.fftfreq(len(t)) / dt
```

## Common Analysis Patterns

### Shielding Effectiveness (SGBC)
```python
# After running simulation with front/back probes
front = Probe(front_files[0])
back = Probe(back_files[0])

INC = np.fft.fft(front['field'])
S transmitted = np.fft.fft(back['field']) / INC
fq = np.fft.fftfreq(len(t)) / dt

# Compare with analytical transmission line model
from skrf.media import Freespace
fq = np.fft.fftfreq(len(t)) / dt
f = fq[idx_min:idx_max]
plt.plot(f, 20*np.log10(np.abs(fdtd_s21)), label='FDTD')
plt.plot(f, 20*np.log10(np.abs(slab.s[:,0,1])), label='Analytical')
plt.ylabel('Shielding Effectiveness (dB)')
```

### RCS (Radar Cross Section) via Far-Field Probes
```python
# Far-field probes produce _df.dat with angular data
farfield = Probe(farfield_files[0])
# Columns: frequency, magnitude, phase (integrated over observation sphere)
```

### MTLN Circuit Analysis
```python
# MTLN probes include voltage and current at each conductor
mtln_probe = Probe(mtln_probe_files[0])
# Check probe.domain_type for "time" or "frequency"
# Columns include voltage/current per conductor
```

### Comparing MTLN vs Non-MTLN Builds
```python
# Run with both builds to compare transmission line effects
solver_mtln = FDTD(fn, SEMBA_MTLN_EXE, run_in_folder=mtln_folder)
solver_mtln.run()
solver_nomtln = FDTD(fn, SEMBA_NOMTLN_EXE, run_in_folder=nomtln_folder)
solver_nomtln.run()

for pf in probe_files_mtln:
    p = Probe(pf)
    plt.plot(p['time'], p['current_0'], label=f'MTLN - {p.name}')
for pf in probe_files_nomtln:
    p = Probe(pf)
    plt.plot(p['time'], p['current'], linestyle='--', label=f'No MTLN - {p.name}')
plt.legend()
```

## Material Types Reference

| Material | Use Case |
|----------|----------|
| `pec` | Perfect electric conductor |
| `pmc` | Perfect magnetic conductor |
| `isotropic` | Dielectric materials (epsilon_r, mu_r, sigma) |
| `wire` | Wire grid materials |
| `shieldedMultiwire` | Multi-conductor cables with PUL parameters |
| `unshieldedMultiwire` | Bundles of parallel wires |
| `terminal` | Port/termination definitions |
| `lumped` | R/L/C lumped elements |
| `multilayeredSurface` | SGBC multi-layer surfaces |
| `thinSlot` | Thin slot models |
| `connector` | SPICE connector definitions |

## Common Gotchas

1. **CFL condition** — The solver automatically adjusts `dt` if the user-specified value would cause instability. Check the log for the adjusted value.
2. **Probe file naming** — Encodes probe name, field component, and cell position: `case_name.fdtd_probeName_Ex_5_10_3.dat`
3. **Domain types** — Probes can output `time`, `frequency`, or `timeFrequency` data. The `_df.dat` and `_tr.dat` files are auto-generated during postprocessing.
4. **Transfer functions** — Specify `magnitudeFile` in a probe's domain to compute transfer function normalization automatically.
5. **MPI output** — When running with MPI, each rank writes its own probe files. Use `solver.getSolvedProbeFilenames()` to find all matching files.
6. **VTK maps** — Use `-mapvtk` flag to generate geometry maps. Open with ParaView and use the `_tag_paraviewfilters.txt` guide for filtering.
7. **Resuming simulations** — Use `-r` flag and `.old` checkpoint files. Ensure the input file hasn't changed since the checkpoint was created.
8. **Mesh origin** — Default is `[0,0,0]`. If your geometry is offset, specify `grid.origin` in the mesh definition.
