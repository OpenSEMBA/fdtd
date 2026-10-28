# Output Plotting

`plot_outputs.py` is the single entry point for plotting the probe outputs of a
finished simulation. Run it from the repository root with the packages from
`requirements.txt` installed; the 3D commands also need ParaView's `pvbatch`
on `PATH`.

```bash
python3 tools/plot_outputs/plot_outputs.py <command> <source> [-o <output>]
```

`<source>` is the solver output directory (the one holding the `.dat` files and
the `.fdtd.json` input). `-o` sets the directory for the PNG files and defaults
to `<source>/plots`.

| Command | What it draws |
| --- | --- |
| `probe` | One PNG per probe `.dat` file. `<source>` may also be a single file. |
| `group` | Probes drawn together in one figure, selected by probe name. |
| `3d-probes` | Probe locations in 3D, colored by peak time-domain value (ParaView). |
| `3d-scene` | The same probes plus the problem geometry from the `__MAP` folder (ParaView). |

## Individual plots

```bash
python3 tools/plot_outputs/plot_outputs.py probe testing/mpi-single-probe/run
```

## Grouped plots

Each `-g TITLE=NAME[,NAME...]` draws the listed probe names in one figure, with
one PNG per domain (`_tm` for time, `_fq` for frequency). Names are the probe
names from the input file, not the file names.

```bash
python3 tools/plot_outputs/plot_outputs.py group testing/mpi-single-probe/run \
  -g "Center=center" \
  -g "Corners=corner_000,corner_N00"
```

Without `-g`, every probe name gets its own figure.

## 3D views

```bash
python3 tools/plot_outputs/plot_outputs.py 3d-probes testing/mpi-single-probe/run
python3 tools/plot_outputs/plot_outputs.py 3d-scene testing/pec_cavity/run
```

Probe positions come from the cell indices in the file names and the grid
(`mesh.grid`) of the `.fdtd.json` input. `3d-scene` needs the `__MAP` folder,
so run the solver with `-mapvtk`. Both commands also write `probes.csv` next
to the image, with the position and peak value of every probe.

## Tests

```bash
.venv/bin/python -m pytest test/tools/test_plot_outputs.py
```

## Limits

- Only point and bulk-current probe files are plotted. Wire, line, far-field,
  energy and MTLN outputs are not recognised and are skipped.
- 3D views use time-domain probes only, since frequency and time peaks have
  different units.
