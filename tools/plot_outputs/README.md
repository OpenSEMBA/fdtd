# Output Plotting

`plot_probe.py` plots two-column time-series output written by scalar probes.
It accepts a flat `.dat` file or a coordinated probe's JSON descriptor.

Run it from the repository root after installing the packages in
`requirements.txt`:

```bash
python3 tools/plot_outputs/plot_probe.py \
  Testing/nodal-source-with-movie.fdtd_electric_field_near_source_Ex_5_15_5_tm.dat
```

Use `--output` to save an image without opening a window:

```bash
python3 tools/plot_outputs/plot_probe.py path/to/probe_tm.dat \
  --output plots/probe.png
```

The descriptor form selects its canonical text artifact and uses the output
quantity as the y-axis label.

To render every supported two-column flat or coordinated text-series probe
below an output directory:

```bash
MPLBACKEND=Agg python3 tools/plot_outputs/plot_all_probes.py Testing Testing/plots
```
