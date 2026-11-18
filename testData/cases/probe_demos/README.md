# Probe demo cases

Small, independent simulations that exercise every probe type and domain the
solver supports. Each `.fdtd.json` runs on its own. Run the solver from this
folder, because `magnitudeFile` is resolved against the working directory.

| Case | Probes | Domains |
| --- | --- | --- |
| `point.fdtd.json` | `point` (E and H, E in frequency and timeFrequency) | time, frequency, timeFrequency |
| `movie.fdtd.json` | `movie` (magnitude, component x) | time, frequency |
| `line.fdtd.json` | `line` | time |
| `bulk_current.fdtd.json` | `bulkCurrent` | time |
| `farfield_small.fdtd.json` | `farField` on a 5 cm PEC sphere (bistatic RCS) | frequency |
| `wire.fdtd.json` | `wire` (current, voltage) | time |

Not covered: the `wire` `charge` field and MTLN bundle probes. Charge and
multiwire voltage need a multiwire (MTLN) material association, which these
single-wire cases do not set up.

## Run

```bash
cd testData/cases/probe_demos
../../../build-rls/bin/semba-fdtd -i point.fdtd.json
```

Repeat with any other case. Each run writes its `.dat` files, `_Report.txt`,
`_Warnings.txt` and `_Energy.dat` next to the input. Those outputs are not
tracked (see `.gitignore`). Runtime is under a minute for most cases; the
far-field case takes a few minutes.

## Plot

From the repository root:

```bash
python3 tools/plot_outputs/plot_outputs.py probe testData/cases/probe_demos
python3 tools/plot_outputs/plot_outputs.py farfield testData/cases/probe_demos
```

`probe` draws the point, bulk current, line and wire outputs. `farfield` draws
the RCS pattern and the backscatter curve of the sphere case. Movie outputs are
not plotted (see `tools/plot_outputs/README.md`).

## Regenerating inputs

- `make_gauss.py` writes `gauss.exc` (30 ps step, used by the other cases) and
  `gauss_wide.exc` (7.7 ps step, a wider pulse with the same peak, used by
  `farfield_small`).
- `make_farfield_small.py` writes `farfield_small.fdtd.json`: a 5 cm PEC sphere
  on a 5 mm grid, probed from 0.05 to 0.25 GHz (ka = 0.05 to 0.26).

## Far-field caveat

The backscatter values (theta = 180 deg) run from about -71 dBsm at 0.05 GHz to
-43 dBsm at 0.25 GHz. That matches the Rayleigh estimate 9 pi a^2 (ka)^4 to within
a few dB, but it has not been checked against an exact Mie solution.

The bistatic cut at phi = 90 deg is smooth. The cut at phi = 0 deg is dominated
by numerical noise at -40 to -120 dBsm, which is below any physical RCS of this
sphere. Treat that cut as unreliable.
