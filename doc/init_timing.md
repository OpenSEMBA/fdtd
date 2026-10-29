# Init Timing Report

The `Init Timing...` block is printed once, after the probes are initialised
and right before the time-stepping loop starts
(`src_main_pub/timestepping.F90`, `InitTiming` in `src_main_pub/progress.F90`).
It summarises what the solver is about to do: probe buffer memory, mesh size,
flush policy and reporting cadence.

## Before and after

Before:

```
Init Timing...
Total Mcells:   0.125000000
NO flushing of restarting FIELDS scheduled
Flushing observation DATA every      10000001  minutes and every         1024  steps
Reporting simulation info every             1  minutes
```

After:

```
Init Timing...
Probe buffer memory (estimate, all ranks):
movie                  1 probes          0.41 MB
total                  1 probes          0.41 MB
Total Mcells:   0.125000000
Restart FIELDS flush: none scheduled
Observation DATA flush: every 128 steps (line, movie, frequency slice); time-based: none
Reporting simulation info every 1 minutes
```

## What was wrong

| Old line | What it claimed | What actually happens |
|---|---|---|
| `every 10000001 minutes` | A time-based flush every ~19 years | `TOPCPUTIME` is the "no time-based flush" sentinel; the extra `1` came from a single-precision `seconds/60.0` division. |
| `every 1024 steps` | Buffer flush every `BUFFOBSE` steps | `BUFFOBSE` is legacy. Volumetric and line probes flush every `OUTPUT_TIME_BUFFER_SIZE` (128) steps; point, wire-current, wire-charge and bulk probes flush every `SMALL_PROBE_BUFFER_SIZE` (4096) steps. |
| `NO flushing of observation DATA scheduled` | No probe data is written | The step-based flush always runs; the warning was removed. |

## Lines

| Line | Source | Values |
|---|---|---|
| `Probe buffer memory` | `report_outputs_memory` (`src_output/output.F90`) | Only printed when there are probes. One row per probe type in use, then the total. Counts are per output object, summed over ranks: `point` counts probe x component, and `movie` counts movies x regions (one per MPI rank holding part of a movie), so a single movie on 4 ranks shows `4 probes`. |
| `Total Mcells` | Sweep extent, summed over MPI ranks | Millions of cells. |
| `Restart FIELDS flush` | `-flush <minutes>`, `maxCPUtime` | `every N minutes`, `at the end, after N minutes`, or `none scheduled`. |
| `Observation DATA flush` | `OUTPUT_TIME_BUFFER_SIZE`, `SMALL_PROBE_BUFFER_SIZE`, `-flushdata <minutes>` | `every 128 steps (line, movie, frequency slice)` and/or `every 4096 steps (point, wire, bulk)`, joined by `, ` for the groups present on any rank; `at the end only` when only far-field probes exist; `no probes` when there is nothing to observe. Then `; time-based: none` or `; time-based: every N minutes`. |
| `Reporting simulation info` | `REPORTINGSECONDS` | Minutes between progress reports. |

## How probe data reaches disk

```mermaid
flowchart LR
  S[time step n] --> U[update_outputs: sample into buffer]
  S --> M{n mod 128 == 0}
  M -->|yes| F[flush_outputs: volumetric and line probes]
  F --> SP{n mod 4096 == 0}
  SP -->|yes| PP[point, wire and bulk probes]
  T[Timing, about every 60 s] --> D{elapsed > flushdata minutes}
  D -->|yes, if -flushdata given| ALL[flush all probes]
```

## Memory estimate

```mermaid
flowchart TD
  A[initializeObservation: buffers allocated] --> B[estimate_outputs_memory]
  B -->|bytes and count per outputID| C[MPI_AllReduce SUM over SUBCOMM_MPI]
  C --> D[print one row per probe type plus total]
```

`estimate_outputs_memory` adds `size * storage_size` of every allocated buffer
held by each probe: time samples, frequency accumulators, coordinates and
material tags. Time buffers have a fixed size of `ceil(flushSteps / tstride)` samples
(`time_buffer_capacity`), with `flushSteps` 4096 (`SMALL_PROBE_BUFFER_SIZE`) for point, wire-current, wire-charge and bulk probes and 128 (`OUTPUT_TIME_BUFFER_SIZE`) for the rest, so the figure is valid for the whole run.

Each row is the sum over all probes of that type. Two rows come from outside
the `outputs` array:

- `far field`: the near-to-far-field DFT accumulators held in `FF`
  (`farfield_buffer_bytes` in `src_main_pub/farfield.F90`). All far-field probes
  share this single accumulator, so the row reports it once.
- `MTLN`: wire probes solved by MTLN (`src_mtln/probes.F90`), reached through
  `GetSolverPtr`. They buffer a single frame, so they are usually 0.00 MB.

Limits:

- HDF5, MPI transport and metadata buffers are not included.
- Under MPI, a probe split across ranks is counted once per rank.
