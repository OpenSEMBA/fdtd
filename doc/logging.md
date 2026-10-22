# Solver logging

This page describes how semba-fdtd reports warnings, errors and diagnostic
messages, and how to add new ones. The implementation lives in
`Report_m` (`src_main_pub/errorreport.F90`).

## Summary

- Every diagnostic message goes through one call:
  `call log_msg(LEVEL, 'ORIGIN', message)`.
- Each line carries the MPI rank, the time step (inside the time loop), the
  subsystem that emitted it and its severity:

  ```text
  (    2/    4) [n=1200] [WIRES] WARN: wir0_WARNING: ...
  └── rank ───┘ └ step ┘ └origin┘ └level┘ └─ message
  ```

- `LOG_ERROR` marks the run as failed. Before the time loop starts, the solver
  aborts unless `-ignoreerrors` was given. `LOG_INFO` and `LOG_WARN` never stop
  the run.
- `LOG_FATAL` is the immediate-abort path and **does not return**: the message
  goes to stdout, `_Report.txt` and the warnings file, then every rank stops
  with exit code 1. It replaces the old `StopOnError`.
- `print11` is unchanged. It writes formatted, human-oriented output (banners,
  tables, progress) to stdout and `_Report.txt`. It is not a diagnostics
  channel.

## Levels

| Level | Constant | Use it for | Sets fatal flag | Typical source |
|---|---|---|---|---|
| `INFO` | `LOG_INFO` | Diagnostics, configuration dumps, verbose traces | no | `-verbose` output, media dump, wire topology |
| `WARN` | `LOG_WARN` | Something odd that the solver can work around | no | ignored sources, write failures that are skipped |
| `ERROR` | `LOG_ERROR` | Invalid input or state; the run must not proceed | **yes** | parser validation, unsupported geometry |
| `FATAL` | `LOG_FATAL` | Unrecoverable state; the process exits right away | n/a | missing input file, MPI failures |

```mermaid
flowchart TD
    A[Something worth reporting] --> B{Can the run continue<br/>with correct results?}
    B -- "No, and we must stop now" --> F["log_msg(LOG_FATAL, ...)"]
    B -- "No, but let the rest of the<br/>input be checked first" --> E["log_msg(LOG_ERROR, ...)"]
    B -- Yes --> C{Does the user need<br/>to act on it?}
    C -- Yes --> W["log_msg(LOG_WARN, ...)"]
    C -- "No, diagnostic only" --> I["log_msg(LOG_INFO, ...)"]
```

## Message flow

```mermaid
flowchart LR
    subgraph callers [Solver modules]
        P[smbjson / preprocess]
        S[timestepping]
        M[MTLN / wires / outputs]
    end

    callers -->|"log_msg(level, origin, msg)"| L[log_msg]
    L --> F{"level >= ERROR?"}
    F -- yes --> FE[(fatalerror = .true.)]
    L --> FMT["log_line()<br/>[n=step] [ORIGIN] LEVEL: msg"]
    FMT --> PFX["prefix rank<br/>(rank/nprocs)"]
    PFX --> U17[("unit 17<br/>&lt;root&gt;&lt;rank&gt;_tmpWarnings.txt<br/>(only while open)")]
    PFX --> ERR[(stderr)]

    S -->|"set_log_step(n)"| STEP[(logStep)]
    STEP -.-> FMT

    SOE["log_msg(LOG_FATAL, ...)"] --> OUT11["print11: stdout + _Report.txt"]
    SOE --> U17
    SOE --> X[["MPI_Abort / stop 1"]]
```

## Warnings file lifecycle

Each rank writes its own temporary file. Rank 0 merges them into
`<root>_Warnings.txt` and reduces the fatal flag across ranks with
`MPI_AllReduce(MPI_LOR)`.

```mermaid
sequenceDiagram
    autonumber
    participant Sw as interpret_switches
    participant Rep as Report_m
    participant Pre as parser / preprocess
    participant TS as timestepping
    participant FS as Files

    Sw->>Rep: INITWARNINGFILE
    Rep->>FS: open unit 17 = <root><rank>_tmpWarnings.txt
    Pre->>Rep: log_msg(LOG_ERROR, 'SMBJSON', ...)
    Rep->>FS: "(1/4) [SMBJSON] ERROR: ..." on unit 17 and stderr
    TS->>Rep: CLOSEWARNINGFILE
    Rep->>FS: rank 0 merges every _tmpWarnings.txt into <root>_Warnings.txt
    Rep-->>TS: fatalerror reduced over all ranks
    alt fatal error and no -ignoreerrors
        TS->>Rep: log_msg(LOG_FATAL, 'SOLVER', 'FATAL ERRORS. Revise *Warnings.txt file. ABORTING...')
    else
        loop every time step n
            TS->>Rep: set_log_step(n)
            TS->>Rep: log_msg(...) goes to stderr only, tagged [n=...]
        end
        TS->>Rep: set_log_step(-1)
    end
```

The warnings file is closed before the time loop starts. Messages emitted
during time stepping therefore go to **stderr only**. That is where the
`[n=...]` tag is most useful.

## Output files

| File | Writer | Content |
|---|---|---|
| `<root><rank>_tmpWarnings.txt` | `log_msg` on every rank | Per-rank messages before the time loop. They are left on disk if the run aborts. |
| `<root>_Warnings.txt` | `CLOSEWARNINGFILE` on rank 0 | Merged messages from every rank |
| `<root>_Report.txt` | `print11` on rank 0 | Human-oriented run report: banners, configuration, progress |
| stderr | `log_msg` on every rank | Every leveled message, including those emitted inside the time loop |
| stdout | `print11` | Same content as `_Report.txt` while printing is on (`OnPrint` / `OffPrint`) |

## Origins

The origin tag names the subsystem, not the file. Reuse an existing tag where
one fits:

| Origin | Area |
|---|---|
| `SMBJSON` | `.fdtd.json` parser |
| `PREPROCESS` | geometry preprocessing (`preprocess_geom.F90`) |
| `ROTATE` | NFDE rotation for MPI layouts |
| `WIRES` | Holland / thin-wire models |
| `MTLN` | multiconductor transmission lines and ngspice |
| `SOLVER` | time-stepping setup and loop |
| `MEDIA` | media dump (`-verbose`) |
| `SOURCES`, `PLANEWAVE`, `LUMPED` | excitations and lumped elements |
| `SGBC`, `EDISPERSIVE`, `MDISPERSIVE` | material models |
| `PML`, `CPML`, `MUR` | boundaries |
| `FARFIELD`, `OUTPUT`, `PROGRESS`, `RESUMING` | outputs, progress and restart |
| `MAIN` | `semba_fdtd.F90` launch flow |
| `MPI` | MPI communication setup (`mpicomm.F90`) |

## Writing a new message

```fortran
use Report_m, only: log_msg, LOG_INFO, LOG_WARN, LOG_ERROR, LOG_FATAL

write(buff, *) 'Thin slot width ', width, ' is below the cell size'
call log_msg(LOG_WARN, 'PREPROCESS', buff)

call log_msg(LOG_ERROR, 'SMBJSON', 'Unknown probe type "'//trim(probeType)//'"')

! Runtime choice of severity
call log_msg(merge(LOG_ERROR, LOG_WARN, isFatal), 'WIRES', buff)

! Unrecoverable: does not return
call log_msg(LOG_FATAL, 'OUTPUT', 'Unable to initialise movie HDF5 output')
```

Guidelines:

- Do not put `WARNING:` or `ERROR:` in the message text: the level is already
  printed. Older messages still carry tags such as `wir0_WARNING:`; remove
  them when you edit those lines.
- Use `LOG_ERROR` only when the run has to stop. It sets the fatal flag and
  lets the rest of the input be checked. Use `LOG_FATAL` when it must stop now;
  do not follow a log call with `stop` or `error stop`.
- Do not use `print *` or `write(*,*)` for diagnostics. Those messages bypass
  the warnings file and carry no rank or step.
- The `(rank/nprocs)` prefix comes from `set_log_rank`, which `semba_fdtd`
  calls right after MPI starts and again after a stochastic run splits into
  halves, so each half has its own rank 0 for the signal files. Call it again
  if your code changes `layoutnumber`. Do not pass the rank in the message.
- Inside the time loop, `log_msg` adds the step automatically. Do not format
  `n` into the message.

## Before and after

| Before | After |
|---|---|
| `call WarnErrReport(buff)` | `call log_msg(LOG_WARN, 'ORIGIN', buff)` |
| `call WarnErrReport(buff, .true.)` | `call log_msg(LOG_ERROR, 'ORIGIN', buff)` |
| `call WarnErrReport(buff, flag)` | `call log_msg(merge(LOG_ERROR, LOG_WARN, flag), 'ORIGIN', buff)` |
| `print *, 'Bug in tags'` then `stop` | `call log_msg(LOG_FATAL, 'PREPROCESS', 'Bug in tags')` (no `stop`) |
| `call StopOnError(layoutnumber, num_procs, msg)` | `call log_msg(LOG_FATAL, 'ORIGIN', msg)` |
| `(    1/    1) some text` | `(    1/    1) [n=42] [WIRES] WARN: some text` |
| `(    1/    1) ERROR: msg` (StopOnError, console only) | `(    1/    1) [OUTPUT] FATAL: msg` (console and warnings file) |

## Known limitations

- **No step tag in `_Warnings.txt`.** The warnings file is merged and closed
  before the time loop, so `[n=...]` only shows up on stderr.
- **Levels were assigned mechanically.** During the migration, the level came
  from the old error flag: messages whose text says `INFO`, or that only print
  under `verbose`, became `INFO`, and the rest became `WARN`. Some `WARN`
  messages are really informational. Reclassify them when you touch them.
- **No level filter.** There is no `-loglevel` switch yet. Use `-verbose` to
  get the extra `INFO` output that already depends on it.
