# Startup memory checks

Before the time-stepping loop starts, `semba-fdtd` compares the memory it is
about to allocate with the memory available to each process.  If the run does
not fit, it stops through `log_msg(LOG_FATAL, 'MEMORY', ...)` with a message
showing the required and available amounts.

## Checkpoints

The check is split into the allocation groups that own the corresponding
arrays:

| Group | Owner | When |
|---|---|---|
| media/tag copies, fields, distances | `memoryBudget_m` | before `create_solver` and before field allocation |
| SGBC | `SGBC_nostoch_m` | before `InitSGBCs` |
| Holland wires | `HollandWires_m` | before `InitWires` |
| MTLN bundles | `mtln_preprocess_m` | before `InitWires_mtln` |
| CPML | `BORDERS_CPML_m` | before `InitCPMLBorders` |
| conformal fields | `conformal_timestepping_m` | before `initializeConformalFields` |
| probes | `output_m` | at the start of `init_outputs` |
| HDF5/transient output work | `output_m` | pre-output and before time stepping |

Each checkpoint uses the resident memory already in use plus only the pending
allocation for that group, so groups are not counted twice.

## Available memory

On Linux the probe reads:

- `MemTotal` and `MemAvailable` from `/proc/meminfo`,
- cgroup v1/v2 memory limits and usage,
- `RLIMIT_AS` and `RLIMIT_DATA` headroom,
- `VmRSS` and `VmSize` from `/proc/self/status`.

On Windows it uses `GlobalMemoryStatusEx`, process memory counters and, when
present, the process memory limit of a Job Object.  If a platform cannot report
available memory, the check is skipped with a warning instead of aborting.

In MPI runs, ranks are grouped by node with `MPI_Comm_split_type`.  The node
headroom is divided equally among the local ranks; the check therefore
identifies the rank that does not fit and prints its host name.

## Estimate accuracy

Estimates are exact for the arrays that the owning module allocates (fields,
distances, SGBC, CPML) and conservative for the pre-allocation probe estimate,
Holland wires, conformal maps and HDF5 staging.  The post-allocation probe
estimate is exact for the buffers that actually exist.  HDF5 library-internal
allocations are not directly observable; the probe reserves the configured
1 MiB chunk target per dataset plus a fixed per-file metadata allowance.

The check is a lower-bound/early-warning mechanism.  Other jobs on the same
node can still consume memory after the check, and the operating system may
overcommit virtual memory.

## Example messages

Serial:

```text
(    1/    1) [MEMORY] FATAL: Not enough memory to run the simulation.
  Required:  48.2 GB
  Available: 31.4 GB
```

MPI:

```text
(    4/   24) [MEMORY] FATAL: Not enough memory on rank 3 (node03).
  Required by process:  12.1 GB
  Available to process:  7.8 GB
```
