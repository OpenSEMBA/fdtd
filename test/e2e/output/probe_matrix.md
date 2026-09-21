# Output Probe Coverage Matrix

This matrix defines the probe families covered by the output module.
It is the authoritative test inventory for unit, multi-worker, and end-to-end
suites.

## Support Matrix

The release gate covers these execution modes:

| Execution mode | Environment | Build requirements |
| --- | --- | --- |
| Single worker | Linux | Output module and visualisation output enabled |
| Single worker | Windows | Output module and visualisation output enabled |
| Two workers, collective publication | Linux | MPI and parallel visualisation output enabled |
| Two workers, root aggregation | Linux | MPI and visualisation output enabled |

The matrix covers both configured time-domain and frequency-domain requests.
The supported precision variants are the default and double-precision builds
where the selected environment provides them.

## Terminology

- **Canonical artifact:** the one logical result exposed to consumers.
- **Fragment:** one worker-owned contribution to a volumetric result.
- **Collective publication:** participating workers publish disjoint regions
  as one logical volumetric result.
- **Root aggregation:** the designated root publishes the gathered logical
  result when collective publication is unavailable.

| Probe family | Canonical result | Worker fragments | Unit lifecycle coverage | End-to-end coverage |
| --- | --- | --- | --- | --- |
| Point electric or magnetic field | Flat `.dat` series only | No | Initialise, update, flush, cleanup, flat layout | One point probe with no sidecars |
| Wire current | Flat `.dat` series only | No | Initialise, update, flush, cleanup, absent local segment | One current probe with no sidecars |
| Wire charge | Flat `.dat` series only | No | Initialise, update, flush, cleanup, absent local segment | Not exposed by the JSON parser |
| Bulk current or magnetic circulation | Flat reduced `.dat` series only | No | Initialise, update, flush, cleanup, flat layout | One bulk probe with no sidecars |
| Line integral | Flat `.dat` series only | No | Initialise, update, flush, cleanup, flat layout | One line probe with no sidecars |
| Far field | Flat `.dat` result only | No | Initialise, update, flush, flat layout | One far-field probe with no sidecars |
| Geometry map | One canonical geometry map without JSON metadata | No | Initialise and publish | One geometry-map probe |
| Time-domain volumetric movie | One canonical binary and visualisation result | Yes | Initialise, update, flush, cleanup, finalise | One movie probe crossing a worker boundary |
| Frequency-domain volumetric slice | One canonical binary and visualisation result | Yes | Initialise, update, flush, cleanup, finalise | One frequency-slice probe crossing a worker boundary |

Every end-to-end case runs in single-worker and two-worker modes.
The canonical result must be equivalent in both modes.
Volumetric cases must additionally expose one fragment for each contributing worker.
No probe family may publish a JSON descriptor or run output manifest.
