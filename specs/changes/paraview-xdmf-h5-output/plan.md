# Plan: ParaView-Compatible XDMF/HDF5 Output

## Architecture

The bounded context is probe visualisation output.
It covers logical probe output folders, sparse selected-point exports,
time-domain visualisation exports, frequency-domain visualisation exports, and
MPI merged visualisation artifacts.

The dependency direction is `infrastructure -> application -> domain`.
The domain layer owns the meaning and invariants of visualisation outputs:
probe identity, selected points, coordinate axes, time states, frequency states,
component groups, vector quantities, frequency entry naming, and rank merge
policy.
The application layer coordinates export use cases such as create logical
output, append time state, publish frequency state, merge rank outputs, and
apply rank-local retention policy.
The infrastructure layer writes concrete files, validates generated XDMF/HDF5
references, and performs any ParaView 5 readability checks.

External visualisation behaviour is treated as an external contract.
The application layer depends on ports for visualisation-file writing,
visualisation-file validation, MPI-rank collection, and rank-local cleanup.
Concrete file formats and external tool checks remain behind those ports.

The plan preserves the existing output subsystem boundaries:
probe registration determines which outputs exist, coordinate selection
determines which grid points are exported, and export publishing turns those
values into user-readable visualisation artifacts.

## Data Models

### Logical Probe Output

Identity is the user-facing probe folder path.
Two logical probe outputs are the same only when their folder identity is the
same.

Fields:
- probe folder path
- probe basename
- output domain: time or frequency
- requested quantity: current density, electric field, or magnetic field
- component mode: full vector or single component
- selected coordinate bounds
- selected coordinate list
- MPI participation: serial, rank-local, or merged logical output

Invariants:
- every generated XDMF/HDF5 file for the logical probe is inside the probe
  folder
- the XDMF file references HDF5 data by a path valid from the XDMF file's
  folder
- the logical probe output is the primary user artifact for ParaView 5
- rank-local files are never the primary user artifact after a successful merge

### Sparse Probe Point Set

Structural equality applies.
Two sparse probe point sets are equivalent when their ordered selected points
and coordinate values are equal.

Fields:
- selected point coordinates
- x coordinate values named `coordsX`
- y coordinate values named `coordsY`
- z coordinate values named `coordsZ`
- selected point count

Invariants:
- each coordinate array is one-dimensional
- `coordsX`, `coordsY`, and `coordsZ` have equal length `N`, where `N` is the
  selected point count
- for every selected point index `p`, the point coordinate is
  `(coordsX[p], coordsY[p], coordsZ[p])`
- selected point indexes are ordered by ascending global physical coordinates:
  first x, then y, then z
- dataset names remain stable
- only selected probe points are represented
- unselected points are omitted, not zero-filled, NaN-filled, or masked
- every visualised field state uses the same selected point order for a logical
  probe

### Time Series

Structural equality applies.
Two time series are equivalent when their ordered time values are equal.

Fields:
- ordered time values named `times`
- logical state count

Invariants:
- the number of time values equals the number of visible time states
- repeated flushes append to the logical series
- repeated flushes do not reset prior visible time states
- one logical time-domain probe has one logical XDMF file and one logical HDF5
  file
- scalar time-domain value arrays have dimensions `T x N`, where `T` is the
  number of time states and `N` is the selected point count
- vector time-domain value arrays have dimensions `T x N x 3`
- vector component order is x, y, then z
- scalar value `values[t, p]` belongs to `times[t]` and selected point `p`
- vector value `vectors[t, p, 0..2]` belongs to `times[t]` and selected
  point `p`

### Frequency Series

Structural equality applies.
Two frequency series are equivalent when their ordered frequency values are
equal.

Fields:
- ordered frequency values
- frequency state count
- named frequency entries

Invariants:
- every requested frequency is visible once in the frequency sequence
- every requested frequency is also available as a named frequency entry
- each named frequency entry uses scientific notation in hertz with prefix
  `freq_` and suffix `Hz`, for example `freq_1.000e6Hz`
- frequency names use a dot as decimal separator, exactly three digits after
  the decimal point, lowercase `e`, no plus sign for positive exponents, no
  leading zeroes in the exponent, and a minus sign for negative exponents
- frequency value `1000000` is named `freq_1.000e6Hz`
- frequency value `0.001` is named `freq_1.000e-3Hz`
- scalar frequency-domain value arrays have dimensions `F x N`, where `F` is
  the requested frequency count and `N` is the selected point count
- vector frequency-domain value arrays have dimensions `F x N x 3`
- vector component order is x, y, then z
- frequency value index `f` belongs to the `f`th requested frequency and
  selected point index `p` follows the sparse point-set alignment
- frequency values are clearly distinguished from physical simulation time

### Visual Attribute

Identity is the logical probe output plus attribute name.

Fields:
- attribute name
- centre: node-centred
- attribute kind: scalar or vector
- quantity: current density, electric field, or magnetic field
- component: x, y, z, or full vector
- data dimensions over selected points and states

Invariants:
- scalar attributes contain exactly one component
- vector attributes contain exactly three components
- full-vector requests publish x, y, and z scalar attributes and one vector
  attribute
- single-component requests publish only the requested scalar component
- attributes contain values only for selected probe points
- unselected points do not appear in scalar or vector attributes

### MPI Merge Policy

Identity is the logical probe output plus execution run.

Fields:
- rank-local output locations
- merged logical output location
- retention setting: retain rank-local files or remove rank-local files
- merge status: not started, complete, or failed

Invariants:
- merged output is created after rank-local output has been collected
- merged coordinate and value data is ordered by ascending global physical
  coordinates: first x, then y, then z
- merged output ordering does not depend on MPI rank number, rank-local file
  discovery order, or output-register entry order
- rank-local XDMF/HDF5 files are removed only after a successful merge and only
  when the retention setting requests removal
- failed merge leaves rank-local XDMF/HDF5 files available for diagnostics

### Readability Check Result

Structural equality applies.
Two readability check results are equivalent when they have the same pass/fail
state and the same reported missing references.

Fields:
- target XDMF path
- ParaView major version
- pass/fail state
- missing file references
- missing dataset references
- manual path edit requirement

Invariants:
- the minimum passing check is that ParaView 5 opens the generated XDMF file
  without missing-file errors, missing-dataset errors, or manual path edits
- a passing result requires major version `5`, target XDMF path equal to the
  generated logical probe XDMF path, no missing file references, no missing
  dataset references, and no manual path edit requirement
- a passing minimum check does not require visual rendering, saved-state reload,
  or verification that every expected array is visible in the user interface

## Interface Contracts

### Use Case: Publish Time-Domain Visualisation Output

Input:
- logical probe output
- sparse probe point set
- ordered time values to append
- scalar component arrays with dimensions `T x N` for requested components
- optional full-vector arrays with dimensions `T x N x 3` for full-vector
  requests

Output:
- one logical XDMF file in the probe folder
- one matching logical HDF5 file in the probe folder
- validation result describing whether the generated output is ParaView 5
  readable

Failure outcomes:
- no successful publication if mandatory selected coordinates, time, or
  requested value arrays are missing
- no successful publication if any value array dimensions do not match the
  sparse point count or time-state count
- no successful publication if XDMF references do not resolve from the XDMF
  folder
- prior logical time states remain intact when an append fails

### Use Case: Publish Frequency-Domain Visualisation Output

Input:
- logical probe output
- sparse probe point set
- ordered frequency values
- scalar component magnitudes with dimensions `F x N` for requested components
- optional full-vector magnitudes with dimensions `F x N x 3` for full-vector
  requests

Output:
- one ParaView-readable frequency sequence with one state per requested
  frequency
- one logical XDMF file in the probe folder
- one matching logical HDF5 file in the probe folder
- separately named frequency entries for every requested frequency using
  `freq_<scientific-value>Hz` names
- validation result describing whether all frequency entries are readable

Failure outcomes:
- no successful publication if any requested frequency lacks required data
- no successful publication if any value array dimensions do not match the
  sparse point count or frequency count
- no successful publication if any named frequency entry does not match the
  required scientific-name contract
- no successful publication if frequency values can be confused with physical
  simulation time in the generated visualisation contract

### Use Case: Merge MPI Visualisation Output

Input:
- logical probe identity
- rank-local output locations
- rank count
- retention setting

Output:
- one merged logical XDMF/HDF5 output per probe
- merged coordinate and value data ordered by ascending global x, then y,
  then z coordinates
- rank-local XDMF/HDF5 files retained or removed according to retention setting
- global output register entry that points users to the merged logical output

Failure outcomes:
- if merge fails, rank-local files are retained
- if merge succeeds and retention requests removal, rank-local files are
  removed only after the merged logical output is readable

### Ports

| Port | Responsibility |
| --- | --- |
| Visualisation file writer | Create logical XDMF/HDF5 files from domain output descriptions |
| Sparse output contract validator | Verify coordinate/value array lengths, dimensions, ordering, and vector packing |
| Visualisation reference validator | Verify that XDMF references resolve to matching HDF5 datasets |
| ParaView readability checker | Verify the minimum ParaView 5 open-cleanly contract |
| MPI rank output collector | Discover rank-local visualisation outputs for a logical probe |
| Rank-local retention manager | Retain or remove rank-local XDMF/HDF5 files after merge |
| Output register publisher | Publish user-facing logical output paths in the output register |

No repository is needed because this change has no database aggregate.
Filesystem persistence is handled by ports and belongs to infrastructure.

### Domain Events

| Event | Meaning |
| --- | --- |
| LogicalVisualisationOutputCreated | The logical probe output files were created in the probe folder |
| TimeStatesAppended | Time-domain states were appended to a logical probe output |
| FrequencyStatesPublished | Frequency-domain states were published for a logical probe output |
| MPIVisualisationOutputMerged | Rank-local visualisation output was merged into the logical output |
| RankLocalFilesRetained | Rank-local files remain available after merge |
| RankLocalFilesRemoved | Rank-local files were removed after successful merge |

## Implementation Phases

1. Output contract model — define the logical probe output, sparse probe point
   set, time series, frequency series, visual attribute, readability check
   result, and MPI merge policy used by all visualisation outputs.
2. Time-domain publishing contract — establish one logical XDMF/HDF5 pair per
   time-domain volumetric probe, including coordinate arrays, time states,
   selected point coordinates, `T x N` scalar arrays, `T x N x 3` vector arrays,
   vector attributes, append semantics, and validation.
3. Frequency-domain publishing contract — establish both frequency-sequence and
   scientific-notation named-frequency visualisation outputs for every
   requested frequency, including `F x N` scalar arrays, `F x N x 3` vector
   arrays, and exact frequency-name formatting.
4. MPI merged-output contract — establish rank-local collection, merged logical
   output publication in ascending global x/y/z order, global register entries,
   and user-selectable retention.
5. Validation and compatibility checks — verify file presence, dataset names,
   coordinate/value dimensions, point ordering, vector packing, XDMF references,
   vector/scalar publication, and the ParaView 5 open-cleanly readability check.
6. Canonical feature spec update — promote the accepted visualisation output
   contract into the relevant canonical feature specs after implementation
   planning and analysis pass.

## Decisions & Rationale

| Decision | Rationale | Alternatives Considered |
| --- | --- | --- |
| Use one logical XDMF/HDF5 pair per time-domain probe | Users need one stable artifact per probe, independent of flush cadence | Multiple files per flush; rank-local files as primary artifacts |
| Keep generated files inside the probe folder | Matches the refined output folder rule and keeps probe artifacts self-contained | Case-root flat files; shared folder for all probes |
| Publish scalar components and vector attributes for full-vector quantities | Supports scripts that need components and ParaView tools that need vectors | Scalars only; vectors only |
| Represent frequency output in both sequence and named-entry forms | Supports browsing by frequency and direct selection of named frequency results | Sequence only; named entries only |
| Name frequency entries with `freq_<scientific-value>Hz` | Gives stable, human-readable frequency identifiers | Index-only names; raw decimal names |
| Use `T x N`, `F x N`, and `* x N x 3` array shapes | Makes selected-point alignment and vector packing testable | Per-component files only; full rectangular arrays |
| Export only selected probe points | Matches the clarified sparse-output contract and avoids invented values outside selection | Zero-fill; NaN-fill; validity mask; full rectangular volume |
| Treat ParaView readability as an external validation port | Keeps domain rules independent from external tooling while preserving user-facing behaviour | Embed tool-specific checks in domain rules; omit validation |
| Define minimum readability as opening cleanly in ParaView 5 | Gives a measurable acceptance check without over-specifying UI rendering | Require visible arrays; require rendered image; require saved-state reload |
| Merge MPI rank-local outputs into one logical probe output ordered by global x/y/z | Users should not manually assemble rank pieces, and ordering must be reproducible | Rank pieces only; direct concurrent writes to one file; rank-number ordering |
| Make rank-local retention user-selectable | Retention is useful for diagnostics but unnecessary for normal visualisation | Always retain; always delete |

## Constraints & Risks

- No `FRAMEWORK.md` file was present, so the plan applies the stated dependency
  rule directly.
- The current canonical specs include compatibility-sensitive names such as
  `CurrenDensityX`, `CurrenDensityY`, and `CurrenDensityZ`; changing those names
  requires an explicit compatibility decision.
- The existing frequency-domain output behaviour includes compatibility quirks
  around binary records and accumulation; the visualisation contract must not
  silently change numerical solver results while improving ParaView readability.
- ParaView 5 readability is an external behaviour and should be verified
  through a validation port rather than assumed from file creation alone.
- The minimum readability check is intentionally limited to opening cleanly;
  it does not prove visual rendering or user-interface array visibility.
- Sparse selected-point export must remain readable by ParaView 5 without a
  full rectangular volume.
- Scientific-notation frequency names must be formatted consistently enough to
  be stable across platforms and runs.
- The plan intentionally requires selected-point arrays rather than full
  rectangular arrays; any ParaView 5 limitation around sparse point rendering
  must be caught by the open-cleanly validation port.
- Rank-local cleanup must be failure-safe: deletion can only occur after a
  successful and readable merged logical output exists.
