# Output System Specification Index

## Scope

These specifications cover the simulation output subsystem only.
They describe how configured observations become output files, in-memory output
state, geometry exports, and external file contracts.

## Capability Partition

| Capability | Canonical scope | Evidence |
| --- | --- | --- |
| Output lifecycle | registration, update, flush, close, domain normalisation | VERIFIED |
| Scalar probes | point electric and magnetic scalar output | VERIFIED |
| Wire probes | wire current and wire charge output | VERIFIED |
| Block probes | bounded current and magnetic-current integral output | VERIFIED |
| Volumetric time output | time-domain volumetric field and current-density output | VERIFIED |
| Volumetric frequency output | frequency-domain volumetric accumulation and export | VERIFIED |
| Far-field output | far-field accumulation and flush gateway | VERIFIED |
| Geometry map output | material and face-tag map export | VERIFIED |
| Structured export contracts | VTK, HDF5, and XDMF file contracts | VERIFIED |

## Concept Map

| Concept | Definition | Identity | Evidence |
| --- | --- | --- | --- |
| Observation | A configured request for one or more output quantities | output name plus observable points | VERIFIED |
| Observable point | A requested quantity at one point or bounded region | coordinate range plus request value | VERIFIED |
| Domain | The time or frequency range over which an output is active | domain type and range fields | VERIFIED |
| Registered output | Runtime output state produced from one observable point | registration order and output type | VERIFIED |
| Point sample | Scalar field value at one grid coordinate and time or frequency | probe folder plus coordinate | VERIFIED |
| Wire sample | Current or charge value bound to a wire segment | segment node plus coordinate suffix | VERIFIED |
| Block sample | One integrated scalar over a rectangular boundary | bounds plus component | VERIFIED |
| Volumetric selection | Coordinates selected for a volumetric output | ordered coordinate list | VERIFIED |
| Geometry map | Exported material and face-tag topology | map bounds plus case/output prefix | VERIFIED |
| Far-field result | External far-field accumulator output | near-field bounds plus angular ranges | VERIFIED |
| Export artifact | Text, binary, HDF5, XDMF, or VTK-compatible file | exact path inside owning probe folder | VERIFIED |

## Use Case Catalog

| Use Case | Area | Actor | Trigger type |
| --- | --- | --- | --- |
| Register outputs | Output lifecycle | simulation runtime | system trigger |
| Normalise output domain | Output lifecycle | simulation runtime | system trigger |
| Register output files | Output lifecycle | simulation runtime | side effect |
| Update registered outputs | Output lifecycle | simulation runtime | time-step trigger |
| Flush registered outputs | Output lifecycle | simulation runtime | persistence trigger |
| Close frequency outputs | Output lifecycle | simulation runtime | shutdown trigger |
| Record point scalar sample | Scalar probes | simulation runtime | time-step trigger |
| Flush point scalar data | Scalar probes | simulation runtime | persistence trigger |
| Bind wire current probe | Wire probes | simulation runtime | registration trigger |
| Record wire current sample | Wire probes | simulation runtime | time-step trigger |
| Bind wire charge probe | Wire probes | simulation runtime | registration trigger |
| Record wire charge sample | Wire probes | simulation runtime | time-step trigger |
| Record block integral | Block probes | simulation runtime | time-step trigger |
| Flush block integral data | Block probes | simulation runtime | persistence trigger |
| Select volumetric coordinates | Volumetric outputs | simulation runtime | registration trigger |
| Record volumetric time data | Volumetric time output | simulation runtime | time-step trigger |
| Flush volumetric time data | Volumetric time output | simulation runtime | persistence trigger |
| Accumulate frequency slices | Volumetric frequency output | simulation runtime | time-step trigger |
| Flush frequency-slice binary | Volumetric frequency output | simulation runtime | persistence trigger |
| Write frequency magnitudes | Volumetric frequency output | simulation runtime | shutdown trigger |
| Initialise far-field output | Far-field output | simulation runtime | registration trigger |
| Update far-field output | Far-field output | simulation runtime | time-step trigger |
| Flush far-field output | Far-field output | simulation runtime | conditional persistence trigger |
| Produce geometry map | Geometry map output | simulation runtime | registration side effect |
| Write structured grid export | Structured export contracts | simulation runtime | file export trigger |
| Write unstructured grid export | Structured export contracts | simulation runtime | file export trigger |
| Write HDF5/XDMF export | Structured export contracts | simulation runtime | file export trigger |

## Cross-Cutting Rules

| Rule | Behaviour | Evidence |
| --- | --- | --- |
| Authentication | No user authentication applies to this subsystem | VERIFIED |
| Authorisation | No role or ownership checks apply to this subsystem | VERIFIED |
| Database | No database table is read or written by this subsystem | VERIFIED |
| Idempotency | Repeated flushes append or replace according to output type | VERIFIED |
| Folder ownership | Every probe-generated file is written inside its probe folder | VERIFIED |
| Ordering | Registered outputs are processed in registration order | VERIFIED |
| Time source | Output updates use the simulation time array by index | VERIFIED |
| Configuration | MPI, transmission-line, and wire-flavour behaviours depend on compile-time options | VERIFIED |
| Failure model | Several invalid states stop execution with fatal messages | VERIFIED |

## Coverage Matrix

| Area | Status | Canonical coverage |
| --- | --- | --- |
| Happy path | Covered | every feature contains at least one happy-path scenario |
| Input contract | Covered with risks | paths, domains, request values, bounds, and sampled arrays |
| Output contract | Covered | text rows, binary records, HDF5 names, XDMF names, VTK topology |
| Persistence writes | Covered | filesystem writes only; no database writes |
| Persistence reads | Not applicable | no database-backed persistence is read |
| Database enforcement | Not applicable | no database contract exists for this subsystem |
| Authorisation | Not applicable | subsystem has no authenticated actors |
| State rules | Covered | registered, pending, flushed, and closed output state |
| Failure modes | Covered with risks | fatal stops and printed no-data conditions |
| Side effects | Covered | file creation, append, replace, and external far-field calls |
| Concurrency | Risk | no locking or concurrent writer contract was verified |
| Configuration | Covered with risks | compile-time behaviour variants are specified |
| Time behaviour | Covered | simulation-indexed update and flush time rules |
| Compatibility quirks | Covered | misspelled dataset name and binary repetition are preserved |
| Evidence | Covered | each feature requirement and scenario is evidence-tagged |
