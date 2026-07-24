# Plan: Robust Output Layer

## Architecture

The bounded context is solver-result publication.
It begins when configured probes are initialised and ends when their artifacts and run manifest are finalised.
It does not change field sampling, numerical values, or input interpretation.

The dependency direction is `infrastructure -> application -> domain`.
The domain layer defines output contracts, artifact descriptions, lifecycle states, and partition ownership rules without importing filesystem, distributed-execution, or visualisation-writer APIs.
The application layer coordinates probe lifecycle operations, selects one writer for scalar outputs, and invokes ports for artifact persistence and distributed aggregation.
The infrastructure layer adapts native filesystems, the external visualisation writer, and distributed communication to those application ports.

`output_m` becomes the sole application entry point for output lifecycle coordination.
Existing probe-specific modules retain responsibility for collecting their physical samples, but delegate artifact publication and metadata updates to the application layer.
The legacy observation implementation is removed after equivalent new-layer coverage exists.

Cross-context communication is limited to application-service calls and output lifecycle events.
No framework, filesystem, or distributed-execution type crosses into the domain layer.
No repositories are required because output publication is an append/finalise workflow rather than aggregate persistence.

## Data Models

`probe_output_t` is the aggregate for one configured probe output.
Its identity is the immutable probe identifier within a simulation run.
It owns the probe definition reference, artifact set, publication state, writer ownership, and, for volumetric data, the output partition.

`output_artifact_t` is a value object.
It has structural equality over artifact kind, relative path, representation, byte layout, and required/optional status.
Artifact kinds include text result, binary result, visualisation metadata, visualisation heavy data, geometry result, and probe metadata.

`probe_metadata_t` is the externally published representation of a `probe_output_t`.
It carries a schema version, probe identity, sampled quantity, bounds where applicable, domain, artifacts, precision and layout information, sample coverage, and lifecycle state.
Its artifact references must be relative to the probe metadata location.

`output_partition_t` is a value object describing global bounds, local owned bounds, local shape, and zero-based global offset.
It has structural equality over all coordinates and shapes.
Its invariants are disjoint ownership between participants, complete coverage of the requested volume, and no duplicate shared boundary plane.

`output_lifecycle_t` is a state value with the states declared, active, finalising, complete, and failed.
Only an active output may accept samples.
Only a complete output may advertise all required artifacts as complete.
A failed output must retain diagnostic context and must not advertise incomplete artifacts as complete.

`run_output_manifest_t` is the aggregate for all probe outputs in one simulation run.
Its identity is the simulation run identifier.
It is published by the designated run writer only after every configured probe has reached a terminal state.

## Interface Contracts

The application layer exposes these use cases:

- `declare_probe_output(probe_definition, execution_context) -> probe_output_t`
- `record_probe_sample(probe_output_t, sample_batch)`
- `flush_probe_output(probe_output_t)`
- `finalise_probe_output(probe_output_t)`
- `finalise_run_outputs(run_output_manifest_t)`

The domain-facing ports are:

- `artifact_store`: create, replace atomically, append, and inspect output artifacts using relative artifact identities.
- `metadata_publisher`: publish an initial descriptor and replace it with a later lifecycle state.
- `scalar_result_writer`: publish owner-selected scalar, wire, bulk, far-field, and geometry payloads.
- `volumetric_result_writer`: define a dataset, publish local or complete volume data, and finalise the visualisation artifact pair.
- `output_collective`: determine scalar ownership, create a probe participant group, exchange partition data, and aggregate data when collective publication is unavailable.

The external visualisation library and native filesystem sit behind `volumetric_result_writer` and `artifact_store` respectively.
Distributed execution sits behind `output_collective`.
The application layer receives status values from each port and converts failures into contextual solver errors.

The application emits `ProbeOutputDeclared`, `ProbeOutputFlushed`, `ProbeOutputFinalised`, and `ProbeOutputFailed` lifecycle events.
The run-manifest application service consumes terminal probe events and publishes the run manifest only when all declared probes are terminal.

## Implementation Phases

1. Output Contract Foundation — Define probe-output, artifact, metadata, lifecycle, and manifest models; establish the versioned metadata schema; specify required artifact sets per probe family; and add contract tests for empty, active, complete, and failed outputs.
2. Portable Artifact Infrastructure — Replace shell-dependent directory and file handling with native portable adapters; support nested paths and spaces; provide atomic metadata replacement; and make all path references relative and platform-neutral.
3. Lifecycle Coordinator and Scalar Ownership — Refactor `output_m` into the application coordinator; declare metadata for every probe during initialisation; select exactly one owner for scalar, wire, bulk, far-field, and geometry outputs; and publish their existing payloads through the common artifact contract.
4. Volumetric Serial Publication — Adapt movie and frequency probes to the external visualisation writer through the volumetric port; create the binary payload and complete visualisation pair for both domains; publish complex frequency data without discarding phase information; and correct the binary record contract.
5. Distributed Volumetric Publication — Apply the existing partition model to all volumetric components; use disjoint local regions for collective publication; make all participants execute lifecycle operations in the same order; and provide a root-aggregation publication mode when collective file publication is unavailable.
6. Run Discovery and Legacy Retirement — Publish an accurate root-owned run manifest from declared artifacts; remove the old observation dispatch path and obsolete direct writers; and preserve unrelated solver facilities only where they are outside the output bounded context.
7. Verification Matrix — Add unit and integration coverage for every probe family, zero-sample probes, nested and spaced paths, single-process execution, distributed collective execution, and distributed aggregation execution; validate artifact existence, metadata consistency, spatial coverage, and absence of duplicate boundary data.

## Decisions & Rationale

| Decision | Rationale | Alternatives Considered |
|----------|-----------|------------------------|
| Publish one versioned metadata descriptor per probe | Gives tools a stable, complete source of truth independent of naming conventions and preserves discovery for zero-sample probes | Per-run text list only; naming conventions only |
| Retain existing human-readable payloads for applicable probes | Preserves established user workflows while adding reliable machine discovery | Replace all payloads with one representation |
| Require binary and visualisation artifacts for volumetric probes | Meets volumetric portability and visualisation requirements without imposing binary payloads on scalar workflows | Binary payloads for every probe; visualisation artifacts for every probe |
| Give one owner responsibility for scalar artifacts | Prevents concurrent writes and matches the single logical result requirement | Rank-suffixed files; concurrent append operations |
| Use disjoint partition publication for volumetric artifacts | Avoids duplicate process-boundary samples and scales without reconstructing complete volumes on every participant | Full-volume reduction on every participant; independent per-rank visualisation files |
| Prefer collective publication with root aggregation fallback | Supports efficient capable environments while retaining correct results where collective publication is not available | Require collective publication; always aggregate to root |
| Isolate external writer, filesystem, and distributed execution behind ports | Protects output semantics from platform and library behaviour and preserves inward dependency direction | Call external APIs directly from probe modules |
| Retire the legacy output implementation after parity tests | Removes conflicting output semantics and eliminates a second unsupported distributed path | Maintain both implementations indefinitely |

## Constraints & Risks

The external visualisation writer requires all distributed participants to perform lifecycle calls in matching order.
The application coordinator must therefore make probe membership and finalisation decisions consistently across participants.

Root aggregation is a correctness fallback, not the preferred large-volume path.
It can require root memory proportional to a complete requested volume and must report insufficient-resource failures clearly.

Shared field planes have component-specific ownership rules.
Partition validation must cover every supported component and rank boundary before distributed publication is enabled.

Portable binary artifacts require an explicit versioned layout, numeric representation, byte order, and complex-value convention in metadata.
The visualisation data pair remains the authoritative portable representation when a consumer cannot read the binary layout.

Atomic replacement semantics and path handling differ by operating system.
The artifact-store adapter must define failure behaviour for partial finalisation and must never mark metadata complete before all required artifacts are durable.

The current codebase contains direct output calls and legacy utilities outside the new dispatcher.
Retirement must be preceded by a dependency audit so that non-output snapshot or reporting capabilities are not unintentionally removed.

No `FRAMEWORK.md` was found in the repository.
This plan nevertheless applies the supplied non-negotiable dependency rule and keeps external dependencies behind ports.
