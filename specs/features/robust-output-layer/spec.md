# Robust Output Layer

## Purpose

The solver MUST publish complete, discoverable probe results in serial and
distributed runs.

## Probe Metadata

WHEN any probe output is published
THEN the solver MUST NOT create a JSON probe descriptor or run output
manifest.

### Empty Probe

WHEN a configured probe records zero samples
THEN its applicable output artifacts MUST remain discoverable.

### Failed Publication

WHEN required output publication fails
THEN the solver MUST report the failure and MUST NOT report successful
completion.

## Artifact Sets

### Scalar Probes

WHEN a point, wire, bulk, line, or far-field probe is published
THEN the solver MUST publish only its applicable formatted `.dat` file.

Non-MTLN `.dat` files MUST be located beside the originating `.fdtd.json`
input.
The solver MUST NOT create a probe-specific directory for those files.

### Geometry Probes

WHEN a geometry probe is published
THEN the solver MUST publish its applicable geometry artifacts without a text
sidecar.

### Volumetric Probes

WHEN a volumetric time-domain or frequency-domain probe is published
THEN the solver MUST publish a binary artifact, visualisation metadata,
and visualisation heavy data without a JSON metadata sidecar.

## Distributed Output

WHEN a probe spans distributed partitions
THEN each sample location MUST have exactly one owner.

WHEN a geometry map spans distributed partitions
THEN each participating rank MUST publish its disjoint `.vtu` piece and the
root rank MUST publish one `.pvtu` descriptor referencing all pieces.

WHEN collective publication is available
THEN participating ranks MUST publish their disjoint regions collectively.

WHEN collective publication is unavailable
THEN the designated root MUST publish the gathered logical result.

## Filesystem Behaviour

WHEN an output path contains nested directories or spaces
THEN output creation and removal MUST succeed on supported Windows and Linux
environments.
