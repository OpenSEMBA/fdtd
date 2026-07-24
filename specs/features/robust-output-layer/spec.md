# Robust Output Layer

## Purpose

The solver MUST publish complete, discoverable probe results in serial and distributed runs.

## Probe Metadata

### Declared Probe

WHEN a configured probe is initialised
THEN the solver MUST create a machine-readable descriptor for that probe.

The descriptor MUST identify the probe, sampled quantity, domain, lifecycle state, and all declared artifacts.
Artifact references MUST be relative to the descriptor location.

### Empty Probe

WHEN a configured probe records zero samples
THEN its descriptor and declared artifacts MUST remain discoverable.

### Failed Publication

WHEN required output publication fails
THEN the descriptor MUST report failure and MUST NOT report completion.

## Artifact Sets

### Scalar And Geometry Probes

WHEN a scalar, wire, bulk, far-field, or geometry probe is published
THEN the solver MUST preserve its applicable human-readable artifact and reference it from probe metadata.

### Volumetric Probes

WHEN a volumetric time-domain or frequency-domain probe is published
THEN the solver MUST publish a binary artifact, visualisation metadata, visualisation heavy data, and probe metadata.

The binary descriptor MUST declare byte order, numeric representation, record size, and complex-value convention.

## Distributed Output

WHEN a probe spans distributed partitions
THEN each sample location MUST have exactly one owner.

WHEN collective publication is available
THEN participating ranks MUST publish their disjoint regions collectively.

WHEN collective publication is unavailable
THEN the designated root MUST publish the gathered logical result.

## Filesystem Behaviour

WHEN an output path contains nested directories or spaces
THEN output creation and removal MUST succeed on supported Windows and Linux environments.
