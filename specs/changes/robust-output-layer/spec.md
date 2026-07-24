# Spec: Robust Output Layer

## Problem

Simulation output is not consistently reliable across supported operating systems and execution modes.
Distributed runs can allow multiple workers to target the same artifact, which can leave incomplete or corrupted results.
Probe types also expose inconsistent artifact sets, so users and downstream tools cannot reliably discover, interpret, or validate generated results.

## Goals

- Ensure every configured probe has a complete, machine-readable description of its generated artifacts.
- Ensure distributed and single-process simulations produce one logically complete result for each configured probe.
- Ensure volumetric probes consistently provide both a portable binary payload and a visualisation-ready data pair.
- Preserve the existing human-readable outputs for probe types that provide them.
- Make output creation, replacement, and finalisation reliable on supported Windows and Linux environments.
- Ensure output artifacts remain discoverable when a probe records no samples.

## Non-Goals

- Changing the physical quantities, sampling cadence, or numerical values produced by probes.
- Adding new probe request types.
- Changing input-file semantics for existing simulations.
- Supporting operating systems outside the currently supported Windows and Linux environments.
- Providing a graphical result browser.

## User Stories

- As a simulation user, I want every requested probe to publish a machine-readable descriptor so that I can reliably discover its result files and interpret their contents. [P1]
- As a simulation user, I want a distributed run to produce the same logical probe results as an equivalent single-process run so that I can scale simulations without corrupting or duplicating output. [P1]
- As a visualisation user, I want every volumetric probe to publish a complete data and metadata pair so that I can open the result without reconstructing missing files. [P1]
- As an automation author, I want artifact metadata to identify the probe, sampled quantity, bounds, domain, representation, and completion state so that post-processing can validate results without relying on file-name conventions. [P1]
- As a simulation user, I want existing text-oriented probe results to remain available so that established manual and scripted workflows continue to work. [P2]
- As a cross-platform user, I want output paths containing nested directories or spaces to work consistently so that result location does not depend on the operating system. [P2]

## Success Criteria

- Each configured probe produces exactly one readable metadata descriptor before simulation output is finalised, including when it records zero samples.
- Each descriptor identifies its probe request, sampled quantity, spatial extent where applicable, sampling domain, generated artifacts, data representation, and whether output is complete.
- Each scalar, wire, bulk, far-field, and geometry probe retains its applicable human-readable result artifact and has a descriptor that references it.
- Each volumetric time-domain and frequency-domain probe produces a binary result artifact, a visualisation-ready metadata-and-heavy-data pair, and a descriptor that references all three.
- Equivalent single-process and distributed executions expose the same logical artifact set and sample coverage for every probe; no artifact contains duplicate ownership at a process boundary.
- A distributed execution creates no conflicting concurrent writes to a shared probe artifact.
- Output creation and finalisation succeed for nested output paths, including paths with spaces, on supported Windows and Linux environments.
- The output subsystem reports a contextual failure when an artifact cannot be created, written, or finalised; it does not silently leave an artifact declared complete when it is not.

## Open Questions

- None.
