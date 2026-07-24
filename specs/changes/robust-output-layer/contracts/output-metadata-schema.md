# Output Metadata Schema

## Version

Schema version `1` describes one probe output.

## Required Fields

Each descriptor identifies the schema version, probe identifier, sampled quantity, domain, lifecycle state, and artifacts.
Spatial probes also identify lower and upper bounds.
Each artifact identifies its kind, relative path, and whether it is required.

Artifact paths are relative to the descriptor file.
Descriptors must not contain absolute paths.

## Lifecycle

An output is `declared` after its descriptor and required artifact set are known.
It is `active` while samples may be recorded.
It is `finalising` while required artifacts are being made durable.
It is `complete` only after every required artifact has been finalised.
It is `failed` when publication cannot complete.

Failed descriptors retain diagnostic context and must not report completion.

## Binary Artifacts

Binary artifacts declare their byte order, numeric representation, record size, and complex-value representation.
Version `1` uses little-endian records.
Complex values are stored as adjacent real and imaginary values.
Consumers must reject unsupported schema versions or binary layouts.
