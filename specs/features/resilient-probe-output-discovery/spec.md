# Resilient Probe Output Discovery

## Run Isolation

WHEN a caller discovers artifacts through a solver instance
THEN discovery MUST be rooted in that solver's simulation run directory.

WHEN the process working directory changes after solver construction
THEN previously constructed solver instances MUST continue to discover their
own artifacts.

Discovered artifact paths MUST be absolute and directly usable by callers.

## Probe Identity

WHEN a requested name matches a configured probe
THEN discovery MUST treat the name as literal text rather than as a regular
expression.

WHEN a legacy caller supplies a generated artifact prefix instead of a
configured probe name
THEN discovery SHOULD preserve literal prefix discovery for compatible
artifacts.

Consumers MUST NOT assign semantic meaning to artifact list position.
A consumer expecting one artifact MUST verify uniqueness.
A consumer of a multi-artifact probe MUST select the required artifact by
identity, marker, or suffix and verify uniqueness.

Legacy artifact names MAY be ambiguous when a probe name, output tag, or
transmission-line bundle name shares the same prefix.
Discovery MUST return matching candidates rather than guess artifact ownership.

## Publication Layouts

WHEN coordinated probe artifacts are published in a probe-specific directory
THEN discovery MUST find their supported payload artifacts recursively.

WHEN transmission-line artifacts are published directly in the run directory
THEN discovery MUST find those flat artifacts through the same interface.

WHEN a far-field text artifact has no filename extension
THEN discovery MUST include it without including probe descriptors or unrelated
file types.

## Artifact Filtering

WHEN binary-inclusive discovery is not requested
THEN raw binary artifacts MUST be excluded.

WHEN binary-inclusive discovery is requested
THEN applicable raw binary artifacts MUST be included.

Supported human-readable, visualisation, and heavy-data artifacts MUST remain
discoverable by default.

## Filename Interpretation

WHEN a probe filename contains a current time-domain or frequency-domain marker
THEN the marker MUST NOT be interpreted as part of coordinates, bounds, or wire
segment identity.

WHEN a probe filename contains the legacy frequency-domain marker
THEN its established interpretation MUST remain supported.

WHEN a probe path contains generated directory names
THEN filename semantics MUST be derived from the final basename rather than
from parent-directory text.
