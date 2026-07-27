# Binary Probe Output

## Binary Artifacts

WHEN a supported sampled probe is initialised
THEN the solver MUST declare a binary artifact for every applicable sampling
series.

WHEN a probe records no samples
THEN its declared binary artifacts MUST remain discoverable as empty artifacts.

WHEN a probe publishes time and frequency series
THEN it MUST publish one binary artifact for each series.

## Record Contract

WHEN a binary artifact is published
THEN its descriptor MUST identify byte order, numeric representation, record
size, component order, and complex-value convention.

Binary records MUST use the highest native precision required by their values.
Binary records MUST use one documented byte order on supported platforms.

WHEN a binary record contains complex values
THEN its descriptor MUST identify their representation and order.

## Result Equivalence

WHEN a probe publishes text and binary results
THEN both MUST contain the same samples in the same order.
Values compared with text output MUST be within the precision exposed by that
text format.

## Distributed Publication

WHEN a probe result is distributed across workers
THEN every logical sample MUST have exactly one owner.
The published binary result MUST NOT contain duplicated or missing samples.
