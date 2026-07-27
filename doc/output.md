# Binary Probe Output

Every sampled probe publishes a binary artifact beside its existing text result.
Artifact metadata declares the artifact path, byte order, numeric representation,
record size, component order, and complex-value representation.

Binary values use little-endian IEEE real64 records.
Each record uses the highest native precision required by its coordinate and
measured values.
Empty binary artifacts are created when their probe is initialised.

Time scalar records contain `time,value`.
Frequency scalar records contain `frequency,value.real,value.imag`.
Bulk and line records use the time scalar layout.
Wire-current records contain time followed by current and voltage values.
Wire-charge records contain `time,charge`.

Volumetric time records contain time, three coordinates, and three values.
Volumetric frequency records contain frequency, three coordinates, and real and
imaginary values for each vector component.
Far-field records contain frequency, theta, phi, field magnitudes and phases,
and arithmetic and geometric radar-cross-section values.

Binary and text results describe the same samples in the same order.
Comparisons against text results use a tolerance matching the text format's
published precision.
