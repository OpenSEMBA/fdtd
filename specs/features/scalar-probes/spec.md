### Requirement: Point probes SHALL record scalar field samples
The system MUST record scalar electric or magnetic field values at a single
grid coordinate for point probes. Evidence: VERIFIED.

#### Scenario: Point probe is initialised for time output
- **WHEN** a point probe has component `electric-x`, coordinate `(4, 4, 4)`,
  case root `case`, output name `pointProbe`, and a time domain
- **THEN** the probe folder is `case_pointProbe_Ex_4_4_4`
- **AND** an empty time data file exists at
  `case_pointProbe_Ex_4_4_4/case_pointProbe_Ex_4_4_4_tm.dat`
- **AND** the time buffer contains `BUFSIZE` slots initialised to `0.0`
- **AND** the value buffer contains `BUFSIZE` slots initialised to `0.0`

#### Scenario: Point probe is initialised for frequency output
- **WHEN** a point probe has component `electric-x`, coordinate `(4, 4, 4)`,
  case root `case`, output name `pointProbe`, and a frequency domain
- **THEN** the probe folder is `case_pointProbe_Ex_4_4_4`
- **AND** an empty frequency data file exists at
  `case_pointProbe_Ex_4_4_4/case_pointProbe_Ex_4_4_4_fq.dat`

#### Scenario: Point probe records two time samples
- **WHEN** a point probe at `(4, 4, 4)` observes `electric-x`, the simulation
  time values are `0.0` and `0.1`, and the sampled field values are `5.0`
  and `-4.0`
- **THEN** the probe stores sample `1` as time `0.0` and value `5.0`
- **AND** the probe stores sample `2` as time `0.1` and value `-4.0`
- **AND** no file write occurs until a flush trigger occurs

### Requirement: Point probes SHALL append time rows when flushed
The system MUST append pending point-probe time samples to the time data file
and clear only the time-domain memory. Evidence: VERIFIED.

#### Scenario: Time point data is flushed
- **WHEN** a point probe has ten pending time samples with times `1` through
  `10` and values `10` through `100` by increments of `10`
- **THEN** the time data file contains ten appended rows
- **AND** each row contains exactly two numeric columns: time and value
- **AND** all in-memory time values are reset to `0.0`
- **AND** all in-memory scalar values are reset to `0.0`
- **AND** the pending time-sample count becomes `0`

#### Scenario: Time point data is flushed twice
- **WHEN** a point probe first flushes ten rows and later flushes ten more
  time rows
- **THEN** the time data file contains twenty rows in flush order
- **AND** rows from the first flush remain present

#### Scenario: Point probe has no time data to flush
- **WHEN** a point probe has zero pending time samples at flush time
- **THEN** the system prints `No data to write.`
- **AND** the time data file receives no new rows

### Requirement: Point probes SHALL replace frequency rows when flushed
The system MUST write point-probe frequency data as the complete current
frequency result, replacing prior frequency file contents. Evidence: VERIFIED.

#### Scenario: Frequency point data is flushed
- **WHEN** a point probe has ten frequency values and ten complex accumulated
  results
- **THEN** the frequency data file contains ten rows
- **AND** each row contains exactly three numeric columns: frequency, real
  component, and imaginary component
- **AND** the frequency count remains unchanged after flush

#### Scenario: Frequency point data is flushed twice
- **WHEN** a point probe flushes frequency rows once and later flushes the
  same frequency slots with different accumulated values
- **THEN** the frequency data file contains only the later frequency rows
- **AND** earlier frequency values in the file are not retained

#### Scenario: Frequency arrays are not allocated
- **WHEN** a point probe is flushed for frequency output without allocated
  frequency values or accumulated complex values
- **THEN** the system prints `Error: arrays not allocated.`
- **AND** no frequency rows are written

### Requirement: Frequency slices SHALL be generated from domain settings
The system MUST generate frequency values from the registered frequency domain.
Evidence: VERIFIED.

#### Scenario: Linear frequency spacing is used
- **WHEN** a frequency domain has start `10.0`, step `9.0`, count `10`, and
  logarithmic spacing false
- **THEN** the generated values are `10.0`, `19.0`, `28.0`, `37.0`, `46.0`,
  `55.0`, `64.0`, `73.0`, `82.0`, and `91.0`

#### Scenario: Logarithmic frequency spacing is used
- **WHEN** a frequency domain has start `1.0`, step `0.5`, count `3`, and
  logarithmic spacing true
- **THEN** the generated values are `10.0`, `31.622776...`, and `100.0`

### Requirement: Field-prefix compatibility SHALL be preserved
The system MUST use the established output prefixes for scalar probe folders.
Evidence: VERIFIED.

#### Scenario: Electric and magnetic scalar prefixes are used
- **WHEN** a scalar probe folder is generated for `electric-x`, `electric-y`,
  `electric-z`, `magnetic-x`, `magnetic-y`, and `magnetic-z`
- **THEN** the corresponding prefixes are `Ex`, `Ey`, `Ez`, `Hx`, `Hy`,
  and `Hz`

#### Scenario: MPI axis rotation is disabled
- **WHEN** MPI axis rotation is not enabled and coordinates are `(1, 2, 3)`
- **THEN** the coordinate suffix is `1_2_3`

#### Scenario: MPI axis rotation uses direction 2
- **WHEN** MPI axis rotation is enabled, the selected direction is `2`, and
  coordinates are `(1, 2, 3)`
- **THEN** the coordinate suffix is `2_3_1`
- **AND** `electric-x`, `electric-y`, and `electric-z` prefixes become `Ez`,
  `Ex`, and `Ey`

#### Scenario: MPI axis rotation uses direction 1
- **WHEN** MPI axis rotation is enabled, the selected direction is `1`, and
  coordinates are `(1, 2, 3)`
- **THEN** the coordinate suffix is `3_1_2`
- **AND** `electric-x`, `electric-y`, and `electric-z` prefixes become `Ey`,
  `Ez`, and `Ex`
