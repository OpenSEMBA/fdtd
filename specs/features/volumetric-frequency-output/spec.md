### Requirement: Volumetric frequency outputs SHALL create frequency folders
The system MUST create a folder and binary file for each volumetric
frequency-domain output. Evidence: VERIFIED.

#### Scenario: Frequency-slice output is initialised
- **WHEN** a volumetric frequency output has case/output prefix
  `case_frequency`, request `current-density`, bounds `(2, 2, 2)` through
  `(5, 5, 5)`, initial frequency `0.0`, final frequency `100.0`, and six
  frequency slots
- **THEN** the output folder path is `case_frequency_BC_2_2_2__5_5_5`
- **AND** the binary file path is
  `case_frequency_BC_2_2_2__5_5_5/case_frequency_BC_2_2_2__5_5_5.bin`
- **AND** six frequency values are allocated
- **AND** one complex x, y, and z value is allocated for every selected point
  at every frequency

### Requirement: Volumetric frequency outputs SHALL accumulate complex values
The system MUST accumulate complex frequency-domain values for selected
coordinates at every update. Evidence: VERIFIED.

#### Scenario: Current-density frequency output is updated
- **WHEN** a current-density frequency output is updated at simulation time `t`
- **THEN** every selected current coordinate accumulates x, y, and z
  curl-derived current-density contributions for every frequency
- **AND** the contribution for each frequency is multiplied by that
  frequency's electric-field exponential raised to time `t`

#### Scenario: Magnetic-field frequency output is updated
- **WHEN** a magnetic-field frequency output is updated at simulation time `t`
- **THEN** every selected magnetic-field coordinate accumulates x, y, and z
  sampled field contributions for every frequency
- **AND** the contribution for each frequency is multiplied by that
  frequency's magnetic-field exponential raised to time `t`

#### Scenario: Component frequency output is updated
- **WHEN** an electric-y frequency output is updated at simulation time `t`
- **THEN** only the y complex value array accumulates sampled values
- **AND** x and z complex value arrays remain unchanged

### Requirement: Volumetric frequency flush SHALL rewrite binary results
The system MUST rewrite the binary frequency-slice file on each flush. Evidence:
VERIFIED.

#### Scenario: Frequency-slice output is flushed
- **WHEN** a volumetric frequency output has `F` frequency slots and `P`
  selected coordinates
- **THEN** the binary file contains `F * P` stream records
- **AND** each record contains frequency, x coordinate, y coordinate,
  z coordinate, x complex value, y complex value, and the x complex value
  again in the third component position
- **AND** previous binary contents are replaced, not appended

### Requirement: Volumetric frequency close SHALL write HDF5 magnitudes
The system MUST write frequency-slice magnitudes to HDF5-compatible datasets
when the output is closed. Evidence: VERIFIED.

#### Scenario: Frequency-slice output is closed
- **WHEN** a volumetric frequency output is closed after accumulated values
  exist
- **THEN** rows are appended to datasets named `xVal`, `yVal`, and `zVal`
  for every frequency
- **AND** each stored value is the magnitude of the corresponding complex
  accumulated value

### Requirement: Frequency-slice compatibility quirks SHALL be preserved
The system MUST preserve existing frequency-slice output quirks unless a
future compatibility decision explicitly changes them. Evidence: VERIFIED.

#### Scenario: Binary third vector component is written
- **WHEN** a frequency-slice binary record is written
- **THEN** the seventh field repeats the x complex value
- **AND** the z complex value is not written in that binary field

#### Scenario: Field accumulation assignment is applied
- **WHEN** a sampled field value contributes to one frequency and coordinate
- **THEN** the current behaviour assigns the calculated value across the
  target complex array rather than only the addressed cell
- **AND** this behaviour is marked compatibility-sensitive
