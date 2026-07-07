### Requirement: Far-field output SHALL require a time-domain observation
The system MUST initialise far-field output only for a time-domain observation.
Evidence: VERIFIED.

#### Scenario: Far-field output has a time domain
- **WHEN** a far-field output is requested with time-domain settings,
  rectangular near-field bounds, angular ranges, frequency range fields,
  a normalisation file name, face-selection settings, and decimation settings
- **THEN** the system registers the far-field output
- **AND** the output folder uses prefix `FF` and the bounded coordinate suffix
- **AND** every far-field file generated for that probe is written inside the
  far-field output folder

#### Scenario: Far-field output has a non-time domain
- **WHEN** a far-field output is requested with any domain type other than
  time-domain
- **THEN** registration stops with fatal message
  `Unexpected domain type for farField probe`
- **AND** no far-field update or flush state is initialised

### Requirement: Far-field output SHALL update from all field components
The system MUST pass the complete electric and magnetic field state to the
far-field accumulator at every update. Evidence: VERIFIED.

#### Scenario: Far-field output is updated
- **WHEN** a far-field output is registered and the update trigger is called
  for time index `n`
- **THEN** the far-field accumulator receives time index `n`, simulation
  bounds, electric x/y/z field arrays, and magnetic x/y/z field arrays

### Requirement: Far-field output SHALL flush only on explicit flush request
The system MUST flush far-field data only when the far-field flush flag is
true. Evidence: VERIFIED.

#### Scenario: Far-field flush flag is false
- **WHEN** far-field output is registered and the flush trigger is called with
  the far-field flush flag false
- **THEN** no far-field flush is performed

#### Scenario: Far-field flush flag is true
- **WHEN** far-field output is registered, the flush trigger is called with
  the far-field flush flag true, and the simulation time at the current index
  is `t`
- **THEN** the far-field writer receives `t`, simulation bounds, electric and
  magnetic grid spacing arrays, and configured far-field face selections
