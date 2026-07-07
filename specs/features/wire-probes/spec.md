### Requirement: Wire current probes SHALL bind to existing wire segments
The system MUST bind each wire current probe to a wire segment matching the
requested original node and orientation data. Evidence: VERIFIED.

#### Scenario: Matching Holland wire segment exists
- **WHEN** a wire current probe requests node `42`, current component `wire-x`,
  coordinates equal to the stored wire-segment coordinates, and wire flavour
  `holland`
- **THEN** the probe binds to the matching segment
- **AND** the probe sign is `-1` if the segment is marked reversed
- **AND** the probe sign is `1` otherwise

#### Scenario: Matching transition wire segment exists
- **WHEN** a wire current probe requests node `42`, current component `wire-x`,
  coordinates equal to the stored wire-segment coordinates, and wire flavour
  `transition`
- **THEN** the probe binds using the same rules as `holland`

#### Scenario: Matching segment is not found
- **WHEN** a wire current probe requests a node and coordinates that do not
  match any available current segment
- **THEN** registration stops with a fatal warning whose text begins
  `ERROR: WIRE probe`
- **AND** no time data for that wire probe is produced

### Requirement: Wire current probes SHALL record wire electrical values
The system MUST record six columns for each wire current sample. Evidence:
VERIFIED.

#### Scenario: Holland wire current is sampled with direct charge values
- **WHEN** a bound Holland wire current probe is updated at time `t`, the wire
  sign is `s`, wire current-past is `I`, wire electric field is `E`, segment
  length is `d`, inductance is `L`, material inverses are `MuInv` and
  `EpsInv`, current plus charge is `Qp`, current minus charge is `Qm`, and
  crank charging is enabled
- **THEN** the stored sample time is `t`
- **AND** current equals `s * I`
- **AND** delta voltage equals `-E * d`
- **AND** plus voltage equals `s * Qp * L * MuInv * EpsInv`
- **AND** minus voltage equals `s * Qm * L * MuInv * EpsInv`
- **AND** voltage difference equals plus voltage minus minus voltage

#### Scenario: Holland wire current is sampled with averaged charge values
- **WHEN** a bound Holland wire current probe is updated and crank charging is
  disabled
- **THEN** plus voltage uses the average of current and previous plus charge
- **AND** minus voltage uses the average of current and previous minus charge
- **AND** current, delta voltage, and voltage difference follow the same rules
  as the direct-charge case

#### Scenario: Wire current data is flushed
- **WHEN** a wire current probe has `N` pending samples
- **THEN** the time data file receives `N` appended rows
- **AND** each row contains time, current, delta voltage, plus voltage, minus
  voltage, and voltage difference
- **AND** the probe clears all pending sample values to `0.0`
- **AND** the pending sample count becomes `0`

### Requirement: Wire charge probes SHALL support Holland-compatible wires only
The system MUST reject charge probes for wire flavours other than `holland` and
`transition`. Evidence: VERIFIED.

#### Scenario: Charge probe uses unsupported wire flavour
- **WHEN** a wire charge probe is requested with a wire flavour other than
  `holland` or `transition`
- **THEN** registration stops with the fatal warning
  `Charge probes only supported for holland wires`
- **AND** no charge time file is produced

#### Scenario: Matching charge segment exists
- **WHEN** a wire charge probe requests a node, charge component, and
  coordinates matching a Holland-compatible current segment
- **THEN** the probe binds to that segment
- **AND** the probe sign is `-1` if the segment is marked reversed
- **AND** the probe sign is `1` otherwise

#### Scenario: Charge segment is not found
- **WHEN** a wire charge probe requests a node and coordinates that do not
  match any available charge segment
- **THEN** registration stops with a fatal warning whose text begins
  `ERROR: CHARGE probe`

### Requirement: Wire charge probes SHALL record charge values
The system MUST record time and minus-side charge for each wire charge sample.
Evidence: VERIFIED.

#### Scenario: Wire charge is sampled
- **WHEN** a bound wire charge probe is updated at time `t` and the bound
  segment has current minus-side charge `Qm`
- **THEN** the stored sample time is `t`
- **AND** the stored charge value is `Qm`

#### Scenario: Wire charge data is flushed
- **WHEN** a wire charge probe has `N` pending samples
- **THEN** the charge data file receives `N` appended rows
- **AND** each row contains time and charge
- **AND** the probe clears all pending sample values to `0.0`
- **AND** the pending sample count becomes `0`

### Requirement: Wire output folder compatibility SHALL be preserved
The system MUST include the wire segment node in wire output folders and MUST
place wire data files inside that folder. Evidence: VERIFIED.

#### Scenario: Wire folder is generated
- **WHEN** a wire probe has case/output prefix `case_wire`, field prefix `Wx`,
  coordinate suffix `1_2_3`, and segment node `42`
- **THEN** the probe folder is `case_wire_Wx_1_2_3_s42`
- **AND** the time data file path is
  `case_wire_Wx_1_2_3_s42/case_wire_Wx_1_2_3_s42_tm.dat`
