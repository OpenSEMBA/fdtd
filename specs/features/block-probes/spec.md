### Requirement: Block probes SHALL register bounded time outputs
The system MUST create time-domain block probes for block electric-current and
magnetic-current requests. Evidence: VERIFIED.

#### Scenario: Block probe is initialised
- **WHEN** a block probe has lower coordinate `(1, 2, 3)`, upper coordinate
  `(4, 5, 6)`, case/output prefix `case_block`, and request `block-current-x`
- **THEN** the probe folder is `case_block_Jx_1_2_3__4_5_6`
- **AND** the time data file is
  `case_block_Jx_1_2_3__4_5_6/case_block_Jx_1_2_3__4_5_6_tm.dat`
- **AND** the time buffer contains `BuffObse` slots initialised to `0.0`
- **AND** the value buffer contains `BuffObse` slots initialised to `0.0`

### Requirement: Block probes SHALL integrate boundary field differences
The system MUST compute one scalar value per update by summing boundary field
differences weighted by grid spacing. Evidence: VERIFIED.

#### Scenario: X-directed electric block current is sampled
- **WHEN** a block-current-x probe spans `i1..i2`, `j1..j2`, `k1..k2`
  and is updated at time `t`
- **THEN** the stored sample time is `t`
- **AND** the stored value equals the sum over `j1..j2` of
  `(Ey(i1,j,k1-1) - Ey(i1,j,k2)) * dy(j)` plus the sum over `k1..k2` of
  `(-Ez(i1,j1-1,k) + Ez(i1,j2,k)) * dz(k)`

#### Scenario: Y-directed electric block current is sampled
- **WHEN** a block-current-y probe spans `i1..i2`, `j1..j2`, `k1..k2`
  and is updated at time `t`
- **THEN** the stored value equals the sum over `k1..k2` of
  `(-Ez(i2,j1,k) + Ez(i1-1,j1,k)) * dz(k)` plus the sum over `i1..i2`
  of `(Ex(i,j1,k2) - Ex(i,j1,k1-1)) * dx(i)`

#### Scenario: Z-directed electric block current is sampled
- **WHEN** a block-current-z probe spans `i1..i2`, `j1..j2`, `k1..k2`
  and is updated at time `t`
- **THEN** the stored value equals the sum over `i1..i2` of
  `(Ex(i,j1-1,k1) - Ex(i,j2,k1)) * dx(i)` plus the sum over `j1..j2`
  of `(-Ey(i1-1,j,k1) + Ey(i2,j,k1)) * dy(j)`

#### Scenario: X-directed magnetic block current is sampled
- **WHEN** a magnetic-block-current-x probe spans `i1..i2`, `j1..j2`,
  `k1..k2` and is updated at time `t`
- **THEN** the stored value equals the sum over `j1..j2` of
  `(-Hy(i1,j,k1) + Hy(i1,j,k2+1)) * dy(j)` plus the sum over `k1..k2`
  of `(Hz(i1,j1,k) - Hz(i1,j2+1,k)) * dz(k)`

#### Scenario: Block probe is flushed
- **WHEN** a block probe has `N` pending samples
- **THEN** the time data file receives `N` appended rows
- **AND** each row contains time and integrated value
- **AND** all pending time and value slots are reset to `0.0`
- **AND** the pending sample count becomes `0`

#### Scenario: Block probe has no data to flush
- **WHEN** a block probe has zero pending samples at flush time
- **THEN** the system prints `No data to write.`
- **AND** the time data file receives no new rows
