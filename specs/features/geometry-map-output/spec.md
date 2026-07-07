### Requirement: Geometry map output SHALL be produced during registration
The system MUST create geometry-map output immediately during output
registration and MUST NOT update it during time stepping. Evidence: VERIFIED.

#### Scenario: Geometry map is requested
- **WHEN** a geometry-map request has bounds `(1, 2, 3)` through `(4, 5, 6)`
  and case/output prefix `case_map`
- **THEN** the output folder path is `case_map_MAP_1_2_3__4_5_6`
- **AND** a VTU file is written inside that folder using the same basename
- **AND** a metadata text file is written inside that folder using the same
  basename and suffix `.txt`

#### Scenario: Geometry map is updated or flushed
- **WHEN** time-step update or flush triggers are called after a geometry map
  has been produced
- **THEN** the geometry-map files remain unchanged by those triggers

### Requirement: Geometry map SHALL select material and face-tag locations
The system MUST include only geometry-map coordinates that represent relevant
material or face-tag information. Evidence: VERIFIED.

#### Scenario: Electric edge is present
- **WHEN** an electric-field coordinate lies on a detected material edge
- **THEN** that coordinate is included in the geometry map
- **AND** its stored current type is the thin-current type for that field axis

#### Scenario: Non-vacuum non-PML magnetic location is present
- **WHEN** a magnetic-field coordinate is inside bounds, belongs to a material
  other than vacuum, and is not PML
- **THEN** that coordinate is included in the geometry map
- **AND** its stored current type is the block-current type for that field axis

#### Scenario: Negative face tag applies outside PML
- **WHEN** a magnetic-field coordinate is inside bounds, has a negative
  face-tag whose bit for the field axis is set, and is not PML
- **THEN** that coordinate is included in the geometry map

### Requirement: Geometry map SHALL write unstructured grid topology
The system MUST encode thin currents as line cells and block currents as quad
cells in geometry-map output. Evidence: VERIFIED.

#### Scenario: Thin x-current coordinate is exported
- **WHEN** a geometry-map coordinate has current type `thin-current-x` at
  coordinate `(i, j, k)`
- **THEN** the unstructured grid contains a line cell from `(i, j, k)` to
  `(i + 1, j, k)`

#### Scenario: Block z-current coordinate is exported
- **WHEN** a geometry-map coordinate has current type `block-current-z` at
  coordinate `(i, j, k)`
- **THEN** the unstructured grid contains a quad with corners `(i, j, k)`,
  `(i + 1, j, k)`, `(i + 1, j + 1, k)`, and `(i, j + 1, k)`

#### Scenario: Real grid coordinates are requested
- **WHEN** map export is configured to use real grid coordinates
- **THEN** node coordinates are written from the x, y, and z grid-position
  arrays rather than raw integer grid indexes

### Requirement: Geometry map metadata SHALL preserve material tag meanings
The system MUST write the established metadata legend for geometry-map output.
Evidence: VERIFIED.

#### Scenario: Metadata is written
- **WHEN** a geometry map is produced
- **THEN** the metadata text contains the legend entries for `PEC=0`,
  `COMPO=3`, `DISPER=1`, `DIEL=2`, `SLOT=4`, `WIRE=7`,
  `WIRE-COLISION=8`, `OTHER=-1`, and border values as base value plus `0.5`
