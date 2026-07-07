### Requirement: Structured grid exports SHALL validate point shape
The system MUST reject structured-grid point arrays whose shape does not match
the declared grid dimensions. Evidence: VERIFIED.

#### Scenario: Structured grid receives valid points
- **WHEN** a structured grid declares dimensions `nx`, `ny`, and `nz`, and
  receives a point array with first dimension `3` and second dimension
  `nx * ny * nz`
- **THEN** those points replace any previously stored points

#### Scenario: Structured grid receives wrong first dimension
- **WHEN** a structured grid receives a point array whose first dimension is
  not `3`
- **THEN** export stops with fatal message
  `add_points_structured: first dim must be 3`

#### Scenario: Structured grid receives wrong point count
- **WHEN** a structured grid declares dimensions `2`, `3`, and `4`, and
  receives a point array whose second dimension is not `24`
- **THEN** export stops with fatal message
  `add_points_structured: wrong number of points`

### Requirement: Unstructured grid exports SHALL preserve explicit topology
The system MUST export unstructured grids using the provided points,
connectivity, offsets, and cell types. Evidence: VERIFIED.

#### Scenario: Unstructured grid receives points
- **WHEN** an unstructured grid receives a point array with first dimension
  `3` and second dimension `P`
- **THEN** the stored point count is `P`
- **AND** those points replace any previously stored points

#### Scenario: Unstructured grid receives topology
- **WHEN** an unstructured grid receives connectivity, offsets, and cell-type
  arrays
- **THEN** those arrays are stored unchanged
- **AND** the stored cell count equals the number of offsets

#### Scenario: Unstructured grid receives wrong first point dimension
- **WHEN** an unstructured grid receives a point array whose first dimension
  is not `3`
- **THEN** export stops with fatal message
  `add_points_unstructured: first dim must be 3`

### Requirement: Grid exports SHALL support point and cell data arrays
The system MUST allow scalar and vector data arrays on points and cells.
Evidence: VERIFIED.

#### Scenario: Scalar data is added
- **WHEN** a scalar data array named `temperature` is added with `N` numeric
  values
- **THEN** the stored data array has name `temperature`, component count `1`,
  and exactly the provided `N` values

#### Scenario: Vector data is added
- **WHEN** a vector data array named `field` is added with numeric values
- **THEN** the stored data array has name `field`, component count `3`, and
  exactly the provided values

### Requirement: HDF5-backed exports SHALL preserve dataset names and ranks
The system MUST write HDF5 datasets and XDMF references using established names
and dimensions. Evidence: VERIFIED.

#### Scenario: Rectilinear coordinates are created
- **WHEN** an HDF5-backed rectilinear export is created with x, y, and z grid
  coordinate arrays
- **THEN** datasets named `coordsX`, `coordsY`, and `coordsZ` exist
- **AND** each dataset stores the corresponding coordinate array as
  double-precision values

#### Scenario: Time dataset is created
- **WHEN** an HDF5-backed time-dependent export is initialised with chunk size
  `BUFSIZE`
- **THEN** an extendable one-dimensional dataset named `times` exists

#### Scenario: XDMF scalar attribute is written
- **WHEN** an XDMF scalar attribute named `ElectricFieldX` references an HDF5
  file and four-dimensional dataset shape `(t, z, y, x)`
- **THEN** the XDMF attribute has name `ElectricFieldX`, centre `Node`, type
  `Scalar`, and an HDF data item pointing to `/ElectricFieldX`

### Requirement: Export API tests SHALL define compatibility coverage
The system MUST preserve behaviours covered by direct export tests. Evidence:
VERIFIED.

#### Scenario: VTK export compatibility is checked
- **WHEN** the direct export tests are enabled
- **THEN** tests cover point allocation, point scalar data, point vector data,
  cell scalar data, cell vector data, structured-grid file creation,
  unstructured-grid file creation, unstructured-grid cell data, and basic
  structured and unstructured file contents

#### Scenario: HDF5 and XDMF compatibility is checked
- **WHEN** the direct export tests are enabled
- **THEN** tests cover HDF5 file creation, one-dimensional dataset writing,
  two-dimensional dataset writing, three-dimensional dataset writing, XDMF file
  creation, and XDMF output that references HDF5 data
