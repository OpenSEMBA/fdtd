### Requirement: Volumetric time outputs SHALL select valid grid coordinates
The system MUST select only coordinates valid for the requested volumetric
quantity. Evidence: VERIFIED.

#### Scenario: Field volume selects bounded points
- **WHEN** a volumetric electric-x field output spans `(1, 1, 1)` through
  `(3, 3, 3)` and all 27 coordinates are within field bounds
- **THEN** the selected point count is `27`
- **AND** selected coordinates are stored as x, y, z columns in traversal order
  z first, then y, then x

#### Scenario: Current volume in empty space selects no points
- **WHEN** a volumetric current output checks a coordinate that is in bounds
  but whose material is neither perfectly conducting nor a thin wire
- **THEN** the coordinate is not selected

#### Scenario: Field coordinate outside bounds is rejected
- **WHEN** a field output checks coordinate `(-1, 5, 5)` against bounds whose
  x lower limit is `0`
- **THEN** the coordinate is not selected

### Requirement: Volumetric time outputs SHALL create folder-based artifacts
The system MUST create a folder and associated binary, HDF5, and XDMF artifacts
for each volumetric time output. Evidence: VERIFIED.

#### Scenario: Current-density movie output is initialised
- **WHEN** a volumetric time output has case/output prefix `case_movie`,
  request `current-density`, and bounds `(2, 2, 2)` through `(5, 5, 5)`
- **THEN** the output folder path is `case_movie_BC_2_2_2__5_5_5`
- **AND** the binary stream path is
  `case_movie_BC_2_2_2__5_5_5/case_movie_BC_2_2_2__5_5_5.bin`
- **AND** the HDF5 path is
  `case_movie_BC_2_2_2__5_5_5/case_movie_BC_2_2_2__5_5_5.h5`
- **AND** the XDMF path is
  `case_movie_BC_2_2_2__5_5_5/case_movie_BC_2_2_2__5_5_5.xdmf`

#### Scenario: HDF5 movie datasets are initialised
- **WHEN** a volumetric time output is initialised
- **THEN** the HDF5 file contains rectilinear coordinate datasets named
  `coordsX`, `coordsY`, and `coordsZ`
- **AND** the HDF5 file contains an extendable time dataset named `times`
- **AND** current-density outputs create datasets named `CurrenDensityX`,
  `CurrenDensityY`, and `CurrenDensityZ` according to requested components
- **AND** electric-field outputs create datasets named `ElectricFieldX`,
  `ElectricFieldY`, and `ElectricFieldZ` according to requested components
- **AND** magnetic-field outputs create datasets named `MagneticFieldX`,
  `MagneticFieldY`, and `MagneticFieldZ` according to requested components

### Requirement: Volumetric time outputs SHALL record requested components
The system MUST store values for all or selected vector components according
to the requested volumetric quantity. Evidence: VERIFIED.

#### Scenario: Full current-density output is updated
- **WHEN** a current-density output is updated at time `0.1`
- **THEN** the pending time count increases by `1`
- **AND** the sample time is `0.1`
- **AND** x, y, and z current-density arrays receive curl-derived values at
  every selected current coordinate

#### Scenario: Electric x-component output is updated
- **WHEN** an electric-x volumetric output is updated at time `0.1`
- **THEN** only the x value array receives sampled field values
- **AND** y and z value arrays remain unchanged from their prior values

#### Scenario: Unsupported volumetric request is updated
- **WHEN** a volumetric time output has a request outside current-density,
  electric-field, and magnetic-field volumetric request values
- **THEN** the system stops with fatal message `Volumic measure not supported`

### Requirement: Volumetric time flush SHALL append HDF5 rows and rewrite XDMF
The system MUST persist all pending volumetric time samples and clear memory
after each flush. Evidence: VERIFIED.

#### Scenario: Movie output is flushed with pending samples
- **WHEN** a volumetric time output has two pending samples and four selected
  coordinates
- **THEN** the binary stream receives eight records
- **AND** each binary record contains time, x coordinate, y coordinate,
  z coordinate, x value, y value, and z value
- **AND** the HDF5 `times` dataset receives the two sample times
- **AND** selected value datasets receive appended rows for the requested
  components
- **AND** the XDMF file is rewritten to reference all flushed times
- **AND** the flushed-sample counter increases by `2`
- **AND** the pending time count becomes `0`

#### Scenario: Movie output is flushed without pending samples
- **WHEN** a volumetric time output has zero pending samples
- **THEN** no binary, HDF5, or XDMF data is written
- **AND** in-memory arrays are still cleared according to the requested
  component group
