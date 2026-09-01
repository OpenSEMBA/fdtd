module xdmf_hdf5_m
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
#ifdef XDMF_HDF5_WITH_MPI
  use mpi
#endif
  use xdmf_model_m, only: xdmf_status_t, xdmf_options_t, &
    xdmf_collection_id_t, xdmf_grid_id_t, xdmf_attribute_id_t, &
    collection_record_t, grid_record_t, attribute_record_t, &
    XDMF_SUCCESS, XDMF_ERROR_ARGUMENT, XDMF_ERROR_STATE, &
    XDMF_ERROR_IO, XDMF_ERROR_HDF5, XDMF_ERROR_CONSISTENCY, &
    XDMF_SERIES_NONE, XDMF_SERIES_TIME, XDMF_SERIES_FREQUENCY, &
    XDMF_SERIES_PARAMETER, XDMF_GEOMETRY_UNIFORM, &
    XDMF_GEOMETRY_RECTILINEAR, XDMF_GEOMETRY_CURVILINEAR, &
    XDMF_GEOMETRY_UNSTRUCTURED, XDMF_TOPOLOGY_POLYVERTEX, &
    XDMF_TOPOLOGY_POLYLINE, XDMF_TOPOLOGY_POLYGON, &
    XDMF_TOPOLOGY_TRIANGLE, XDMF_TOPOLOGY_QUADRILATERAL, &
    XDMF_TOPOLOGY_TETRAHEDRON, XDMF_TOPOLOGY_PYRAMID, &
    XDMF_TOPOLOGY_WEDGE, XDMF_TOPOLOGY_HEXAHEDRON, &
    XDMF_TOPOLOGY_POLYHEDRON, XDMF_TOPOLOGY_EDGE_3, &
    XDMF_TOPOLOGY_QUADRILATERAL_9, XDMF_TOPOLOGY_TRIANGLE_6, &
    XDMF_TOPOLOGY_QUADRILATERAL_8, XDMF_TOPOLOGY_TETRAHEDRON_10, &
    XDMF_TOPOLOGY_PYRAMID_13, XDMF_TOPOLOGY_WEDGE_15, &
    XDMF_TOPOLOGY_WEDGE_18, XDMF_TOPOLOGY_HEXAHEDRON_20, &
    XDMF_TOPOLOGY_HEXAHEDRON_24, XDMF_TOPOLOGY_HEXAHEDRON_27, &
    XDMF_TOPOLOGY_MIXED, XDMF_TOPOLOGY_2D_SMESH, &
    XDMF_TOPOLOGY_2D_RECTMESH, XDMF_TOPOLOGY_2D_CORECTMESH, &
    XDMF_TOPOLOGY_3D_SMESH, XDMF_TOPOLOGY_3D_RECTMESH, &
    XDMF_TOPOLOGY_3D_CORECTMESH, XDMF_CENTER_NODE, &
    XDMF_CENTER_EDGE, XDMF_CENTER_FACE, XDMF_CENTER_CELL, &
    XDMF_CENTER_GRID, XDMF_ATTRIBUTE_SCALAR, XDMF_ATTRIBUTE_VECTOR, &
    XDMF_ATTRIBUTE_TENSOR, XDMF_ATTRIBUTE_TENSOR6, &
    XDMF_ATTRIBUTE_MATRIX, XDMF_ATTRIBUTE_GLOBAL_ID, &
    XDMF_NUMERIC_REAL32, XDMF_NUMERIC_REAL64, XDMF_NUMERIC_INT32, &
    XDMF_NUMERIC_INT64, set_status_success, set_status_error, &
    make_collection_id, make_grid_id, make_attribute_id, &
    collection_id_value, grid_id_value, attribute_id_value, &
    collection_id_owner, grid_id_owner, attribute_id_owner, &
    topology_name, topology_nodes_per_element, topology_is_supported, &
    center_name, attribute_type_name, numeric_type_name, product_int64
  use xdmf_hdf5_backend_m, only: hdf5_file_t, hdf_create_file, &
    hdf_file_is_open, hdf_close_file, hdf_flush_file, hdf_create_group, &
    hdf_write_dataset, &
    hdf_create_series_dataset, hdf_append_series, &
    hdf_append_series_hyperslab, hdf_truncate_series
  use xdmf_xml_m, only: write_xdmf_document

  implicit none

  private

  public :: XDMF_SUCCESS, XDMF_ERROR_ARGUMENT, XDMF_ERROR_STATE
  public :: XDMF_ERROR_IO, XDMF_ERROR_HDF5
  public :: XDMF_ERROR_CONSISTENCY
  public :: XDMF_SERIES_NONE, XDMF_SERIES_TIME
  public :: XDMF_SERIES_FREQUENCY, XDMF_SERIES_PARAMETER
  public :: XDMF_GEOMETRY_UNIFORM, XDMF_GEOMETRY_RECTILINEAR
  public :: XDMF_GEOMETRY_CURVILINEAR, XDMF_GEOMETRY_UNSTRUCTURED
  public :: XDMF_TOPOLOGY_POLYVERTEX, XDMF_TOPOLOGY_POLYLINE
  public :: XDMF_TOPOLOGY_POLYGON, XDMF_TOPOLOGY_TRIANGLE
  public :: XDMF_TOPOLOGY_QUADRILATERAL, XDMF_TOPOLOGY_TETRAHEDRON
  public :: XDMF_TOPOLOGY_PYRAMID, XDMF_TOPOLOGY_WEDGE
  public :: XDMF_TOPOLOGY_HEXAHEDRON, XDMF_TOPOLOGY_POLYHEDRON
  public :: XDMF_TOPOLOGY_EDGE_3, XDMF_TOPOLOGY_QUADRILATERAL_9
  public :: XDMF_TOPOLOGY_TRIANGLE_6, XDMF_TOPOLOGY_QUADRILATERAL_8
  public :: XDMF_TOPOLOGY_TETRAHEDRON_10, XDMF_TOPOLOGY_PYRAMID_13
  public :: XDMF_TOPOLOGY_WEDGE_15, XDMF_TOPOLOGY_WEDGE_18
  public :: XDMF_TOPOLOGY_HEXAHEDRON_20, XDMF_TOPOLOGY_HEXAHEDRON_24
  public :: XDMF_TOPOLOGY_HEXAHEDRON_27, XDMF_TOPOLOGY_MIXED
  public :: XDMF_TOPOLOGY_2D_SMESH, XDMF_TOPOLOGY_2D_RECTMESH
  public :: XDMF_TOPOLOGY_2D_CORECTMESH, XDMF_TOPOLOGY_3D_SMESH
  public :: XDMF_TOPOLOGY_3D_RECTMESH, XDMF_TOPOLOGY_3D_CORECTMESH
  public :: XDMF_CENTER_NODE, XDMF_CENTER_EDGE, XDMF_CENTER_FACE
  public :: XDMF_CENTER_CELL, XDMF_CENTER_GRID
  public :: XDMF_ATTRIBUTE_SCALAR, XDMF_ATTRIBUTE_VECTOR
  public :: XDMF_ATTRIBUTE_TENSOR, XDMF_ATTRIBUTE_TENSOR6
  public :: XDMF_ATTRIBUTE_MATRIX, XDMF_ATTRIBUTE_GLOBAL_ID
  public :: XDMF_NUMERIC_REAL32, XDMF_NUMERIC_REAL64
  public :: XDMF_NUMERIC_INT32, XDMF_NUMERIC_INT64
  public :: xdmf_status_t, xdmf_options_t, xdmf_collection_id_t
  public :: xdmf_grid_id_t, xdmf_attribute_id_t

  character(len=*), parameter :: SERIES_VALUES_PATH = '/series/values'
  integer(int64), save :: next_writer_token = 1_int64

  type, public :: xdmf_writer_t
    private
    type(hdf5_file_t) :: hdf5_file
    type(xdmf_options_t) :: options
    type(collection_record_t), allocatable :: collections(:)
    type(grid_record_t), allocatable :: grids(:)
    type(attribute_record_t), allocatable :: attributes(:)
    real(real64), allocatable :: series_values(:)
    character(len=:), allocatable :: xdmf_path
    character(len=:), allocatable :: hdf5_path
    character(len=:), allocatable :: hdf5_name
    integer :: next_collection_id = 1
    integer :: next_grid_id = 1
    integer :: next_attribute_id = 1
    integer :: committed_steps = 0
    integer :: communicator = 0
    integer :: rank = 0
    integer :: root_rank = 0
    integer(int64) :: owner_token = 0_int64
    real(real64) :: active_step_value = 0.0_real64
    logical :: is_open = .false.
    logical :: definitions_locked = .false.
    logical :: step_is_active = .false.
    logical :: is_poisoned = .false.
    logical :: is_collective = .false.
  contains
    procedure, private :: writer_create_with_options
    procedure, private :: writer_create_default
    generic, public :: create => writer_create_with_options, &
      writer_create_default
    procedure, public :: define_collection => writer_define_collection
    procedure, private :: writer_define_uniform_grid_r4
    procedure, private :: writer_define_uniform_grid_r8
    generic, public :: define_uniform_grid => &
      writer_define_uniform_grid_r4, writer_define_uniform_grid_r8
    procedure, private :: writer_define_rectilinear_grid_2d_r4
    procedure, private :: writer_define_rectilinear_grid_2d_r8
    procedure, private :: writer_define_rectilinear_grid_3d_r4
    procedure, private :: writer_define_rectilinear_grid_3d_r8
    generic, public :: define_rectilinear_grid => &
      writer_define_rectilinear_grid_2d_r4, &
      writer_define_rectilinear_grid_2d_r8, &
      writer_define_rectilinear_grid_3d_r4, &
      writer_define_rectilinear_grid_3d_r8
    procedure, private :: writer_define_curvilinear_grid_r4
    procedure, private :: writer_define_curvilinear_grid_r8
    generic, public :: define_curvilinear_grid => &
      writer_define_curvilinear_grid_r4, &
      writer_define_curvilinear_grid_r8
    procedure, private :: writer_define_unstructured_grid_r4
    procedure, private :: writer_define_unstructured_grid_r8
    generic, public :: define_unstructured_grid => &
      writer_define_unstructured_grid_r4, &
      writer_define_unstructured_grid_r8
    procedure, private :: writer_define_mixed_grid_r4
    procedure, private :: writer_define_mixed_grid_r8
    generic, public :: define_mixed_grid => writer_define_mixed_grid_r4, &
      writer_define_mixed_grid_r8
    procedure, private :: writer_define_attribute_with_series
    procedure, private :: writer_define_static_attribute
    generic, public :: define_attribute => &
      writer_define_attribute_with_series, writer_define_static_attribute
    procedure, private :: writer_write_attribute_r4
    procedure, private :: writer_write_attribute_r8
    procedure, private :: writer_write_attribute_i4
    procedure, private :: writer_write_attribute_i8
    generic, public :: write_attribute => writer_write_attribute_r4, &
      writer_write_attribute_r8, writer_write_attribute_i4, &
      writer_write_attribute_i8
    procedure, private :: writer_write_attribute_hyperslab_r4
    procedure, private :: writer_write_attribute_hyperslab_r8
    procedure, private :: writer_write_attribute_hyperslab_i4
    procedure, private :: writer_write_attribute_hyperslab_i8
    generic, public :: write_attribute_hyperslab => &
      writer_write_attribute_hyperslab_r4, &
      writer_write_attribute_hyperslab_r8, &
      writer_write_attribute_hyperslab_i4, &
      writer_write_attribute_hyperslab_i8
    procedure, public :: begin_step => writer_begin_step
    procedure, public :: end_step => writer_end_step
    procedure, public :: flush => writer_flush
    procedure, public :: close => writer_close
    final :: writer_finalize
  end type xdmf_writer_t

contains

  subroutine writer_create_with_options(this, path, options, status)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    type(xdmf_options_t), intent(in) :: options
    type(xdmf_status_t), intent(out) :: status

    call writer_create_impl(this, path, options, status)
  end subroutine writer_create_with_options

  subroutine writer_create_default(this, path, status)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    type(xdmf_status_t), intent(out) :: status

    type(xdmf_options_t) :: options

    call writer_create_impl(this, path, options, status)
  end subroutine writer_create_default

  subroutine writer_create_impl(this, path, options, status)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: path
    type(xdmf_options_t), intent(in) :: options
    type(xdmf_status_t), intent(out) :: status

    integer(int64) :: scalar_shape(0)
    integer :: mpi_error
    logical :: xdmf_exists

    call set_status_success(status)
    if (this%is_open) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'The XDMF writer is already open')
      return
    end if
    if (len_trim(path) == 0) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'The output path must not be empty')
      return
    end if
    call validate_options(options, status)
    if (status%is_error()) return

    call reset_writer_metadata(this)
    this%options = options
    this%is_collective = options%collective_io
    this%communicator = options%communicator
    this%root_rank = options%root_rank
#ifdef XDMF_HDF5_WITH_MPI
    if (this%is_collective) then
      call MPI_Comm_rank(this%communicator, this%rank, mpi_error)
      if (mpi_error /= MPI_SUCCESS) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Could not query the collective-output MPI communicator')
        return
      end if
    end if
#endif
    call derive_output_paths(path, this%xdmf_path, this%hdf5_path, &
      this%hdf5_name)
    if (.not. options%overwrite) then
      inquire(file=this%xdmf_path, exist=xdmf_exists)
      if (xdmf_exists) then
        call set_status_error(status, XDMF_ERROR_IO, &
          'Output file already exists: '//this%xdmf_path)
        return
      end if
    end if

    call hdf_create_file(this%hdf5_file, this%hdf5_path, options, status)
    call synchronize_collective_status(this, status, &
      'Collective HDF5 file creation failed')
    if (status%is_error()) then
      if (hdf_file_is_open(this%hdf5_file)) then
        this%is_open = .true.
        this%is_poisoned = .true.
      end if
      return
    end if
    this%is_open = .true.

    call hdf_create_group(this%hdf5_file, '/grids', status)
    if (status%is_error()) then
      call close_after_create_failure(this, status)
      return
    end if
    call hdf_create_group(this%hdf5_file, '/attributes', status)
    if (status%is_error()) then
      call close_after_create_failure(this, status)
      return
    end if
    call hdf_create_group(this%hdf5_file, '/series', status)
    if (status%is_error()) then
      call close_after_create_failure(this, status)
      return
    end if
    if (options%series_kind /= XDMF_SERIES_NONE) then
      call hdf_create_series_dataset(this%hdf5_file, SERIES_VALUES_PATH, &
        XDMF_NUMERIC_REAL64, scalar_shape, options, status)
      if (status%is_error()) call close_after_create_failure(this, status)
    end if
  end subroutine writer_create_impl

  subroutine writer_define_collection(this, name, collection_id, status)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(xdmf_collection_id_t), intent(out) :: collection_id
    type(xdmf_status_t), intent(out) :: status

    type(collection_record_t) :: record
    integer :: index

    collection_id = make_collection_id(0)
    call check_definition_state(this, name, status)
    if (status%is_error()) return
    do index = 1, size(this%collections)
      if (this%collections(index)%name == trim(name)) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Collection names must be unique')
        return
      end if
    end do

    record%id = this%next_collection_id
    record%name = trim(name)
    this%next_collection_id = this%next_collection_id + 1
    call append_collection(this, record, status)
    if (.not. status%is_error()) then
      collection_id = make_collection_id(record%id, this%owner_token)
    end if
  end subroutine writer_define_collection

  subroutine writer_define_uniform_grid_r4(this, name, dimensions, origin, &
      spacing, grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: dimensions(:)
    real(real32), intent(in) :: origin(:), spacing(:)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record

    grid_id = make_grid_id(0)
    call prepare_uniform_grid(this, name, dimensions, origin_size=size(origin), &
      spacing_size=size(spacing), record=record, status=status, &
      collection_id=collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL32

    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%origin_path, origin, &
        [int(record%dimension, int64)], status)
    end if
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%spacing_path, spacing, &
        [int(record%dimension, int64)], status)
    end if
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_uniform_grid_r4

  subroutine writer_define_uniform_grid_r8(this, name, dimensions, origin, &
      spacing, grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: dimensions(:)
    real(real64), intent(in) :: origin(:), spacing(:)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record

    grid_id = make_grid_id(0)
    call prepare_uniform_grid(this, name, dimensions, origin_size=size(origin), &
      spacing_size=size(spacing), record=record, status=status, &
      collection_id=collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL64

    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%origin_path, origin, &
        [int(record%dimension, int64)], status)
    end if
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%spacing_path, spacing, &
        [int(record%dimension, int64)], status)
    end if
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_uniform_grid_r8

  subroutine prepare_uniform_grid(this, name, dimensions, origin_size, &
      spacing_size, record, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: dimensions(:)
    integer, intent(in) :: origin_size, spacing_size
    type(grid_record_t), intent(out) :: record
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    call prepare_structured_grid(this, name, dimensions, record, status, &
      collection_id)
    if (status%is_error()) return
    if (origin_size /= record%dimension .or. &
        spacing_size /= record%dimension) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Uniform-grid origin and spacing must match the grid dimension')
      return
    end if

    record%geometry_type = XDMF_GEOMETRY_UNIFORM
    if (record%dimension == 2) then
      record%topology_type = XDMF_TOPOLOGY_2D_CORECTMESH
    else
      record%topology_type = XDMF_TOPOLOGY_3D_CORECTMESH
    end if
    record%origin_path = record%group_path//'/origin'
    record%spacing_path = record%group_path//'/spacing'
  end subroutine prepare_uniform_grid

  subroutine writer_define_rectilinear_grid_2d_r4(this, name, x, y, &
      grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: x(:), y(:)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record

    grid_id = make_grid_id(0)
    call prepare_rectilinear_grid(this, name, &
      [int(size(x), int64), int(size(y), int64)], record, status, &
      collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL32
    call write_rectilinear_group_r4(this, record, x, y, status=status)
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_rectilinear_grid_2d_r4

  subroutine writer_define_rectilinear_grid_2d_r8(this, name, x, y, &
      grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: x(:), y(:)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record

    grid_id = make_grid_id(0)
    call prepare_rectilinear_grid(this, name, &
      [int(size(x), int64), int(size(y), int64)], record, status, &
      collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL64
    call write_rectilinear_group_r8(this, record, x, y, status=status)
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_rectilinear_grid_2d_r8

  subroutine writer_define_rectilinear_grid_3d_r4(this, name, x, y, z, &
      grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: x(:), y(:), z(:)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record

    grid_id = make_grid_id(0)
    call prepare_rectilinear_grid(this, name, &
      [int(size(x), int64), int(size(y), int64), int(size(z), int64)], &
      record, status, collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL32
    call write_rectilinear_group_r4(this, record, x, y, z, status)
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_rectilinear_grid_3d_r4

  subroutine writer_define_rectilinear_grid_3d_r8(this, name, x, y, z, &
      grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: x(:), y(:), z(:)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record

    grid_id = make_grid_id(0)
    call prepare_rectilinear_grid(this, name, &
      [int(size(x), int64), int(size(y), int64), int(size(z), int64)], &
      record, status, collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL64
    call write_rectilinear_group_r8(this, record, x, y, z, status)
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_rectilinear_grid_3d_r8

  subroutine prepare_rectilinear_grid(this, name, dimensions, record, &
      status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: dimensions(:)
    type(grid_record_t), intent(out) :: record
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    integer :: path_length

    call prepare_structured_grid(this, name, dimensions, record, status, &
      collection_id)
    if (status%is_error()) return
    record%geometry_type = XDMF_GEOMETRY_RECTILINEAR
    if (record%dimension == 2) then
      record%topology_type = XDMF_TOPOLOGY_2D_RECTMESH
    else
      record%topology_type = XDMF_TOPOLOGY_3D_RECTMESH
    end if
    record%axis_sizes = dimensions
    path_length = len(record%group_path) + len('/axis_x')
    allocate(character(len=path_length) :: &
      record%axis_paths(record%dimension))
    record%axis_paths(1) = record%group_path//'/axis_x'
    record%axis_paths(2) = record%group_path//'/axis_y'
    if (record%dimension == 3) then
      record%axis_paths(3) = record%group_path//'/axis_z'
    end if
  end subroutine prepare_rectilinear_grid

  subroutine write_rectilinear_group_r4(this, record, x, y, z, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(grid_record_t), intent(in) :: record
    real(real32), intent(in) :: x(:), y(:)
    real(real32), intent(in), optional :: z(:)
    type(xdmf_status_t), intent(out) :: status

    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%axis_paths(1), x, &
        [int(size(x), int64)], status)
    end if
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%axis_paths(2), y, &
        [int(size(y), int64)], status)
    end if
    if (present(z) .and. .not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%axis_paths(3), z, &
        [int(size(z), int64)], status)
    end if
  end subroutine write_rectilinear_group_r4

  subroutine write_rectilinear_group_r8(this, record, x, y, z, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(grid_record_t), intent(in) :: record
    real(real64), intent(in) :: x(:), y(:)
    real(real64), intent(in), optional :: z(:)
    type(xdmf_status_t), intent(out) :: status

    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%axis_paths(1), x, &
        [int(size(x), int64)], status)
    end if
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%axis_paths(2), y, &
        [int(size(y), int64)], status)
    end if
    if (present(z) .and. .not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%axis_paths(3), z, &
        [int(size(z), int64)], status)
    end if
  end subroutine write_rectilinear_group_r8

  subroutine writer_define_curvilinear_grid_r4(this, name, dimensions, &
      points, grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: dimensions(:)
    real(real32), intent(in) :: points(:, :)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record
    integer(int64), allocatable :: shape(:)

    grid_id = make_grid_id(0)
    call prepare_curvilinear_grid(this, name, dimensions, size(points, 1), &
      size(points, 2), record, status, collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL32
    shape = [int(record%dimension, int64), dimensions]
    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%points_path, &
        reshape(points, [size(points)]), shape, status)
    end if
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_curvilinear_grid_r4

  subroutine writer_define_curvilinear_grid_r8(this, name, dimensions, &
      points, grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: dimensions(:)
    real(real64), intent(in) :: points(:, :)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record
    integer(int64), allocatable :: shape(:)

    grid_id = make_grid_id(0)
    call prepare_curvilinear_grid(this, name, dimensions, size(points, 1), &
      size(points, 2), record, status, collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL64
    shape = [int(record%dimension, int64), dimensions]
    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%points_path, &
        reshape(points, [size(points)]), shape, status)
    end if
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_curvilinear_grid_r8

  subroutine prepare_curvilinear_grid(this, name, dimensions, point_dimension, &
      point_count, record, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: dimensions(:)
    integer, intent(in) :: point_dimension, point_count
    type(grid_record_t), intent(out) :: record
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    call prepare_structured_grid(this, name, dimensions, record, status, &
      collection_id)
    if (status%is_error()) return
    if (point_dimension /= record%dimension .or. &
        int(point_count, int64) /= record%number_of_points) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Curvilinear points must have shape (dimension, number_of_points)')
      return
    end if

    record%geometry_type = XDMF_GEOMETRY_CURVILINEAR
    if (record%dimension == 2) then
      record%topology_type = XDMF_TOPOLOGY_2D_SMESH
    else
      record%topology_type = XDMF_TOPOLOGY_3D_SMESH
    end if
    record%points_path = record%group_path//'/points'
  end subroutine prepare_curvilinear_grid

  subroutine writer_define_unstructured_grid_r4(this, name, topology, &
      points, connectivity, grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: topology
    real(real32), intent(in) :: points(:, :)
    integer(int64), intent(in) :: connectivity(:, :)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record
    integer(int64), allocatable :: zero_connectivity(:, :)

    grid_id = make_grid_id(0)
    call prepare_unstructured_grid(this, name, topology, size(points, 1), &
      size(points, 2), connectivity, record, zero_connectivity, status, &
      collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL32
    call write_unstructured_group_r4(this, record, points, &
      zero_connectivity, status)
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_unstructured_grid_r4

  subroutine writer_define_unstructured_grid_r8(this, name, topology, &
      points, connectivity, grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: topology
    real(real64), intent(in) :: points(:, :)
    integer(int64), intent(in) :: connectivity(:, :)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record
    integer(int64), allocatable :: zero_connectivity(:, :)

    grid_id = make_grid_id(0)
    call prepare_unstructured_grid(this, name, topology, size(points, 1), &
      size(points, 2), connectivity, record, zero_connectivity, status, &
      collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL64
    call write_unstructured_group_r8(this, record, points, &
      zero_connectivity, status)
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_unstructured_grid_r8

  subroutine prepare_unstructured_grid(this, name, topology, point_dimension, &
      point_count, connectivity, record, zero_connectivity, status, &
      collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: topology, point_dimension, point_count
    integer(int64), intent(in) :: connectivity(:, :)
    type(grid_record_t), intent(out) :: record
    integer(int64), allocatable, intent(out) :: zero_connectivity(:, :)
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    integer :: expected_nodes

    call prepare_grid_record(this, name, record, status, collection_id)
    if (status%is_error()) return
    if (point_dimension /= 2 .and. point_dimension /= 3) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Unstructured points must have dimension 2 or 3')
      return
    end if
    if (point_count <= 0 .or. size(connectivity, 1) <= 0 .or. &
        size(connectivity, 2) <= 0) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Unstructured points and connectivity must not be empty')
      return
    end if
    if (.not. topology_is_supported(topology) .or. &
        topology == XDMF_TOPOLOGY_MIXED .or. &
        topology == XDMF_TOPOLOGY_POLYHEDRON .or. &
        is_structured_topology(topology)) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Unsupported fixed unstructured topology: '//topology_name(topology))
      return
    end if

    expected_nodes = topology_nodes_per_element(topology)
    if (expected_nodes > 0 .and. size(connectivity, 1) /= expected_nodes) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Connectivity node count does not match the topology')
      return
    end if
    if (topology == XDMF_TOPOLOGY_POLYLINE .and. &
        size(connectivity, 1) < 2) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'A polyline requires at least two nodes')
      return
    end if
    if (topology == XDMF_TOPOLOGY_POLYGON .and. &
        size(connectivity, 1) < 3) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'A polygon requires at least three nodes')
      return
    end if
    call validate_node_indices(connectivity, int(point_count, int64), status)
    if (status%is_error()) return

    zero_connectivity = connectivity - 1_int64
    record%geometry_type = XDMF_GEOMETRY_UNSTRUCTURED
    record%topology_type = topology
    record%dimension = point_dimension
    record%nodes_per_element = size(connectivity, 1)
    record%number_of_points = int(point_count, int64)
    record%number_of_elements = int(size(connectivity, 2), int64)
    record%points_path = record%group_path//'/points'
    record%connectivity_path = record%group_path//'/connectivity'
  end subroutine prepare_unstructured_grid

  subroutine write_unstructured_group_r4(this, record, points, connectivity, &
      status)
    class(xdmf_writer_t), intent(inout) :: this
    type(grid_record_t), intent(in) :: record
    real(real32), intent(in) :: points(:, :)
    integer(int64), intent(in) :: connectivity(:, :)
    type(xdmf_status_t), intent(out) :: status

    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%points_path, &
        reshape(points, [size(points)]), &
        [int(size(points, 1), int64), int(size(points, 2), int64)], status)
    end if
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%connectivity_path, &
        reshape(connectivity, [size(connectivity)]), &
        [int(size(connectivity, 1), int64), &
         int(size(connectivity, 2), int64)], status)
    end if
  end subroutine write_unstructured_group_r4

  subroutine write_unstructured_group_r8(this, record, points, connectivity, &
      status)
    class(xdmf_writer_t), intent(inout) :: this
    type(grid_record_t), intent(in) :: record
    real(real64), intent(in) :: points(:, :)
    integer(int64), intent(in) :: connectivity(:, :)
    type(xdmf_status_t), intent(out) :: status

    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%points_path, &
        reshape(points, [size(points)]), &
        [int(size(points, 1), int64), int(size(points, 2), int64)], status)
    end if
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%connectivity_path, &
        reshape(connectivity, [size(connectivity)]), &
        [int(size(connectivity, 1), int64), &
         int(size(connectivity, 2), int64)], status)
    end if
  end subroutine write_unstructured_group_r8

  subroutine writer_define_mixed_grid_r4(this, name, points, connectivity, &
      grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(real32), intent(in) :: points(:, :)
    integer(int64), intent(in) :: connectivity(:)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record
    integer(int64), allocatable :: zero_connectivity(:)

    grid_id = make_grid_id(0)
    call prepare_mixed_grid(this, name, size(points, 1), size(points, 2), &
      connectivity, record, zero_connectivity, status, collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL32
    call write_mixed_group_r4(this, record, points, zero_connectivity, status)
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_mixed_grid_r4

  subroutine writer_define_mixed_grid_r8(this, name, points, connectivity, &
      grid_id, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    real(real64), intent(in) :: points(:, :)
    integer(int64), intent(in) :: connectivity(:)
    type(xdmf_grid_id_t), intent(out) :: grid_id
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    type(grid_record_t) :: record
    integer(int64), allocatable :: zero_connectivity(:)

    grid_id = make_grid_id(0)
    call prepare_mixed_grid(this, name, size(points, 1), size(points, 2), &
      connectivity, record, zero_connectivity, status, collection_id)
    if (status%is_error()) return
    record%geometry_numeric_type = XDMF_NUMERIC_REAL64
    call write_mixed_group_r8(this, record, points, zero_connectivity, status)
    call finish_grid_definition(this, record, grid_id, status)
  end subroutine writer_define_mixed_grid_r8

  subroutine prepare_mixed_grid(this, name, point_dimension, point_count, &
      connectivity, record, zero_connectivity, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer, intent(in) :: point_dimension, point_count
    integer(int64), intent(in) :: connectivity(:)
    type(grid_record_t), intent(out) :: record
    integer(int64), allocatable, intent(out) :: zero_connectivity(:)
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    integer(int64) :: element_count

    call prepare_grid_record(this, name, record, status, collection_id)
    if (status%is_error()) return
    if ((point_dimension /= 2 .and. point_dimension /= 3) .or. &
        point_count <= 0) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Mixed-grid points must be a nonempty dimension-by-points array')
      return
    end if

    call convert_mixed_connectivity(connectivity, int(point_count, int64), &
      zero_connectivity, element_count, status)
    if (status%is_error()) return

    record%geometry_type = XDMF_GEOMETRY_UNSTRUCTURED
    record%topology_type = XDMF_TOPOLOGY_MIXED
    record%dimension = point_dimension
    record%number_of_points = int(point_count, int64)
    record%number_of_elements = element_count
    record%mixed_connectivity_size = int(size(connectivity), int64)
    record%points_path = record%group_path//'/points'
    record%connectivity_path = record%group_path//'/connectivity'
  end subroutine prepare_mixed_grid

  subroutine write_mixed_group_r4(this, record, points, connectivity, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(grid_record_t), intent(in) :: record
    real(real32), intent(in) :: points(:, :)
    integer(int64), intent(in) :: connectivity(:)
    type(xdmf_status_t), intent(out) :: status

    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%points_path, &
        reshape(points, [size(points)]), &
        [int(size(points, 1), int64), int(size(points, 2), int64)], status)
    end if
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%connectivity_path, &
        connectivity, [int(size(connectivity), int64)], status)
    end if
  end subroutine write_mixed_group_r4

  subroutine write_mixed_group_r8(this, record, points, connectivity, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(grid_record_t), intent(in) :: record
    real(real64), intent(in) :: points(:, :)
    integer(int64), intent(in) :: connectivity(:)
    type(xdmf_status_t), intent(out) :: status

    call hdf_create_group(this%hdf5_file, record%group_path, status)
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%points_path, &
        reshape(points, [size(points)]), &
        [int(size(points, 1), int64), int(size(points, 2), int64)], status)
    end if
    if (.not. status%is_error()) then
      call hdf_write_dataset(this%hdf5_file, record%connectivity_path, &
        connectivity, [int(size(connectivity), int64)], status)
    end if
  end subroutine write_mixed_group_r8

  subroutine writer_define_attribute_with_series(this, grid_id, name, &
      center, attribute_type, numeric_type, is_series, attribute_id, status, &
      entity_count, component_shape)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_grid_id_t), intent(in) :: grid_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: center, attribute_type, numeric_type
    logical, intent(in) :: is_series
    type(xdmf_attribute_id_t), intent(out) :: attribute_id
    type(xdmf_status_t), intent(out) :: status
    integer(int64), intent(in), optional :: entity_count
    integer(int64), intent(in), optional :: component_shape(:)

    call define_attribute_impl(this, grid_id, name, center, attribute_type, &
      numeric_type, is_series, attribute_id, status, entity_count, &
      component_shape)
  end subroutine writer_define_attribute_with_series

  subroutine writer_define_static_attribute(this, grid_id, name, center, &
      attribute_type, numeric_type, attribute_id, status, entity_count, &
      component_shape)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_grid_id_t), intent(in) :: grid_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: center, attribute_type, numeric_type
    type(xdmf_attribute_id_t), intent(out) :: attribute_id
    type(xdmf_status_t), intent(out) :: status
    integer(int64), intent(in), optional :: entity_count
    integer(int64), intent(in), optional :: component_shape(:)

    call define_attribute_impl(this, grid_id, name, center, attribute_type, &
      numeric_type, .false., attribute_id, status, entity_count, &
      component_shape)
  end subroutine writer_define_static_attribute

  subroutine define_attribute_impl(this, grid_id, name, center, &
      attribute_type, numeric_type, is_series, attribute_id, status, &
      entity_count, component_shape)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_grid_id_t), intent(in) :: grid_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: center, attribute_type, numeric_type
    logical, intent(in) :: is_series
    type(xdmf_attribute_id_t), intent(out) :: attribute_id
    type(xdmf_status_t), intent(out) :: status
    integer(int64), intent(in), optional :: entity_count
    integer(int64), intent(in), optional :: component_shape(:)

    type(attribute_record_t) :: record
    character(len=:), allocatable :: group_path
    integer :: grid_index, index

    attribute_id = make_attribute_id(0)
    call check_definition_state(this, name, status)
    if (status%is_error()) return
    grid_index = find_grid(this, grid_id)
    if (grid_index == 0) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'The attribute references an unknown grid')
      return
    end if
    if (len(center_name(center)) == 0 .or. &
        len(attribute_type_name(attribute_type)) == 0 .or. &
        len(numeric_type_name(numeric_type)) == 0) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Invalid attribute center, type, or numeric type')
      return
    end if
    if (is_series .and. this%options%series_kind == XDMF_SERIES_NONE) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Series attributes require a series-enabled writer')
      return
    end if
    if (attribute_type == XDMF_ATTRIBUTE_GLOBAL_ID .and. &
        numeric_type /= XDMF_NUMERIC_INT32 .and. &
        numeric_type /= XDMF_NUMERIC_INT64) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'GlobalID attributes require an integer numeric type')
      return
    end if
    do index = 1, size(this%attributes)
      if (this%attributes(index)%grid_id == grid_id_value(grid_id) .and. &
          this%attributes(index)%name == trim(name)) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Attribute names must be unique within a grid')
        return
      end if
    end do

    record%id = this%next_attribute_id
    record%grid_id = grid_id_value(grid_id)
    record%name = trim(name)
    record%center = center
    record%attribute_type = attribute_type
    record%numeric_type = numeric_type
    record%is_series = is_series
    record%last_step = this%committed_steps
    call infer_attribute_shape(this%grids(grid_index), center, attribute_type, &
      entity_count, component_shape, record, status)
    if (status%is_error()) return

    group_path = indexed_path('/attributes/a', record%id)
    record%dataset_path = group_path//'/values'
    this%next_attribute_id = this%next_attribute_id + 1
    call hdf_create_group(this%hdf5_file, group_path, status)
    if (status%is_error()) return
    if (is_series) then
      call hdf_create_series_dataset(this%hdf5_file, record%dataset_path, &
        numeric_type, record%storage_shape, this%options, status)
      if (status%is_error()) return
    end if

    call append_attribute(this, record, status)
    if (.not. status%is_error()) then
      attribute_id = make_attribute_id(record%id, this%owner_token)
    end if
  end subroutine define_attribute_impl

  subroutine writer_write_attribute_r4(this, attribute_id, values, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    real(real32), intent(in) :: values(:)
    type(xdmf_status_t), intent(out) :: status

    integer :: index

    call prepare_attribute_write(this, attribute_id, XDMF_NUMERIC_REAL32, &
      size(values, kind=int64), index, status)
    if (status%is_error()) return
    call write_attribute_data_r4(this, index, values, status)
    call finish_attribute_write(this, index, status)
  end subroutine writer_write_attribute_r4

  subroutine writer_write_attribute_r8(this, attribute_id, values, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    real(real64), intent(in) :: values(:)
    type(xdmf_status_t), intent(out) :: status

    integer :: index

    call prepare_attribute_write(this, attribute_id, XDMF_NUMERIC_REAL64, &
      size(values, kind=int64), index, status)
    if (status%is_error()) return
    call write_attribute_data_r8(this, index, values, status)
    call finish_attribute_write(this, index, status)
  end subroutine writer_write_attribute_r8

  subroutine writer_write_attribute_i4(this, attribute_id, values, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    integer(int32), intent(in) :: values(:)
    type(xdmf_status_t), intent(out) :: status

    integer :: index

    call prepare_attribute_write(this, attribute_id, XDMF_NUMERIC_INT32, &
      size(values, kind=int64), index, status)
    if (status%is_error()) return
    call write_attribute_data_i4(this, index, values, status)
    call finish_attribute_write(this, index, status)
  end subroutine writer_write_attribute_i4

  subroutine writer_write_attribute_i8(this, attribute_id, values, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    integer(int64), intent(in) :: values(:)
    type(xdmf_status_t), intent(out) :: status

    integer :: index

    call prepare_attribute_write(this, attribute_id, XDMF_NUMERIC_INT64, &
      size(values, kind=int64), index, status)
    if (status%is_error()) return
    call write_attribute_data_i8(this, index, values, status)
    call finish_attribute_write(this, index, status)
  end subroutine writer_write_attribute_i8

  subroutine writer_write_attribute_hyperslab_r4(this, attribute_id, values, &
      spatial_offset, spatial_count, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    real(real32), intent(in) :: values(:)
    integer(int64), intent(in) :: spatial_offset(:), spatial_count(:)
    type(xdmf_status_t), intent(out) :: status

    integer(int64), allocatable :: storage_offset(:), storage_count(:)
    integer :: index

    call prepare_attribute_hyperslab_write(this, attribute_id, &
      XDMF_NUMERIC_REAL32, size(values, kind=int64), spatial_offset, &
      spatial_count, index, storage_offset, storage_count, status)
    if (status%is_error()) return
    if (this%attributes(index)%is_series) then
    call hdf_append_series_hyperslab(this%hdf5_file, &
      this%attributes(index)%dataset_path, values, &
      this%attributes(index)%storage_shape, storage_offset, storage_count, &
      this%committed_steps, status)
    else
      call hdf_write_dataset(this%hdf5_file, &
                             this%attributes(index)%dataset_path, values, &
                             this%attributes(index)%storage_shape, status, storage_offset, &
                             storage_count)
    end if
    call synchronize_collective_status(this, status, &
      'Collective attribute hyperslab write failed')
    call finish_attribute_write(this, index, status)
  end subroutine writer_write_attribute_hyperslab_r4

  subroutine writer_write_attribute_hyperslab_r8(this, attribute_id, values, &
      spatial_offset, spatial_count, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    real(real64), intent(in) :: values(:)
    integer(int64), intent(in) :: spatial_offset(:), spatial_count(:)
    type(xdmf_status_t), intent(out) :: status

    integer(int64), allocatable :: storage_offset(:), storage_count(:)
    integer :: index

    call prepare_attribute_hyperslab_write(this, attribute_id, &
      XDMF_NUMERIC_REAL64, size(values, kind=int64), spatial_offset, &
      spatial_count, index, storage_offset, storage_count, status)
    if (status%is_error()) return
    if (this%attributes(index)%is_series) then
    call hdf_append_series_hyperslab(this%hdf5_file, &
      this%attributes(index)%dataset_path, values, &
      this%attributes(index)%storage_shape, storage_offset, storage_count, &
      this%committed_steps, status)
    else
      call hdf_write_dataset(this%hdf5_file, &
                             this%attributes(index)%dataset_path, values, &
                             this%attributes(index)%storage_shape, status, storage_offset, &
                             storage_count)
    end if
    call synchronize_collective_status(this, status, &
      'Collective attribute hyperslab write failed')
    call finish_attribute_write(this, index, status)
  end subroutine writer_write_attribute_hyperslab_r8

  subroutine writer_write_attribute_hyperslab_i4(this, attribute_id, values, &
      spatial_offset, spatial_count, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    integer(int32), intent(in) :: values(:)
    integer(int64), intent(in) :: spatial_offset(:), spatial_count(:)
    type(xdmf_status_t), intent(out) :: status

    integer(int64), allocatable :: storage_offset(:), storage_count(:)
    integer :: index

    call prepare_attribute_hyperslab_write(this, attribute_id, &
      XDMF_NUMERIC_INT32, size(values, kind=int64), spatial_offset, &
      spatial_count, index, storage_offset, storage_count, status)
    if (status%is_error()) return
    if (this%attributes(index)%is_series) then
    call hdf_append_series_hyperslab(this%hdf5_file, &
      this%attributes(index)%dataset_path, values, &
      this%attributes(index)%storage_shape, storage_offset, storage_count, &
      this%committed_steps, status)
    else
      call hdf_write_dataset(this%hdf5_file, &
                             this%attributes(index)%dataset_path, values, &
                             this%attributes(index)%storage_shape, status, storage_offset, &
                             storage_count)
    end if
    call synchronize_collective_status(this, status, &
      'Collective attribute hyperslab write failed')
    call finish_attribute_write(this, index, status)
  end subroutine writer_write_attribute_hyperslab_i4

  subroutine writer_write_attribute_hyperslab_i8(this, attribute_id, values, &
      spatial_offset, spatial_count, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    integer(int64), intent(in) :: values(:)
    integer(int64), intent(in) :: spatial_offset(:), spatial_count(:)
    type(xdmf_status_t), intent(out) :: status

    integer(int64), allocatable :: storage_offset(:), storage_count(:)
    integer :: index

    call prepare_attribute_hyperslab_write(this, attribute_id, &
      XDMF_NUMERIC_INT64, size(values, kind=int64), spatial_offset, &
      spatial_count, index, storage_offset, storage_count, status)
    if (status%is_error()) return
    if (this%attributes(index)%is_series) then
    call hdf_append_series_hyperslab(this%hdf5_file, &
      this%attributes(index)%dataset_path, values, &
      this%attributes(index)%storage_shape, storage_offset, storage_count, &
      this%committed_steps, status)
    else
      call hdf_write_dataset(this%hdf5_file, &
                             this%attributes(index)%dataset_path, values, &
                             this%attributes(index)%storage_shape, status, storage_offset, &
                             storage_count)
    end if
    call synchronize_collective_status(this, status, &
      'Collective attribute hyperslab write failed')
    call finish_attribute_write(this, index, status)
  end subroutine writer_write_attribute_hyperslab_i8

  subroutine write_attribute_data_r4(this, index, values, status)
    class(xdmf_writer_t), intent(inout) :: this
    integer, intent(in) :: index
    real(real32), intent(in) :: values(:)
    type(xdmf_status_t), intent(out) :: status

    real(real32), allocatable :: packed(:)

    if (size(this%attributes(index)%component_shape) == 2) then
      call pack_components_r4(values, this%attributes(index), packed, status)
      if (status%is_error()) return
      if (this%attributes(index)%is_series) then
        call hdf_append_series(this%hdf5_file, &
          this%attributes(index)%dataset_path, packed, &
          this%attributes(index)%storage_shape, this%committed_steps, status)
      else
        call hdf_write_dataset(this%hdf5_file, &
          this%attributes(index)%dataset_path, packed, &
          this%attributes(index)%storage_shape, status)
      end if
    else if (this%attributes(index)%is_series) then
      call hdf_append_series(this%hdf5_file, &
        this%attributes(index)%dataset_path, values, &
        this%attributes(index)%storage_shape, this%committed_steps, status)
    else
      call hdf_write_dataset(this%hdf5_file, &
        this%attributes(index)%dataset_path, values, &
        this%attributes(index)%storage_shape, status)
    end if
  end subroutine write_attribute_data_r4

  subroutine write_attribute_data_r8(this, index, values, status)
    class(xdmf_writer_t), intent(inout) :: this
    integer, intent(in) :: index
    real(real64), intent(in) :: values(:)
    type(xdmf_status_t), intent(out) :: status

    real(real64), allocatable :: packed(:)

    if (size(this%attributes(index)%component_shape) == 2) then
      call pack_components_r8(values, this%attributes(index), packed, status)
      if (status%is_error()) return
      if (this%attributes(index)%is_series) then
        call hdf_append_series(this%hdf5_file, &
          this%attributes(index)%dataset_path, packed, &
          this%attributes(index)%storage_shape, this%committed_steps, status)
      else
        call hdf_write_dataset(this%hdf5_file, &
          this%attributes(index)%dataset_path, packed, &
          this%attributes(index)%storage_shape, status)
      end if
    else if (this%attributes(index)%is_series) then
      call hdf_append_series(this%hdf5_file, &
        this%attributes(index)%dataset_path, values, &
        this%attributes(index)%storage_shape, this%committed_steps, status)
    else
      call hdf_write_dataset(this%hdf5_file, &
        this%attributes(index)%dataset_path, values, &
        this%attributes(index)%storage_shape, status)
    end if
  end subroutine write_attribute_data_r8

  subroutine write_attribute_data_i4(this, index, values, status)
    class(xdmf_writer_t), intent(inout) :: this
    integer, intent(in) :: index
    integer(int32), intent(in) :: values(:)
    type(xdmf_status_t), intent(out) :: status

    integer(int32), allocatable :: packed(:)

    if (size(this%attributes(index)%component_shape) == 2) then
      call pack_components_i4(values, this%attributes(index), packed, status)
      if (status%is_error()) return
      if (this%attributes(index)%is_series) then
        call hdf_append_series(this%hdf5_file, &
          this%attributes(index)%dataset_path, packed, &
          this%attributes(index)%storage_shape, this%committed_steps, status)
      else
        call hdf_write_dataset(this%hdf5_file, &
          this%attributes(index)%dataset_path, packed, &
          this%attributes(index)%storage_shape, status)
      end if
    else if (this%attributes(index)%is_series) then
      call hdf_append_series(this%hdf5_file, &
        this%attributes(index)%dataset_path, values, &
        this%attributes(index)%storage_shape, this%committed_steps, status)
    else
      call hdf_write_dataset(this%hdf5_file, &
        this%attributes(index)%dataset_path, values, &
        this%attributes(index)%storage_shape, status)
    end if
  end subroutine write_attribute_data_i4

  subroutine write_attribute_data_i8(this, index, values, status)
    class(xdmf_writer_t), intent(inout) :: this
    integer, intent(in) :: index
    integer(int64), intent(in) :: values(:)
    type(xdmf_status_t), intent(out) :: status

    integer(int64), allocatable :: packed(:)

    if (size(this%attributes(index)%component_shape) == 2) then
      call pack_components_i8(values, this%attributes(index), packed, status)
      if (status%is_error()) return
      if (this%attributes(index)%is_series) then
        call hdf_append_series(this%hdf5_file, &
          this%attributes(index)%dataset_path, packed, &
          this%attributes(index)%storage_shape, this%committed_steps, status)
      else
        call hdf_write_dataset(this%hdf5_file, &
          this%attributes(index)%dataset_path, packed, &
          this%attributes(index)%storage_shape, status)
      end if
    else if (this%attributes(index)%is_series) then
      call hdf_append_series(this%hdf5_file, &
        this%attributes(index)%dataset_path, values, &
        this%attributes(index)%storage_shape, this%committed_steps, status)
    else
      call hdf_write_dataset(this%hdf5_file, &
        this%attributes(index)%dataset_path, values, &
        this%attributes(index)%storage_shape, status)
    end if
  end subroutine write_attribute_data_i8

  subroutine pack_components_r4(values, attribute, packed, status)
    real(real32), intent(in) :: values(:)
    type(attribute_record_t), intent(in) :: attribute
    real(real32), allocatable, intent(out) :: packed(:)
    type(xdmf_status_t), intent(out) :: status

    integer :: allocation_status, column, entity, input_index, output_index
    integer :: row, rows, columns, values_per_entity

    allocate(packed(size(values)), stat=allocation_status)
    if (allocation_status /= 0) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Could not allocate attribute packing storage')
      return
    end if
    if (size(attribute%component_shape) /= 2) then
      packed = values
      call set_status_success(status)
      return
    end if

    rows = int(attribute%component_shape(1))
    columns = int(attribute%component_shape(2))
    values_per_entity = rows * columns
    do entity = 0, int(attribute%entity_count) - 1
      do column = 1, columns
        do row = 1, rows
          input_index = entity * values_per_entity + row + (column - 1) * rows
          output_index = entity * values_per_entity + column + (row - 1) * columns
          packed(output_index) = values(input_index)
        end do
      end do
    end do
    call set_status_success(status)
  end subroutine pack_components_r4

  subroutine pack_components_r8(values, attribute, packed, status)
    real(real64), intent(in) :: values(:)
    type(attribute_record_t), intent(in) :: attribute
    real(real64), allocatable, intent(out) :: packed(:)
    type(xdmf_status_t), intent(out) :: status

    integer :: allocation_status, column, entity, input_index, output_index
    integer :: row, rows, columns, values_per_entity

    allocate(packed(size(values)), stat=allocation_status)
    if (allocation_status /= 0) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Could not allocate attribute packing storage')
      return
    end if
    if (size(attribute%component_shape) /= 2) then
      packed = values
      call set_status_success(status)
      return
    end if

    rows = int(attribute%component_shape(1))
    columns = int(attribute%component_shape(2))
    values_per_entity = rows * columns
    do entity = 0, int(attribute%entity_count) - 1
      do column = 1, columns
        do row = 1, rows
          input_index = entity * values_per_entity + row + (column - 1) * rows
          output_index = entity * values_per_entity + column + (row - 1) * columns
          packed(output_index) = values(input_index)
        end do
      end do
    end do
    call set_status_success(status)
  end subroutine pack_components_r8

  subroutine pack_components_i4(values, attribute, packed, status)
    integer(int32), intent(in) :: values(:)
    type(attribute_record_t), intent(in) :: attribute
    integer(int32), allocatable, intent(out) :: packed(:)
    type(xdmf_status_t), intent(out) :: status

    integer :: allocation_status, column, entity, input_index, output_index
    integer :: row, rows, columns, values_per_entity

    allocate(packed(size(values)), stat=allocation_status)
    if (allocation_status /= 0) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Could not allocate attribute packing storage')
      return
    end if
    if (size(attribute%component_shape) /= 2) then
      packed = values
      call set_status_success(status)
      return
    end if

    rows = int(attribute%component_shape(1))
    columns = int(attribute%component_shape(2))
    values_per_entity = rows * columns
    do entity = 0, int(attribute%entity_count) - 1
      do column = 1, columns
        do row = 1, rows
          input_index = entity * values_per_entity + row + (column - 1) * rows
          output_index = entity * values_per_entity + column + (row - 1) * columns
          packed(output_index) = values(input_index)
        end do
      end do
    end do
    call set_status_success(status)
  end subroutine pack_components_i4

  subroutine pack_components_i8(values, attribute, packed, status)
    integer(int64), intent(in) :: values(:)
    type(attribute_record_t), intent(in) :: attribute
    integer(int64), allocatable, intent(out) :: packed(:)
    type(xdmf_status_t), intent(out) :: status

    integer :: allocation_status, column, entity, input_index, output_index
    integer :: row, rows, columns, values_per_entity

    allocate(packed(size(values)), stat=allocation_status)
    if (allocation_status /= 0) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Could not allocate attribute packing storage')
      return
    end if
    if (size(attribute%component_shape) /= 2) then
      packed = values
      call set_status_success(status)
      return
    end if

    rows = int(attribute%component_shape(1))
    columns = int(attribute%component_shape(2))
    values_per_entity = rows * columns
    do entity = 0, int(attribute%entity_count) - 1
      do column = 1, columns
        do row = 1, rows
          input_index = entity * values_per_entity + row + (column - 1) * rows
          output_index = entity * values_per_entity + column + (row - 1) * columns
          packed(output_index) = values(input_index)
        end do
      end do
    end do
    call set_status_success(status)
  end subroutine pack_components_i8

  subroutine writer_begin_step(this, value, status)
    class(xdmf_writer_t), intent(inout) :: this
    real(real64), intent(in) :: value
    type(xdmf_status_t), intent(out) :: status

    call check_open(this, status)
    if (status%is_error()) return
    if (this%options%series_kind == XDMF_SERIES_NONE) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'This writer was not created for series output')
      return
    end if
    if (this%step_is_active) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'A series step is already active')
      return
    end if

    this%definitions_locked = .true.
    this%step_is_active = .true.
    this%active_step_value = value
    call set_status_success(status)
  end subroutine writer_begin_step

  subroutine writer_end_step(this, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_status_t), intent(out) :: status

    type(xdmf_status_t) :: original_status, rollback_status
    real(real64), allocatable :: new_values(:)
    real(real64) :: step_value(1)
    integer(int64) :: scalar_shape(0)
    integer :: allocation_status, index, next_step

    call check_open(this, status)
    if (status%is_error()) return
    if (.not. this%step_is_active) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'No series step is active')
      return
    end if

    next_step = this%committed_steps + 1
    do index = 1, size(this%attributes)
      if (this%attributes(index)%is_series .and. &
          this%attributes(index)%last_step /= next_step) then
        call set_status_error(status, XDMF_ERROR_STATE, &
          'Every series attribute must be written exactly once per step')
        call rollback_active_step(this, rollback_status)
        if (rollback_status%is_error()) status = rollback_status
        return
      end if
    end do

    allocate(new_values(this%committed_steps + 1), stat=allocation_status)
    if (allocation_status /= 0) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Could not extend the in-memory series metadata')
      call rollback_active_step(this, rollback_status)
      if (rollback_status%is_error()) status = rollback_status
      return
    end if
    if (this%committed_steps > 0) then
      new_values(:this%committed_steps) = this%series_values
    end if
    new_values(this%committed_steps + 1) = this%active_step_value

    step_value(1) = this%active_step_value
    call hdf_append_series(this%hdf5_file, SERIES_VALUES_PATH, step_value, &
      scalar_shape, this%committed_steps, status)
    call synchronize_collective_status(this, status, &
      'Collective series-coordinate write failed')
    if (status%is_error()) then
      original_status = status
      if (status%error_code() == XDMF_ERROR_CONSISTENCY) then
        this%is_poisoned = .true.
      end if
      call rollback_active_step(this, rollback_status)
      if (rollback_status%is_error()) then
        status = rollback_status
      else
        status = original_status
      end if
      return
    end if

    call move_alloc(new_values, this%series_values)
    this%committed_steps = next_step
    this%step_is_active = .false.
    call set_status_success(status)
  end subroutine writer_end_step

  subroutine writer_flush(this, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_status_t), intent(out) :: status

    integer :: index

    call check_open(this, status)
    if (status%is_error()) return
    if (this%step_is_active) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Cannot flush while a series step is active')
      return
    end if
    do index = 1, size(this%attributes)
      if (.not. this%attributes(index)%is_series .and. &
          .not. this%attributes(index)%is_written) then
        call set_status_error(status, XDMF_ERROR_STATE, &
          'Every static attribute must be written before flushing')
        return
      end if
    end do

    call hdf_flush_file(this%hdf5_file, status)
    call synchronize_collective_status(this, status, &
      'Collective HDF5 flush failed')
    if (status%is_error()) return
    if (.not. this%is_collective .or. this%rank == this%root_rank) then
      call write_xdmf_document(this%xdmf_path, this%hdf5_name, &
        this%collections, this%grids, this%attributes, &
        this%options%series_kind, this%series_values, status)
    else
      call set_status_success(status)
    end if
    call synchronize_collective_status(this, status, &
      'The collective-output root could not publish XDMF metadata')
  end subroutine writer_flush

  subroutine writer_close(this, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_status_t), intent(out) :: status

    type(xdmf_status_t) :: operation_status, close_status

    call set_status_success(status)
    if (.not. this%is_open) return

    if (this%is_poisoned) then
      call hdf_close_file(this%hdf5_file, close_status)
      call synchronize_collective_status(this, close_status, &
        'Collective HDF5 close failed')
      if (close_status%is_error()) then
        status = close_status
      else
        this%is_open = .false.
        this%step_is_active = .false.
        call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
          'Writer consistency was lost; HDF5 closed without writing XDMF')
      end if
      return
    end if

    if (this%step_is_active) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Cannot close while a series step is active')
      call rollback_active_step(this, operation_status)
      if (operation_status%is_error()) status = operation_status
      return
    end if

    call writer_flush(this, operation_status)
    if (operation_status%is_error()) then
      status = operation_status
      return
    end if

    call hdf_close_file(this%hdf5_file, close_status)
    call synchronize_collective_status(this, close_status, &
      'Collective HDF5 close failed')
    if (close_status%is_error()) then
      status = close_status
    else
      this%is_open = .false.
      this%step_is_active = .false.
    end if
  end subroutine writer_close

  impure elemental subroutine writer_finalize(this)
    type(xdmf_writer_t), intent(inout) :: this

    type(xdmf_status_t) :: status

    if (this%is_open) call hdf_close_file(this%hdf5_file, status)
    if (.not. status%is_error()) this%is_open = .false.
  end subroutine writer_finalize

  subroutine prepare_structured_grid(this, name, dimensions, record, status, &
      collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    integer(int64), intent(in) :: dimensions(:)
    type(grid_record_t), intent(out) :: record
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    call prepare_grid_record(this, name, record, status, collection_id)
    if (status%is_error()) return
    if ((size(dimensions) /= 2 .and. size(dimensions) /= 3) .or. &
        any(dimensions <= 0_int64)) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Structured dimensions must contain positive I/J or I/J/K sizes')
      return
    end if

    record%dimension = size(dimensions)
    record%dimensions = dimensions
    record%number_of_points = product_int64(dimensions)
    if (record%number_of_points < 0_int64) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Structured grid dimensions overflow int64')
      return
    end if
    if (all(dimensions > 1_int64)) then
      record%number_of_elements = product_int64(dimensions - 1_int64)
      if (record%number_of_elements < 0_int64) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Structured cell dimensions overflow int64')
        return
      end if
    else
      record%number_of_elements = 0_int64
    end if
  end subroutine prepare_structured_grid

  subroutine prepare_grid_record(this, name, record, status, collection_id)
    class(xdmf_writer_t), intent(inout) :: this
    character(len=*), intent(in) :: name
    type(grid_record_t), intent(out) :: record
    type(xdmf_status_t), intent(out) :: status
    type(xdmf_collection_id_t), intent(in), optional :: collection_id

    integer :: collection_value, index

    call check_definition_state(this, name, status)
    if (status%is_error()) return
    collection_value = 0
    if (present(collection_id)) then
      if (collection_id_owner(collection_id) /= this%owner_token) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'The collection identifier belongs to another writer')
        return
      end if
      collection_value = collection_id_value(collection_id)
      if (collection_value /= 0 .and. &
          .not. collection_exists(this, collection_value)) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'The grid references an unknown collection')
        return
      end if
    end if
    do index = 1, size(this%grids)
      if (this%grids(index)%name == trim(name)) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Grid names must be unique')
        return
      end if
    end do

    record%id = this%next_grid_id
    record%collection_id = collection_value
    record%name = trim(name)
    record%group_path = indexed_path('/grids/g', record%id)
    this%next_grid_id = this%next_grid_id + 1
    call set_status_success(status)
  end subroutine prepare_grid_record

  subroutine finish_grid_definition(this, record, grid_id, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(grid_record_t), intent(in) :: record
    type(xdmf_grid_id_t), intent(inout) :: grid_id
    type(xdmf_status_t), intent(inout) :: status

    if (status%is_error()) return
    call append_grid(this, record, status)
    if (.not. status%is_error()) then
      grid_id = make_grid_id(record%id, this%owner_token)
    end if
  end subroutine finish_grid_definition

  subroutine infer_attribute_shape(grid, center, attribute_type, entity_count, &
      requested_components, attribute, status)
    type(grid_record_t), intent(in) :: grid
    integer, intent(in) :: center, attribute_type
    integer(int64), intent(in), optional :: entity_count
    integer(int64), intent(in), optional :: requested_components(:)
    type(attribute_record_t), intent(inout) :: attribute
    type(xdmf_status_t), intent(out) :: status

    integer(int64), allocatable :: spatial_shape(:), components(:)
    integer :: component_rank

    call set_status_success(status)
    select case (center)
    case (XDMF_CENTER_NODE)
      if (allocated(grid%dimensions)) then
        spatial_shape = grid%dimensions
      else
        spatial_shape = [grid%number_of_points]
      end if
    case (XDMF_CENTER_CELL)
      if (allocated(grid%dimensions)) then
        if (any(grid%dimensions <= 1_int64)) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'A cell-centred attribute requires at least one structured cell')
          return
        end if
        spatial_shape = grid%dimensions - 1_int64
      else
        spatial_shape = [grid%number_of_elements]
      end if
    case (XDMF_CENTER_GRID)
      spatial_shape = [1_int64]
    case (XDMF_CENTER_EDGE, XDMF_CENTER_FACE)
      if (.not. present(entity_count)) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Edge- and face-centred attributes require entity_count')
        return
      end if
      if (entity_count <= 0_int64) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Edge- and face-centred attributes require positive entity_count')
        return
      end if
      spatial_shape = [entity_count]
    end select
    if (any(spatial_shape <= 0_int64)) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'The attribute center has no entities on this grid')
      return
    end if
    if (present(entity_count)) then
      if (center /= XDMF_CENTER_EDGE .and. center /= XDMF_CENTER_FACE) then
        if (entity_count /= product_int64(spatial_shape)) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'entity_count conflicts with the inferred attribute shape')
          return
        end if
      end if
    end if

    if (present(requested_components)) then
      if (any(requested_components <= 0_int64)) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Attribute component dimensions must be positive')
        return
      end if
      select case (attribute_type)
      case (XDMF_ATTRIBUTE_SCALAR, XDMF_ATTRIBUTE_GLOBAL_ID)
        if (size(requested_components) /= 0) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'Scalar and GlobalID attributes cannot have component dimensions')
          return
        end if
      case (XDMF_ATTRIBUTE_VECTOR)
        if (size(requested_components) /= 1) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'Vector attributes require one component dimension')
          return
        end if
      case (XDMF_ATTRIBUTE_TENSOR)
        if (size(requested_components) /= 2) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'Tensor attributes require component shape [3, 3]')
          return
        end if
        if (any(requested_components /= [3_int64, 3_int64])) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'Tensor attributes require component shape [3, 3]')
          return
        end if
      case (XDMF_ATTRIBUTE_TENSOR6)
        if (size(requested_components) /= 1) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'Tensor6 attributes require component shape [6]')
          return
        end if
        if (requested_components(1) /= 6_int64) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'Tensor6 attributes require component shape [6]')
          return
        end if
      case (XDMF_ATTRIBUTE_MATRIX)
        if (size(requested_components) /= 2) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'Matrix attributes require two component dimensions')
          return
        end if
      end select
      components = requested_components
    else
      select case (attribute_type)
      case (XDMF_ATTRIBUTE_SCALAR, XDMF_ATTRIBUTE_GLOBAL_ID)
        allocate(components(0))
      case (XDMF_ATTRIBUTE_VECTOR)
        allocate(components(1))
        components(1) = int(grid%dimension, int64)
      case (XDMF_ATTRIBUTE_TENSOR)
        allocate(components(2))
        components = 3_int64
      case (XDMF_ATTRIBUTE_TENSOR6)
        allocate(components(1))
        components(1) = 6_int64
      case (XDMF_ATTRIBUTE_MATRIX)
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Matrix attributes require component_shape')
        return
      end select
    end if

    component_rank = size(components)
    allocate(attribute%component_shape(component_rank))
    if (component_rank > 0) attribute%component_shape = components
    allocate(attribute%storage_shape(component_rank + size(spatial_shape)))
    if (component_rank > 0) then
      attribute%storage_shape(:component_rank) = &
        components(component_rank:1:-1)
    end if
    attribute%storage_shape(component_rank + 1:) = spatial_shape
    attribute%entity_count = product_int64(spatial_shape)
    if (attribute%entity_count < 0_int64 .or. &
        product_int64(attribute%storage_shape) < 0_int64) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Attribute dimensions overflow int64')
      return
    end if
  end subroutine infer_attribute_shape

  subroutine prepare_attribute_write(this, attribute_id, numeric_type, &
      value_count, index, status)
    class(xdmf_writer_t), intent(in) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    integer, intent(in) :: numeric_type
    integer(int64), intent(in) :: value_count
    integer, intent(out) :: index
    type(xdmf_status_t), intent(out) :: status

    call check_open(this, status)
    if (status%is_error()) return
    index = find_attribute(this, attribute_id)
    if (index == 0) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Unknown attribute identifier')
      return
    end if
    if (this%attributes(index)%numeric_type /= numeric_type) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Attribute values have the wrong numeric type')
      return
    end if
    if (value_count /= product_int64(this%attributes(index)%storage_shape)) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Attribute value count does not match its defined shape')
      return
    end if

    if (this%attributes(index)%is_series) then
      if (.not. this%step_is_active) then
        call set_status_error(status, XDMF_ERROR_STATE, &
          'A series attribute can only be written inside an active step')
        return
      end if
      if (this%attributes(index)%last_step == this%committed_steps + 1) then
        call set_status_error(status, XDMF_ERROR_STATE, &
          'A series attribute can only be written once per step')
        return
      end if
    else if (this%attributes(index)%is_written) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'A static attribute can only be written once')
      return
    end if
  end subroutine prepare_attribute_write

  subroutine prepare_attribute_hyperslab_write(this, attribute_id, &
      numeric_type, value_count, spatial_offset, spatial_count, index, &
      storage_offset, storage_count, status)
    class(xdmf_writer_t), intent(in) :: this
    type(xdmf_attribute_id_t), intent(in) :: attribute_id
    integer, intent(in) :: numeric_type
    integer(int64), intent(in) :: value_count
    integer(int64), intent(in) :: spatial_offset(:), spatial_count(:)
    integer, intent(out) :: index
    integer(int64), allocatable, intent(out) :: storage_offset(:)
    integer(int64), allocatable, intent(out) :: storage_count(:)
    type(xdmf_status_t), intent(out) :: status

    integer(int64) :: expected_count
    logical :: empty_selection

    index = 0
    call check_open(this, status)
    if (.not. status%is_error() .and. .not. this%is_collective) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Attribute hyperslabs require collective HDF5 output')
    end if
    if (.not. status%is_error()) then
      index = find_attribute(this, attribute_id)
      if (index == 0) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Unknown attribute identifier')
      end if
    end if
    if (.not. status%is_error()) then
      if (this%attributes(index)%numeric_type /= numeric_type) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Attribute values have the wrong numeric type')
      else if (size(this%attributes(index)%component_shape) /= 0) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Collective hyperslabs currently support scalar attributes only')
      end if
    end if
    if (.not. status%is_error()) then
      if (size(spatial_offset) /= size(this%attributes(index)%storage_shape) .or. &
          size(spatial_count) /= size(this%attributes(index)%storage_shape)) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Hyperslab offset and count must match the spatial rank')
      else if (any(spatial_offset < 0_int64) .or. &
          any(spatial_count < 0_int64)) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Hyperslab offsets and counts must not be negative')
      end if
    end if
    if (.not. status%is_error()) then
      empty_selection = all(spatial_count == 0_int64)
      if (.not. empty_selection .and. any(spatial_count == 0_int64)) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'A hyperslab must be entirely empty or positive in every dimension')
      else if (.not. empty_selection) then
        if (any(spatial_count > this%attributes(index)%storage_shape) .or. &
            any(spatial_offset > &
              this%attributes(index)%storage_shape - spatial_count)) then
          call set_status_error(status, XDMF_ERROR_ARGUMENT, &
            'Attribute hyperslab lies outside the defined shape')
        end if
      end if
    end if
    if (.not. status%is_error()) then
      if (all(spatial_count == 0_int64)) then
        expected_count = 0_int64
      else
        expected_count = product_int64(spatial_count)
      end if
      if (expected_count < 0_int64 .or. value_count /= expected_count) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Hyperslab value count does not match its local shape')
      else if (this%attributes(index)%is_series) then
        if (.not. this%step_is_active) then
        call set_status_error(status, XDMF_ERROR_STATE, &
          'A series attribute can only be written inside an active step')
      else if (this%attributes(index)%last_step == &
          this%committed_steps + 1) then
        call set_status_error(status, XDMF_ERROR_STATE, &
          'A series attribute can only be written once per step')
      end if
      else if (this%attributes(index)%is_written) then
        call set_status_error(status, XDMF_ERROR_STATE, &
                              'A static attribute can only be written once')
      end if
    end if

    call synchronize_collective_status(this, status, &
      'Collective attribute hyperslab validation failed')
    if (status%is_error()) return
    storage_offset = spatial_offset
    storage_count = spatial_count
  end subroutine prepare_attribute_hyperslab_write

  subroutine finish_attribute_write(this, index, status)
    class(xdmf_writer_t), intent(inout) :: this
    integer, intent(in) :: index
    type(xdmf_status_t), intent(inout) :: status

    type(xdmf_status_t) :: original_status, rollback_status

    if (status%is_error()) then
      if (this%attributes(index)%is_series) then
        original_status = status
        call rollback_active_step(this, rollback_status)
        if (rollback_status%is_error()) then
          this%is_poisoned = .true.
          status = rollback_status
        else
          status = original_status
          if (status%error_code() == XDMF_ERROR_CONSISTENCY) then
            this%is_poisoned = .true.
          end if
        end if
      else if (status%error_code() == XDMF_ERROR_CONSISTENCY) then
        this%is_poisoned = .true.
      end if
      return
    end if
    if (this%attributes(index)%is_series) then
      this%attributes(index)%last_step = this%committed_steps + 1
    else
      this%attributes(index)%is_written = .true.
    end if
  end subroutine finish_attribute_write

  subroutine rollback_active_step(this, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_status_t), intent(out) :: status

    type(xdmf_status_t) :: truncate_status
    integer :: index

    call set_status_success(status)
    do index = 1, size(this%attributes)
      if (this%attributes(index)%is_series .and. &
          this%attributes(index)%last_step > this%committed_steps) then
        call hdf_truncate_series(this%hdf5_file, &
          this%attributes(index)%dataset_path, &
          this%attributes(index)%storage_shape, this%committed_steps, &
          truncate_status)
        if (truncate_status%is_error() .and. .not. status%is_error()) then
          status = truncate_status
        end if
        if (truncate_status%is_error()) then
          this%is_poisoned = .true.
        else
          this%attributes(index)%last_step = this%committed_steps
        end if
      end if
    end do
    this%step_is_active = .false.
  end subroutine rollback_active_step

  subroutine convert_mixed_connectivity(input, number_of_points, output, &
      number_of_elements, status)
    integer(int64), intent(in) :: input(:)
    integer(int64), intent(in) :: number_of_points
    integer(int64), allocatable, intent(out) :: output(:)
    integer(int64), intent(out) :: number_of_elements
    type(xdmf_status_t), intent(out) :: status

    integer(int64) :: token, count_value
    integer :: position, topology, node_count, face, face_count

    call set_status_success(status)
    number_of_elements = 0_int64
    if (size(input) == 0) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Mixed connectivity must not be empty')
      return
    end if
    output = input
    position = 1
    do while (position <= size(input))
      token = input(position)
      if (token < int(-huge(topology), int64) .or. &
          token > int(huge(topology), int64)) then
        call mixed_encoding_error(status)
        return
      end if
      topology = int(token)
      position = position + 1

      select case (topology)
      case (XDMF_TOPOLOGY_POLYVERTEX, XDMF_TOPOLOGY_POLYLINE, &
          XDMF_TOPOLOGY_POLYGON)
        if (position > size(input)) then
          call mixed_encoding_error(status)
          return
        end if
        count_value = input(position)
        if (count_value > int(huge(node_count), int64) .or. &
            count_value < 1_int64) then
          call mixed_encoding_error(status)
          return
        end if
        node_count = int(count_value)
        if ((topology == XDMF_TOPOLOGY_POLYLINE .and. node_count < 2) .or. &
            (topology == XDMF_TOPOLOGY_POLYGON .and. node_count < 3)) then
          call mixed_encoding_error(status)
          return
        end if
        position = position + 1
        call convert_mixed_nodes(input, output, position, node_count, &
          number_of_points, status)
        if (status%is_error()) return

      case (XDMF_TOPOLOGY_POLYHEDRON)
        if (position > size(input)) then
          call mixed_encoding_error(status)
          return
        end if
        if (input(position) < 1_int64 .or. &
            input(position) > int(huge(face_count), int64)) then
          call mixed_encoding_error(status)
          return
        end if
        face_count = int(input(position))
        position = position + 1
        do face = 1, face_count
          if (position > size(input)) then
            call mixed_encoding_error(status)
            return
          end if
          if (input(position) < 3_int64 .or. &
              input(position) > int(huge(node_count), int64)) then
            call mixed_encoding_error(status)
            return
          end if
          node_count = int(input(position))
          position = position + 1
          call convert_mixed_nodes(input, output, position, node_count, &
            number_of_points, status)
          if (status%is_error()) return
        end do

      case default
        if (.not. topology_is_supported(topology) .or. &
            is_structured_topology(topology) .or. &
            topology == XDMF_TOPOLOGY_MIXED) then
          call mixed_encoding_error(status)
          return
        end if
        node_count = topology_nodes_per_element(topology)
        if (node_count <= 0) then
          call mixed_encoding_error(status)
          return
        end if
        call convert_mixed_nodes(input, output, position, node_count, &
          number_of_points, status)
        if (status%is_error()) return
      end select
      number_of_elements = number_of_elements + 1_int64
    end do
  end subroutine convert_mixed_connectivity

  subroutine convert_mixed_nodes(input, output, position, node_count, &
      number_of_points, status)
    integer(int64), intent(in) :: input(:)
    integer(int64), intent(inout) :: output(:)
    integer, intent(inout) :: position
    integer, intent(in) :: node_count
    integer(int64), intent(in) :: number_of_points
    type(xdmf_status_t), intent(out) :: status

    integer :: last

    call set_status_success(status)
    if (node_count > size(input) - position + 1) then
      call mixed_encoding_error(status)
      return
    end if
    last = position + node_count - 1
    if (any(input(position:last) < 1_int64) .or. &
        any(input(position:last) > number_of_points)) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Mixed connectivity contains a node index outside the point range')
      return
    end if
    output(position:last) = input(position:last) - 1_int64
    position = last + 1
  end subroutine convert_mixed_nodes

  subroutine mixed_encoding_error(status)
    type(xdmf_status_t), intent(out) :: status

    call set_status_error(status, XDMF_ERROR_ARGUMENT, &
      'Invalid or truncated standard XDMF mixed connectivity encoding')
  end subroutine mixed_encoding_error

  subroutine validate_node_indices(connectivity, number_of_points, status)
    integer(int64), intent(in) :: connectivity(:, :)
    integer(int64), intent(in) :: number_of_points
    type(xdmf_status_t), intent(out) :: status

    if (any(connectivity < 1_int64) .or. &
        any(connectivity > number_of_points)) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Connectivity contains a node index outside the point range')
    else
      call set_status_success(status)
    end if
  end subroutine validate_node_indices

  subroutine validate_options(options, status)
    type(xdmf_options_t), intent(in) :: options
    type(xdmf_status_t), intent(out) :: status

    integer :: mpi_error, rank_count
    logical :: mpi_is_initialized

    select case (options%series_kind)
    case (XDMF_SERIES_NONE, XDMF_SERIES_TIME, XDMF_SERIES_FREQUENCY, &
        XDMF_SERIES_PARAMETER)
    case default
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Invalid series kind')
      return
    end select
    if (options%compression_level < 0 .or. options%compression_level > 9) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'The HDF5 compression level must be between 0 and 9')
      return
    end if
    if (options%chunk_target_bytes <= 0_int64) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'The HDF5 chunk target must be positive')
      return
    end if
    if (options%collective_io .and. options%compression_level /= 0) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Collective HDF5 output does not currently support compression')
      return
    end if
    if (options%collective_io) then
#ifdef XDMF_HDF5_WITH_MPI
      call MPI_Initialized(mpi_is_initialized, mpi_error)
      if (mpi_error /= MPI_SUCCESS .or. .not. mpi_is_initialized .or. &
          options%communicator == MPI_COMM_NULL) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'Collective HDF5 output requires an initialized MPI communicator')
        return
      end if
      call MPI_Comm_size(options%communicator, rank_count, mpi_error)
      if (mpi_error /= MPI_SUCCESS .or. options%root_rank < 0 .or. &
          options%root_rank >= rank_count) then
        call set_status_error(status, XDMF_ERROR_ARGUMENT, &
          'The collective-output root rank is outside the communicator')
        return
      end if
#ifdef XDMF_HDF5_WITH_PARALLEL_HDF5
#else
      call set_status_error(status, XDMF_ERROR_STATE, &
        'This library was built with MPI but the selected HDF5 is serial')
      return
#endif
#else
      call set_status_error(status, XDMF_ERROR_STATE, &
        'This library was built without MPI support')
      return
#endif
    end if
    call set_status_success(status)
  end subroutine validate_options

  subroutine check_definition_state(this, name, status)
    class(xdmf_writer_t), intent(in) :: this
    character(len=*), intent(in) :: name
    type(xdmf_status_t), intent(out) :: status

    call check_open(this, status)
    if (status%is_error()) return
    if (this%definitions_locked) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Definitions are locked after the first series step begins')
      return
    end if
    if (len_trim(name) == 0) then
      call set_status_error(status, XDMF_ERROR_ARGUMENT, &
        'Definition names must not be empty')
    end if
  end subroutine check_definition_state

  subroutine check_open(this, status)
    class(xdmf_writer_t), intent(in) :: this
    type(xdmf_status_t), intent(out) :: status

    if (this%is_open .and. .not. this%is_poisoned) then
      call set_status_success(status)
    else if (this%is_poisoned) then
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'The XDMF writer is poisoned and must be closed')
    else
      call set_status_error(status, XDMF_ERROR_STATE, &
        'The XDMF writer is not open')
    end if
  end subroutine check_open

  subroutine synchronize_collective_status(this, status, context)
    class(xdmf_writer_t), intent(in) :: this
    type(xdmf_status_t), intent(inout) :: status
    character(len=*), intent(in) :: context

#ifdef XDMF_HDF5_WITH_MPI
    integer :: local_code, global_code, mpi_error
#endif

    if (.not. this%is_collective) return
#ifdef XDMF_HDF5_WITH_MPI
    local_code = status%error_code()
    call MPI_Allreduce(local_code, global_code, 1, MPI_INTEGER, MPI_MAX, &
      this%communicator, mpi_error)
    if (mpi_error /= MPI_SUCCESS) then
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'Could not synchronize collective writer status')
    else if (global_code /= XDMF_SUCCESS .and. local_code == XDMF_SUCCESS) then
      call set_status_error(status, global_code, context)
    end if
#else
    call set_status_error(status, XDMF_ERROR_STATE, &
      'This library was built without parallel HDF5 support')
#endif
  end subroutine synchronize_collective_status

  logical function collection_exists(this, id)
    class(xdmf_writer_t), intent(in) :: this
    integer, intent(in) :: id

    integer :: index

    collection_exists = .false.
    do index = 1, size(this%collections)
      if (this%collections(index)%id == id) then
        collection_exists = .true.
        return
      end if
    end do
  end function collection_exists

  integer function find_grid(this, id)
    class(xdmf_writer_t), intent(in) :: this
    type(xdmf_grid_id_t), intent(in) :: id

    integer :: index, value

    find_grid = 0
    if (grid_id_owner(id) /= this%owner_token) return
    value = grid_id_value(id)
    do index = 1, size(this%grids)
      if (this%grids(index)%id == value) then
        find_grid = index
        return
      end if
    end do
  end function find_grid

  integer function find_attribute(this, id)
    class(xdmf_writer_t), intent(in) :: this
    type(xdmf_attribute_id_t), intent(in) :: id

    integer :: index, value

    find_attribute = 0
    if (attribute_id_owner(id) /= this%owner_token) return
    value = attribute_id_value(id)
    do index = 1, size(this%attributes)
      if (this%attributes(index)%id == value) then
        find_attribute = index
        return
      end if
    end do
  end function find_attribute

  logical function is_structured_topology(topology)
    integer, intent(in) :: topology

    select case (topology)
    case (XDMF_TOPOLOGY_2D_SMESH, XDMF_TOPOLOGY_2D_RECTMESH, &
        XDMF_TOPOLOGY_2D_CORECTMESH, XDMF_TOPOLOGY_3D_SMESH, &
        XDMF_TOPOLOGY_3D_RECTMESH, XDMF_TOPOLOGY_3D_CORECTMESH)
      is_structured_topology = .true.
    case default
      is_structured_topology = .false.
    end select
  end function is_structured_topology

  subroutine append_collection(this, record, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(collection_record_t), intent(in) :: record
    type(xdmf_status_t), intent(out) :: status

    type(collection_record_t), allocatable :: records(:)
    integer :: allocation_status, count

    count = size(this%collections)
    allocate(records(count + 1), stat=allocation_status)
    if (allocation_status /= 0) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Could not store collection metadata')
      return
    end if
    if (count > 0) records(:count) = this%collections
    records(count + 1) = record
    call move_alloc(records, this%collections)
    call set_status_success(status)
  end subroutine append_collection

  subroutine append_grid(this, record, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(grid_record_t), intent(in) :: record
    type(xdmf_status_t), intent(out) :: status

    type(grid_record_t), allocatable :: records(:)
    integer :: allocation_status, count

    count = size(this%grids)
    allocate(records(count + 1), stat=allocation_status)
    if (allocation_status /= 0) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Could not store grid metadata')
      return
    end if
    if (count > 0) records(:count) = this%grids
    records(count + 1) = record
    call move_alloc(records, this%grids)
    call set_status_success(status)
  end subroutine append_grid

  subroutine append_attribute(this, record, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(attribute_record_t), intent(in) :: record
    type(xdmf_status_t), intent(out) :: status

    type(attribute_record_t), allocatable :: records(:)
    integer :: allocation_status, count

    count = size(this%attributes)
    allocate(records(count + 1), stat=allocation_status)
    if (allocation_status /= 0) then
      call set_status_error(status, XDMF_ERROR_STATE, &
        'Could not store attribute metadata')
      return
    end if
    if (count > 0) records(:count) = this%attributes
    records(count + 1) = record
    call move_alloc(records, this%attributes)
    call set_status_success(status)
  end subroutine append_attribute

  subroutine reset_writer_metadata(this)
    class(xdmf_writer_t), intent(inout) :: this

    if (allocated(this%collections)) deallocate(this%collections)
    if (allocated(this%grids)) deallocate(this%grids)
    if (allocated(this%attributes)) deallocate(this%attributes)
    if (allocated(this%series_values)) deallocate(this%series_values)
    allocate(this%collections(0), this%grids(0), this%attributes(0))
    allocate(this%series_values(0))
    this%next_collection_id = 1
    this%next_grid_id = 1
    this%next_attribute_id = 1
    this%committed_steps = 0
    this%owner_token = next_writer_token
    if (next_writer_token == huge(next_writer_token)) then
      next_writer_token = 1_int64
    else
      next_writer_token = next_writer_token + 1_int64
    end if
    this%active_step_value = 0.0_real64
    this%communicator = 0
    this%rank = 0
    this%root_rank = 0
    this%definitions_locked = .false.
    this%step_is_active = .false.
    this%is_poisoned = .false.
    this%is_collective = .false.
  end subroutine reset_writer_metadata

  subroutine close_after_create_failure(this, status)
    class(xdmf_writer_t), intent(inout) :: this
    type(xdmf_status_t), intent(inout) :: status

    type(xdmf_status_t) :: original_status, close_status

    original_status = status
    call hdf_close_file(this%hdf5_file, close_status)
    if (close_status%is_error()) then
      this%is_open = .true.
      this%is_poisoned = .true.
      status = close_status
    else
      this%is_open = .false.
      status = original_status
    end if
  end subroutine close_after_create_failure

  subroutine derive_output_paths(path, xdmf_path, hdf5_path, hdf5_name)
    character(len=*), intent(in) :: path
    character(len=:), allocatable, intent(out) :: xdmf_path
    character(len=:), allocatable, intent(out) :: hdf5_path
    character(len=:), allocatable, intent(out) :: hdf5_name

    character(len=:), allocatable :: trimmed, stem
    integer :: separator

    trimmed = trim(path)
    if (ends_with(trimmed, '.xdmf')) then
      xdmf_path = trimmed
      stem = trimmed(:len(trimmed) - len('.xdmf'))
      hdf5_path = stem//'.h5'
    else if (ends_with(trimmed, '.xmf')) then
      xdmf_path = trimmed
      stem = trimmed(:len(trimmed) - len('.xmf'))
      hdf5_path = stem//'.h5'
    else if (ends_with(trimmed, '.hdf5')) then
      hdf5_path = trimmed
      stem = trimmed(:len(trimmed) - len('.hdf5'))
      xdmf_path = stem//'.xdmf'
    else if (ends_with(trimmed, '.h5')) then
      hdf5_path = trimmed
      stem = trimmed(:len(trimmed) - len('.h5'))
      xdmf_path = stem//'.xdmf'
    else
      xdmf_path = trimmed//'.xdmf'
      hdf5_path = trimmed//'.h5'
    end if

    separator = max(index(hdf5_path, '/', back=.true.), &
      index(hdf5_path, achar(92), back=.true.))
    hdf5_name = hdf5_path(separator + 1:)
  end subroutine derive_output_paths

  logical function ends_with(value, suffix)
    character(len=*), intent(in) :: value, suffix

    if (len(value) < len(suffix)) then
      ends_with = .false.
    else
      ends_with = value(len(value) - len(suffix) + 1:) == suffix
    end if
  end function ends_with

  function indexed_path(prefix, id) result(path)
    character(len=*), intent(in) :: prefix
    integer, intent(in) :: id
    character(len=:), allocatable :: path

    character(len=32) :: number

    if (id < 10000) then
      write(number, '(I4.4)') id
    else
      write(number, '(I0)') id
    end if
    path = prefix//trim(number)
  end function indexed_path

end module xdmf_hdf5_m
