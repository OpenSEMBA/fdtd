module xdmf_model_m
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64

  implicit none

  private

  integer, parameter, public :: XDMF_SUCCESS = 0
  integer, parameter, public :: XDMF_ERROR_ARGUMENT = 1
  integer, parameter, public :: XDMF_ERROR_STATE = 2
  integer, parameter, public :: XDMF_ERROR_IO = 3
  integer, parameter, public :: XDMF_ERROR_HDF5 = 4
  integer, parameter, public :: XDMF_ERROR_CONSISTENCY = 5

  integer, parameter, public :: XDMF_SERIES_NONE = 0
  integer, parameter, public :: XDMF_SERIES_TIME = 1
  integer, parameter, public :: XDMF_SERIES_FREQUENCY = 2
  integer, parameter, public :: XDMF_SERIES_PARAMETER = 3

  integer, parameter, public :: XDMF_GEOMETRY_UNIFORM = 1
  integer, parameter, public :: XDMF_GEOMETRY_RECTILINEAR = 2
  integer, parameter, public :: XDMF_GEOMETRY_CURVILINEAR = 3
  integer, parameter, public :: XDMF_GEOMETRY_UNSTRUCTURED = 4

  integer, parameter, public :: XDMF_TOPOLOGY_POLYVERTEX = 1
  integer, parameter, public :: XDMF_TOPOLOGY_POLYLINE = 2
  integer, parameter, public :: XDMF_TOPOLOGY_POLYGON = 3
  integer, parameter, public :: XDMF_TOPOLOGY_TRIANGLE = 4
  integer, parameter, public :: XDMF_TOPOLOGY_QUADRILATERAL = 5
  integer, parameter, public :: XDMF_TOPOLOGY_TETRAHEDRON = 6
  integer, parameter, public :: XDMF_TOPOLOGY_PYRAMID = 7
  integer, parameter, public :: XDMF_TOPOLOGY_WEDGE = 8
  integer, parameter, public :: XDMF_TOPOLOGY_HEXAHEDRON = 9
  integer, parameter, public :: XDMF_TOPOLOGY_POLYHEDRON = 16
  integer, parameter, public :: XDMF_TOPOLOGY_EDGE_3 = 34
  integer, parameter, public :: XDMF_TOPOLOGY_QUADRILATERAL_9 = 35
  integer, parameter, public :: XDMF_TOPOLOGY_TRIANGLE_6 = 36
  integer, parameter, public :: XDMF_TOPOLOGY_QUADRILATERAL_8 = 37
  integer, parameter, public :: XDMF_TOPOLOGY_TETRAHEDRON_10 = 38
  integer, parameter, public :: XDMF_TOPOLOGY_PYRAMID_13 = 39
  integer, parameter, public :: XDMF_TOPOLOGY_WEDGE_15 = 40
  integer, parameter, public :: XDMF_TOPOLOGY_WEDGE_18 = 41
  integer, parameter, public :: XDMF_TOPOLOGY_HEXAHEDRON_20 = 48
  integer, parameter, public :: XDMF_TOPOLOGY_HEXAHEDRON_24 = 49
  integer, parameter, public :: XDMF_TOPOLOGY_HEXAHEDRON_27 = 50
  integer, parameter, public :: XDMF_TOPOLOGY_MIXED = 100
  integer, parameter, public :: XDMF_TOPOLOGY_2D_SMESH = 201
  integer, parameter, public :: XDMF_TOPOLOGY_2D_RECTMESH = 202
  integer, parameter, public :: XDMF_TOPOLOGY_2D_CORECTMESH = 203
  integer, parameter, public :: XDMF_TOPOLOGY_3D_SMESH = 204
  integer, parameter, public :: XDMF_TOPOLOGY_3D_RECTMESH = 205
  integer, parameter, public :: XDMF_TOPOLOGY_3D_CORECTMESH = 206

  integer, parameter, public :: XDMF_CENTER_NODE = 1
  integer, parameter, public :: XDMF_CENTER_EDGE = 2
  integer, parameter, public :: XDMF_CENTER_FACE = 3
  integer, parameter, public :: XDMF_CENTER_CELL = 4
  integer, parameter, public :: XDMF_CENTER_GRID = 5

  integer, parameter, public :: XDMF_ATTRIBUTE_SCALAR = 1
  integer, parameter, public :: XDMF_ATTRIBUTE_VECTOR = 2
  integer, parameter, public :: XDMF_ATTRIBUTE_TENSOR = 3
  integer, parameter, public :: XDMF_ATTRIBUTE_TENSOR6 = 4
  integer, parameter, public :: XDMF_ATTRIBUTE_MATRIX = 5
  integer, parameter, public :: XDMF_ATTRIBUTE_GLOBAL_ID = 6

  integer, parameter, public :: XDMF_NUMERIC_REAL32 = 1
  integer, parameter, public :: XDMF_NUMERIC_REAL64 = 2
  integer, parameter, public :: XDMF_NUMERIC_INT32 = 3
  integer, parameter, public :: XDMF_NUMERIC_INT64 = 4

  type, public :: xdmf_status_t
    private
    integer :: code = XDMF_SUCCESS
    character(len=:), allocatable :: detail
  contains
    procedure, public :: is_error => status_is_error
    procedure, public :: error_code => status_error_code
    procedure, public :: message => status_message
  end type xdmf_status_t

  type, public :: xdmf_options_t
    logical :: overwrite = .false.
    integer :: series_kind = XDMF_SERIES_NONE
    integer :: compression_level = 0
    integer(int64) :: chunk_target_bytes = 1048576_int64
    logical :: collective_io = .false.
    integer :: communicator = 0
    integer :: root_rank = 0
  end type xdmf_options_t

  type, public :: xdmf_collection_id_t
    private
    integer :: value = 0
    integer(int64) :: owner = 0_int64
  end type xdmf_collection_id_t

  type, public :: xdmf_grid_id_t
    private
    integer :: value = 0
    integer(int64) :: owner = 0_int64
  end type xdmf_grid_id_t

  type, public :: xdmf_attribute_id_t
    private
    integer :: value = 0
    integer(int64) :: owner = 0_int64
  end type xdmf_attribute_id_t

  type, public :: collection_record_t
    integer :: id = 0
    character(len=:), allocatable :: name
  end type collection_record_t

  type, public :: grid_record_t
    integer :: id = 0
    integer :: collection_id = 0
    integer :: geometry_type = 0
    integer :: topology_type = 0
    integer :: dimension = 0
    integer :: geometry_numeric_type = XDMF_NUMERIC_REAL64
    integer :: topology_numeric_type = XDMF_NUMERIC_INT64
    integer :: nodes_per_element = 0
    integer(int64) :: number_of_points = 0_int64
    integer(int64) :: number_of_elements = 0_int64
    integer(int64) :: mixed_connectivity_size = 0_int64
    integer(int64), allocatable :: dimensions(:)
    integer(int64), allocatable :: axis_sizes(:)
    character(len=:), allocatable :: name
    character(len=:), allocatable :: group_path
    character(len=:), allocatable :: origin_path
    character(len=:), allocatable :: spacing_path
    character(len=:), allocatable :: points_path
    character(len=:), allocatable :: connectivity_path
    character(len=:), allocatable :: axis_paths(:)
  end type grid_record_t

  type, public :: attribute_record_t
    integer :: id = 0
    integer :: grid_id = 0
    integer :: center = 0
    integer :: attribute_type = 0
    integer :: numeric_type = 0
    integer :: last_step = 0
    logical :: is_series = .false.
    logical :: is_written = .false.
    integer(int64) :: entity_count = 0_int64
    integer(int64), allocatable :: component_shape(:)
    integer(int64), allocatable :: storage_shape(:)
    character(len=:), allocatable :: name
    character(len=:), allocatable :: dataset_path
  end type attribute_record_t

  public :: set_status_success
  public :: set_status_error
  public :: make_collection_id
  public :: make_grid_id
  public :: make_attribute_id
  public :: collection_id_value
  public :: grid_id_value
  public :: attribute_id_value
  public :: collection_id_owner
  public :: grid_id_owner
  public :: attribute_id_owner
  public :: topology_name
  public :: topology_nodes_per_element
  public :: topology_is_supported
  public :: center_name
  public :: attribute_type_name
  public :: numeric_type_name
  public :: numeric_type_precision
  public :: numeric_type_size
  public :: product_int64

contains

  logical function status_is_error(this)
    class(xdmf_status_t), intent(in) :: this

    status_is_error = this%code /= XDMF_SUCCESS
  end function status_is_error

  integer function status_error_code(this)
    class(xdmf_status_t), intent(in) :: this

    status_error_code = this%code
  end function status_error_code

  function status_message(this) result(message)
    class(xdmf_status_t), intent(in) :: this
    character(len=:), allocatable :: message

    if (allocated(this%detail)) then
      message = this%detail
    else
      message = ''
    end if
  end function status_message

  subroutine set_status_success(status)
    type(xdmf_status_t), intent(out) :: status

    status%code = XDMF_SUCCESS
    status%detail = ''
  end subroutine set_status_success

  subroutine set_status_error(status, code, message)
    type(xdmf_status_t), intent(out) :: status
    integer, intent(in) :: code
    character(len=*), intent(in) :: message

    status%code = code
    status%detail = trim(message)
  end subroutine set_status_error

  function make_collection_id(value, owner) result(id)
    integer, intent(in) :: value
    integer(int64), intent(in), optional :: owner
    type(xdmf_collection_id_t) :: id

    id%value = value
    if (present(owner)) id%owner = owner
  end function make_collection_id

  function make_grid_id(value, owner) result(id)
    integer, intent(in) :: value
    integer(int64), intent(in), optional :: owner
    type(xdmf_grid_id_t) :: id

    id%value = value
    if (present(owner)) id%owner = owner
  end function make_grid_id

  function make_attribute_id(value, owner) result(id)
    integer, intent(in) :: value
    integer(int64), intent(in), optional :: owner
    type(xdmf_attribute_id_t) :: id

    id%value = value
    if (present(owner)) id%owner = owner
  end function make_attribute_id

  integer function collection_id_value(id)
    type(xdmf_collection_id_t), intent(in) :: id

    collection_id_value = id%value
  end function collection_id_value

  integer function grid_id_value(id)
    type(xdmf_grid_id_t), intent(in) :: id

    grid_id_value = id%value
  end function grid_id_value

  integer function attribute_id_value(id)
    type(xdmf_attribute_id_t), intent(in) :: id

    attribute_id_value = id%value
  end function attribute_id_value

  integer(int64) function collection_id_owner(id)
    type(xdmf_collection_id_t), intent(in) :: id

    collection_id_owner = id%owner
  end function collection_id_owner

  integer(int64) function grid_id_owner(id)
    type(xdmf_grid_id_t), intent(in) :: id

    grid_id_owner = id%owner
  end function grid_id_owner

  integer(int64) function attribute_id_owner(id)
    type(xdmf_attribute_id_t), intent(in) :: id

    attribute_id_owner = id%owner
  end function attribute_id_owner

  logical function topology_is_supported(topology)
    integer, intent(in) :: topology

    topology_is_supported = len(topology_name(topology)) > 0
  end function topology_is_supported

  function topology_name(topology) result(name)
    integer, intent(in) :: topology
    character(len=:), allocatable :: name

    select case (topology)
    case (XDMF_TOPOLOGY_POLYVERTEX); name = 'Polyvertex'
    case (XDMF_TOPOLOGY_POLYLINE); name = 'Polyline'
    case (XDMF_TOPOLOGY_POLYGON); name = 'Polygon'
    case (XDMF_TOPOLOGY_TRIANGLE); name = 'Triangle'
    case (XDMF_TOPOLOGY_QUADRILATERAL); name = 'Quadrilateral'
    case (XDMF_TOPOLOGY_TETRAHEDRON); name = 'Tetrahedron'
    case (XDMF_TOPOLOGY_PYRAMID); name = 'Pyramid'
    case (XDMF_TOPOLOGY_WEDGE); name = 'Wedge'
    case (XDMF_TOPOLOGY_HEXAHEDRON); name = 'Hexahedron'
    case (XDMF_TOPOLOGY_POLYHEDRON); name = 'Polyhedron'
    case (XDMF_TOPOLOGY_EDGE_3); name = 'Edge_3'
    case (XDMF_TOPOLOGY_QUADRILATERAL_9); name = 'Quadrilateral_9'
    case (XDMF_TOPOLOGY_TRIANGLE_6); name = 'Triangle_6'
    case (XDMF_TOPOLOGY_QUADRILATERAL_8); name = 'Quadrilateral_8'
    case (XDMF_TOPOLOGY_TETRAHEDRON_10); name = 'Tetrahedron_10'
    case (XDMF_TOPOLOGY_PYRAMID_13); name = 'Pyramid_13'
    case (XDMF_TOPOLOGY_WEDGE_15); name = 'Wedge_15'
    case (XDMF_TOPOLOGY_WEDGE_18); name = 'Wedge_18'
    case (XDMF_TOPOLOGY_HEXAHEDRON_20); name = 'Hexahedron_20'
    case (XDMF_TOPOLOGY_HEXAHEDRON_24); name = 'Hexahedron_24'
    case (XDMF_TOPOLOGY_HEXAHEDRON_27); name = 'Hexahedron_27'
    case (XDMF_TOPOLOGY_MIXED); name = 'Mixed'
    case (XDMF_TOPOLOGY_2D_SMESH); name = '2DSMesh'
    case (XDMF_TOPOLOGY_2D_RECTMESH); name = '2DRectMesh'
    case (XDMF_TOPOLOGY_2D_CORECTMESH); name = '2DCoRectMesh'
    case (XDMF_TOPOLOGY_3D_SMESH); name = '3DSMesh'
    case (XDMF_TOPOLOGY_3D_RECTMESH); name = '3DRectMesh'
    case (XDMF_TOPOLOGY_3D_CORECTMESH); name = '3DCoRectMesh'
    case default; name = ''
    end select
  end function topology_name

  integer function topology_nodes_per_element(topology)
    integer, intent(in) :: topology

    select case (topology)
    case (XDMF_TOPOLOGY_POLYVERTEX); topology_nodes_per_element = 1
    case (XDMF_TOPOLOGY_TRIANGLE); topology_nodes_per_element = 3
    case (XDMF_TOPOLOGY_QUADRILATERAL); topology_nodes_per_element = 4
    case (XDMF_TOPOLOGY_TETRAHEDRON); topology_nodes_per_element = 4
    case (XDMF_TOPOLOGY_PYRAMID); topology_nodes_per_element = 5
    case (XDMF_TOPOLOGY_WEDGE); topology_nodes_per_element = 6
    case (XDMF_TOPOLOGY_HEXAHEDRON); topology_nodes_per_element = 8
    case (XDMF_TOPOLOGY_EDGE_3); topology_nodes_per_element = 3
    case (XDMF_TOPOLOGY_QUADRILATERAL_9); topology_nodes_per_element = 9
    case (XDMF_TOPOLOGY_TRIANGLE_6); topology_nodes_per_element = 6
    case (XDMF_TOPOLOGY_QUADRILATERAL_8); topology_nodes_per_element = 8
    case (XDMF_TOPOLOGY_TETRAHEDRON_10); topology_nodes_per_element = 10
    case (XDMF_TOPOLOGY_PYRAMID_13); topology_nodes_per_element = 13
    case (XDMF_TOPOLOGY_WEDGE_15); topology_nodes_per_element = 15
    case (XDMF_TOPOLOGY_WEDGE_18); topology_nodes_per_element = 18
    case (XDMF_TOPOLOGY_HEXAHEDRON_20); topology_nodes_per_element = 20
    case (XDMF_TOPOLOGY_HEXAHEDRON_24); topology_nodes_per_element = 24
    case (XDMF_TOPOLOGY_HEXAHEDRON_27); topology_nodes_per_element = 27
    case default; topology_nodes_per_element = 0
    end select
  end function topology_nodes_per_element

  function center_name(center) result(name)
    integer, intent(in) :: center
    character(len=:), allocatable :: name

    select case (center)
    case (XDMF_CENTER_NODE); name = 'Node'
    case (XDMF_CENTER_EDGE); name = 'Edge'
    case (XDMF_CENTER_FACE); name = 'Face'
    case (XDMF_CENTER_CELL); name = 'Cell'
    case (XDMF_CENTER_GRID); name = 'Grid'
    case default; name = ''
    end select
  end function center_name

  function attribute_type_name(attribute_type) result(name)
    integer, intent(in) :: attribute_type
    character(len=:), allocatable :: name

    select case (attribute_type)
    case (XDMF_ATTRIBUTE_SCALAR); name = 'Scalar'
    case (XDMF_ATTRIBUTE_VECTOR); name = 'Vector'
    case (XDMF_ATTRIBUTE_TENSOR); name = 'Tensor'
    case (XDMF_ATTRIBUTE_TENSOR6); name = 'Tensor6'
    case (XDMF_ATTRIBUTE_MATRIX); name = 'Matrix'
    case (XDMF_ATTRIBUTE_GLOBAL_ID); name = 'GlobalID'
    case default; name = ''
    end select
  end function attribute_type_name

  function numeric_type_name(numeric_type) result(name)
    integer, intent(in) :: numeric_type
    character(len=:), allocatable :: name

    select case (numeric_type)
    case (XDMF_NUMERIC_REAL32, XDMF_NUMERIC_REAL64); name = 'Float'
    case (XDMF_NUMERIC_INT32, XDMF_NUMERIC_INT64); name = 'Int'
    case default; name = ''
    end select
  end function numeric_type_name

  integer function numeric_type_precision(numeric_type)
    integer, intent(in) :: numeric_type

    select case (numeric_type)
    case (XDMF_NUMERIC_REAL32, XDMF_NUMERIC_INT32); numeric_type_precision = 4
    case (XDMF_NUMERIC_REAL64, XDMF_NUMERIC_INT64); numeric_type_precision = 8
    case default; numeric_type_precision = 0
    end select
  end function numeric_type_precision

  integer(int64) function numeric_type_size(numeric_type)
    integer, intent(in) :: numeric_type

    numeric_type_size = int(numeric_type_precision(numeric_type), int64)
  end function numeric_type_size

  integer(int64) function product_int64(values)
    integer(int64), intent(in) :: values(:)
    integer :: i

    product_int64 = 1_int64
    do i = 1, size(values)
      if (values(i) < 0_int64) then
        product_int64 = -1_int64
        return
      end if
      if (values(i) > 0_int64) then
        if (product_int64 > huge(product_int64) / values(i)) then
          product_int64 = -1_int64
          return
        end if
      end if
      product_int64 = product_int64 * values(i)
    end do
  end function product_int64

end module xdmf_model_m
