module xdmf_xml_m
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use xdmf_model_m, only: xdmf_status_t, collection_record_t, &
    grid_record_t, attribute_record_t, XDMF_ERROR_IO, &
    XDMF_SERIES_NONE, XDMF_SERIES_TIME, XDMF_SERIES_FREQUENCY, &
    XDMF_GEOMETRY_UNIFORM, XDMF_GEOMETRY_RECTILINEAR, &
    XDMF_GEOMETRY_CURVILINEAR, XDMF_GEOMETRY_UNSTRUCTURED, &
    XDMF_TOPOLOGY_MIXED, XDMF_TOPOLOGY_2D_SMESH, &
    XDMF_TOPOLOGY_2D_RECTMESH, XDMF_TOPOLOGY_2D_CORECTMESH, &
    XDMF_TOPOLOGY_3D_SMESH, XDMF_TOPOLOGY_3D_RECTMESH, &
    XDMF_TOPOLOGY_3D_CORECTMESH, topology_name, center_name, &
    attribute_type_name, numeric_type_name, numeric_type_precision, &
    set_status_success, set_status_error

  implicit none

  private

  public :: write_xdmf_document

contains

  subroutine write_xdmf_document(path, hdf5_name, collections, grids, &
      attributes, series_kind, series_values, status)
    character(len=*), intent(in) :: path, hdf5_name
    type(collection_record_t), intent(in) :: collections(:)
    type(grid_record_t), intent(in) :: grids(:)
    type(attribute_record_t), intent(in) :: attributes(:)
    integer, intent(in) :: series_kind
    real(real64), intent(in) :: series_values(:)
    type(xdmf_status_t), intent(out) :: status

    integer :: unit, error, step

    call set_status_success(status)
    open(newunit=unit, file=trim(path), status='replace', action='write', &
      iostat=error)
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_IO, &
        'Could not create XDMF file: '//trim(path))
      return
    end if

    call put_line(unit, 0, '<?xml version="1.0" encoding="UTF-8"?>', status)
    call put_line(unit, 0, '<Xdmf Version="3.0">', status)
    call put_line(unit, 1, '<Domain Name="OpenSEMBA XDMF-HDF5">', status)
    call put_line(unit, 2, &
      '<Information Name="Schema" Value="OpenSEMBA XDMF-HDF5 1.0"/>', &
      status)

    if (series_kind == XDMF_SERIES_NONE) then
      call write_spatial_contents(unit, 2, hdf5_name, collections, grids, &
        attributes, 0, 0, status)
    else
      if (series_kind == XDMF_SERIES_TIME) then
        call put_line(unit, 2, &
          '<Grid Name="Time Series" GridType="Collection" '// &
          'CollectionType="Temporal">', status)
      else
        call put_line(unit, 2, &
          '<Grid Name="Parameter Series" GridType="Collection" '// &
          'CollectionType="Spatial">', status)
      end if

      do step = 1, size(series_values)
        call write_series_step_header(unit, 3, series_kind, step, &
          series_values(step), status)
        call write_spatial_contents(unit, 4, hdf5_name, collections, grids, &
          attributes, step, size(series_values), status)
        call put_line(unit, 3, '</Grid>', status)
      end do
      call put_line(unit, 2, '</Grid>', status)
    end if

    call put_line(unit, 1, '</Domain>', status)
    call put_line(unit, 0, '</Xdmf>', status)
    close(unit, iostat=error)
    if (error /= 0 .and. .not. status%is_error()) then
      call set_status_error(status, XDMF_ERROR_IO, &
        'Could not close XDMF file: '//trim(path))
    end if
  end subroutine write_xdmf_document

  subroutine write_series_step_header(unit, indent, series_kind, step, &
      value, status)
    integer, intent(in) :: unit, indent, series_kind, step
    real(real64), intent(in) :: value
    type(xdmf_status_t), intent(inout) :: status

    character(len=:), allocatable :: step_name, value_text

    step_name = 'Step '//integer_string(int(step, int64))
    value_text = real_string(value)
    call put_line(unit, indent, '<Grid Name="'//xml_escape(step_name)// &
      '" GridType="Collection" CollectionType="Spatial">', status)
    if (series_kind == XDMF_SERIES_TIME) then
      call put_line(unit, indent + 1, '<Time Value="'//value_text//'"/>', status)
    else if (series_kind == XDMF_SERIES_FREQUENCY) then
      call put_line(unit, indent + 1, &
        '<Information Name="Frequency" Value="'//value_text//'"/>', status)
    else
      call put_line(unit, indent + 1, &
        '<Information Name="Parameter" Value="'//value_text//'"/>', status)
    end if
  end subroutine write_series_step_header

  subroutine write_spatial_contents(unit, indent, hdf5_name, collections, &
      grids, attributes, step, number_of_steps, status)
    integer, intent(in) :: unit, indent, step, number_of_steps
    character(len=*), intent(in) :: hdf5_name
    type(collection_record_t), intent(in) :: collections(:)
    type(grid_record_t), intent(in) :: grids(:)
    type(attribute_record_t), intent(in) :: attributes(:)
    type(xdmf_status_t), intent(inout) :: status

    integer :: collection_index, grid_index

    do collection_index = 1, size(collections)
      call put_line(unit, indent, '<Grid Name="'// &
        xml_escape(collections(collection_index)%name)// &
        '" GridType="Collection" CollectionType="Spatial">', status)
      do grid_index = 1, size(grids)
        if (grids(grid_index)%collection_id == &
            collections(collection_index)%id) then
          call write_grid(unit, indent + 1, hdf5_name, grids(grid_index), &
            attributes, step, number_of_steps, status)
        end if
      end do
      call put_line(unit, indent, '</Grid>', status)
    end do

    do grid_index = 1, size(grids)
      if (grids(grid_index)%collection_id == 0) then
        call write_grid(unit, indent, hdf5_name, grids(grid_index), &
          attributes, step, number_of_steps, status)
      end if
    end do
  end subroutine write_spatial_contents

  subroutine write_grid(unit, indent, hdf5_name, grid, attributes, step, &
      number_of_steps, status)
    integer, intent(in) :: unit, indent, step, number_of_steps
    character(len=*), intent(in) :: hdf5_name
    type(grid_record_t), intent(in) :: grid
    type(attribute_record_t), intent(in) :: attributes(:)
    type(xdmf_status_t), intent(inout) :: status

    integer :: attribute_index

    call put_line(unit, indent, '<Grid Name="'//xml_escape(grid%name)// &
      '" GridType="Uniform">', status)
    call write_topology(unit, indent + 1, hdf5_name, grid, status)
    call write_geometry(unit, indent + 1, hdf5_name, grid, status)
    do attribute_index = 1, size(attributes)
      if (attributes(attribute_index)%grid_id == grid%id) then
        if (.not. attributes(attribute_index)%is_series .or. step > 0) then
          call write_attribute(unit, indent + 1, hdf5_name, &
            attributes(attribute_index), step, number_of_steps, status)
        end if
      end if
    end do
    call put_line(unit, indent, '</Grid>', status)
  end subroutine write_grid

  subroutine write_topology(unit, indent, hdf5_name, grid, status)
    integer, intent(in) :: unit, indent
    character(len=*), intent(in) :: hdf5_name
    type(grid_record_t), intent(in) :: grid
    type(xdmf_status_t), intent(inout) :: status

    character(len=:), allocatable :: line, dimensions
    integer(int64), allocatable :: connectivity_shape(:)

    if (is_structured_topology(grid%topology_type)) then
      dimensions = dimensions_string(grid%dimensions, .true.)
      line = '<Topology TopologyType="'//topology_name(grid%topology_type)// &
        '" Dimensions="'//dimensions//'"/>'
      call put_line(unit, indent, line, status)
      return
    end if

    line = '<Topology TopologyType="'//topology_name(grid%topology_type)// &
      '" NumberOfElements="'//integer_string(grid%number_of_elements)//'"'
    if (grid%nodes_per_element > 0 .and. &
        grid%topology_type /= XDMF_TOPOLOGY_MIXED) then
      line = line//' NodesPerElement="'// &
        integer_string(int(grid%nodes_per_element, int64))//'"'
    end if
    call put_line(unit, indent, line//'>', status)

    if (grid%topology_type == XDMF_TOPOLOGY_MIXED) then
      connectivity_shape = [grid%mixed_connectivity_size]
    else
      connectivity_shape = [int(grid%nodes_per_element, int64), &
        grid%number_of_elements]
    end if
    call write_hdf_data_item(unit, indent + 1, hdf5_name, &
      grid%connectivity_path, connectivity_shape, &
      grid%topology_numeric_type, status)
    call put_line(unit, indent, '</Topology>', status)
  end subroutine write_topology

  subroutine write_geometry(unit, indent, hdf5_name, grid, status)
    integer, intent(in) :: unit, indent
    character(len=*), intent(in) :: hdf5_name
    type(grid_record_t), intent(in) :: grid
    type(xdmf_status_t), intent(inout) :: status

    integer :: axis
    integer(int64), allocatable :: shape(:)
    character(len=:), allocatable :: geometry_name

    select case (grid%geometry_type)
    case (XDMF_GEOMETRY_UNIFORM)
      if (grid%dimension == 2) then
        geometry_name = 'ORIGIN_DXDY'
      else
        geometry_name = 'ORIGIN_DXDYDZ'
      end if
      call put_line(unit, indent, '<Geometry GeometryType="'// &
        geometry_name//'">', status)
      shape = [int(grid%dimension, int64)]
      call write_hdf_data_item(unit, indent + 1, hdf5_name, &
        grid%origin_path, shape, grid%geometry_numeric_type, status)
      call write_hdf_data_item(unit, indent + 1, hdf5_name, &
        grid%spacing_path, shape, grid%geometry_numeric_type, status)
      call put_line(unit, indent, '</Geometry>', status)

    case (XDMF_GEOMETRY_RECTILINEAR)
      if (grid%dimension == 2) then
        geometry_name = 'VXVY'
      else
        geometry_name = 'VXVYVZ'
      end if
      call put_line(unit, indent, '<Geometry GeometryType="'// &
        geometry_name//'">', status)
      do axis = 1, grid%dimension
        shape = [grid%axis_sizes(axis)]
        call write_hdf_data_item(unit, indent + 1, hdf5_name, &
          grid%axis_paths(axis), shape, grid%geometry_numeric_type, status)
      end do
      call put_line(unit, indent, '</Geometry>', status)

    case (XDMF_GEOMETRY_CURVILINEAR)
      if (grid%dimension == 2) then
        geometry_name = 'XY'
      else
        geometry_name = 'XYZ'
      end if
      call put_line(unit, indent, '<Geometry GeometryType="'// &
        geometry_name//'">', status)
      allocate(shape(size(grid%dimensions) + 1))
      shape(1) = int(grid%dimension, int64)
      shape(2:) = grid%dimensions
      call write_hdf_data_item(unit, indent + 1, hdf5_name, &
        grid%points_path, shape, grid%geometry_numeric_type, status)
      call put_line(unit, indent, '</Geometry>', status)

    case (XDMF_GEOMETRY_UNSTRUCTURED)
      if (grid%dimension == 2) then
        geometry_name = 'XY'
      else
        geometry_name = 'XYZ'
      end if
      call put_line(unit, indent, '<Geometry GeometryType="'// &
        geometry_name//'">', status)
      shape = [int(grid%dimension, int64), grid%number_of_points]
      call write_hdf_data_item(unit, indent + 1, hdf5_name, &
        grid%points_path, shape, grid%geometry_numeric_type, status)
      call put_line(unit, indent, '</Geometry>', status)
    end select
  end subroutine write_geometry

  subroutine write_attribute(unit, indent, hdf5_name, attribute, step, &
      number_of_steps, status)
    integer, intent(in) :: unit, indent, step, number_of_steps
    character(len=*), intent(in) :: hdf5_name
    type(attribute_record_t), intent(in) :: attribute
    type(xdmf_status_t), intent(inout) :: status

    character(len=:), allocatable :: line

    line = '<Attribute Name="'//xml_escape(attribute%name)// &
      '" Center="'//center_name(attribute%center)// &
      '" AttributeType="'//attribute_type_name(attribute%attribute_type)//'">'
    call put_line(unit, indent, line, status)
    if (attribute%is_series) then
      call write_hyperslab_data_item(unit, indent + 1, hdf5_name, &
        attribute, step, number_of_steps, status)
    else
      call write_hdf_data_item(unit, indent + 1, hdf5_name, &
        attribute%dataset_path, attribute%storage_shape, &
        attribute%numeric_type, status)
    end if
    call put_line(unit, indent, '</Attribute>', status)
  end subroutine write_attribute

  subroutine write_hyperslab_data_item(unit, indent, hdf5_name, attribute, &
      step, number_of_steps, status)
    integer, intent(in) :: unit, indent, step, number_of_steps
    character(len=*), intent(in) :: hdf5_name
    type(attribute_record_t), intent(in) :: attribute
    type(xdmf_status_t), intent(inout) :: status

    integer :: rank
    integer(int64), allocatable :: external_shape(:), source_shape(:)
    integer(int64), allocatable :: start(:), stride(:), count(:)

    external_shape = reverse_int64(attribute%storage_shape)
    rank = size(external_shape) + 1
    allocate(source_shape(rank), start(rank), stride(rank), count(rank))
    source_shape(1) = int(number_of_steps, int64)
    source_shape(2:) = external_shape
    start = 0_int64
    start(1) = int(step - 1, int64)
    stride = 1_int64
    count(1) = 1_int64
    count(2:) = external_shape

    call put_line(unit, indent, '<DataItem ItemType="HyperSlab" '// &
      'Dimensions="'//dimensions_string(external_shape, .false.)//'">', status)
    call put_line(unit, indent + 1, '<DataItem Dimensions="3 '// &
      integer_string(int(rank, int64))//'" Format="XML">', status)
    call put_line(unit, indent + 2, &
      dimensions_string(start, .false.), status)
    call put_line(unit, indent + 2, &
      dimensions_string(stride, .false.), status)
    call put_line(unit, indent + 2, &
      dimensions_string(count, .false.), status)
    call put_line(unit, indent + 1, '</DataItem>', status)
    call write_hdf_data_item_external_shape(unit, indent + 1, hdf5_name, &
      attribute%dataset_path, source_shape, attribute%numeric_type, status)
    call put_line(unit, indent, '</DataItem>', status)
  end subroutine write_hyperslab_data_item

  subroutine write_hdf_data_item(unit, indent, hdf5_name, dataset_path, &
      fortran_shape, numeric_type, status)
    integer, intent(in) :: unit, indent, numeric_type
    character(len=*), intent(in) :: hdf5_name, dataset_path
    integer(int64), intent(in) :: fortran_shape(:)
    type(xdmf_status_t), intent(inout) :: status

    call write_hdf_data_item_external_shape(unit, indent, hdf5_name, &
      dataset_path, reverse_int64(fortran_shape), numeric_type, status)
  end subroutine write_hdf_data_item

  subroutine write_hdf_data_item_external_shape(unit, indent, hdf5_name, &
      dataset_path, external_shape, numeric_type, status)
    integer, intent(in) :: unit, indent, numeric_type
    character(len=*), intent(in) :: hdf5_name, dataset_path
    integer(int64), intent(in) :: external_shape(:)
    type(xdmf_status_t), intent(inout) :: status

    character(len=:), allocatable :: line

    line = '<DataItem Dimensions="'// &
      dimensions_string(external_shape, .false.)//'" NumberType="'// &
      numeric_type_name(numeric_type)//'" Precision="'// &
      integer_string(int(numeric_type_precision(numeric_type), int64))// &
      '" Format="HDF">'
    call put_line(unit, indent, line, status)
    call put_line(unit, indent + 1, xml_escape(trim(hdf5_name))//':'// &
      trim(dataset_path), status)
    call put_line(unit, indent, '</DataItem>', status)
  end subroutine write_hdf_data_item_external_shape

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

  function reverse_int64(values) result(reversed)
    integer(int64), intent(in) :: values(:)
    integer(int64), allocatable :: reversed(:)
    integer :: i

    allocate(reversed(size(values)))
    do i = 1, size(values)
      reversed(i) = values(size(values) - i + 1)
    end do
  end function reverse_int64

  function dimensions_string(values, reverse) result(text)
    integer(int64), intent(in) :: values(:)
    logical, intent(in) :: reverse
    character(len=:), allocatable :: text
    integer :: i, index

    text = ''
    do i = 1, size(values)
      if (reverse) then
        index = size(values) - i + 1
      else
        index = i
      end if
      if (i > 1) text = text//' '
      text = text//integer_string(values(index))
    end do
  end function dimensions_string

  function integer_string(value) result(text)
    integer(int64), intent(in) :: value
    character(len=:), allocatable :: text
    character(len=32) :: buffer

    write(buffer, '(I0)') value
    text = trim(buffer)
  end function integer_string

  function real_string(value) result(text)
    real(real64), intent(in) :: value
    character(len=:), allocatable :: text
    character(len=64) :: buffer

    write(buffer, '(G0.17)') value
    text = trim(adjustl(buffer))
  end function real_string

  function xml_escape(value) result(escaped)
    character(len=*), intent(in) :: value
    character(len=:), allocatable :: escaped
    integer :: i

    escaped = ''
    do i = 1, len_trim(value)
      select case (value(i:i))
      case ('&'); escaped = escaped//'&amp;'
      case ('<'); escaped = escaped//'&lt;'
      case ('>'); escaped = escaped//'&gt;'
      case ('"'); escaped = escaped//'&quot;'
      case ("'"); escaped = escaped//'&apos;'
      case default; escaped = escaped//value(i:i)
      end select
    end do
  end function xml_escape

  subroutine put_line(unit, indent, line, status)
    integer, intent(in) :: unit, indent
    character(len=*), intent(in) :: line
    type(xdmf_status_t), intent(inout) :: status

    integer :: error

    if (status%is_error()) return
    write(unit, '(A)', iostat=error) repeat('  ', indent)//trim(line)
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_IO, &
        'Could not write XDMF metadata')
    end if
  end subroutine put_line

end module xdmf_xml_m
