program generate_cases
  use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use xdmf_hdf5_m

  implicit none

  character(len=1024) :: output_dir
  character(len=32) :: option
  character(kind=c_char), allocatable :: c_output_dir(:)
  integer :: argument_index, failures, output_dir_length
  integer(c_int) :: directory_error
  logical :: examples_only, replace_existing

  interface
    function xdmf_create_directory(directory) bind(C) result(error)
      import :: c_char, c_int
      character(kind=c_char), intent(in) :: directory(*)
      integer(c_int) :: error
    end function xdmf_create_directory
  end interface

  failures = 0
  examples_only = .false.
  replace_existing = .false.
  if (command_argument_count() == 0) then
    call print_help()
    stop
  end if
  call get_command_argument(1, output_dir)
  if (trim(output_dir) == '--help') then
    call print_help()
    stop
  end if
  if (len_trim(output_dir) == 0) then
    call print_help()
    stop 1
  end if
  output_dir_length = len_trim(output_dir)
  allocate(c_output_dir(output_dir_length + 1))
  do argument_index = 1, output_dir_length
    c_output_dir(argument_index) = output_dir(argument_index:argument_index)
  end do
  c_output_dir(output_dir_length + 1) = c_null_char
  directory_error = xdmf_create_directory(c_output_dir)
  if (directory_error /= 0) then
    write(*, '(A,I0,A)') 'Error: could not create output directory (system error ', &
      directory_error, '): '//trim(output_dir)
    stop 1
  end if
  do argument_index = 2, command_argument_count()
    call get_command_argument(argument_index, option)
    select case (trim(option))
    case ('--examples')
      examples_only = .true.
    case ('--replace')
      replace_existing = .true.
    case default
      write(*, '(A)') 'Error: unknown option: '//trim(option)
      call print_help()
      stop 1
    end select
  end do
  if (.not. replace_existing) then
    call require_available_output_paths(trim(output_dir), examples_only)
  end if

  call generate_uniform(trim(output_dir), failures)
  call generate_volume(trim(output_dir), failures)
  call generate_curvilinear(trim(output_dir), failures)
  call generate_unstructured(trim(output_dir), failures)
  call generate_mixed(trim(output_dir), failures)
  call generate_time_series(trim(output_dir), failures)
  if (.not. examples_only) then
    call generate_rectilinear(trim(output_dir), failures)
    call generate_frequency_series(trim(output_dir), failures)
    call generate_escaped_path(trim(output_dir), failures)
    call exercise_validation(trim(output_dir), failures)
  end if

  if (failures /= 0) then
    write(*, '(A,I0)') 'XDMF/HDF5 generator failures: ', failures
    stop 1
  end if

contains

  subroutine generate_uniform(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_grid_id_t) :: grid
    type(xdmf_attribute_id_t) :: field
    real(real32) :: values(6)
    integer :: i

    options%overwrite = replace_existing
    call writer%create(pair_path(directory, 'uniform'), options, status)
    call require_ok(status, 'create uniform', failures)
    call writer%define_uniform_grid('uniform-2d', [3_int64, 2_int64], &
      [0.0_real32, -1.0_real32], [0.5_real32, 2.0_real32], &
      grid, status)
    call require_ok(status, 'define uniform grid', failures)
    call writer%define_attribute(grid, 'temperature', XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL32, field, status)
    call require_ok(status, 'define uniform attribute', failures)
    do i = 1, size(values)
      values(i) = real(i, real32) * 0.25_real32
    end do
    call writer%write_attribute(field, values, status)
    call require_ok(status, 'write uniform attribute', failures)
    call writer%close(status)
    call require_ok(status, 'close uniform', failures)
  end subroutine generate_uniform

  subroutine generate_volume(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    integer, parameter :: nx = 21, ny = 21, nz = 21
    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_grid_id_t) :: grid
    type(xdmf_attribute_id_t) :: field
    real(real64), allocatable :: values(:)
    real(real64) :: x, y, z
    integer :: i, j, k, value_index

    options%overwrite = replace_existing
    call writer%create(pair_path(directory, 'volume'), options, status)
    call require_ok(status, 'create volume', failures)
    call writer%define_uniform_grid('electric-field-volume', &
      [int(nx, int64), int(ny, int64), int(nz, int64)], &
      [-1.0_real64, -1.0_real64, -1.0_real64], &
      [0.1_real64, 0.1_real64, 0.1_real64], grid, status)
    call require_ok(status, 'define volume grid', failures)
    call writer%define_attribute(grid, 'electric-field-magnitude', &
      XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, &
      field, status)
    call require_ok(status, 'define volume field', failures)

    allocate(values(nx * ny * nz))
    do k = 1, nz
      z = -1.0_real64 + 0.1_real64 * real(k - 1, real64)
      do j = 1, ny
        y = -1.0_real64 + 0.1_real64 * real(j - 1, real64)
        do i = 1, nx
          x = -1.0_real64 + 0.1_real64 * real(i - 1, real64)
          value_index = i + nx * ((j - 1) + ny * (k - 1))
          values(value_index) = exp(-8.0_real64 * (x*x + y*y + z*z)) &
            + 0.35_real64 * exp(-20.0_real64 * &
              ((x - 0.45_real64)**2 + (y + 0.25_real64)**2 &
                + (z - 0.2_real64)**2))
        end do
      end do
    end do

    call writer%write_attribute(field, values, status)
    call require_ok(status, 'write volume field', failures)
    call writer%close(status)
    call require_ok(status, 'close volume', failures)
  end subroutine generate_volume

  subroutine generate_rectilinear(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_collection_id_t) :: collection
    type(xdmf_grid_id_t) :: grid
    type(xdmf_attribute_id_t) :: velocity, material, identifier
    type(xdmf_attribute_id_t) :: tensor, tensor6, matrix
    real(real64) :: x(2), y(3), z(4), vector_values(72)
    real(real64) :: tensor_values(216), tensor6_values(36)
    real(real32) :: matrix_values(6)
    integer(int32) :: material_values(6)
    integer(int64) :: identifier_value(1)
    integer :: i

    options%overwrite = replace_existing
    call writer%create(pair_path(directory, 'rectilinear'), options, status)
    call require_ok(status, 'create rectilinear', failures)
    call writer%define_collection('collection & <escaped>', collection, status)
    call require_ok(status, 'define escaped collection', failures)
    x = [-2.0_real64, 0.5_real64]
    y = [0.0_real64, 1.0_real64, 4.0_real64]
    z = [10.0_real64, 11.0_real64, 13.0_real64, 17.0_real64]
    call writer%define_rectilinear_grid('grid & "quoted"', x, y, z, &
      grid, status, collection)
    call require_ok(status, 'define rectilinear grid', failures)
    call writer%define_attribute(grid, 'velocity <xyz>', XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_VECTOR, XDMF_NUMERIC_REAL64, velocity, status)
    call require_ok(status, 'define vector attribute', failures)
    call writer%define_attribute(grid, 'material', XDMF_CENTER_CELL, &
      XDMF_ATTRIBUTE_GLOBAL_ID, XDMF_NUMERIC_INT32, material, status)
    call require_ok(status, 'define cell IDs', failures)
    call writer%define_attribute(grid, 'grid-id', XDMF_CENTER_GRID, &
      XDMF_ATTRIBUTE_GLOBAL_ID, XDMF_NUMERIC_INT64, identifier, status)
    call require_ok(status, 'define grid ID', failures)
    call writer%define_attribute(grid, 'stress', XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_TENSOR, XDMF_NUMERIC_REAL64, tensor, status)
    call require_ok(status, 'define tensor', failures)
    call writer%define_attribute(grid, 'symmetric-stress', XDMF_CENTER_CELL, &
      XDMF_ATTRIBUTE_TENSOR6, XDMF_NUMERIC_REAL64, tensor6, status)
    call require_ok(status, 'define tensor6', failures)
    call writer%define_attribute(grid, 'transform', XDMF_CENTER_GRID, &
      XDMF_ATTRIBUTE_MATRIX, XDMF_NUMERIC_REAL32, matrix, status, &
      component_shape=[2_int64, 3_int64])
    call require_ok(status, 'define matrix', failures)
    do i = 1, size(vector_values)
      vector_values(i) = real(i, real64)
    end do
    material_values = [11_int32, 12_int32, 13_int32, 14_int32, &
      15_int32, 16_int32]
    identifier_value = [9876543210_int64]
    do i = 1, size(tensor_values)
      tensor_values(i) = real(2000 + i, real64)
    end do
    do i = 1, size(tensor6_values)
      tensor6_values(i) = real(3000 + i, real64)
    end do
    do i = 1, size(matrix_values)
      matrix_values(i) = real(4000 + i, real32)
    end do
    call writer%write_attribute(velocity, vector_values, status)
    call require_ok(status, 'write vector attribute', failures)
    call writer%write_attribute(material, material_values, status)
    call require_ok(status, 'write cell IDs', failures)
    call writer%write_attribute(identifier, identifier_value, status)
    call require_ok(status, 'write grid ID', failures)
    call writer%write_attribute(tensor, tensor_values, status)
    call require_ok(status, 'write tensor', failures)
    call writer%write_attribute(tensor6, tensor6_values, status)
    call require_ok(status, 'write tensor6', failures)
    call writer%write_attribute(matrix, matrix_values, status)
    call require_ok(status, 'write matrix', failures)
    call writer%close(status)
    call require_ok(status, 'close rectilinear', failures)
  end subroutine generate_rectilinear

  subroutine generate_curvilinear(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_grid_id_t) :: grid
    real(real32) :: points(2, 6)
    integer :: i

    do i = 1, 6
      points(1, i) = real(mod(i - 1, 2), real32) + &
        0.1_real32 * real((i - 1) / 2, real32)
      points(2, i) = real((i - 1) / 2, real32)
    end do
    options%overwrite = replace_existing
    call writer%create(pair_path(directory, 'curvilinear'), options, status)
    call require_ok(status, 'create curvilinear', failures)
    call writer%define_curvilinear_grid('curved-surface', &
      [2_int64, 3_int64], points, grid, status)
    call require_ok(status, 'define curvilinear grid', failures)
    call writer%close(status)
    call require_ok(status, 'close curvilinear', failures)
  end subroutine generate_curvilinear

  subroutine generate_unstructured(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_collection_id_t) :: collection
    real(real64) :: points(3, 27)
    integer :: i

    do i = 1, 27
      points(:, i) = [real(i - 1, real64), &
        real(mod(i - 1, 5), real64), real(mod(i - 1, 3), real64)]
    end do
    options%overwrite = replace_existing
    call writer%create(pair_path(directory, 'unstructured'), options, status)
    call require_ok(status, 'create unstructured', failures)
    call writer%define_collection('topology-suite', collection, status)
    call require_ok(status, 'define topology collection', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_POLYVERTEX, &
      1, 'polyvertex', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_POLYLINE, &
      3, 'polyline', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_POLYGON, &
      5, 'polygon', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_TRIANGLE, &
      3, 'triangle', failures)
    call add_topology(writer, collection, points, &
      XDMF_TOPOLOGY_QUADRILATERAL, 4, 'quadrilateral', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_TETRAHEDRON, &
      4, 'tetrahedron', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_PYRAMID, &
      5, 'pyramid', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_WEDGE, &
      6, 'wedge', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_HEXAHEDRON, &
      8, 'hexahedron', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_EDGE_3, &
      3, 'edge-3', failures)
    call add_topology(writer, collection, points, &
      XDMF_TOPOLOGY_QUADRILATERAL_9, 9, 'quadrilateral-9', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_TRIANGLE_6, &
      6, 'triangle-6', failures)
    call add_topology(writer, collection, points, &
      XDMF_TOPOLOGY_QUADRILATERAL_8, 8, 'quadrilateral-8', failures)
    call add_topology(writer, collection, points, &
      XDMF_TOPOLOGY_TETRAHEDRON_10, 10, 'tetrahedron-10', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_PYRAMID_13, &
      13, 'pyramid-13', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_WEDGE_15, &
      15, 'wedge-15', failures)
    call add_topology(writer, collection, points, XDMF_TOPOLOGY_WEDGE_18, &
      18, 'wedge-18', failures)
    call add_topology(writer, collection, points, &
      XDMF_TOPOLOGY_HEXAHEDRON_20, 20, 'hexahedron-20', failures)
    call add_topology(writer, collection, points, &
      XDMF_TOPOLOGY_HEXAHEDRON_24, 24, 'hexahedron-24', failures)
    call add_topology(writer, collection, points, &
      XDMF_TOPOLOGY_HEXAHEDRON_27, 27, 'hexahedron-27', failures)
    call writer%close(status)
    call require_ok(status, 'close unstructured', failures)
  end subroutine generate_unstructured

  subroutine add_topology(writer, collection, points, topology, node_count, &
      name, failures)
    type(xdmf_writer_t), intent(inout) :: writer
    type(xdmf_collection_id_t), intent(in) :: collection
    real(real64), intent(in) :: points(:, :)
    integer, intent(in) :: topology, node_count
    character(len=*), intent(in) :: name
    integer, intent(inout) :: failures

    type(xdmf_grid_id_t) :: grid
    type(xdmf_status_t) :: status
    integer(int64), allocatable :: connectivity(:, :)
    integer :: i

    allocate(connectivity(node_count, 1))
    do i = 1, node_count
      connectivity(i, 1) = int(i, int64)
    end do
    call writer%define_unstructured_grid(name, topology, points, &
      connectivity, grid, status, collection)
    call require_ok(status, 'define topology '//trim(name), failures)
  end subroutine add_topology

  subroutine generate_mixed(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_grid_id_t) :: grid
    real(real64) :: points(3, 8)
    integer(int64) :: connectivity(55)
    integer :: i

    do i = 1, 8
      points(:, i) = [real(mod(i - 1, 2), real64), &
        real(mod((i - 1) / 2, 2), real64), &
        real((i - 1) / 4, real64)]
    end do
    connectivity = [ &
      int(XDMF_TOPOLOGY_TRIANGLE, int64), 1_int64, 2_int64, 3_int64, &
      int(XDMF_TOPOLOGY_QUADRILATERAL, int64), &
        1_int64, 2_int64, 4_int64, 3_int64, &
      int(XDMF_TOPOLOGY_POLYGON, int64), 5_int64, &
        1_int64, 2_int64, 6_int64, 7_int64, 3_int64, &
      int(XDMF_TOPOLOGY_TETRAHEDRON, int64), &
        1_int64, 2_int64, 3_int64, 5_int64, &
      int(XDMF_TOPOLOGY_POLYVERTEX, int64), 2_int64, &
        7_int64, 8_int64, &
      int(XDMF_TOPOLOGY_POLYLINE, int64), 3_int64, &
        1_int64, 2_int64, 4_int64, &
      int(XDMF_TOPOLOGY_POLYHEDRON, int64), 4_int64, &
        3_int64, 1_int64, 2_int64, 3_int64, &
        3_int64, 1_int64, 2_int64, 5_int64, &
        3_int64, 2_int64, 3_int64, 5_int64, &
        3_int64, 3_int64, 1_int64, 5_int64, &
      int(XDMF_TOPOLOGY_TRIANGLE_6, int64), &
        1_int64, 2_int64, 3_int64, 4_int64, 5_int64, 6_int64]
    options%overwrite = replace_existing
    call writer%create(pair_path(directory, 'mixed'), options, status)
    call require_ok(status, 'create mixed', failures)
    call writer%define_mixed_grid('mixed-elements', points, connectivity, &
      grid, status)
    call require_ok(status, 'define mixed grid', failures)
    call writer%close(status)
    call require_ok(status, 'close mixed', failures)
  end subroutine generate_mixed

  subroutine generate_time_series(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    integer, parameter :: nx = 25, ny = 25, nz = 25
    integer, parameter :: npoints = nx * ny * nz, nsteps = 20
    real(real32), parameter :: origin = -1.2_real32
    real(real32), parameter :: spacing = 0.1_real32
    real(real32), parameter :: two_pi = 6.2831853_real32
    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_grid_id_t) :: grid
    type(xdmf_attribute_id_t) :: magnitude, electric_field
    real(real32), allocatable :: magnitude_values(:), electric_values(:)
    real(real32) :: x, y, z, progress, pulse_center, envelope, ey
    real(real64) :: time
    integer :: i, j, k, point_index, step

    options%overwrite = replace_existing
    options%series_kind = XDMF_SERIES_TIME
    options%chunk_target_bytes = 1048576_int64
    call writer%create(pair_path(directory, 'time-series'), options, status)
    call require_ok(status, 'create time series', failures)
    call writer%define_uniform_grid('travelling-wave-pulse', &
      [int(nx, int64), int(ny, int64), int(nz, int64)], &
      [origin, origin, origin], [spacing, spacing, spacing], grid, status)
    call require_ok(status, 'define time grid', failures)
    call writer%define_attribute(grid, 'electric-field-magnitude', &
      XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL32, &
      .true., magnitude, status)
    call require_ok(status, 'define time magnitude', failures)
    call writer%define_attribute(grid, 'electric-field', XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_VECTOR, XDMF_NUMERIC_REAL32, .true., &
      electric_field, status)
    call require_ok(status, 'define time electric field', failures)

    allocate(magnitude_values(npoints), electric_values(3 * npoints))
    do step = 1, nsteps
      progress = real(step - 1, real32) / real(nsteps - 1, real32)
      time = real(progress, real64)
      pulse_center = -0.8_real32 + 1.6_real32 * progress
      electric_values = 0.0_real32
      do k = 1, nz
        z = origin + spacing * real(k - 1, real32)
        do j = 1, ny
          y = origin + spacing * real(j - 1, real32)
          do i = 1, nx
            x = origin + spacing * real(i - 1, real32)
            point_index = i + nx * ((j - 1) + ny * (k - 1))
            envelope = exp(-((x - pulse_center) / 0.32_real32)**2 &
              - (y / 0.55_real32)**2 - (z / 0.55_real32)**2)
            ey = envelope * cos(two_pi * (x - pulse_center) / 0.55_real32)
            magnitude_values(point_index) = abs(ey)
            electric_values(3 * point_index - 1) = ey
          end do
        end do
      end do

      call writer%begin_step(time, status)
      call require_ok(status, 'begin time step', failures)
      call writer%write_attribute(magnitude, magnitude_values, status)
      call require_ok(status, 'write time magnitude', failures)
      call writer%write_attribute(electric_field, electric_values, status)
      call require_ok(status, 'write time electric field', failures)
      call writer%end_step(status)
      call require_ok(status, 'end time step', failures)
      if (step == nsteps / 2) then
        call writer%flush(status)
        call require_ok(status, 'flush time series', failures)
      end if
    end do
    call writer%close(status)
    call require_ok(status, 'close time series', failures)
  end subroutine generate_time_series

  subroutine generate_frequency_series(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_grid_id_t) :: grid
    type(xdmf_attribute_id_t) :: real_part, imaginary_part
    real(real64) :: real_values(4), imaginary_values(4)
    integer :: step

    options%overwrite = replace_existing
    options%series_kind = XDMF_SERIES_FREQUENCY
    call writer%create(pair_path(directory, 'frequency-series'), &
      options, status)
    call require_ok(status, 'create frequency series', failures)
    call writer%define_uniform_grid('frequency-plane', &
      [2_int64, 2_int64], [0.0_real64, 0.0_real64], &
      [0.5_real64, 0.5_real64], grid, status)
    call require_ok(status, 'define frequency grid', failures)
    call writer%define_attribute(grid, 'field.real', XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., &
      real_part, status)
    call require_ok(status, 'define real part', failures)
    call writer%define_attribute(grid, 'field.imaginary', XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., &
      imaginary_part, status)
    call require_ok(status, 'define imaginary part', failures)
    do step = 1, 2
      call writer%begin_step(real(step, real64) * 1.0e6_real64, status)
      call require_ok(status, 'begin frequency step', failures)
      real_values = real(step, real64) * [1.0_real64, 2.0_real64, &
        3.0_real64, 4.0_real64]
      imaginary_values = -real_values
      call writer%write_attribute(real_part, real_values, status)
      call require_ok(status, 'write real part', failures)
      call writer%write_attribute(imaginary_part, imaginary_values, status)
      call require_ok(status, 'write imaginary part', failures)
      call writer%end_step(status)
      call require_ok(status, 'end frequency step', failures)
    end do
    call writer%close(status)
    call require_ok(status, 'close frequency series', failures)
  end subroutine generate_frequency_series

  subroutine generate_escaped_path(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_grid_id_t) :: grid

    options%overwrite = replace_existing
    call writer%create(pair_path(directory, 'escaped-&-pair'), &
      options, status)
    call require_ok(status, 'create escaped path', failures)
    call writer%define_uniform_grid('escaped-path-grid', &
      [2_int64, 2_int64], [0.0_real64, 0.0_real64], &
      [1.0_real64, 1.0_real64], grid, status)
    call require_ok(status, 'define escaped path grid', failures)
    call writer%close(status)
    call require_ok(status, 'close escaped path', failures)
  end subroutine generate_escaped_path

  subroutine exercise_validation(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_grid_id_t) :: grid, unused_grid
    type(xdmf_attribute_id_t) :: field
    real(real64) :: points(3, 3), correct_values(4)
    integer(int64) :: bad_connectivity(3, 1)

    options%overwrite = replace_existing
    call writer%create(pair_path(directory, 'validation-static'), &
      options, status)
    call require_ok(status, 'create validation static', failures)
    call writer%define_uniform_grid('validation-grid', &
      [2_int64, 2_int64], [0.0_real64, 0.0_real64], &
      [1.0_real64, 1.0_real64], grid, status)
    call require_ok(status, 'define validation grid', failures)
    call writer%define_uniform_grid('validation-grid', &
      [2_int64, 2_int64], [0.0_real64, 0.0_real64], &
      [1.0_real64, 1.0_real64], unused_grid, status)
    call require_error(status, 'duplicate grid name', failures)
    points = reshape([0.0_real64, 0.0_real64, 0.0_real64, &
      1.0_real64, 0.0_real64, 0.0_real64, &
      0.0_real64, 1.0_real64, 0.0_real64], shape(points))
    bad_connectivity(:, 1) = [1_int64, 2_int64, 99_int64]
    call writer%define_unstructured_grid('bad-connectivity', &
      XDMF_TOPOLOGY_TRIANGLE, points, bad_connectivity, unused_grid, status)
    call require_error(status, 'out-of-range connectivity', failures)
    call writer%define_attribute(grid, 'required-static', XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, field, status)
    call require_ok(status, 'define required static', failures)
    call writer%write_attribute(field, [1.0_real64], status)
    call require_error(status, 'wrong value count', failures)
    call writer%write_attribute(field, &
      [1.0_real32, 2.0_real32, 3.0_real32, 4.0_real32], status)
    call require_error(status, 'wrong value type', failures)
    call writer%close(status)
    call require_error(status, 'retryable close', failures)
    correct_values = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    call writer%write_attribute(field, correct_values, status)
    call require_ok(status, 'write after failed close', failures)
    call writer%close(status)
    call require_ok(status, 'close validation static', failures)

    call exercise_series_validation(directory, failures)
    call exercise_writer_ownership(directory, failures)
    call exercise_pair_protection(directory, failures)
  end subroutine exercise_validation

  subroutine exercise_series_validation(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_collection_id_t) :: collection
    type(xdmf_grid_id_t) :: grid
    type(xdmf_attribute_id_t) :: first, second
    real(real64) :: values(4)

    options%overwrite = replace_existing
    options%series_kind = XDMF_SERIES_TIME
    call writer%create(pair_path(directory, 'validation-series'), &
      options, status)
    call require_ok(status, 'create validation series', failures)
    call writer%define_uniform_grid('validation-series-grid', &
      [2_int64, 2_int64], [0.0_real64, 0.0_real64], &
      [1.0_real64, 1.0_real64], grid, status)
    call require_ok(status, 'define validation series grid', failures)
    call writer%define_attribute(grid, 'first', XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., first, status)
    call require_ok(status, 'define first series field', failures)
    call writer%define_attribute(grid, 'second', XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., second, status)
    call require_ok(status, 'define second series field', failures)
    values = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    call writer%begin_step(0.0_real64, status)
    call require_ok(status, 'begin incomplete step', failures)
    call writer%define_collection('too-late', collection, status)
    call require_error(status, 'definition after begin step', failures)
    call writer%write_attribute(first, values, status)
    call require_ok(status, 'write incomplete step', failures)
    call writer%end_step(status)
    call require_error(status, 'incomplete step rollback', failures)
    call writer%begin_step(1.0_real64, status)
    call require_ok(status, 'begin recovered step', failures)
    call writer%write_attribute(first, values, status)
    call require_ok(status, 'write recovered first', failures)
    call writer%write_attribute(second, -values, status)
    call require_ok(status, 'write recovered second', failures)
    call writer%end_step(status)
    call require_ok(status, 'commit recovered step', failures)
    call writer%close(status)
    call require_ok(status, 'close validation series', failures)
  end subroutine exercise_series_validation

  subroutine exercise_writer_ownership(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: first_writer, second_writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    type(xdmf_grid_id_t) :: first_grid, second_grid
    type(xdmf_attribute_id_t) :: field
    real(real64) :: values(4)

    options%overwrite = replace_existing
    call first_writer%create(pair_path(directory, 'owner-first'), &
      options, status)
    call require_ok(status, 'create first owner', failures)
    call first_writer%define_uniform_grid('first-grid', &
      [2_int64, 2_int64], [0.0_real64, 0.0_real64], &
      [1.0_real64, 1.0_real64], first_grid, status)
    call require_ok(status, 'define first owner grid', failures)

    call second_writer%create(pair_path(directory, 'owner-second'), &
      options, status)
    call require_ok(status, 'create second owner', failures)
    call second_writer%define_uniform_grid('second-grid', &
      [2_int64, 2_int64], [0.0_real64, 0.0_real64], &
      [1.0_real64, 1.0_real64], second_grid, status)
    call require_ok(status, 'define second owner grid', failures)
    call second_writer%define_attribute(first_grid, 'foreign-grid', &
      XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, &
      field, status)
    call require_error(status, 'cross-writer grid ID', failures)

    call first_writer%close(status)
    call require_ok(status, 'close first concurrent writer', failures)
    call second_writer%define_attribute(second_grid, 'still-open', &
      XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, &
      field, status)
    call require_ok(status, 'use second writer after first close', failures)
    values = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
    call second_writer%write_attribute(field, values, status)
    call require_ok(status, 'write second concurrent writer', failures)
    call second_writer%close(status)
    call require_ok(status, 'close second concurrent writer', failures)

    ! This conformance check intentionally reopens the pair created above.
    options%overwrite = .true.
    call first_writer%create(pair_path(directory, 'owner-first'), &
      options, status)
    call require_ok(status, 'reopen first owner', failures)
    call first_writer%define_attribute(first_grid, 'stale-grid', &
      XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, &
      field, status)
    call require_error(status, 'stale grid ID', failures)
    call first_writer%close(status)
    call require_ok(status, 'close reopened owner', failures)
  end subroutine exercise_writer_ownership

  subroutine exercise_pair_protection(directory, failures)
    character(len=*), intent(in) :: directory
    integer, intent(inout) :: failures

    type(xdmf_writer_t) :: writer
    type(xdmf_options_t) :: options
    type(xdmf_status_t) :: status
    character(len=:), allocatable :: path
    integer :: unit

    path = pair_path(directory, 'protected')
    open(newunit=unit, file=path//'.xdmf', status='replace', action='write')
    write(unit, '(A)') 'sentinel'
    close(unit)
    options%overwrite = .false.
    call writer%create(path, options, status)
    call require_error(status, 'existing XDMF protection', failures)
  end subroutine exercise_pair_protection

  function pair_path(directory, name) result(path)
    character(len=*), intent(in) :: directory, name
    character(len=:), allocatable :: path

    path = trim(directory)//'/'//trim(name)
  end function pair_path

  subroutine require_available_output_paths(directory, examples_only)
    character(len=*), intent(in) :: directory
    logical, intent(in) :: examples_only

    character(len=20), parameter :: example_names(6) = [character(len=20) :: &
      'uniform', 'volume', 'curvilinear', 'unstructured', 'mixed', &
      'time-series']
    character(len=20), parameter :: all_names(14) = [character(len=20) :: &
      example_names, 'rectilinear', 'frequency-series', 'escaped-&-pair', &
      'validation-static', 'validation-series', 'owner-first', 'owner-second', &
      'protected']
    integer :: i, name_count
    logical :: hdf5_exists, xdmf_exists

    if (examples_only) then
      name_count = size(example_names)
    else
      name_count = size(all_names)
    end if
    do i = 1, name_count
      if (examples_only) then
        inquire(file=pair_path(directory, example_names(i))//'.h5', &
          exist=hdf5_exists)
        inquire(file=pair_path(directory, example_names(i))//'.xdmf', &
          exist=xdmf_exists)
      else
        inquire(file=pair_path(directory, all_names(i))//'.h5', &
          exist=hdf5_exists)
        inquire(file=pair_path(directory, all_names(i))//'.xdmf', &
          exist=xdmf_exists)
      end if
      if (hdf5_exists .or. xdmf_exists) then
        write(*, '(A)') 'Error: output already exists: '// &
          merge(pair_path(directory, example_names(i)), &
            pair_path(directory, all_names(i)), examples_only)
        write(*, '(A)') 'Use --replace to replace generated output pairs.'
        stop 1
      end if
    end do
  end subroutine require_available_output_paths

  subroutine require_ok(status, context, failures)
    type(xdmf_status_t), intent(in) :: status
    character(len=*), intent(in) :: context
    integer, intent(inout) :: failures

    if (status%is_error()) then
      failures = failures + 1
      write(*, '(A)') trim(context)//': '//status%message()
    end if
  end subroutine require_ok

  subroutine require_error(status, context, failures)
    type(xdmf_status_t), intent(in) :: status
    character(len=*), intent(in) :: context
    integer, intent(inout) :: failures

    if (.not. status%is_error()) then
      failures = failures + 1
      write(*, '(A)') trim(context)//': expected an error status'
    end if
  end subroutine require_error

  subroutine print_help()
    write(*, '(A)') 'Usage'
    write(*, '(A)') '  xdmf_hdf5_generate_cases <output-directory> [options]'
    write(*, '()')
    write(*, '(A)') 'Generate XDMF/HDF5 example files in <output-directory>.'
    write(*, '()')
    write(*, '(A)') 'Options'
    write(*, '(A)') '  --examples  Generate only the committed examples.'
    write(*, '(A)') '  --help      Print this help information and exit.'
    write(*, '(A)') '  --replace   Replace generated files that already exist.'
    write(*, '()')
    write(*, '(A)') 'The output directory is created when necessary. Existing generated files'
    write(*, '(A)') 'are preserved unless --replace is provided. Relative paths are resolved'
    write(*, '(A)') 'from the current working directory.'
  end subroutine print_help

end program generate_cases
