program generate_parallel_case
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use mpi
  use xdmf_hdf5_m

  implicit none

  type(xdmf_writer_t) :: writer
  type(xdmf_options_t) :: options
  type(xdmf_status_t) :: status
  type(xdmf_grid_id_t) :: grid
  type(xdmf_attribute_id_t) :: field, static_field
  character(len=1024) :: output_directory
  integer(int64) :: local_offset(3), local_count(3)
  real(real64), allocatable :: values(:)
  integer :: error, i, j, k, local_k, rank, rank_count, step, value_index

  call MPI_Init(error)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, error)
  call MPI_Comm_size(MPI_COMM_WORLD, rank_count, error)
  if (rank_count /= 3) call abort_test('This test requires exactly three ranks')

  call get_command_argument(1, output_directory)
  if (len_trim(output_directory) == 0) then
    call abort_test('An output directory argument is required')
  end if

  options%overwrite = .true.
  options%series_kind = XDMF_SERIES_TIME
  options%collective_io = .true.
  options%communicator = MPI_COMM_WORLD
  options%root_rank = 0
  call writer%create(trim(output_directory)//'/parallel-hyperslab', &
    options, status)
  call require_ok(status, 'create collective writer')

  call writer%define_uniform_grid('global-volume', &
    [2_int64, 2_int64, 4_int64], &
    [0.0_real64, 0.0_real64, 0.0_real64], &
    [1.0_real64, 1.0_real64, 1.0_real64], grid, status)
  call require_ok(status, 'define global grid')
  call writer%define_attribute(grid, 'field', XDMF_CENTER_NODE, &
    XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., field, status)
  call require_ok(status, 'define distributed field')
  call writer%define_attribute(grid, 'static-field', XDMF_CENTER_NODE, &
    XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, static_field, status)
  call require_ok(status, 'define distributed static field')

  if (rank < 2) then
    local_offset = [0_int64, 0_int64, int(2 * rank, int64)]
    local_count = [2_int64, 2_int64, 2_int64]
    allocate(values(8))
  else
    local_offset = 0_int64
    local_count = 0_int64
    allocate(values(0))
  end if

  if (rank < 2) then
    value_index = 0
    do local_k = 1, 2
      k = 2 * rank + local_k
      do j = 1, 2
        do i = 1, 2
          value_index = value_index + 1
          values(value_index) = real(100 * k + 10 * j + i, real64)
        end do
      end do
    end do
  end if
  call writer%write_attribute_hyperslab(static_field, values, local_offset, &
    local_count, status)
  call require_ok(status, 'write collective static hyperslab')

  do step = 1, 2
    if (rank < 2) then
      value_index = 0
      do local_k = 1, 2
        k = 2 * rank + local_k
        do j = 1, 2
          do i = 1, 2
            value_index = value_index + 1
            values(value_index) = real(1000 * step + 100 * k + 10 * j + i, &
              real64)
          end do
        end do
      end do
    end if

    call writer%begin_step(0.5_real64 * real(step - 1, real64), status)
    call require_ok(status, 'begin collective step')
    call writer%write_attribute_hyperslab(field, values, local_offset, &
      local_count, status)
    call require_ok(status, 'write collective hyperslab')
    call writer%end_step(status)
    call require_ok(status, 'end collective step')
  end do

  call writer%close(status)
  call require_ok(status, 'close collective writer')
  call MPI_Finalize(error)

contains

  subroutine require_ok(result, operation)
    type(xdmf_status_t), intent(in) :: result
    character(len=*), intent(in) :: operation

    if (result%is_error()) then
      call abort_test(trim(operation)//': '//result%message())
    end if
  end subroutine require_ok

  subroutine abort_test(message)
    character(len=*), intent(in) :: message
    integer :: abort_error

    if (rank == 0) write(*, '(A)') trim(message)
    call MPI_Abort(MPI_COMM_WORLD, 1, abort_error)
  end subroutine abort_test

end program generate_parallel_case
