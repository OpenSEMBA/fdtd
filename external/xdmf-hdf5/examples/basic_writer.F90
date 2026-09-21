program basic_writer
  use, intrinsic :: iso_fortran_env, only: int64, real64
  use xdmf_hdf5_m

  implicit none

  type(xdmf_writer_t) :: writer
  type(xdmf_options_t) :: options
  type(xdmf_status_t) :: status
  type(xdmf_grid_id_t) :: grid
  type(xdmf_attribute_id_t) :: pressure
  real(real64) :: values(24)
  integer :: index

  options%overwrite = .true.
  options%series_kind = XDMF_SERIES_TIME
  call writer%create('example', options, status)
  call check(status)

  call writer%define_uniform_grid('volume', &
    [2_int64, 3_int64, 4_int64], &
    [0.0_real64, 0.0_real64, 0.0_real64], &
    [1.0_real64, 1.0_real64, 1.0_real64], grid, status)
  call check(status)
  call writer%define_attribute(grid, 'pressure', XDMF_CENTER_NODE, &
    XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., pressure, status)
  call check(status)

  do index = 1, size(values)
    values(index) = real(index, real64)
  end do
  call writer%begin_step(0.0_real64, status)
  call check(status)
  call writer%write_attribute(pressure, values, status)
  call check(status)
  call writer%end_step(status)
  call check(status)
  call writer%close(status)
  call check(status)

contains

  subroutine check(result)
    type(xdmf_status_t), intent(in) :: result

    if (result%is_error()) then
      write(*, '(A)') result%message()
      error stop 1
    end if
  end subroutine check

end program basic_writer
