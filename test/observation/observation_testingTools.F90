module observation_testingTools
  use FDETYPES
  use OBSERVA
  implicit none
  public check_shape_real, check_shape_complex, check_size
  public approx_equal
contains
  subroutine check_shape_real(arr, n_expected, test_err, name)
    use Observa
    real(kind=8), intent(in), dimension(:, :) :: arr
    integer, intent(in) :: n_expected
    integer, intent(inout) :: test_err
    character(len=*), intent(in), optional :: name

    integer :: rank_arr
    integer, allocatable :: shp(:)
    character(len=:), allocatable :: nm

    nm = merge(name, "array", present(name))

    rank_arr = rank(arr)
    shp = shape(arr)

    ! Expect exactly two dimensions [1, n_expected]
    if (rank_arr /= 2 .or. any(shp /= [1, n_expected])) then
      test_err = test_err + 1
    end if
  end subroutine check_shape_real

  subroutine check_shape_complex(arr, n_expected, test_err, name)
    use Observa
    complex(kind=8), intent(in), dimension(:, :) :: arr
    integer, intent(in) :: n_expected
    integer, intent(inout) :: test_err
    character(len=*), intent(in), optional :: name

    integer :: rank_arr
    integer, allocatable :: shp(:)
    character(len=:), allocatable :: nm

    nm = merge(name, "array", present(name))

    rank_arr = rank(arr)
    shp = shape(arr)

    ! Expect exactly two dimensions [1, n_expected]
    if (rank_arr /= 2 .or. any(shp /= [1, n_expected])) then
      test_err = test_err + 1
    end if
  end subroutine check_shape_complex

  subroutine check_size(arr, n_expected, test_err, name)
    use Observa
    integer, intent(in), dimension(:) :: arr
    integer, intent(in) :: n_expected
    integer, intent(inout) :: test_err
    character(len=*), intent(in), optional :: name

    integer :: rank_arr
    integer :: siz
    character(len=:), allocatable :: nm

    nm = merge(name, "array", present(name))

    rank_arr = rank(arr)
    siz = size(arr)

    if (rank_arr /= 1 .or. siz /= n_expected) then
      test_err = test_err + 1
    end if
  end subroutine check_size

  logical function approx_equal(a, b, tol)
    real(kind=RKIND), intent(in) :: a, b, tol
    approx_equal = abs(a - b) <= tol
  end function approx_equal

  function create_time_array(array_size, interval) result(arr)
    integer, intent(in) :: array_size
    real(kind=RKIND_tiempo) :: interval

    real(kind=RKIND_tiempo), dimension(:), allocatable :: arr
    allocate (arr(array_size))

    DO i = 1, array_size
      arr(i) = (i - 1)*interval
    END DO
  end function

end module
