module observation_testingTools
   use FDETYPES
   type :: dummyFields_t
      real(kind=RKIND),allocatable, dimension(:,:,:) :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=RKIND),allocatable, dimension(:) :: dxe, dye, dze, dxh, dyh, dzh
      contains
      procedure, public :: createDummyFields => create_dummy_fields
   end type dummyFields_t

   public dummyFields
   public check_shape_real, check_shape_complex, check_size
   public approx_equal
contains
   subroutine create_dummy_fields(this, lower, upper, delta)
      class(dummyFields_t), intent(inout) :: this
      integer, intent(in) :: lower, upper
      real(kind=rkind), intent(in) :: delta
      allocate(&
         this%Ex(lower:upper, lower:upper, lower:upper),&
         this%Ey(lower:upper, lower:upper, lower:upper),&
         this%Ez(lower:upper, lower:upper, lower:upper),&
         this%Hx(lower:upper, lower:upper, lower:upper),&
         this%Hy(lower:upper, lower:upper, lower:upper),&
         this%Hz(lower:upper, lower:upper, lower:upper)&
      )

      this%Ex = 0.0_RKIND
      this%Ey = 0.0_RKIND
      this%Ez = 0.0_RKIND
      this%Hx = 0.0_RKIND
      this%Hy = 0.0_RKIND
      this%Hz = 0.0_RKIND

      allocate(&
         this%dxh(lower:upper), &
         this%dyh(lower:upper), &
         this%dzh(lower:upper), &
         this%dxe(lower:upper), &
         this%dye(lower:upper), &
         this%dze(lower:upper)&
      )
      this%dxh = delta
      this%dyh = delta
      this%dzh = delta
      this%dxe = delta
      this%dye = delta
      this%dze = delta
   end subroutine create_dummy_fields

   subroutine check_shape_real(arr, n_expected, test_err, name)
      use Observa
      use FDETYPES
      real(kind=RKIND), intent(in), dimension(:, :) :: arr
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
      use FDETYPES
      complex(kind=CKIND), intent(in), dimension(:, :) :: arr
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
      use FDETYPES
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

   logical function approx_equal(a, b, tol) result(equal)
      use FDETYPES
      real(kind=RKIND), intent(in) :: a, b, tol
      equal = abs(a - b) <= tol
   end function approx_equal
end module
