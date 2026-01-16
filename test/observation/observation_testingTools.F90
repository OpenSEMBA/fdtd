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
      real(kind=RKIND), intent(in), dimension(:, :) :: arr
      integer, intent(in) :: n_expected
      integer, intent(inout) :: test_err
      character(len=*), intent(in), optional :: name

      integer :: rank_arr
      integer, allocatable :: shp(:)
      character(len=:), allocatable :: nm

      if (present(name)) then
          nm = trim(adjustl(name))
      else
          nm = "array"
      end if

      rank_arr = rank(arr)
      shp = shape(arr)

      ! Expect exactly two dimensions [1, n_expected]
      if (rank_arr /= 2 .or. any(shp /= [1, n_expected])) then
         test_err = test_err + 1
      end if
   end subroutine check_shape_real

   subroutine check_shape_complex(arr, n_expected, test_err, name)
      use Observa
      complex(kind=CKIND), intent(in), dimension(:, :) :: arr
      integer, intent(in) :: n_expected
      integer, intent(inout) :: test_err
      character(len=*), intent(in), optional :: name

      integer :: rank_arr
      integer, allocatable :: shp(:)
      character(len=:), allocatable :: nm

      if (present(name)) then
          nm = trim(adjustl(name))
      else
          nm = "array"
      end if

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

      if (present(name)) then
          nm = trim(adjustl(name))
      else
          nm = "array"
      end if

      rank_arr = rank(arr)
      siz = size(arr)

      if (rank_arr /= 1 .or. siz /= n_expected) then
         test_err = test_err + 1
      end if
   end subroutine check_size

   logical function approx_equal(a, b, tol) result(equal)
      real(kind=RKIND), intent(in) :: a, b, tol
      equal = abs(a - b) <= tol
   end function approx_equal

   function create_xyz_limit_array(XI,YI,ZI,XE,YE,ZE) result(arr)
      type(XYZlimit_t), dimension(1:6) :: arr
      integer (kind=4), intent(in) :: XI,YI,ZI,XE,YE,ZE
      integer :: i
      do i = 1, 6
         arr(i)%XI = XI
         arr(i)%XE = XE
         arr(i)%YI = YI
         arr(i)%YE = YE
         arr(i)%ZI = ZI
         arr(i)%ZE = ZE
      end do
   end function create_xyz_limit_array

   function create_basic_media () result(media)
      type(MediaData_t) :: media
   end function create_basic_media
end module
