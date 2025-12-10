module mod_testOutputUtils
   use FDETYPES
   use FDETYPES_TOOLS

   implicit none
   type :: dummyFields_t
      real(kind=RKIND), allocatable, dimension(:, :, :) :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=RKIND), allocatable, dimension(:) :: dxe, dye, dze, dxh, dyh, dzh
   contains
      procedure, public :: createDummyFields => create_dummy_fields
   end type dummyFields_t
contains
   function create_point_probe_observable() result(obs)
      type(Obses_t) :: obs

      type(observable_t), dimension(:), allocatable :: P
      type(observation_domain_t) :: domain

      call initialize_time_domain(domain, 0.0_RKIND, 10.0_RKIND, 0.1_RKIND)
      allocate (P(1))
      P(1) = create_observable(4, 4, 4, 6, 6, 6, iEx)
      call set_observable(obs, P, 'poinProbe', domain, 'DummyFileNormalize')

   end function
   
   function create_volumic_probe_observable() result(obs)
      type(Obses_t) :: obs

      type(observable_t), dimension(:), allocatable :: P
      type(observation_domain_t) :: domain

      call initialize_time_domain(domain, 0.0_RKIND, 10.0_RKIND, 0.1_RKIND)
      call initialize_frequency_domain(domain, 0.0_RKIND, 1000.0_RKIND, 50.0_RKIND)
      allocate (P(1))
      P(1) = create_observable(4, 4, 4, 6, 6, 6, iCurX)
      call set_observable(obs, P, 'volumicProbe', domain, 'DummyFileNormalize')
   end function create_volumic_probe_observable

   subroutine create_dummy_fields(this, lower, upper, delta)
      class(dummyFields_t), intent(inout) :: this
      integer, intent(in) :: lower, upper
      real(kind=rkind), intent(in) :: delta
      allocate ( &
         this%Ex(lower:upper, lower:upper, lower:upper), &
         this%Ey(lower:upper, lower:upper, lower:upper), &
         this%Ez(lower:upper, lower:upper, lower:upper), &
         this%Hx(lower:upper, lower:upper, lower:upper), &
         this%Hy(lower:upper, lower:upper, lower:upper), &
         this%Hz(lower:upper, lower:upper, lower:upper) &
         )

      this%Ex = 0.0_RKIND
      this%Ey = 0.0_RKIND
      this%Ez = 0.0_RKIND
      this%Hx = 0.0_RKIND
      this%Hy = 0.0_RKIND
      this%Hz = 0.0_RKIND

      allocate ( &
         this%dxh(lower:upper), &
         this%dyh(lower:upper), &
         this%dzh(lower:upper), &
         this%dxe(lower:upper), &
         this%dye(lower:upper), &
         this%dze(lower:upper) &
         )
      this%dxh = delta
      this%dyh = delta
      this%dzh = delta
      this%dxe = delta
      this%dye = delta
      this%dze = delta
   end subroutine create_dummy_fields

   function assert_integer_equal(val, expected, errorMessage) result(err)

      integer, intent(in) :: val
      integer, intent(in) :: expected
      character(*), intent(in) :: errorMessage
      integer :: err

      if (val == expected) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED: ', trim(errorMessage)
         print *, "  Value: ", val, ". Expected: ", expected
      end if
   end function assert_integer_equal

   function assert_real_equal(val, expected, tolerance, errorMessage) result(err)

      real(kind=rkind), intent(in) :: val
      real(kind=rkind), intent(in) :: expected
      real(kind=rkind), intent(in) :: tolerance
      character(*), intent(in) :: errorMessage
      integer :: err

      if (abs(val - expected) <= tolerance) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED: ', trim(errorMessage)
         print *, "  Value: ", val, ". Expected: ", expected, ". Tolerance: ", tolerance
      end if
   end function assert_real_equal

   function assert_real_time_equal(val, expected, tolerance, errorMessage) result(err)

      real(kind=RKIND_tiempo), intent(in) :: val
      real(kind=RKIND_tiempo), intent(in) :: expected
      real(kind=RKIND_tiempo), intent(in) :: tolerance
      character(*), intent(in) :: errorMessage
      integer :: err

      if (abs(val - expected) <= tolerance) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED: ', trim(errorMessage)
         print *, "  Value: ", val, ". Expected: ", expected, ". Tolerance: ", tolerance
      end if
   end function assert_real_time_equal

   function assert_string_equal(val, expected, errorMessage) result(err)

      character(*), intent(in) :: val
      character(*), intent(in) :: expected
      character(*), intent(in) :: errorMessage
      integer :: err

      if (trim(val) == trim(expected)) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED: ', trim(errorMessage)
         print *, '  Value: "', trim(val), '". Expected: "', trim(expected), '"'
      end if
   end function assert_string_equal

   integer function assert_written_output_file(filename) result(code)
      implicit none
      character(len=*), intent(in) :: filename
      logical :: ex
      integer :: filesize

      code = 0

      inquire (file=filename, exist=ex, size=filesize)

      if (.not. ex) then
         print *, "ERROR: Output file not created:", trim(filename)
         code = 1
      else if (filesize <= 0) then
         print *, "ERROR: Output file is empty:", trim(filename)
         code = 2
      end if
   end function assert_written_output_file

   integer function assert_file_content(unit, expectedValues, nRows, nCols, headers) result(flag)
      implicit none
      integer(kind=SINGLE), intent(in) :: unit
      real(kind=RKIND), intent(in) :: expectedValues(:, :)
      integer(kind=SINGLE), intent(in) :: nRows, nCols
      character(len=*), intent(in), optional :: headers(:)
      integer(kind=SINGLE) :: i, j, ios
      real(kind=RKIND), dimension(nCols) :: val
      character(len=BUFSIZE) :: line
      flag = 0

      if (present(headers)) then
         read (unit, '(F12.6,1X,F12.6)', iostat=ios) line
         if (ios /= 0) return
      end if

      do i = 1, nRows
         read (unit, *, iostat=ios) val
         if (ios /= 0) then
            flag = flag + 1
            return
         end if
         do j = 1, nCols
            if (abs(val(j) - expectedValues(i, j)) > 1d-6) then
               flag = flag + 1
            end if
         end do
      end do
   end function assert_file_content

end module mod_testOutputUtils
