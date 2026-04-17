module assertionTools_m
   use FDETYPES_m
   use arrayAssertionTools_m
   use iso_fortran_env, only: real32, real64
   implicit none

   private :: assert_real_equal_impl
   private :: assert_real_RKIND_equal_impl
#ifndef CompileWithReal8
   private :: assert_real_time_equal_impl
#endif
   ! Generic interface for real assertions
   interface assert_real_equal
      module procedure assert_real_equal_impl
      module procedure assert_real_RKIND_equal_impl
#ifndef CompileWithReal8
      module procedure assert_real_time_equal_impl
#endif
   end interface

contains

   !---------------------------------------
   ! Logical assertion
   !---------------------------------------
   function assert_true(boolean, errorMessage) result(err)
      logical, intent(in) :: boolean
      character(*), intent(in) :: errorMessage
      integer :: err
      if (boolean) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   !---------------------------------------
   ! Integer equality
   !---------------------------------------
   function assert_integer_equal(val, expected, errorMessage) result(err)
      integer, intent(in) :: val, expected
      character(*), intent(in) :: errorMessage
      integer :: err
      if (val == expected) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED: ', trim(errorMessage)
         print *, "  Value: ", val, ". Expected: ", expected
      end if
   end function

   !---------------------------------------
   ! Real equality implementations
   !---------------------------------------
   function assert_real_equal_impl(val, expected, tolerance, errorMessage) result(err)
      real(real32), intent(in) :: val, expected, tolerance
      character(*), intent(in) :: errorMessage
      integer :: err
      if (abs(val - expected) <= tolerance) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED (real32): ', trim(errorMessage)
         print *, '  Value: ', val, '. Expected: ', expected, '. Tolerance: ', tolerance
      end if
   end function

   function assert_real_RKIND_equal_impl(val, expected, tolerance, errorMessage) result(err)
      real(RKIND), intent(in) :: val, expected, tolerance
      character(*), intent(in) :: errorMessage
      integer :: err
      if (abs(val - expected) <= tolerance) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED (RKIND): ', trim(errorMessage)
         print *, '  Value: ', val, '. Expected: ', expected, '. Tolerance: ', tolerance
      end if
   end function

   function assert_real_time_equal_impl(val, expected, tolerance, errorMessage) result(err)
      real(kind=RKIND_tiempo), intent(in) :: val, expected, tolerance
      character(*), intent(in) :: errorMessage
      integer :: err
      if (abs(val - expected) <= tolerance) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED (time): ', trim(errorMessage)
         print *, '  Value: ', val, '. Expected: ', expected, '. Tolerance: ', tolerance
      end if
   end function

   !---------------------------------------
   ! Complex equality
   !---------------------------------------
   function assert_complex_equal(val, expected, tolerance, errorMessage) result(err)
      complex(kind=CKIND), intent(in) :: val, expected
      real(kind=RKIND), intent(in) :: tolerance
      character(len=*), intent(in) :: errorMessage
      integer :: err
      if (abs(val - expected) <= tolerance) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED: ', trim(errorMessage)
         print *, '  Value:    ', val
         print *, '  Expected: ', expected
         print *, '  Delta:    ', abs(val - expected)
         print *, '  Tolerance:', tolerance
      end if
   end function

   !---------------------------------------
   ! String equality
   !---------------------------------------
   function assert_string_equal(val, expected, errorMessage) result(err)
      character(*), intent(in) :: val, expected
      character(*), intent(in) :: errorMessage
      integer :: err
      if (trim(val) == trim(expected)) then
         err = 0
      else
         err = 1
         print *, 'ASSERTION FAILED: ', trim(errorMessage)
         print *, '  Value: "', trim(val), '". Expected: "', trim(expected), '"'
      end if
   end function

   !---------------------------------------
   ! Check if file was written
   !---------------------------------------
   integer function assert_written_output_file(filename) result(code)
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
   end function

   !---------------------------------------
   ! Check file content
   !---------------------------------------
   integer function assert_file_content(unit, expectedValues, nRows, nCols, tolerance, headers) result(flag)
      integer(kind=SINGLE), intent(in) :: unit
      real(kind=RKIND), intent(in) :: expectedValues(:, :)
      integer(kind=SINGLE), intent(in) :: nRows, nCols
      real(kind=RKIND), intent(in) :: tolerance
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
            if (abs(val(j) - expectedValues(i, j)) > tolerance) then
               flag = flag + 1
            end if
         end do
      end do
   end function

   !---------------------------------------
   ! Check file exists
   !---------------------------------------
   integer function assert_file_exists(fileName) result(err)
      character(len=*), intent(in) :: filename
      integer :: unit, ios
      err = 0
      open (newunit=unit, file=filename, status='old', iostat=ios)
      close (unit)
      if (ios /= 0) err = 1
   end function

end module assertionTools_m
