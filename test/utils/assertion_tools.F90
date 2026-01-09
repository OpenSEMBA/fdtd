module mod_assertionTools
  use FDETYPES
  implicit none
  
contains
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

function assert_complex_equal(val, expected, tolerance, errorMessage) result(err)
    complex(kind=CKIND), intent(in) :: val, expected
    real   (kind=RKIND), intent(in) :: tolerance
    character(len=*), intent(in)    :: errorMessage
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
end function assert_complex_equal

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

   integer function assert_file_exists(fileName) result(err)
      character(len=*), intent(in) :: filename
      integer :: unit, ios
      err = 0
      open(newunit=unit, file=filename, status='old', iostat=ios)
      close(unit)
      if (ios/=0) err = 1
   end function
end module mod_assertionTools