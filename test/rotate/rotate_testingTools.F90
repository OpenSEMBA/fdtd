module rotate_testingTools
   use NFDETypes
   implicit none

   character(len=*), parameter :: PATH_TO_TEST_DATA = 'testData/'
   character(len=*), parameter :: INPUT_EXAMPLES='input_examples/'
contains
    subroutine expect_eq_int(test_err, expected, actual, message) 
      integer, intent(inout) :: test_err
      integer, intent(in) :: expected, actual
      character (len=*), intent(in), optional :: message
      if (expected /= actual) then
          test_err = test_err + 1
          print *, "Error: ", message
          print *, "Expected: ", expected
          print *, "Actual:   ", actual
      end if
   end subroutine
   subroutine expect_eq_real(test_err, expected, actual, message) 
      integer , intent(inout) :: test_err
      real (KIND=RK), intent(in) :: expected, actual
      character (len=*), intent(in), optional :: message
      if (abs(expected - actual) > 1.0d-6) then
          test_err = test_err + 1
          print *, "Error: ", message
          print *, "Expected: ", expected
          print *, "Actual:   ", actual
      end if
   end subroutine

   subroutine expect_eq_real_vect(err, ex, pr, msg) 
      integer , intent(inout) :: err
      real (KIND=RK), DIMENSION (:), intent(in) :: ex, pr
      character (len=*), intent(in), optional :: msg
      if (all(ex/=pr)) call testFails(err, msg)
   end subroutine

   subroutine expect_eq_double(test_err, expected, actual, message)
      integer, intent(inout) :: test_err
      real(8), intent(in) :: expected, actual
      character(len=*), intent(in) :: message
      
      if (abs(expected - actual) > 1.0d-7) then
          test_err = test_err + 1
          print *, "Error: ", message
          print *, "Expected: ", expected
          print *, "Actual:   ", actual
      end if
  end subroutine expect_eq_double

  subroutine expect_eq_complx(test_err, expected, actual, message)
      integer, intent(inout) :: test_err
      complex, intent(in) :: expected, actual
      character(len=*), intent(in) :: message
      
      if (abs(expected - actual) > 1.0d-7) then
          test_err = test_err + 1
          print *, "Error: ", message
          print *, "Expected: ", expected
          print *, "Actual:   ", actual
      end if
  end subroutine expect_eq_complx

   subroutine testFails(err, msg)
      integer, intent(inout) :: err
      character(len=*), intent(in), optional :: msg
      err = err + 1
      if (present(msg)) write(*, *)  "FAIL: "// msg
   end subroutine

end module