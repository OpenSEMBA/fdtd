module rotate_testingTools
   implicit none

   character(len=*), parameter :: PATH_TO_TEST_DATA = 'testData/'
   character(len=*), parameter :: INPUT_EXAMPLES='input_examples/'
contains
    subroutine expect_eq_int(err, ex, pr, msg) 
      integer, intent(inout) :: err
      integer, intent(in) :: ex, pr
      character (len=*), intent(in), optional :: msg
      if (ex /= pr) call testFails(err, msg)
   end subroutine
   subroutine expect_eq_float(err, ex, pr, msg) 
      integer , intent(inout) :: err
      real, intent(in) :: ex, pr
      character (len=*), intent(in), optional :: msg
      if (ex /= pr) call testFails(err, msg)
   end subroutine
end module