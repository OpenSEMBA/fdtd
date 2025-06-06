module rotate_testingTools
   use NFDETypes
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
   subroutine expect_eq_real(err, ex, pr, msg) 
      integer , intent(inout) :: err
      real (KIND=RK), intent(in) :: ex, pr
      character (len=*), intent(in), optional :: msg
      if (ex /= pr) call testFails(err, msg)
   end subroutine

   subroutine expect_eq_real_vect(err, ex, pr, msg) 
      integer , intent(inout) :: err
      real (KIND=RK), DIMENSION (:), intent(in) :: ex, pr
      character (len=*), intent(in), optional :: msg
      if (all(ex/=pr)) call testFails(err, msg)
   end subroutine

   subroutine testFails(err, msg)
      integer, intent(inout) :: err
      character(len=*), intent(in), optional :: msg
      err = err + 1
      if (present(msg)) write(*, *)  "FAIL: "// msg
   end subroutine
  
   subroutine createDemoDesplazamiento(des)
      type(Desplazamiento) :: des
      des%desX = 1.0_RKIND
      des%desY = 2.0_RKIND
      des%desZ = 3.0_RKIND
      des%mx1 = 1
      des%mx2 = 4
      des%my1 = 2
      des%my2 = 5
      des%mz1 = 3
      des%mz2 = 6
      des%nX = 10
      des%nY = 20 
      des%nZ = 30
      des%originx = 1.0_RKIND
      des%originy = 2.0_RKIND
      des%originz = 3.0_RKIND
   end subroutine

end module