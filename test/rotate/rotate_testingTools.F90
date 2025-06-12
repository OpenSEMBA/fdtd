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
      if (expected /= actual) then
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
      
      if (abs(expected - actual) > 1.0d-10) then
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
      
      if (abs(expected - actual) > 1.0d-10) then
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
  
   subroutine createDemoDesplazamiento(des)
      type(Desplazamiento) :: des
      REAL(kind=RKIND), dimension(1), target :: x=(/1.0/), y=(/2.0/), z=(/3.0/)
      des%desX => x
      des%desY => y
      des%desZ => z
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

   subroutine createDemoMatrizMedios(matriz)
      type(MatrizMedios) :: matriz
      matriz%totalX = 1
      matriz%totalY = 2
      matriz%totalZ = 3
   end subroutine

end module