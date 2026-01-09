module mod_arrayAssertionTools
   use FDETYPES
   implicit none
   real(RKIND), parameter :: tol = 1.0e-12_RKIND
   private
   !-----------------------------
   ! Public assertion procedures
   !-----------------------------
   public :: assert_arrays_equal
   public :: assert_array_value

   !---------------------------------------
   ! GENERIC INTERFACES
   !---------------------------------------
   interface assert_arrays_equal
      module procedure &
         assert_arrays_equal_int1, assert_arrays_equal_int2, assert_arrays_equal_int3, &
         assert_arrays_equal_real1, assert_arrays_equal_real2, assert_arrays_equal_real3, &
         assert_arrays_equal_complex1, assert_arrays_equal_complex2, assert_arrays_equal_complex3
   end interface

   interface assert_array_value
      module procedure &
         assert_array_value_int1, assert_array_value_int2, assert_array_value_int3, &
         assert_array_value_real1, assert_array_value_real2, assert_array_value_real3, &
         assert_array_value_complex1, assert_array_value_complex2, assert_array_value_complex3
   end interface

contains

   !---------------------------------------
   ! 1D Integer arrays
   !---------------------------------------
   integer function assert_arrays_equal_int1(A, B, errorMessage)
      integer, intent(in) :: A(:), B(:)
      character(*), intent(in), optional :: errorMessage

      if (any(shape(A) /= shape(B))) then
         assert_arrays_equal_int1 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
         return
      end if

      if (all(A == B)) then
         assert_arrays_equal_int1 = 0
      else
         assert_arrays_equal_int1 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   integer function assert_array_value_int1(A, val, errorMessage)
      integer, intent(in) :: A(:)
      integer, intent(in) :: val
      character(*), intent(in), optional :: errorMessage

      if (all(A == val)) then
         assert_array_value_int1 = 0
      else
         assert_array_value_int1 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   !---------------------------------------
   ! 2D Integer arrays
   !---------------------------------------
   integer function assert_arrays_equal_int2(A, B, errorMessage)
      integer, intent(in) :: A(:, :), B(:, :)
      character(*), intent(in), optional :: errorMessage

      if (any(shape(A) /= shape(B))) then
         assert_arrays_equal_int2 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
         return
      end if

      if (all(A == B)) then
         assert_arrays_equal_int2 = 0
      else
         assert_arrays_equal_int2 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   integer function assert_array_value_int2(A, val, errorMessage)
      integer, intent(in) :: A(:, :)
      integer, intent(in) :: val
      character(*), intent(in), optional :: errorMessage

      if (all(A == val)) then
         assert_array_value_int2 = 0
      else
         assert_array_value_int2 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   !---------------------------------------
   ! 3D Integer arrays
   !---------------------------------------
   integer function assert_arrays_equal_int3(A, B, errorMessage)
      integer, intent(in) :: A(:, :, :), B(:, :, :)
      character(*), intent(in), optional :: errorMessage

      if (any(shape(A) /= shape(B))) then
         assert_arrays_equal_int3 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
         return
      end if

      if (all(A == B)) then
         assert_arrays_equal_int3 = 0
      else
         assert_arrays_equal_int3 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   integer function assert_array_value_int3(A, val, errorMessage)
      integer, intent(in) :: A(:, :, :)
      integer, intent(in) :: val
      character(*), intent(in), optional :: errorMessage

      if (all(A == val)) then
         assert_array_value_int3 = 0
      else
         assert_array_value_int3 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   !---------------------------------------
   ! REAL arrays (1D, 2D, 3D)
   !---------------------------------------
   integer function assert_arrays_equal_real1(A, B, errorMessage)
      real(RKIND), intent(in) :: A(:), B(:)
      character(*), intent(in), optional :: errorMessage

      if (any(shape(A) /= shape(B))) then
         assert_arrays_equal_real1 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
         return
      end if

      if (all(abs(A - B) < tol)) then
         assert_arrays_equal_real1 = 0
      else
         assert_arrays_equal_real1 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   integer function assert_array_value_real1(A, val, errorMessage)
      real(RKIND), intent(in) :: A(:)
      real(RKIND), intent(in) :: val
      character(*), intent(in), optional :: errorMessage

      if (all(abs(A - val) < tol)) then
         assert_array_value_real1 = 0
      else
         assert_array_value_real1 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   !---------------------------------------
   ! REAL 2D
   !---------------------------------------
   integer function assert_arrays_equal_real2(A, B, errorMessage)
      real(RKIND), intent(in) :: A(:, :), B(:, :)
      character(*), intent(in), optional :: errorMessage

      if (any(shape(A) /= shape(B))) then
         assert_arrays_equal_real2 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
         return
      end if

      if (all(abs(A - B) < tol)) then
         assert_arrays_equal_real2 = 0
      else
         assert_arrays_equal_real2 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   integer function assert_array_value_real2(A, val, errorMessage)
      real(RKIND), intent(in) :: A(:, :)
      real(RKIND), intent(in) :: val
      character(*), intent(in), optional :: errorMessage

      if (all(abs(A - val) < tol)) then
         assert_array_value_real2 = 0
      else
         assert_array_value_real2 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   !---------------------------------------
   ! REAL 3D
   !---------------------------------------
   integer function assert_arrays_equal_real3(A, B, errorMessage)
      real(RKIND), intent(in) :: A(:, :, :), B(:, :, :)
      character(*), intent(in), optional :: errorMessage

      if (any(shape(A) /= shape(B))) then
         assert_arrays_equal_real3 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
         return
      end if

      if (all(abs(A - B) < tol)) then
         assert_arrays_equal_real3 = 0
      else
         assert_arrays_equal_real3 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   integer function assert_array_value_real3(A, val, errorMessage)
      real(RKIND), intent(in) :: A(:, :, :)
      real(RKIND), intent(in) :: val
      character(*), intent(in), optional :: errorMessage

      if (all(abs(A - val) < tol)) then
         assert_array_value_real3 = 0
      else
         assert_array_value_real3 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   !---------------------------------------
   ! COMPLEX 1D arrays
   !---------------------------------------
   integer function assert_arrays_equal_complex1(A, B, errorMessage)
      complex(CKIND), intent(in) :: A(:), B(:)
      character(*), intent(in), optional :: errorMessage

      if (any(shape(A) /= shape(B))) then
         assert_arrays_equal_complex1 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
         return
      end if

      if (all(abs(A - B) < tol)) then
         assert_arrays_equal_complex1 = 0
      else
         assert_arrays_equal_complex1 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   integer function assert_array_value_complex1(A, val, errorMessage)
      complex(CKIND), intent(in) :: A(:)
      complex(CKIND), intent(in) :: val
      character(*), intent(in), optional :: errorMessage

      if (all(abs(A - val) < tol)) then
         assert_array_value_complex1 = 0
      else
         assert_array_value_complex1 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

!---------------------------------------
! COMPLEX 2D arrays
!---------------------------------------
   integer function assert_arrays_equal_complex2(A, B, errorMessage)
      complex(CKIND), intent(in) :: A(:, :), B(:, :)
      character(*), intent(in), optional :: errorMessage

      if (any(shape(A) /= shape(B))) then
         assert_arrays_equal_complex2 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
         return
      end if

      if (all(abs(A - B) < tol)) then
         assert_arrays_equal_complex2 = 0
      else
         assert_arrays_equal_complex2 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   integer function assert_array_value_complex2(A, val, errorMessage)
      complex(CKIND), intent(in) :: A(:, :)
      complex(CKIND), intent(in) :: val
      character(*), intent(in), optional :: errorMessage

      if (all(abs(A - val) < tol)) then
         assert_array_value_complex2 = 0
      else
         assert_array_value_complex2 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

!---------------------------------------
! COMPLEX 3D arrays
!---------------------------------------
   integer function assert_arrays_equal_complex3(A, B, errorMessage)
      complex(CKIND), intent(in) :: A(:, :, :), B(:, :, :)
      character(*), intent(in), optional :: errorMessage

      if (any(shape(A) /= shape(B))) then
         assert_arrays_equal_complex3 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
         return
      end if

      if (all(abs(A - B) < tol)) then
         assert_arrays_equal_complex3 = 0
      else
         assert_arrays_equal_complex3 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

   integer function assert_array_value_complex3(A, val, errorMessage)
      complex(CKIND), intent(in) :: A(:, :, :)
      complex(CKIND), intent(in) :: val
      character(*), intent(in), optional :: errorMessage

      if (all(abs(A - val) < tol)) then
         assert_array_value_complex3 = 0
      else
         assert_array_value_complex3 = 1
         if (present(errorMessage)) print *, 'ASSERTION FAILED: ', trim(errorMessage)
      end if
   end function

end module mod_arrayAssertionTools
