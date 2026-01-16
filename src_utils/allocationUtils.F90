module mod_allocationUtils
   use FDETYPES, only: RKIND, CKIND, SINGLE, RKIND_tiempo, IKINDMTAG, INTEGERSIZEOFMEDIAMATRICES
   implicit none
   private
   public :: alloc_and_init

   interface alloc_and_init
      procedure alloc_and_init_int_1D
      procedure alloc_and_init_int_2D
      procedure alloc_and_init_int_3D
      procedure alloc_and_init_real_1D
      procedure alloc_and_init_real_2D
      procedure alloc_and_init_real_3D
      procedure alloc_and_init_complex_1D
      procedure alloc_and_init_complex_2D
      procedure alloc_and_init_complex_3D
      procedure alloc_and_init_int_3D_tag
      procedure alloc_and_init_int_3D_med
   end interface
contains

   subroutine alloc_and_init_int_1D(array, n1, initVal)
      integer(SINGLE), allocatable, intent(inout) :: array(:)
      integer, intent(IN) :: n1
      integer(SINGLE), intent(IN) :: initVal

      allocate (array(n1))
      array = initVal
   END subroutine alloc_and_init_int_1D

   subroutine alloc_and_init_int_2D(array, n1, n2, initVal)
      integer(SINGLE), allocatable, intent(inout) :: array(:, :)
      integer, intent(IN) :: n1, n2
      integer(SINGLE), intent(IN) :: initVal

      allocate (array(n1, n2))
      array = initVal
   END subroutine alloc_and_init_int_2D

   subroutine alloc_and_init_int_3D(array, n1, n2, n3, initVal)
      integer(SINGLE), allocatable, intent(inout) :: array(:, :, :)
      integer, intent(IN) :: n1, n2, n3
      integer(SINGLE), intent(IN) :: initVal

      allocate (array(n1, n2, n3))
      array = initVal
   END subroutine alloc_and_init_int_3D

   ! Allocate array of kind=IKINDMTAG
   subroutine alloc_and_init_int_3D_tag(array, n1_min, n1_max, n2_min, n2_max, n3_min, n3_max, initVal)
      integer(kind=IKINDMTAG), allocatable, intent(inout) :: array(:, :, :)
      integer, intent(in) :: n1_min, n1_max, n2_min, n2_max, n3_min, n3_max
      integer(kind=IKINDMTAG), intent(in) :: initVal

      if (allocated(array)) deallocate (array)
      allocate (array(n1_min:n1_max, n2_min:n2_max, n3_min:n3_max))
      array = initVal
   end subroutine

   ! Allocate array of kind=INTEGERSIZEOFMEDIAMATRICES
   subroutine alloc_and_init_int_3D_med(array, n1_min, n1_max, n2_min, n2_max, n3_min, n3_max, initVal)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES), allocatable, intent(inout) :: array(:, :, :)
      integer, intent(in) :: n1_min, n1_max, n2_min, n2_max, n3_min, n3_max
      integer(kind=INTEGERSIZEOFMEDIAMATRICES), intent(in) :: initVal

      if (allocated(array)) deallocate (array)
      allocate (array(n1_min:n1_max, n2_min:n2_max, n3_min:n3_max))
      array = initVal
   end subroutine

   subroutine alloc_and_init_real_1D(array, n1, initVal)
      REAL(RKIND), allocatable, intent(inout) :: array(:)
      integer, intent(IN) :: n1
      REAL(RKIND), intent(IN) :: initVal

      allocate (array(n1))
      array = initVal
   END subroutine alloc_and_init_real_1D

   subroutine alloc_and_init_real_2D(array, n1, n2, initVal)
      REAL(RKIND), allocatable, intent(inout) :: array(:, :)
      integer, intent(IN) :: n1, n2
      REAL(RKIND), intent(IN) :: initVal

      allocate (array(n1, n2))
      array = initVal
   END subroutine alloc_and_init_real_2D

   subroutine alloc_and_init_real_3D(array, n1, n2, n3, initVal)
      REAL(RKIND), allocatable, intent(inout) :: array(:, :, :)
      integer, intent(IN) :: n1, n2, n3
      REAL(RKIND), intent(IN) :: initVal

      allocate (array(n1, n2, n3))
      array = initVal
   END subroutine alloc_and_init_real_3D

   subroutine alloc_and_init_complex_1D(array, n1, initVal)
      COMPLEX(CKIND), allocatable, intent(inout) :: array(:)
      integer, intent(IN) :: n1
      COMPLEX(CKIND), intent(IN) :: initVal

      allocate (array(n1))
      array = initVal
   END subroutine alloc_and_init_complex_1D

   subroutine alloc_and_init_complex_2D(array, n1, n2, initVal)
      COMPLEX(CKIND), allocatable, intent(inout) :: array(:, :)
      integer, intent(IN) :: n1, n2
      COMPLEX(CKIND), intent(IN) :: initVal

      allocate (array(n1, n2))
      array = initVal
   END subroutine alloc_and_init_complex_2D

   subroutine alloc_and_init_complex_3D(array, n1, n2, n3, initVal)
      COMPLEX(CKIND), allocatable, intent(inout) :: array(:, :, :)
      integer, intent(IN) :: n1, n2, n3
      COMPLEX(CKIND), intent(IN) :: initVal

      allocate (array(n1, n2, n3))
      array = initVal
   END subroutine alloc_and_init_complex_3D
end module mod_allocationUtils
