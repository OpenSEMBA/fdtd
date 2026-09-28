module allocationUtils_m
   use FDETYPES_m, only: RKIND, CKIND, SINGLE, RKIND_TIME, IKINDMTAG, INTEGERSIZEOFMEDIAMATRICES
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
#ifndef CompileWithReal8
      procedure alloc_and_init_real_time_1D
#endif
   end interface
contains
#ifndef CompileWithReal8
   subroutine alloc_and_init_real_time_1D(array, n1, initVal)
      real(RKIND_TIME), allocatable, intent(inout) :: array(:)
      integer, intent(in) :: n1
      real(RKIND_TIME), intent(in) :: initVal

      allocate (array(n1))
      array = initVal
   end subroutine alloc_and_init_real_time_1D
#endif
   subroutine alloc_and_init_int_1D(array, n1, initVal)
      integer(SINGLE), allocatable, intent(inout) :: array(:)
      integer, intent(in) :: n1
      integer(SINGLE), intent(in) :: initVal

      allocate (array(n1))
      array = initVal
   end subroutine alloc_and_init_int_1D

   subroutine alloc_and_init_int_2D(array, n1, n2, initVal)
      integer(SINGLE), allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: n1, n2
      integer(SINGLE), intent(in) :: initVal

      allocate (array(n1, n2))
      array = initVal
   end subroutine alloc_and_init_int_2D

   subroutine alloc_and_init_int_3D(array, n1, n2, n3, initVal)
      integer(SINGLE), allocatable, intent(inout) :: array(:, :, :)
      integer, intent(in) :: n1, n2, n3
      integer(SINGLE), intent(in) :: initVal

      allocate (array(n1, n2, n3))
      array = initVal
   end subroutine alloc_and_init_int_3D

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
      real(RKIND), allocatable, intent(inout) :: array(:)
      integer, intent(in) :: n1
      real(RKIND), intent(in) :: initVal

      allocate (array(n1))
      array = initVal
   end subroutine alloc_and_init_real_1D

   subroutine alloc_and_init_real_2D(array, n1, n2, initVal)
      real(RKIND), allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: n1, n2
      real(RKIND), intent(in) :: initVal

      allocate (array(n1, n2))
      array = initVal
   end subroutine alloc_and_init_real_2D

   subroutine alloc_and_init_real_3D(array, n1, n2, n3, initVal)
      real(RKIND), allocatable, intent(inout) :: array(:, :, :)
      integer, intent(in) :: n1, n2, n3
      real(RKIND), intent(in) :: initVal

      allocate (array(n1, n2, n3))
      array = initVal
   end subroutine alloc_and_init_real_3D

   subroutine alloc_and_init_complex_1D(array, n1, initVal)
      complex(CKIND), allocatable, intent(inout) :: array(:)
      integer, intent(in) :: n1
      complex(CKIND), intent(in) :: initVal

      allocate (array(n1))
      array = initVal
   end subroutine alloc_and_init_complex_1D

   subroutine alloc_and_init_complex_2D(array, n1, n2, initVal)
      complex(CKIND), allocatable, intent(inout) :: array(:, :)
      integer, intent(in) :: n1, n2
      complex(CKIND), intent(in) :: initVal

      allocate (array(n1, n2))
      array = initVal
   end subroutine alloc_and_init_complex_2D

   subroutine alloc_and_init_complex_3D(array, n1, n2, n3, initVal)
      complex(CKIND), allocatable, intent(inout) :: array(:, :, :)
      integer, intent(in) :: n1, n2, n3
      complex(CKIND), intent(in) :: initVal

      allocate (array(n1, n2, n3))
      array = initVal
   end subroutine alloc_and_init_complex_3D
end module allocationUtils_m
