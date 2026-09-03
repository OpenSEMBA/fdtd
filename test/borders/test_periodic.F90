module periodic_border_test_m
   use FDETYPES_m
   use BORDERS_other_m, only: CloneMagneticPeriodic
   implicit none

contains

   subroutine initialize_bounds(alloc, sweep)
      type(XYZlimit_t), dimension(1:6), intent(out) :: alloc, sweep
      integer :: field

      do field = 1, 6
         alloc(field)%XI = 0
         alloc(field)%XE = 3
         alloc(field)%YI = 0
         alloc(field)%YE = 3
         alloc(field)%ZI = 0
         alloc(field)%ZE = 3
         sweep(field)%XI = 1
         sweep(field)%XE = 2
         sweep(field)%YI = 1
         sweep(field)%YE = 2
         sweep(field)%ZI = 1
         sweep(field)%ZE = 2
      end do
   end subroutine initialize_bounds

   subroutine initialize_border(border)
      type(Border_t), intent(out) :: border

      border%IsBackPeriodic = .false.
      border%IsFrontPeriodic = .false.
      border%IsLeftPeriodic = .false.
      border%IsRightPeriodic = .false.
      border%IsDownPeriodic = .false.
      border%IsUpPeriodic = .false.
   end subroutine initialize_border

   subroutine initialize_fields(Hx, Hy, Hz)
      real(kind=RKIND), dimension(0:3,0:3,0:3), intent(out) :: Hx, Hy, Hz
      integer :: i, j, k

      do k = 0, 3
         do j = 0, 3
            do i = 0, 3
               Hx(i,j,k) = real(1000 + 100*i + 10*j + k, RKIND)
               Hy(i,j,k) = real(2000 + 100*i + 10*j + k, RKIND)
               Hz(i,j,k) = real(3000 + 100*i + 10*j + k, RKIND)
            end do
         end do
      end do
   end subroutine initialize_fields

end module periodic_border_test_m

integer function test_periodic_x_ghost_planes() bind(C) result(err)
   use periodic_border_test_m
   implicit none
   type(XYZlimit_t), dimension(1:6) :: alloc, sweep
   type(Border_t) :: border
   real(kind=RKIND), dimension(0:3,0:3,0:3) :: Hx, Hy, Hz, Hx_before

   err = 0
   call initialize_bounds(alloc, sweep)
   call initialize_border(border)
   call initialize_fields(Hx, Hy, Hz)
   Hx_before = Hx
   border%IsBackPeriodic = .true.
   border%IsFrontPeriodic = .true.

   call CloneMagneticPeriodic(alloc, border, Hx, Hy, Hz, sweep, 0, 1)

   if (any(Hy(0,:,:) /= Hy(2,:,:))) err = err + 1
   if (any(Hy(3,:,:) /= Hy(1,:,:))) err = err + 1
   if (any(Hz(0,:,:) /= Hz(2,:,:))) err = err + 1
   if (any(Hz(3,:,:) /= Hz(1,:,:))) err = err + 1
   if (any(Hx /= Hx_before)) err = err + 1
end function test_periodic_x_ghost_planes

integer function test_periodic_y_ghost_planes() bind(C) result(err)
   use periodic_border_test_m
   implicit none
   type(XYZlimit_t), dimension(1:6) :: alloc, sweep
   type(Border_t) :: border
   real(kind=RKIND), dimension(0:3,0:3,0:3) :: Hx, Hy, Hz, Hy_before

   err = 0
   call initialize_bounds(alloc, sweep)
   call initialize_border(border)
   call initialize_fields(Hx, Hy, Hz)
   Hy_before = Hy
   border%IsLeftPeriodic = .true.
   border%IsRightPeriodic = .true.

   call CloneMagneticPeriodic(alloc, border, Hx, Hy, Hz, sweep, 0, 1)

   if (any(Hx(:,0,:) /= Hx(:,2,:))) err = err + 1
   if (any(Hx(:,3,:) /= Hx(:,1,:))) err = err + 1
   if (any(Hz(:,0,:) /= Hz(:,2,:))) err = err + 1
   if (any(Hz(:,3,:) /= Hz(:,1,:))) err = err + 1
   if (any(Hy /= Hy_before)) err = err + 1
end function test_periodic_y_ghost_planes

integer function test_periodic_z_ghost_planes() bind(C) result(err)
   use periodic_border_test_m
   implicit none
   type(XYZlimit_t), dimension(1:6) :: alloc, sweep
   type(Border_t) :: border
   real(kind=RKIND), dimension(0:3,0:3,0:3) :: Hx, Hy, Hz, Hz_before

   err = 0
   call initialize_bounds(alloc, sweep)
   call initialize_border(border)
   call initialize_fields(Hx, Hy, Hz)
   Hz_before = Hz
   border%IsDownPeriodic = .true.
   border%IsUpPeriodic = .true.

   ! Rank zero owns only the lower periodic ghost plane.
   call CloneMagneticPeriodic(alloc, border, Hx, Hy, Hz, sweep, 0, 2)
   if (any(Hx(:,:,0) /= Hx(:,:,2))) err = err + 1
   if (any(Hy(:,:,0) /= Hy(:,:,2))) err = err + 1
   if (any(Hx(:,:,3) == Hx(:,:,1))) err = err + 1
   if (any(Hy(:,:,3) == Hy(:,:,1))) err = err + 1

   call initialize_fields(Hx, Hy, Hz)
   ! The last rank owns only the upper periodic ghost plane.
   call CloneMagneticPeriodic(alloc, border, Hx, Hy, Hz, sweep, 1, 2)
   if (any(Hx(:,:,3) /= Hx(:,:,1))) err = err + 1
   if (any(Hy(:,:,3) /= Hy(:,:,1))) err = err + 1
   if (any(Hx(:,:,0) == Hx(:,:,2))) err = err + 1
   if (any(Hy(:,:,0) == Hy(:,:,2))) err = err + 1
   if (any(Hz /= Hz_before)) err = err + 1
end function test_periodic_z_ghost_planes
