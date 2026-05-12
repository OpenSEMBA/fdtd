!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU YEE KERNELS MODULE - CUDA Fortran (CUF) accelerated YEE kernels
!  Fields stay on device between timesteps. Only kernel wrappers here.
!  Split from gpu_kernels_cuf.F90 to avoid NVHPC compiler file-size limit.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_yee_m

   use FDETYPES_m
   use Report_m
   use cudafor
   use gpu_core_m

   implicit none

contains

   !--------------------------------------------------------------------------------
   ! Advance Ex field - GPU accelerated YEE kernel (fields on device, no copy)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceEx(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      if (.not. this%initialized) return

      call gpu_advanceEx_kernel(this%Ex_d, this%Hy_d, this%Hz_d, this%sggMiEx_d, &
                                this%Idyh_d, this%Idzh_d, this%g1_d, this%g2_d, &
                                bounds%sweepEx%NX, bounds%sweepEx%NY, bounds%sweepEx%NZ)

   end subroutine gpu_advanceEx

   subroutine gpu_advanceEx_kernel(Ex_d, Hy_d, Hz_d, sggMiEx_d, Idyh_d, Idzh_d, g1_d, g2_d, nx, ny, nz)
      integer(kind=4), intent(in) :: nx, ny, nz
      real(kind=rkind), device, dimension(:,:,:) :: Ex_d, Hy_d, Hz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEx_d
      real(kind=rkind), device, dimension(:) :: Idyh_d, Idzh_d, g1_d, g2_d

      integer(kind=4) :: i, j, k
      real(kind=rkind) :: Idzhk, Idyhj
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=1,nz
         do j=1,ny
            do i=1,nx
               Idzhk = Idzh_d(k)
               Idyhj = Idyh_d(j)
               medio = sggMiEx_d(i,j,k)
               Ex_d(i,j,k) = g1_d(medio)*Ex_d(i,j,k) + g2_d(medio)* &
                  ((Hz_d(i,j,k) - Hz_d(i,j-1,k))*Idyhj - (Hy_d(i,j,k) - Hy_d(i,j,k-1))*Idzhk)
            end do
         end do
      end do
   end subroutine gpu_advanceEx_kernel

   !--------------------------------------------------------------------------------
   ! Advance Ey field - GPU accelerated YEE kernel (fields on device, no copy)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceEy(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      if (.not. this%initialized) return

      call gpu_advanceEy_kernel(this%Ey_d, this%Hz_d, this%Hx_d, this%sggMiEy_d, &
                                this%Idzh_d, this%Idxh_d, this%g1_d, this%g2_d, &
                                bounds%sweepEy%NX, bounds%sweepEy%NY, bounds%sweepEy%NZ)

   end subroutine gpu_advanceEy

   subroutine gpu_advanceEy_kernel(Ey_d, Hz_d, Hx_d, sggMiEy_d, Idzh_d, Idxh_d, g1_d, g2_d, nx, ny, nz)
      integer(kind=4), intent(in) :: nx, ny, nz
      real(kind=rkind), device, dimension(:,:,:) :: Ey_d, Hz_d, Hx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEy_d
      real(kind=rkind), device, dimension(:) :: Idzh_d, Idxh_d, g1_d, g2_d

      integer(kind=4) :: i, j, k
      real(kind=rkind) :: Idzhk
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=1,nz
         do j=1,ny
            do i=1,nx
               Idzhk = Idzh_d(k)
               medio = sggMiEy_d(i,j,k)
               Ey_d(i,j,k) = g1_d(medio)*Ey_d(i,j,k) + g2_d(medio)* &
                  ((Hx_d(i,j,k)-Hx_d(i,j,k-1))*Idzhk - (Hz_d(i,j,k)-Hz_d(i-1,j,k))*Idxh_d(i))
            end do
         end do
      end do
   end subroutine gpu_advanceEy_kernel

   !--------------------------------------------------------------------------------
   ! Advance Ez field - GPU accelerated YEE kernel (fields on device, no copy)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceEz(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      if (.not. this%initialized) return

      call gpu_advanceEz_kernel(this%Ez_d, this%Hx_d, this%Hy_d, this%sggMiEz_d, &
                                this%Idyh_d, this%Idxh_d, this%g1_d, this%g2_d, &
                                bounds%sweepEz%NX, bounds%sweepEz%NY, bounds%sweepEz%NZ)

   end subroutine gpu_advanceEz

   subroutine gpu_advanceEz_kernel(Ez_d, Hx_d, Hy_d, sggMiEz_d, Idyh_d, Idxh_d, g1_d, g2_d, nx, ny, nz)
      integer(kind=4), intent(in) :: nx, ny, nz
      real(kind=rkind), device, dimension(:,:,:) :: Ez_d, Hx_d, Hy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEz_d
      real(kind=rkind), device, dimension(:) :: Idyh_d, Idxh_d, g1_d, g2_d

      integer(kind=4) :: i, j, k
      real(kind=rkind) :: Idyhj
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=1,nz
         do j=1,ny
            do i=1,nx
               Idyhj = Idyh_d(j)
               medio = sggMiEz_d(i,j,k)
               Ez_d(i,j,k) = g1_d(medio)*Ez_d(i,j,k) + g2_d(medio)* &
                  ((Hy_d(i,j,k)-Hy_d(i-1,j,k))*Idxh_d(i) - (Hx_d(i,j,k)-Hx_d(i,j-1,k))*Idyhj)
            end do
         end do
      end do
   end subroutine gpu_advanceEz_kernel

   !--------------------------------------------------------------------------------
   ! Advance Hx field - GPU accelerated YEE kernel (fields on device, no copy)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceHx(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      if (.not. this%initialized) return

      call gpu_advanceHx_kernel(this%Hx_d, this%Ey_d, this%Ez_d, this%sggMiHx_d, &
                                this%Idye_d, this%Idze_d, this%gm1_d, this%gm2_d, &
                                bounds%sweepHx%NX, bounds%sweepHx%NY, bounds%sweepHx%NZ)

   end subroutine gpu_advanceHx

   subroutine gpu_advanceHx_kernel(Hx_d, Ey_d, Ez_d, sggMiHx_d, Idye_d, Idze_d, gm1_d, gm2_d, nx, ny, nz)
      integer(kind=4), intent(in) :: nx, ny, nz
      real(kind=rkind), device, dimension(:,:,:) :: Hx_d, Ey_d, Ez_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), device, dimension(:) :: Idye_d, Idze_d, gm1_d, gm2_d

      integer(kind=4) :: i, j, k
      real(kind=rkind) :: Idzek, Idyej
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=1,nz
         do j=1,ny
            do i=1,nx
               Idzek = Idze_d(k)
               Idyej = Idye_d(j)
               medio = sggMiHx_d(i,j,k)
               Hx_d(i,j,k) = gm1_d(medio)*Hx_d(i,j,k) + gm2_d(medio)* &
                  ((Ey_d(i,j,k+1)-Ey_d(i,j,k))*Idzek - (Ez_d(i,j+1,k)-Ez_d(i,j,k))*Idyej)
            end do
         end do
      end do
   end subroutine gpu_advanceHx_kernel

   !--------------------------------------------------------------------------------
   ! Advance Hy field - GPU accelerated YEE kernel (fields on device, no copy)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceHy(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      if (.not. this%initialized) return

      call gpu_advanceHy_kernel(this%Hy_d, this%Ez_d, this%Ex_d, this%sggMiHy_d, &
                                this%Idze_d, this%Idxe_d, this%gm1_d, this%gm2_d, &
                                bounds%sweepHy%NX, bounds%sweepHy%NY, bounds%sweepHy%NZ)

   end subroutine gpu_advanceHy

   subroutine gpu_advanceHy_kernel(Hy_d, Ez_d, Ex_d, sggMiHy_d, Idze_d, Idxe_d, gm1_d, gm2_d, nx, ny, nz)
      integer(kind=4), intent(in) :: nx, ny, nz
      real(kind=rkind), device, dimension(:,:,:) :: Hy_d, Ez_d, Ex_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), device, dimension(:) :: Idze_d, Idxe_d, gm1_d, gm2_d

      integer(kind=4) :: i, j, k
      real(kind=rkind) :: Idzek
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=1,nz
         do j=1,ny
            do i=1,nx
               Idzek = Idze_d(k)
               medio = sggMiHy_d(i,j,k)
               Hy_d(i,j,k) = gm1_d(medio)*Hy_d(i,j,k) + gm2_d(medio)* &
                  ((Ez_d(i+1,j,k)-Ez_d(i,j,k))*Idxe_d(i) - (Ex_d(i,j,k+1)-Ex_d(i,j,k))*Idzek)
            end do
         end do
      end do
   end subroutine gpu_advanceHy_kernel

   !--------------------------------------------------------------------------------
   ! Advance Hz field - GPU accelerated YEE kernel (fields on device, no copy)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceHz(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      if (.not. this%initialized) return

      call gpu_advanceHz_kernel(this%Hz_d, this%Ex_d, this%Ey_d, this%sggMiHz_d, &
                                this%Idye_d, this%Idxe_d, this%gm1_d, this%gm2_d, &
                                bounds%sweepHz%NX, bounds%sweepHz%NY, bounds%sweepHz%NZ)

   end subroutine gpu_advanceHz

   subroutine gpu_advanceHz_kernel(Hz_d, Ex_d, Ey_d, sggMiHz_d, Idye_d, Idxe_d, gm1_d, gm2_d, nx, ny, nz)
      integer(kind=4), intent(in) :: nx, ny, nz
      real(kind=rkind), device, dimension(:,:,:) :: Hz_d, Ex_d, Ey_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), device, dimension(:) :: Idye_d, Idxe_d, gm1_d, gm2_d

      integer(kind=4) :: i, j, k
      real(kind=rkind) :: Idyej
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=1,nz
         do j=1,ny
            do i=1,nx
               Idyej = Idye_d(j)
               medio = sggMiHz_d(i,j,k)
               Hz_d(i,j,k) = gm1_d(medio)*Hz_d(i,j,k) + gm2_d(medio)* &
                  ((Ex_d(i,j+1,k)-Ex_d(i,j,k))*Idyej - (Ey_d(i+1,j,k)-Ey_d(i,j,k))*Idxe_d(i))
            end do
         end do
      end do
   end subroutine gpu_advanceHz_kernel

end module gpu_yee_m
