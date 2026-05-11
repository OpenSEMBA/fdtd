!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU KERNELS MODULE - CUDA Fortran (CUF) accelerated YEE kernels
!  Implements the same API as gpu_kernels_m used by timestepping.F90.
!  Current implementation uses local device allocations per kernel call.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_kernels_m

   use FDETYPES_m
   use Report_m
   use cudafor

   implicit none

   type gpu_state_t
      real(kind=rkind), pointer, dimension(:,:,:), contiguous :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=rkind), pointer, dimension(:), contiguous :: Idxe, Idye, Idze, Idxh, Idyh, Idzh, dxe, dye, dze, dxh, dyh, dzh
      integer(kind=integersizeofmediamatrices), pointer, dimension(:,:,:), contiguous :: sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz
      real(kind=rkind), pointer, dimension(:), contiguous :: g1, g2, gm1, gm2

      logical :: initialized = .false.
   end type

contains

   subroutine gpu_init(this, Ex, Ey, Ez, Hx, Hy, Hz, sggMiEx, sggMiEy, sggMiEz, &
                       sggMiHx, sggMiHy, sggMiHz, g1, g2, gm1, gm2, &
                       Idxe_in, Idye_in, Idze_in, Idxh_in, Idyh_in, Idzh_in, dxe_in, dye_in, dze_in, dxh_in, dyh_in, dzh_in)
      class(gpu_state_t), intent(inout) :: this
      real(kind=rkind), dimension(:,:,:), pointer, contiguous, intent(in) :: Ex, Ey, Ez, Hx, Hy, Hz
      integer(kind=integersizeofmediamatrices), dimension(:,:,:), pointer, contiguous, intent(in) :: sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz
      real(kind=rkind), dimension(:), pointer, intent(in) :: g1, g2, gm1, gm2
      real(kind=rkind), dimension(:), pointer, intent(in) :: Idxe_in, Idye_in, Idze_in, Idxh_in, Idyh_in, Idzh_in, dxe_in, dye_in, dze_in, dxh_in, dyh_in, dzh_in

      integer(kind=4) :: ndev
      integer(kind=4) :: cuda_status
      integer(kind=4) :: env_status
      character(len=16) :: enable_cuf

      call get_environment_variable("SEMBA_FDTD_ENABLE_CUF_RUNTIME", enable_cuf, status=env_status)
      if (env_status /= 0 .or. trim(enable_cuf) /= "1") then
         this%initialized = .false.
         return
      endif

      cuda_status = cudaGetDeviceCount(ndev)
      if (cuda_status /= cudaSuccess .or. ndev <= 0) then
         this%initialized = .false.
         return
      endif

      this%Ex => Ex
      this%Ey => Ey
      this%Ez => Ez
      this%Hx => Hx
      this%Hy => Hy
      this%Hz => Hz

      this%sggMiEx => sggMiEx
      this%sggMiEy => sggMiEy
      this%sggMiEz => sggMiEz
      this%sggMiHx => sggMiHx
      this%sggMiHy => sggMiHy
      this%sggMiHz => sggMiHz

      this%g1 => g1
      this%g2 => g2
      this%gm1 => gm1
      this%gm2 => gm2

      this%Idxe => Idxe_in
      this%Idye => Idye_in
      this%Idze => Idze_in
      this%Idxh => Idxh_in
      this%Idyh => Idyh_in
      this%Idzh => Idzh_in
      this%dxe => dxe_in
      this%dye => dye_in
      this%dze => dze_in
      this%dxh => dxh_in
      this%dyh => dyh_in
      this%dzh => dzh_in

      this%initialized = .true.

   end subroutine gpu_init

   subroutine gpu_destroy(this)
      class(gpu_state_t), intent(inout) :: this

      nullify(this%Ex)
      nullify(this%Ey)
      nullify(this%Ez)
      nullify(this%Hx)
      nullify(this%Hy)
      nullify(this%Hz)
      nullify(this%sggMiEx)
      nullify(this%sggMiEy)
      nullify(this%sggMiEz)
      nullify(this%sggMiHx)
      nullify(this%sggMiHy)
      nullify(this%sggMiHz)
      nullify(this%g1)
      nullify(this%g2)
      nullify(this%gm1)
      nullify(this%gm2)
      nullify(this%Idxe)
      nullify(this%Idye)
      nullify(this%Idze)
      nullify(this%Idxh)
      nullify(this%Idyh)
      nullify(this%Idzh)
      nullify(this%dxe)
      nullify(this%dye)
      nullify(this%dze)
      nullify(this%dxh)
      nullify(this%dyh)
      nullify(this%dzh)

      this%initialized = .false.

   end subroutine gpu_destroy

   subroutine gpu_advanceEx(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds

      real(kind=rkind), allocatable, device, dimension(:,:,:) :: Ex_d, Hy_d, Hz_d
      integer(kind=integersizeofmediamatrices), allocatable, device, dimension(:,:,:) :: sggMiEx_d
      real(kind=rkind), allocatable, device, dimension(:) :: Idyh_d, Idzh_d, g1_d, g2_d

      allocate(Ex_d(lbound(this%Ex,1):ubound(this%Ex,1), lbound(this%Ex,2):ubound(this%Ex,2), lbound(this%Ex,3):ubound(this%Ex,3)))
      allocate(Hy_d(lbound(this%Hy,1):ubound(this%Hy,1), lbound(this%Hy,2):ubound(this%Hy,2), lbound(this%Hy,3):ubound(this%Hy,3)))
      allocate(Hz_d(lbound(this%Hz,1):ubound(this%Hz,1), lbound(this%Hz,2):ubound(this%Hz,2), lbound(this%Hz,3):ubound(this%Hz,3)))
      allocate(sggMiEx_d(lbound(this%sggMiEx,1):ubound(this%sggMiEx,1), lbound(this%sggMiEx,2):ubound(this%sggMiEx,2), lbound(this%sggMiEx,3):ubound(this%sggMiEx,3)))
      allocate(Idyh_d(lbound(this%Idyh,1):ubound(this%Idyh,1)))
      allocate(Idzh_d(lbound(this%Idzh,1):ubound(this%Idzh,1)))
      allocate(g1_d(lbound(this%g1,1):ubound(this%g1,1)))
      allocate(g2_d(lbound(this%g2,1):ubound(this%g2,1)))

      Ex_d = this%Ex
      Hy_d = this%Hy
      Hz_d = this%Hz
      sggMiEx_d = this%sggMiEx
      Idyh_d = this%Idyh
      Idzh_d = this%Idzh
      g1_d = this%g1
      g2_d = this%g2

      call gpu_advanceEx_kernel(Ex_d, Hy_d, Hz_d, sggMiEx_d, &
                                Idyh_d, Idzh_d, g1_d, g2_d, &
                                bounds%sweepEx%NX, bounds%sweepEx%NY, bounds%sweepEx%NZ)

      this%Ex = Ex_d

      deallocate(Ex_d, Hy_d, Hz_d, sggMiEx_d, Idyh_d, Idzh_d, g1_d, g2_d)

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

   subroutine gpu_advanceEy(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds
      real(kind=rkind), allocatable, device, dimension(:,:,:) :: Ey_d, Hz_d, Hx_d
      integer(kind=integersizeofmediamatrices), allocatable, device, dimension(:,:,:) :: sggMiEy_d
      real(kind=rkind), allocatable, device, dimension(:) :: Idzh_d, Idxh_d, g1_d, g2_d

      allocate(Ey_d(lbound(this%Ey,1):ubound(this%Ey,1), lbound(this%Ey,2):ubound(this%Ey,2), lbound(this%Ey,3):ubound(this%Ey,3)))
      allocate(Hz_d(lbound(this%Hz,1):ubound(this%Hz,1), lbound(this%Hz,2):ubound(this%Hz,2), lbound(this%Hz,3):ubound(this%Hz,3)))
      allocate(Hx_d(lbound(this%Hx,1):ubound(this%Hx,1), lbound(this%Hx,2):ubound(this%Hx,2), lbound(this%Hx,3):ubound(this%Hx,3)))
      allocate(sggMiEy_d(lbound(this%sggMiEy,1):ubound(this%sggMiEy,1), lbound(this%sggMiEy,2):ubound(this%sggMiEy,2), lbound(this%sggMiEy,3):ubound(this%sggMiEy,3)))
      allocate(Idzh_d(lbound(this%Idzh,1):ubound(this%Idzh,1)))
      allocate(Idxh_d(lbound(this%Idxh,1):ubound(this%Idxh,1)))
      allocate(g1_d(lbound(this%g1,1):ubound(this%g1,1)))
      allocate(g2_d(lbound(this%g2,1):ubound(this%g2,1)))

      Ey_d = this%Ey
      Hz_d = this%Hz
      Hx_d = this%Hx
      sggMiEy_d = this%sggMiEy
      Idzh_d = this%Idzh
      Idxh_d = this%Idxh
      g1_d = this%g1
      g2_d = this%g2

      call gpu_advanceEy_kernel(Ey_d, Hz_d, Hx_d, sggMiEy_d, Idzh_d, Idxh_d, g1_d, g2_d, &
                                bounds%sweepEy%NX, bounds%sweepEy%NY, bounds%sweepEy%NZ)

      this%Ey = Ey_d

      deallocate(Ey_d, Hz_d, Hx_d, sggMiEy_d, Idzh_d, Idxh_d, g1_d, g2_d)
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

   subroutine gpu_advanceEz(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds
      real(kind=rkind), allocatable, device, dimension(:,:,:) :: Ez_d, Hx_d, Hy_d
      integer(kind=integersizeofmediamatrices), allocatable, device, dimension(:,:,:) :: sggMiEz_d
      real(kind=rkind), allocatable, device, dimension(:) :: Idyh_d, Idxh_d, g1_d, g2_d

      allocate(Ez_d(lbound(this%Ez,1):ubound(this%Ez,1), lbound(this%Ez,2):ubound(this%Ez,2), lbound(this%Ez,3):ubound(this%Ez,3)))
      allocate(Hx_d(lbound(this%Hx,1):ubound(this%Hx,1), lbound(this%Hx,2):ubound(this%Hx,2), lbound(this%Hx,3):ubound(this%Hx,3)))
      allocate(Hy_d(lbound(this%Hy,1):ubound(this%Hy,1), lbound(this%Hy,2):ubound(this%Hy,2), lbound(this%Hy,3):ubound(this%Hy,3)))
      allocate(sggMiEz_d(lbound(this%sggMiEz,1):ubound(this%sggMiEz,1), lbound(this%sggMiEz,2):ubound(this%sggMiEz,2), lbound(this%sggMiEz,3):ubound(this%sggMiEz,3)))
      allocate(Idyh_d(lbound(this%Idyh,1):ubound(this%Idyh,1)))
      allocate(Idxh_d(lbound(this%Idxh,1):ubound(this%Idxh,1)))
      allocate(g1_d(lbound(this%g1,1):ubound(this%g1,1)))
      allocate(g2_d(lbound(this%g2,1):ubound(this%g2,1)))

      Ez_d = this%Ez
      Hx_d = this%Hx
      Hy_d = this%Hy
      sggMiEz_d = this%sggMiEz
      Idyh_d = this%Idyh
      Idxh_d = this%Idxh
      g1_d = this%g1
      g2_d = this%g2

      call gpu_advanceEz_kernel(Ez_d, Hx_d, Hy_d, sggMiEz_d, Idyh_d, Idxh_d, g1_d, g2_d, &
                                bounds%sweepEz%NX, bounds%sweepEz%NY, bounds%sweepEz%NZ)

      this%Ez = Ez_d

      deallocate(Ez_d, Hx_d, Hy_d, sggMiEz_d, Idyh_d, Idxh_d, g1_d, g2_d)
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

   subroutine gpu_advanceHx(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds
      real(kind=rkind), allocatable, device, dimension(:,:,:) :: Hx_d, Ey_d, Ez_d
      integer(kind=integersizeofmediamatrices), allocatable, device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), allocatable, device, dimension(:) :: Idye_d, Idze_d, gm1_d, gm2_d

      allocate(Hx_d(lbound(this%Hx,1):ubound(this%Hx,1), lbound(this%Hx,2):ubound(this%Hx,2), lbound(this%Hx,3):ubound(this%Hx,3)))
      allocate(Ey_d(lbound(this%Ey,1):ubound(this%Ey,1), lbound(this%Ey,2):ubound(this%Ey,2), lbound(this%Ey,3):ubound(this%Ey,3)))
      allocate(Ez_d(lbound(this%Ez,1):ubound(this%Ez,1), lbound(this%Ez,2):ubound(this%Ez,2), lbound(this%Ez,3):ubound(this%Ez,3)))
      allocate(sggMiHx_d(lbound(this%sggMiHx,1):ubound(this%sggMiHx,1), lbound(this%sggMiHx,2):ubound(this%sggMiHx,2), lbound(this%sggMiHx,3):ubound(this%sggMiHx,3)))
      allocate(Idye_d(lbound(this%Idye,1):ubound(this%Idye,1)))
      allocate(Idze_d(lbound(this%Idze,1):ubound(this%Idze,1)))
      allocate(gm1_d(lbound(this%gm1,1):ubound(this%gm1,1)))
      allocate(gm2_d(lbound(this%gm2,1):ubound(this%gm2,1)))

      Hx_d = this%Hx
      Ey_d = this%Ey
      Ez_d = this%Ez
      sggMiHx_d = this%sggMiHx
      Idye_d = this%Idye
      Idze_d = this%Idze
      gm1_d = this%gm1
      gm2_d = this%gm2

      call gpu_advanceHx_kernel(Hx_d, Ey_d, Ez_d, sggMiHx_d, Idye_d, Idze_d, gm1_d, gm2_d, &
                                bounds%sweepHx%NX, bounds%sweepHx%NY, bounds%sweepHx%NZ)

      this%Hx = Hx_d

      deallocate(Hx_d, Ey_d, Ez_d, sggMiHx_d, Idye_d, Idze_d, gm1_d, gm2_d)
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

   subroutine gpu_advanceHy(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds
      real(kind=rkind), allocatable, device, dimension(:,:,:) :: Hy_d, Ez_d, Ex_d
      integer(kind=integersizeofmediamatrices), allocatable, device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), allocatable, device, dimension(:) :: Idze_d, Idxe_d, gm1_d, gm2_d

      allocate(Hy_d(lbound(this%Hy,1):ubound(this%Hy,1), lbound(this%Hy,2):ubound(this%Hy,2), lbound(this%Hy,3):ubound(this%Hy,3)))
      allocate(Ez_d(lbound(this%Ez,1):ubound(this%Ez,1), lbound(this%Ez,2):ubound(this%Ez,2), lbound(this%Ez,3):ubound(this%Ez,3)))
      allocate(Ex_d(lbound(this%Ex,1):ubound(this%Ex,1), lbound(this%Ex,2):ubound(this%Ex,2), lbound(this%Ex,3):ubound(this%Ex,3)))
      allocate(sggMiHy_d(lbound(this%sggMiHy,1):ubound(this%sggMiHy,1), lbound(this%sggMiHy,2):ubound(this%sggMiHy,2), lbound(this%sggMiHy,3):ubound(this%sggMiHy,3)))
      allocate(Idze_d(lbound(this%Idze,1):ubound(this%Idze,1)))
      allocate(Idxe_d(lbound(this%Idxe,1):ubound(this%Idxe,1)))
      allocate(gm1_d(lbound(this%gm1,1):ubound(this%gm1,1)))
      allocate(gm2_d(lbound(this%gm2,1):ubound(this%gm2,1)))

      Hy_d = this%Hy
      Ez_d = this%Ez
      Ex_d = this%Ex
      sggMiHy_d = this%sggMiHy
      Idze_d = this%Idze
      Idxe_d = this%Idxe
      gm1_d = this%gm1
      gm2_d = this%gm2

      call gpu_advanceHy_kernel(Hy_d, Ez_d, Ex_d, sggMiHy_d, Idze_d, Idxe_d, gm1_d, gm2_d, &
                                bounds%sweepHy%NX, bounds%sweepHy%NY, bounds%sweepHy%NZ)

      this%Hy = Hy_d

      deallocate(Hy_d, Ez_d, Ex_d, sggMiHy_d, Idze_d, Idxe_d, gm1_d, gm2_d)
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

   subroutine gpu_advanceHz(this, bounds)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: bounds
      real(kind=rkind), allocatable, device, dimension(:,:,:) :: Hz_d, Ex_d, Ey_d
      integer(kind=integersizeofmediamatrices), allocatable, device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), allocatable, device, dimension(:) :: Idye_d, Idxe_d, gm1_d, gm2_d

      allocate(Hz_d(lbound(this%Hz,1):ubound(this%Hz,1), lbound(this%Hz,2):ubound(this%Hz,2), lbound(this%Hz,3):ubound(this%Hz,3)))
      allocate(Ex_d(lbound(this%Ex,1):ubound(this%Ex,1), lbound(this%Ex,2):ubound(this%Ex,2), lbound(this%Ex,3):ubound(this%Ex,3)))
      allocate(Ey_d(lbound(this%Ey,1):ubound(this%Ey,1), lbound(this%Ey,2):ubound(this%Ey,2), lbound(this%Ey,3):ubound(this%Ey,3)))
      allocate(sggMiHz_d(lbound(this%sggMiHz,1):ubound(this%sggMiHz,1), lbound(this%sggMiHz,2):ubound(this%sggMiHz,2), lbound(this%sggMiHz,3):ubound(this%sggMiHz,3)))
      allocate(Idye_d(lbound(this%Idye,1):ubound(this%Idye,1)))
      allocate(Idxe_d(lbound(this%Idxe,1):ubound(this%Idxe,1)))
      allocate(gm1_d(lbound(this%gm1,1):ubound(this%gm1,1)))
      allocate(gm2_d(lbound(this%gm2,1):ubound(this%gm2,1)))

      Hz_d = this%Hz
      Ex_d = this%Ex
      Ey_d = this%Ey
      sggMiHz_d = this%sggMiHz
      Idye_d = this%Idye
      Idxe_d = this%Idxe
      gm1_d = this%gm1
      gm2_d = this%gm2

      call gpu_advanceHz_kernel(Hz_d, Ex_d, Ey_d, sggMiHz_d, Idye_d, Idxe_d, gm1_d, gm2_d, &
                                bounds%sweepHz%NX, bounds%sweepHz%NY, bounds%sweepHz%NZ)

      this%Hz = Hz_d

      deallocate(Hz_d, Ex_d, Ey_d, sggMiHz_d, Idye_d, Idxe_d, gm1_d, gm2_d)
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

end module gpu_kernels_m
