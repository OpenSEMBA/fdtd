!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU KERNELS MODULE - CUDA Fortran (CUF) accelerated YEE kernels
!  Fields stay on device between timesteps. Only upload at init,
!  only download when output/probes are needed.
!  Eliminates per-timestep cudaMemcpy overhead (~4-6s saved).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_kernels_m

   use FDETYPES_m
   use Report_m
   use cudafor

   implicit none

   type gpu_state_t
      ! Host-side pointers (same as original)
      real(kind=rkind), pointer, dimension(:,:,:), contiguous :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=rkind), pointer, dimension(:), contiguous :: Idxe, Idye, Idze, Idxh, Idyh, Idzh, dxe, dye, dze, dxh, dyh, dzh
      integer(kind=integersizeofmediamatrices), pointer, dimension(:,:,:), contiguous :: sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz
      real(kind=rkind), pointer, dimension(:), contiguous :: g1, g2, gm1, gm2

      ! Persistent device memory buffers
      real(kind=rkind), pointer, device, dimension(:,:,:) :: Ex_d, Ey_d, Ez_d, Hx_d, Hy_d, Hz_d
      real(kind=rkind), pointer, device, dimension(:) :: Idxe_d, Idye_d, Idze_d, Idxh_d, Idyh_d, Idzh_d, dxe_d, dye_d, dze_d, dxh_d, dyh_d, dzh_d
      integer(kind=integersizeofmediamatrices), pointer, device, dimension(:,:,:) :: sggMiEx_d, sggMiEy_d, sggMiEz_d, sggMiHx_d, sggMiHy_d, sggMiHz_d
      real(kind=rkind), pointer, device, dimension(:) :: g1_d, g2_d, gm1_d, gm2_d

      ! Dimensions for device memory
      integer(kind=4) :: Ex_nx, Ex_ny, Ex_nz, Ey_nx, Ey_ny, Ey_nz, Ez_nx, Ez_ny, Ez_nz
      integer(kind=4) :: Hx_nx, Hx_ny, Hx_nz, Hy_nx, Hy_ny, Hy_nz, Hz_nx, Hz_ny, Hz_nz
      integer(kind=4) :: Idxe_n, Idye_n, Idze_n, Idxh_n, Idyh_n, Idzh_n
      integer(kind=4) :: dxe_n, dye_n, dze_n, dxh_n, dyh_n, dzh_n
      integer(kind=4) :: g1_n, g2_n, gm1_n, gm2_n
      integer(kind=4) :: sggMi_nx, sggMi_ny, sggMi_nz

      logical :: initialized = .false.
      logical :: fields_on_device = .false.
   end type

contains

   !--------------------------------------------------------------------------------
   ! Initialize GPU data regions - called once at simulation start
   !--------------------------------------------------------------------------------
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

      ! Store host pointers
      this%Ex => Ex; this%Ey => Ey; this%Ez => Ez
      this%Hx => Hx; this%Hy => Hy; this%Hz => Hz
      this%sggMiEx => sggMiEx; this%sggMiEy => sggMiEy; this%sggMiEz => sggMiEz
      this%sggMiHx => sggMiHx; this%sggMiHy => sggMiHy; this%sggMiHz => sggMiHz
      this%g1 => g1; this%g2 => g2; this%gm1 => gm1; this%gm2 => gm2
      this%Idxe => Idxe_in; this%Idye => Idye_in; this%Idze => Idze_in
      this%Idxh => Idxh_in; this%Idyh => Idyh_in; this%Idzh => Idzh_in
      this%dxe => dxe_in; this%dye => dye_in; this%dze => dze_in
      this%dxh => dxh_in; this%dyh => dyh_in; this%dzh => dzh_in

      ! Store dimensions
      this%Ex_nx = ubound(Ex,1) - lbound(Ex,1) + 1
      this%Ex_ny = ubound(Ex,2) - lbound(Ex,2) + 1
      this%Ex_nz = ubound(Ex,3) - lbound(Ex,3) + 1
      this%Ey_nx = ubound(Ey,1) - lbound(Ey,1) + 1
      this%Ey_ny = ubound(Ey,2) - lbound(Ey,2) + 1
      this%Ey_nz = ubound(Ey,3) - lbound(Ey,3) + 1
      this%Ez_nx = ubound(Ez,1) - lbound(Ez,1) + 1
      this%Ez_ny = ubound(Ez,2) - lbound(Ez,2) + 1
      this%Ez_nz = ubound(Ez,3) - lbound(Ez,3) + 1
      this%Hx_nx = ubound(Hx,1) - lbound(Hx,1) + 1
      this%Hx_ny = ubound(Hx,2) - lbound(Hx,2) + 1
      this%Hx_nz = ubound(Hx,3) - lbound(Hx,3) + 1
      this%Hy_nx = ubound(Hy,1) - lbound(Hy,1) + 1
      this%Hy_ny = ubound(Hy,2) - lbound(Hy,2) + 1
      this%Hy_nz = ubound(Hy,3) - lbound(Hy,3) + 1
      this%Hz_nx = ubound(Hz,1) - lbound(Hz,1) + 1
      this%Hz_ny = ubound(Hz,2) - lbound(Hz,2) + 1
      this%Hz_nz = ubound(Hz,3) - lbound(Hz,3) + 1

      this%Idxe_n = ubound(Idxe_in,1) - lbound(Idxe_in,1) + 1
      this%Idye_n = ubound(Idye_in,1) - lbound(Idye_in,1) + 1
      this%Idze_n = ubound(Idze_in,1) - lbound(Idze_in,1) + 1
      this%Idxh_n = ubound(Idxh_in,1) - lbound(Idxh_in,1) + 1
      this%Idyh_n = ubound(Idyh_in,1) - lbound(Idyh_in,1) + 1
      this%Idzh_n = ubound(Idzh_in,1) - lbound(Idzh_in,1) + 1
      this%dxe_n = ubound(dxe_in,1) - lbound(dxe_in,1) + 1
      this%dye_n = ubound(dye_in,1) - lbound(dye_in,1) + 1
      this%dze_n = ubound(dze_in,1) - lbound(dze_in,1) + 1
      this%dxh_n = ubound(dxh_in,1) - lbound(dxh_in,1) + 1
      this%dyh_n = ubound(dyh_in,1) - lbound(dyh_in,1) + 1
      this%dzh_n = ubound(dzh_in,1) - lbound(dzh_in,1) + 1
      this%g1_n = ubound(g1,1) - lbound(g1,1) + 1
      this%g2_n = ubound(g2,1) - lbound(g2,1) + 1
      this%gm1_n = ubound(gm1,1) - lbound(gm1,1) + 1
      this%gm2_n = ubound(gm2,1) - lbound(gm2,1) + 1

      this%sggMi_nx = ubound(sggMiEx,1) - lbound(sggMiEx,1) + 1
      this%sggMi_ny = ubound(sggMiEx,2) - lbound(sggMiEx,2) + 1
      this%sggMi_nz = ubound(sggMiEx,3) - lbound(sggMiEx,3) + 1

      ! Allocate persistent device memory
      allocate(this%Ex_d(this%Ex_nx, this%Ex_ny, this%Ex_nz))
      allocate(this%Ey_d(this%Ey_nx, this%Ey_ny, this%Ey_nz))
      allocate(this%Ez_d(this%Ez_nx, this%Ez_ny, this%Ez_nz))
      allocate(this%Hx_d(this%Hx_nx, this%Hx_ny, this%Hx_nz))
      allocate(this%Hy_d(this%Hy_nx, this%Hy_ny, this%Hy_nz))
      allocate(this%Hz_d(this%Hz_nx, this%Hz_ny, this%Hz_nz))

      allocate(this%Idxe_d(this%Idxe_n))
      allocate(this%Idye_d(this%Idye_n))
      allocate(this%Idze_d(this%Idze_n))
      allocate(this%Idxh_d(this%Idxh_n))
      allocate(this%Idyh_d(this%Idyh_n))
      allocate(this%Idzh_d(this%Idzh_n))
      allocate(this%dxe_d(this%dxe_n))
      allocate(this%dye_d(this%dye_n))
      allocate(this%dze_d(this%dze_n))
      allocate(this%dxh_d(this%dxh_n))
      allocate(this%dyh_d(this%dyh_n))
      allocate(this%dzh_d(this%dzh_n))

      allocate(this%sggMiEx_d(this%sggMi_nx, this%sggMi_ny, this%sggMi_nz))
      allocate(this%sggMiEy_d(this%sggMi_nx, this%sggMi_ny, this%sggMi_nz))
      allocate(this%sggMiEz_d(this%sggMi_nx, this%sggMi_ny, this%sggMi_nz))
      allocate(this%sggMiHx_d(this%sggMi_nx, this%sggMi_ny, this%sggMi_nz))
      allocate(this%sggMiHy_d(this%sggMi_nx, this%sggMi_ny, this%sggMi_nz))
      allocate(this%sggMiHz_d(this%sggMi_nx, this%sggMi_ny, this%sggMi_nz))

      allocate(this%g1_d(this%g1_n))
      allocate(this%g2_d(this%g2_n))
      allocate(this%gm1_d(this%gm1_n))
      allocate(this%gm2_d(this%gm2_n))

      ! Initial upload: host -> device (once at startup)
      this%Ex_d = this%Ex
      this%Ey_d = this%Ey
      this%Ez_d = this%Ez
      this%Hx_d = this%Hx
      this%Hy_d = this%Hy
      this%Hz_d = this%Hz
      this%fields_on_device = .true.

      this%initialized = .true.

    end subroutine gpu_init

    subroutine gpu_destroy(this)
       class(gpu_state_t), intent(inout) :: this

       ! Deallocate device memory
       if (associated(this%Ex_d)) deallocate(this%Ex_d)
       if (associated(this%Ey_d)) deallocate(this%Ey_d)
       if (associated(this%Ez_d)) deallocate(this%Ez_d)
       if (associated(this%Hx_d)) deallocate(this%Hx_d)
       if (associated(this%Hy_d)) deallocate(this%Hy_d)
       if (associated(this%Hz_d)) deallocate(this%Hz_d)

       if (associated(this%Idxe_d)) deallocate(this%Idxe_d)
       if (associated(this%Idye_d)) deallocate(this%Idye_d)
       if (associated(this%Idze_d)) deallocate(this%Idze_d)
       if (associated(this%Idxh_d)) deallocate(this%Idxh_d)
       if (associated(this%Idyh_d)) deallocate(this%Idyh_d)
       if (associated(this%Idzh_d)) deallocate(this%Idzh_d)
       if (associated(this%dxe_d)) deallocate(this%dxe_d)
       if (associated(this%dye_d)) deallocate(this%dye_d)
       if (associated(this%dze_d)) deallocate(this%dze_d)
       if (associated(this%dxh_d)) deallocate(this%dxh_d)
       if (associated(this%dyh_d)) deallocate(this%dyh_d)
       if (associated(this%dzh_d)) deallocate(this%dzh_d)

       if (associated(this%sggMiEx_d)) deallocate(this%sggMiEx_d)
       if (associated(this%sggMiEy_d)) deallocate(this%sggMiEy_d)
       if (associated(this%sggMiEz_d)) deallocate(this%sggMiEz_d)
       if (associated(this%sggMiHx_d)) deallocate(this%sggMiHx_d)
       if (associated(this%sggMiHy_d)) deallocate(this%sggMiHy_d)
       if (associated(this%sggMiHz_d)) deallocate(this%sggMiHz_d)

       if (associated(this%g1_d)) deallocate(this%g1_d)
       if (associated(this%g2_d)) deallocate(this%g2_d)
       if (associated(this%gm1_d)) deallocate(this%gm1_d)
       if (associated(this%gm2_d)) deallocate(this%gm2_d)

       ! Nullify host pointers
       nullify(this%Ex); nullify(this%Ey); nullify(this%Ez)
       nullify(this%Hx); nullify(this%Hy); nullify(this%Hz)
       nullify(this%sggMiEx); nullify(this%sggMiEy); nullify(this%sggMiEz)
       nullify(this%sggMiHx); nullify(this%sggMiHy); nullify(this%sggMiHz)
       nullify(this%g1); nullify(this%g2); nullify(this%gm1); nullify(this%gm2)
       nullify(this%Idxe); nullify(this%Idye); nullify(this%Idze)
       nullify(this%Idxh); nullify(this%Idyh); nullify(this%Idzh)
       nullify(this%dxe); nullify(this%dye); nullify(this%dze)
       nullify(this%dxh); nullify(this%dyh); nullify(this%dzh)

       this%initialized = .false.
       this%fields_on_device = .false.

    end subroutine gpu_destroy

    !--------------------------------------------------------------------------------
    ! Download device data to host - called only when output/probes are needed
    ! Fields stay on device between timesteps to avoid cudaMemcpy overhead
    !--------------------------------------------------------------------------------
    subroutine gpu_download(this)
       class(gpu_state_t), intent(inout) :: this

       if (.not. this%initialized) return

       this%Ex = this%Ex_d
       this%Ey = this%Ey_d
       this%Ez = this%Ez_d
       this%Hx = this%Hx_d
       this%Hy = this%Hy_d
       this%Hz = this%Hz_d
       this%fields_on_device = .false.

    end subroutine gpu_download

    !--------------------------------------------------------------------------------
    ! Upload host data to device - called only when fields are modified on host
    !--------------------------------------------------------------------------------
    subroutine gpu_upload(this)
       class(gpu_state_t), intent(inout) :: this

       if (.not. this%initialized) return

       this%Ex_d = this%Ex
       this%Ey_d = this%Ey
       this%Ez_d = this%Ez
       this%Hx_d = this%Hx
       this%Hy_d = this%Hy
       this%Hz_d = this%Hz
       this%fields_on_device = .true.

    end subroutine gpu_upload

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

end module gpu_kernels_m
