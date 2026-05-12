!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU CORE MODULE - CUDA Fortran (CUF)
!  gpu_state_t type definition + init/upload/download/destroy
!  Split from gpu_kernels_cuf.F90 to avoid NVHPC compiler file-size limit.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_core_m

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

       ! Persistent device memory buffers - YEE fields
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

       ! CPML left boundary - persistent device arrays
       ! 4 psi arrays for left boundary (same global bounds as host arrays)
       integer(kind=4) :: pml_left_Exy_nx, pml_left_Exy_ny, pml_left_Exy_nz
       integer(kind=4) :: pml_left_Ezy_nx, pml_left_Ezy_ny, pml_left_Ezy_nz
       integer(kind=4) :: pml_left_Hxy_nx, pml_left_Hxy_ny, pml_left_Hxy_nz
       integer(kind=4) :: pml_left_Hzy_nx, pml_left_Hzy_ny, pml_left_Hzy_nz
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Exy_left
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Ezy_left
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Hxy_left
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Hzy_left

       ! CPML left boundary coefficients (device, updated every step)
       integer(kind=4) :: pml_coeff_y_n
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_be_y_left
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_ce_y_left
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_bm_y_left
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_cm_y_left

       ! Left boundary PML limits (device)
       integer(kind=4) :: pml_left_Ex_ii, pml_left_Ex_ij, pml_left_Ex_ji, pml_left_Ex_jj, pml_left_Ex_ki, pml_left_Ex_kj
       integer(kind=4) :: pml_left_Ez_ii, pml_left_Ez_ij, pml_left_Ez_ji, pml_left_Ez_jj, pml_left_Ez_ki, pml_left_Ez_kj
       integer(kind=4) :: pml_left_Hx_ii, pml_left_Hx_ij, pml_left_Hx_ji, pml_left_Hx_jj, pml_left_Hx_ki, pml_left_Hx_kj
       integer(kind=4) :: pml_left_Hz_ii, pml_left_Hz_ij, pml_left_Hz_ji, pml_left_Hz_jj, pml_left_Hz_ki, pml_left_Hz_kj

       ! CPML right boundary - same structure as left
       integer(kind=4) :: pml_right_Ex_ii, pml_right_Ex_ij, pml_right_Ex_ji, pml_right_Ex_jj, pml_right_Ex_ki, pml_right_Ex_kj
       integer(kind=4) :: pml_right_Ez_ii, pml_right_Ez_ij, pml_right_Ez_ji, pml_right_Ez_jj, pml_right_Ez_ki, pml_right_Ez_kj
       integer(kind=4) :: pml_right_Hx_ii, pml_right_Hx_ij, pml_right_Hx_ji, pml_right_Hx_jj, pml_right_Hx_ki, pml_right_Hx_kj
       integer(kind=4) :: pml_right_Hz_ii, pml_right_Hz_ij, pml_right_Hz_ji, pml_right_Hz_jj, pml_right_Hz_ki, pml_right_Hz_kj

       ! CPML down/up boundary - z-dependent coefficients
       integer(kind=4) :: pml_down_Ex_ii, pml_down_Ex_ij, pml_down_Ex_ji, pml_down_Ex_jj, pml_down_Ex_ki, pml_down_Ex_kj
       integer(kind=4) :: pml_down_Ey_ii, pml_down_Ey_ij, pml_down_Ey_ji, pml_down_Ey_jj, pml_down_Ey_ki, pml_down_Ey_kj
       integer(kind=4) :: pml_down_Hx_ii, pml_down_Hx_ij, pml_down_Hx_ji, pml_down_Hx_jj, pml_down_Hx_ki, pml_down_Hx_kj
       integer(kind=4) :: pml_down_Hy_ii, pml_down_Hy_ij, pml_down_Hy_ji, pml_down_Hy_jj, pml_down_Hy_ki, pml_down_Hy_kj
       integer(kind=4) :: pml_down_Exz_nx, pml_down_Exz_ny, pml_down_Exz_nz
       integer(kind=4) :: pml_down_Eyz_nx, pml_down_Eyz_ny, pml_down_Eyz_nz
       integer(kind=4) :: pml_down_Hxz_nx, pml_down_Hxz_ny, pml_down_Hxz_nz
       integer(kind=4) :: pml_down_Hyz_nx, pml_down_Hyz_ny, pml_down_Hyz_nz
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Exz_down
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Eyz_down
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Hxz_down
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Hyz_down
       integer(kind=4) :: pml_coeff_z_n
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_be_z_down
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_ce_z_down
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_bm_z_down
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_cm_z_down
       integer(kind=4) :: pml_up_Ex_ii, pml_up_Ex_ij, pml_up_Ex_ji, pml_up_Ex_jj, pml_up_Ex_ki, pml_up_Ex_kj
       integer(kind=4) :: pml_up_Ey_ii, pml_up_Ey_ij, pml_up_Ey_ji, pml_up_Ey_jj, pml_up_Ey_ki, pml_up_Ey_kj
       integer(kind=4) :: pml_up_Hx_ii, pml_up_Hx_ij, pml_up_Hx_ji, pml_up_Hx_jj, pml_up_Hx_ki, pml_up_Hx_kj
       integer(kind=4) :: pml_up_Hy_ii, pml_up_Hy_ij, pml_up_Hy_ji, pml_up_Hy_jj, pml_up_Hy_ki, pml_up_Hy_kj

       ! CPML back/front boundary - x-dependent coefficients
       integer(kind=4) :: pml_back_Ez_ii, pml_back_Ez_ij, pml_back_Ez_ji, pml_back_Ez_jj, pml_back_Ez_ki, pml_back_Ez_kj
       integer(kind=4) :: pml_back_Ey_ii, pml_back_Ey_ij, pml_back_Ey_ji, pml_back_Ey_jj, pml_back_Ey_ki, pml_back_Ey_kj
       integer(kind=4) :: pml_back_Hz_ii, pml_back_Hz_ij, pml_back_Hz_ji, pml_back_Hz_jj, pml_back_Hz_ki, pml_back_Hz_kj
       integer(kind=4) :: pml_back_Hy_ii, pml_back_Hy_ij, pml_back_Hy_ji, pml_back_Hy_jj, pml_back_Hy_ki, pml_back_Hy_kj
       integer(kind=4) :: pml_back_Ezx_nx, pml_back_Ezx_ny, pml_back_Ezx_nz
       integer(kind=4) :: pml_back_Eyx_nx, pml_back_Eyx_ny, pml_back_Eyx_nz
       integer(kind=4) :: pml_back_Hzx_nx, pml_back_Hzx_ny, pml_back_Hzx_nz
       integer(kind=4) :: pml_back_Hyx_nx, pml_back_Hyx_ny, pml_back_Hyx_nz
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Ezx_back
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Eyx_back
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Hzx_back
       real(kind=rkind), pointer, device, dimension(:,:,:) :: pml_psi_Hyx_back
       integer(kind=4) :: pml_coeff_x_n
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_be_x_back
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_ce_x_back
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_bm_x_back
       real(kind=rkind), pointer, device, dimension(:) :: pml_P_cm_x_back
       integer(kind=4) :: pml_front_Ez_ii, pml_front_Ez_ij, pml_front_Ez_ji, pml_front_Ez_jj, pml_front_Ez_ki, pml_front_Ez_kj
       integer(kind=4) :: pml_front_Ey_ii, pml_front_Ey_ij, pml_front_Ey_ji, pml_front_Ey_jj, pml_front_Ey_ki, pml_front_Ey_kj
       integer(kind=4) :: pml_front_Hz_ii, pml_front_Hz_ij, pml_front_Hz_ji, pml_front_Hz_jj, pml_front_Hz_ki, pml_front_Hz_kj
       integer(kind=4) :: pml_front_Hy_ii, pml_front_Hy_ij, pml_front_Hy_ji, pml_front_Hy_jj, pml_front_Hy_ki, pml_front_Hy_kj

       ! Flags
       logical :: initialized = .false.
       logical :: fields_on_device = .false.
       logical :: pml_left_initialized = .false.
       logical :: pml_right_initialized = .false.
       logical :: pml_down_initialized = .false.
       logical :: pml_up_initialized = .false.
       logical :: pml_back_initialized = .false.
       logical :: pml_front_initialized = .false.
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

      ! Allocate persistent device memory - YEE fields
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

   !--------------------------------------------------------------------------------
   ! Initialize CPML left boundary on GPU - called after InitCPMLBorders
   !--------------------------------------------------------------------------------
   subroutine gpu_init_pml_left(this, P_be_y, P_ce_y, P_bm_y, P_cm_y, &
                                Ex_iEx_ii, Ex_iEx_ij, Ex_iEx_ji, Ex_iEx_jj, Ex_iEx_ki, Ex_iEx_kj, &
                                Ex_iEz_ii, Ex_iEz_ij, Ex_iEz_ji, Ex_iEz_jj, Ex_iEz_ki, Ex_iEz_kj, &
                                Hx_iHx_ii, Hx_iHx_ij, Hx_iHx_ji, Hx_iHx_jj, Hx_iHx_ki, Hx_iHx_kj, &
                                Hx_iHz_ii, Hx_iHz_ij, Hx_iHz_ji, Hx_iHz_jj, Hx_iHz_ki, Hx_iHz_kj, &
                                Ex_ny, Ex_nz, Ez_ny, Ez_nz, Hx_ny, Hx_nz, Hz_ny, Hz_nz)
      class(gpu_state_t), intent(inout) :: this
      real(kind=rkind), dimension(:), intent(in) :: P_be_y, P_ce_y, P_bm_y, P_cm_y
      integer(kind=4), intent(in) :: Ex_iEx_ii, Ex_iEx_ij, Ex_iEx_ji, Ex_iEx_jj, Ex_iEx_ki, Ex_iEx_kj
      integer(kind=4), intent(in) :: Ex_iEz_ii, Ex_iEz_ij, Ex_iEz_ji, Ex_iEz_jj, Ex_iEz_ki, Ex_iEz_kj
      integer(kind=4), intent(in) :: Hx_iHx_ii, Hx_iHx_ij, Hx_iHx_ji, Hx_iHx_jj, Hx_iHx_ki, Hx_iHx_kj
      integer(kind=4), intent(in) :: Hx_iHz_ii, Hx_iHz_ij, Hx_iHz_ji, Hx_iHz_jj, Hx_iHz_ki, Hx_iHz_kj
      integer(kind=4), intent(in) :: Ex_ny, Ex_nz, Ez_ny, Ez_nz, Hx_ny, Hx_nz, Hz_ny, Hz_nz

      integer(kind=4) :: y_lo, y_hi

      if (.not. this%initialized) return

      ! Store PML limits (convert from global 0-based to GPU 1-based)
      this%pml_left_Ex_ii = Ex_iEx_ii + 1; this%pml_left_Ex_ij = Ex_iEx_ij + 1
      this%pml_left_Ex_ji = Ex_iEx_ji + 1; this%pml_left_Ex_jj = Ex_iEx_jj + 1
      this%pml_left_Ex_ki = Ex_iEx_ki + 1; this%pml_left_Ex_kj = Ex_iEx_kj + 1

      this%pml_left_Ez_ii = Ex_iEz_ii + 1; this%pml_left_Ez_ij = Ex_iEz_ij + 1
      this%pml_left_Ez_ji = Ex_iEz_ji + 1; this%pml_left_Ez_jj = Ex_iEz_jj + 1
      this%pml_left_Ez_ki = Ex_iEz_ki + 1; this%pml_left_Ez_kj = Ex_iEz_kj + 1

      this%pml_left_Hx_ii = Hx_iHx_ii + 1; this%pml_left_Hx_ij = Hx_iHx_ij + 1
      this%pml_left_Hx_ji = Hx_iHx_ji + 1; this%pml_left_Hx_jj = Hx_iHx_jj + 1
      this%pml_left_Hx_ki = Hx_iHx_ki + 1; this%pml_left_Hx_kj = Hx_iHx_kj + 1

      this%pml_left_Hz_ii = Hx_iHz_ii + 1; this%pml_left_Hz_ij = Hx_iHz_ij + 1
      this%pml_left_Hz_ji = Hx_iHz_ji + 1; this%pml_left_Hz_jj = Hx_iHz_jj + 1
      this%pml_left_Hz_ki = Hx_iHz_ki + 1; this%pml_left_Hz_kj = Hx_iHz_kj + 1

      ! Allocate psi arrays with global bounds (same as host)
      this%pml_left_Exy_nx = Ex_iEx_ij - Ex_iEx_ii + 1
      this%pml_left_Exy_ny = Ex_iEx_jj - Ex_iEx_ji + 1
      this%pml_left_Exy_nz = Ex_iEx_kj - Ex_iEx_ki + 1
      allocate(this%pml_psi_Exy_left(this%pml_left_Exy_nx, this%pml_left_Exy_ny, this%pml_left_Exy_nz))

      this%pml_left_Ezy_nx = Ex_iEz_ij - Ex_iEz_ii + 1
      this%pml_left_Ezy_ny = Ex_iEz_jj - Ex_iEz_ji + 1
      this%pml_left_Ezy_nz = Ex_iEz_kj - Ex_iEz_ki + 1
      allocate(this%pml_psi_Ezy_left(this%pml_left_Ezy_nx, this%pml_left_Ezy_ny, this%pml_left_Ezy_nz))

      this%pml_left_Hxy_nx = Hx_iHx_ij - Hx_iHx_ii + 1
      this%pml_left_Hxy_ny = Hx_iHx_jj - Hx_iHx_ji + 1
      this%pml_left_Hxy_nz = Hx_iHx_kj - Hx_iHx_ki + 1
      allocate(this%pml_psi_Hxy_left(this%pml_left_Hxy_nx, this%pml_left_Hxy_ny, this%pml_left_Hxy_nz))

      this%pml_left_Hzy_nx = Hx_iHz_ij - Hx_iHz_ii + 1
      this%pml_left_Hzy_ny = Hx_iHz_jj - Hx_iHz_ji + 1
      this%pml_left_Hzy_nz = Hx_iHz_kj - Hx_iHz_ki + 1
      allocate(this%pml_psi_Hzy_left(this%pml_left_Hzy_nx, this%pml_left_Hzy_ny, this%pml_left_Hzy_nz))

      ! Initialize psi arrays to zero on device
      this%pml_psi_Exy_left = 0.0_rkind
      this%pml_psi_Ezy_left = 0.0_rkind
      this%pml_psi_Hxy_left = 0.0_rkind
      this%pml_psi_Hzy_left = 0.0_rkind

      ! Allocate and copy coefficient arrays
      y_lo = lbound(P_be_y, 1)
      y_hi = ubound(P_be_y, 1)
      this%pml_coeff_y_n = y_hi - y_lo + 1
      allocate(this%pml_P_be_y_left(y_lo:y_hi))
      allocate(this%pml_P_ce_y_left(y_lo:y_hi))
      allocate(this%pml_P_bm_y_left(y_lo:y_hi))
      allocate(this%pml_P_cm_y_left(y_lo:y_hi))

      this%pml_P_be_y_left = P_be_y
      this%pml_P_ce_y_left = P_ce_y
      this%pml_P_bm_y_left = P_bm_y
      this%pml_P_cm_y_left = P_cm_y

      this%pml_left_initialized = .true.

   end subroutine gpu_init_pml_left

   !--------------------------------------------------------------------------------
   ! Initialize CPML right boundary on GPU
   !--------------------------------------------------------------------------------
   subroutine gpu_init_pml_right(this, &
                                 Ex_iEx_ii, Ex_iEx_ij, Ex_iEx_ji, Ex_iEx_jj, Ex_iEx_ki, Ex_iEx_kj, &
                                 Ex_iEz_ii, Ex_iEz_ij, Ex_iEz_ji, Ex_iEz_jj, Ex_iEz_ki, Ex_iEz_kj, &
                                 Hx_iHx_ii, Hx_iHx_ij, Hx_iHx_ji, Hx_iHx_jj, Hx_iHx_ki, Hx_iHx_kj, &
                                 Hx_iHz_ii, Hx_iHz_ij, Hx_iHz_ji, Hx_iHz_jj, Hx_iHz_ki, Hx_iHz_kj)
      class(gpu_state_t), intent(inout) :: this
      integer(kind=4), intent(in) :: Ex_iEx_ii, Ex_iEx_ij, Ex_iEx_ji, Ex_iEx_jj, Ex_iEx_ki, Ex_iEx_kj
      integer(kind=4), intent(in) :: Ex_iEz_ii, Ex_iEz_ij, Ex_iEz_ji, Ex_iEz_jj, Ex_iEz_ki, Ex_iEz_kj
      integer(kind=4), intent(in) :: Hx_iHx_ii, Hx_iHx_ij, Hx_iHx_ji, Hx_iHx_jj, Hx_iHx_ki, Hx_iHx_kj
      integer(kind=4), intent(in) :: Hx_iHz_ii, Hx_iHz_ij, Hx_iHz_ji, Hx_iHz_jj, Hx_iHz_ki, Hx_iHz_kj

      if (.not. this%initialized) return

      this%pml_right_Ex_ii = Ex_iEx_ii + 1; this%pml_right_Ex_ij = Ex_iEx_ij + 1
      this%pml_right_Ex_ji = Ex_iEx_ji + 1; this%pml_right_Ex_jj = Ex_iEx_jj + 1
      this%pml_right_Ex_ki = Ex_iEx_ki + 1; this%pml_right_Ex_kj = Ex_iEx_kj + 1

      this%pml_right_Ez_ii = Ex_iEz_ii + 1; this%pml_right_Ez_ij = Ex_iEz_ij + 1
      this%pml_right_Ez_ji = Ex_iEz_ji + 1; this%pml_right_Ez_jj = Ex_iEz_jj + 1
      this%pml_right_Ez_ki = Ex_iEz_ki + 1; this%pml_right_Ez_kj = Ex_iEz_kj + 1

      this%pml_right_Hx_ii = Hx_iHx_ii + 1; this%pml_right_Hx_ij = Hx_iHx_ij + 1
      this%pml_right_Hx_ji = Hx_iHx_ji + 1; this%pml_right_Hx_jj = Hx_iHx_jj + 1
      this%pml_right_Hx_ki = Hx_iHx_ki + 1; this%pml_right_Hx_kj = Hx_iHx_kj + 1

      this%pml_right_Hz_ii = Hx_iHz_ii + 1; this%pml_right_Hz_ij = Hx_iHz_ij + 1
      this%pml_right_Hz_ji = Hx_iHz_ji + 1; this%pml_right_Hz_jj = Hx_iHz_jj + 1
      this%pml_right_Hz_ki = Hx_iHz_ki + 1; this%pml_right_Hz_kj = Hx_iHz_kj + 1

      this%pml_right_initialized = .true.

   end subroutine gpu_init_pml_right

   !--------------------------------------------------------------------------------
   ! Download device data to host - called only when output/probes are needed
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
   ! Destroy GPU state - called at simulation end
   !--------------------------------------------------------------------------------
   subroutine gpu_destroy(this)
      class(gpu_state_t), intent(inout) :: this

      if (.not. this%initialized) return

      ! Deallocate device memory - YEE fields
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

    ! Deallocate CPML left boundary device memory
       if (associated(this%pml_psi_Exy_left)) deallocate(this%pml_psi_Exy_left)
       if (associated(this%pml_psi_Ezy_left)) deallocate(this%pml_psi_Ezy_left)
       if (associated(this%pml_psi_Hxy_left)) deallocate(this%pml_psi_Hxy_left)
       if (associated(this%pml_psi_Hzy_left)) deallocate(this%pml_psi_Hzy_left)
       if (associated(this%pml_P_be_y_left)) deallocate(this%pml_P_be_y_left)
       if (associated(this%pml_P_ce_y_left)) deallocate(this%pml_P_ce_y_left)
       if (associated(this%pml_P_bm_y_left)) deallocate(this%pml_P_bm_y_left)
       if (associated(this%pml_P_cm_y_left)) deallocate(this%pml_P_cm_y_left)

       ! Deallocate CPML down boundary device memory
       if (associated(this%pml_psi_Exz_down)) deallocate(this%pml_psi_Exz_down)
       if (associated(this%pml_psi_Eyz_down)) deallocate(this%pml_psi_Eyz_down)
       if (associated(this%pml_psi_Hxz_down)) deallocate(this%pml_psi_Hxz_down)
       if (associated(this%pml_psi_Hyz_down)) deallocate(this%pml_psi_Hyz_down)
       if (associated(this%pml_P_be_z_down)) deallocate(this%pml_P_be_z_down)
       if (associated(this%pml_P_ce_z_down)) deallocate(this%pml_P_ce_z_down)
       if (associated(this%pml_P_bm_z_down)) deallocate(this%pml_P_bm_z_down)
       if (associated(this%pml_P_cm_z_down)) deallocate(this%pml_P_cm_z_down)

       ! Deallocate CPML back boundary device memory
       if (associated(this%pml_psi_Ezx_back)) deallocate(this%pml_psi_Ezx_back)
       if (associated(this%pml_psi_Eyx_back)) deallocate(this%pml_psi_Eyx_back)
       if (associated(this%pml_psi_Hzx_back)) deallocate(this%pml_psi_Hzx_back)
       if (associated(this%pml_psi_Hyx_back)) deallocate(this%pml_psi_Hyx_back)
       if (associated(this%pml_P_be_x_back)) deallocate(this%pml_P_be_x_back)
       if (associated(this%pml_P_ce_x_back)) deallocate(this%pml_P_ce_x_back)
       if (associated(this%pml_P_bm_x_back)) deallocate(this%pml_P_bm_x_back)
       if (associated(this%pml_P_cm_x_back)) deallocate(this%pml_P_cm_x_back)

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
       this%pml_left_initialized = .false.
       this%pml_right_initialized = .false.
       this%pml_down_initialized = .false.
       this%pml_up_initialized = .false.
       this%pml_back_initialized = .false.
       this%pml_front_initialized = .false.

    end subroutine gpu_destroy

    !--------------------------------------------------------------------------------
    ! Update CPML down/up boundary coefficients on device - called every step
    !--------------------------------------------------------------------------------
    subroutine gpu_update_pml_down_coeffs(this, P_be_z, P_ce_z, P_bm_z, P_cm_z)
       class(gpu_state_t), intent(inout) :: this
       real(kind=rkind), dimension(:), intent(in) :: P_be_z, P_ce_z, P_bm_z, P_cm_z

       if (.not. this%pml_down_initialized) return

       this%pml_P_be_z_down = P_be_z
       this%pml_P_ce_z_down = P_ce_z
       this%pml_P_bm_z_down = P_bm_z
       this%pml_P_cm_z_down = P_cm_z

    end subroutine gpu_update_pml_down_coeffs

    !--------------------------------------------------------------------------------
    ! Update CPML back/front boundary coefficients on device - called every step
    !--------------------------------------------------------------------------------
    subroutine gpu_update_pml_back_coeffs(this, P_be_x, P_ce_x, P_bm_x, P_cm_x)
       class(gpu_state_t), intent(inout) :: this
       real(kind=rkind), dimension(:), intent(in) :: P_be_x, P_ce_x, P_bm_x, P_cm_x

       if (.not. this%pml_back_initialized) return

       this%pml_P_be_x_back = P_be_x
       this%pml_P_ce_x_back = P_ce_x
       this%pml_P_bm_x_back = P_bm_x
       this%pml_P_cm_x_back = P_cm_x

    end subroutine gpu_update_pml_back_coeffs

    !--------------------------------------------------------------------------------
    ! Initialize CPML down boundary on GPU
    !--------------------------------------------------------------------------------
    subroutine gpu_init_pml_down(this, P_be_z, P_ce_z, P_bm_z, P_cm_z, &
                                 Ex_iEx_ii, Ex_iEx_ij, Ex_iEx_ji, Ex_iEx_jj, Ex_iEx_ki, Ex_iEx_kj, &
                                 Ex_iEy_ii, Ex_iEy_ij, Ex_iEy_ji, Ex_iEy_jj, Ex_iEy_ki, Ex_iEy_kj, &
                                 Hx_iHx_ii, Hx_iHx_ij, Hx_iHx_ji, Hx_iHx_jj, Hx_iHx_ki, Hx_iHx_kj, &
                                 Hx_iHy_ii, Hx_iHy_ij, Hx_iHy_ji, Hx_iHy_jj, Hx_iHy_ki, Hx_iHy_kj)
       class(gpu_state_t), intent(inout) :: this
       real(kind=rkind), dimension(:), intent(in) :: P_be_z, P_ce_z, P_bm_z, P_cm_z
       integer(kind=4), intent(in) :: Ex_iEx_ii, Ex_iEx_ij, Ex_iEx_ji, Ex_iEx_jj, Ex_iEx_ki, Ex_iEx_kj
       integer(kind=4), intent(in) :: Ex_iEy_ii, Ex_iEy_ij, Ex_iEy_ji, Ex_iEy_jj, Ex_iEy_ki, Ex_iEy_kj
       integer(kind=4), intent(in) :: Hx_iHx_ii, Hx_iHx_ij, Hx_iHx_ji, Hx_iHx_jj, Hx_iHx_ki, Hx_iHx_kj
       integer(kind=4), intent(in) :: Hx_iHy_ii, Hx_iHy_ij, Hx_iHy_ji, Hx_iHy_jj, Hx_iHy_ki, Hx_iHy_kj

       integer(kind=4) :: z_lo, z_hi

       if (.not. this%initialized) return

       this%pml_down_Ex_ii = Ex_iEx_ii + 1; this%pml_down_Ex_ij = Ex_iEx_ij + 1
       this%pml_down_Ex_ji = Ex_iEx_ji + 1; this%pml_down_Ex_jj = Ex_iEx_jj + 1
       this%pml_down_Ex_ki = Ex_iEx_ki + 1; this%pml_down_Ex_kj = Ex_iEx_kj + 1

       this%pml_down_Ey_ii = Ex_iEy_ii + 1; this%pml_down_Ey_ij = Ex_iEy_ij + 1
       this%pml_down_Ey_ji = Ex_iEy_ji + 1; this%pml_down_Ey_jj = Ex_iEy_jj + 1
       this%pml_down_Ey_ki = Ex_iEy_ki + 1; this%pml_down_Ey_kj = Ex_iEy_kj + 1

       this%pml_down_Hx_ii = Hx_iHx_ii + 1; this%pml_down_Hx_ij = Hx_iHx_ij + 1
       this%pml_down_Hx_ji = Hx_iHx_ji + 1; this%pml_down_Hx_jj = Hx_iHx_jj + 1
       this%pml_down_Hx_ki = Hx_iHx_ki + 1; this%pml_down_Hx_kj = Hx_iHx_kj + 1

       this%pml_down_Hy_ii = Hx_iHy_ii + 1; this%pml_down_Hy_ij = Hx_iHy_ij + 1
       this%pml_down_Hy_ji = Hx_iHy_ji + 1; this%pml_down_Hy_jj = Hx_iHy_jj + 1
       this%pml_down_Hy_ki = Hx_iHy_ki + 1; this%pml_down_Hy_kj = Hx_iHy_kj + 1

       this%pml_down_Exz_nx = Ex_iEx_ij - Ex_iEx_ii + 1
       this%pml_down_Exz_ny = Ex_iEx_jj - Ex_iEx_ji + 1
       this%pml_down_Exz_nz = Ex_iEx_kj - Ex_iEx_ki + 1
       allocate(this%pml_psi_Exz_down(this%pml_down_Exz_nx, this%pml_down_Exz_ny, this%pml_down_Exz_nz))

       this%pml_down_Eyz_nx = Ex_iEy_ij - Ex_iEy_ii + 1
       this%pml_down_Eyz_ny = Ex_iEy_jj - Ex_iEy_ji + 1
       this%pml_down_Eyz_nz = Ex_iEy_kj - Ex_iEy_ki + 1
       allocate(this%pml_psi_Eyz_down(this%pml_down_Eyz_nx, this%pml_down_Eyz_ny, this%pml_down_Eyz_nz))

       this%pml_down_Hxz_nx = Hx_iHx_ij - Hx_iHx_ii + 1
       this%pml_down_Hxz_ny = Hx_iHx_jj - Hx_iHx_ji + 1
       this%pml_down_Hxz_nz = Hx_iHx_kj - Hx_iHx_ki + 1
       allocate(this%pml_psi_Hxz_down(this%pml_down_Hxz_nx, this%pml_down_Hxz_ny, this%pml_down_Hxz_nz))

       this%pml_down_Hyz_nx = Hx_iHy_ij - Hx_iHy_ii + 1
       this%pml_down_Hyz_ny = Hx_iHy_jj - Hx_iHy_ji + 1
       this%pml_down_Hyz_nz = Hx_iHy_kj - Hx_iHy_ki + 1
       allocate(this%pml_psi_Hyz_down(this%pml_down_Hyz_nx, this%pml_down_Hyz_ny, this%pml_down_Hyz_nz))

       this%pml_psi_Exz_down = 0.0_rkind
       this%pml_psi_Eyz_down = 0.0_rkind
       this%pml_psi_Hxz_down = 0.0_rkind
       this%pml_psi_Hyz_down = 0.0_rkind

       z_lo = lbound(P_be_z, 1)
       z_hi = ubound(P_be_z, 1)
       this%pml_coeff_z_n = z_hi - z_lo + 1
       allocate(this%pml_P_be_z_down(z_lo:z_hi))
       allocate(this%pml_P_ce_z_down(z_lo:z_hi))
       allocate(this%pml_P_bm_z_down(z_lo:z_hi))
       allocate(this%pml_P_cm_z_down(z_lo:z_hi))

       this%pml_P_be_z_down = P_be_z
       this%pml_P_ce_z_down = P_ce_z
       this%pml_P_bm_z_down = P_bm_z
       this%pml_P_cm_z_down = P_cm_z

       this%pml_down_initialized = .true.

    end subroutine gpu_init_pml_down

    !--------------------------------------------------------------------------------
    ! Initialize CPML up boundary on GPU
    !--------------------------------------------------------------------------------
    subroutine gpu_init_pml_up(this, &
                               Ex_iEx_ii, Ex_iEx_ij, Ex_iEx_ji, Ex_iEx_jj, Ex_iEx_ki, Ex_iEx_kj, &
                               Ex_iEy_ii, Ex_iEy_ij, Ex_iEy_ji, Ex_iEy_jj, Ex_iEy_ki, Ex_iEy_kj, &
                               Hx_iHx_ii, Hx_iHx_ij, Hx_iHx_ji, Hx_iHx_jj, Hx_iHx_ki, Hx_iHx_kj, &
                               Hx_iHy_ii, Hx_iHy_ij, Hx_iHy_ji, Hx_iHy_jj, Hx_iHy_ki, Hx_iHy_kj)
       class(gpu_state_t), intent(inout) :: this
       integer(kind=4), intent(in) :: Ex_iEx_ii, Ex_iEx_ij, Ex_iEx_ji, Ex_iEx_jj, Ex_iEx_ki, Ex_iEx_kj
       integer(kind=4), intent(in) :: Ex_iEy_ii, Ex_iEy_ij, Ex_iEy_ji, Ex_iEy_jj, Ex_iEy_ki, Ex_iEy_kj
       integer(kind=4), intent(in) :: Hx_iHx_ii, Hx_iHx_ij, Hx_iHx_ji, Hx_iHx_jj, Hx_iHx_ki, Hx_iHx_kj
       integer(kind=4), intent(in) :: Hx_iHy_ii, Hx_iHy_ij, Hx_iHy_ji, Hx_iHy_jj, Hx_iHy_ki, Hx_iHy_kj

       if (.not. this%initialized) return

       this%pml_up_Ex_ii = Ex_iEx_ii + 1; this%pml_up_Ex_ij = Ex_iEx_ij + 1
       this%pml_up_Ex_ji = Ex_iEx_ji + 1; this%pml_up_Ex_jj = Ex_iEx_jj + 1
       this%pml_up_Ex_ki = Ex_iEx_ki + 1; this%pml_up_Ex_kj = Ex_iEx_kj + 1

       this%pml_up_Ey_ii = Ex_iEy_ii + 1; this%pml_up_Ey_ij = Ex_iEy_ij + 1
       this%pml_up_Ey_ji = Ex_iEy_ji + 1; this%pml_up_Ey_jj = Ex_iEy_jj + 1
       this%pml_up_Ey_ki = Ex_iEy_ki + 1; this%pml_up_Ey_kj = Ex_iEy_kj + 1

       this%pml_up_Hx_ii = Hx_iHx_ii + 1; this%pml_up_Hx_ij = Hx_iHx_ij + 1
       this%pml_up_Hx_ji = Hx_iHx_ji + 1; this%pml_up_Hx_jj = Hx_iHx_jj + 1
       this%pml_up_Hx_ki = Hx_iHx_ki + 1; this%pml_up_Hx_kj = Hx_iHx_kj + 1

       this%pml_up_Hy_ii = Hx_iHy_ii + 1; this%pml_up_Hy_ij = Hx_iHy_ij + 1
       this%pml_up_Hy_ji = Hx_iHy_ji + 1; this%pml_up_Hy_jj = Hx_iHy_jj + 1
       this%pml_up_Hy_ki = Hx_iHy_ki + 1; this%pml_up_Hy_kj = Hx_iHy_kj + 1

       this%pml_up_initialized = .true.

    end subroutine gpu_init_pml_up

    !--------------------------------------------------------------------------------
    ! Initialize CPML back boundary on GPU
    !--------------------------------------------------------------------------------
    subroutine gpu_init_pml_back(this, P_be_x, P_ce_x, P_bm_x, P_cm_x, &
                                 Ez_iEz_ii, Ez_iEz_ij, Ez_iEz_ji, Ez_iEz_jj, Ez_iEz_ki, Ez_iEz_kj, &
                                 Ez_iEy_ii, Ez_iEy_ij, Ez_iEy_ji, Ez_iEy_jj, Ez_iEy_ki, Ez_iEy_kj, &
                                 Hz_iHz_ii, Hz_iHz_ij, Hz_iHz_ji, Hz_iHz_jj, Hz_iHz_ki, Hz_iHz_kj, &
                                 Hz_iHy_ii, Hz_iHy_ij, Hz_iHy_ji, Hz_iHy_jj, Hz_iHy_ki, Hz_iHy_kj)
       class(gpu_state_t), intent(inout) :: this
       real(kind=rkind), dimension(:), intent(in) :: P_be_x, P_ce_x, P_bm_x, P_cm_x
       integer(kind=4), intent(in) :: Ez_iEz_ii, Ez_iEz_ij, Ez_iEz_ji, Ez_iEz_jj, Ez_iEz_ki, Ez_iEz_kj
       integer(kind=4), intent(in) :: Ez_iEy_ii, Ez_iEy_ij, Ez_iEy_ji, Ez_iEy_jj, Ez_iEy_ki, Ez_iEy_kj
       integer(kind=4), intent(in) :: Hz_iHz_ii, Hz_iHz_ij, Hz_iHz_ji, Hz_iHz_jj, Hz_iHz_ki, Hz_iHz_kj
       integer(kind=4), intent(in) :: Hz_iHy_ii, Hz_iHy_ij, Hz_iHy_ji, Hz_iHy_jj, Hz_iHy_ki, Hz_iHy_kj

       integer(kind=4) :: x_lo, x_hi

       if (.not. this%initialized) return

       this%pml_back_Ez_ii = Ez_iEz_ii + 1; this%pml_back_Ez_ij = Ez_iEz_ij + 1
       this%pml_back_Ez_ji = Ez_iEz_ji + 1; this%pml_back_Ez_jj = Ez_iEz_jj + 1
       this%pml_back_Ez_ki = Ez_iEz_ki + 1; this%pml_back_Ez_kj = Ez_iEz_kj + 1

       this%pml_back_Ey_ii = Ez_iEy_ii + 1; this%pml_back_Ey_ij = Ez_iEy_ij + 1
       this%pml_back_Ey_ji = Ez_iEy_ji + 1; this%pml_back_Ey_jj = Ez_iEy_jj + 1
       this%pml_back_Ey_ki = Ez_iEy_ki + 1; this%pml_back_Ey_kj = Ez_iEy_kj + 1

       this%pml_back_Hz_ii = Hz_iHz_ii + 1; this%pml_back_Hz_ij = Hz_iHz_ij + 1
       this%pml_back_Hz_ji = Hz_iHz_ji + 1; this%pml_back_Hz_jj = Hz_iHz_jj + 1
       this%pml_back_Hz_ki = Hz_iHz_ki + 1; this%pml_back_Hz_kj = Hz_iHz_kj + 1

       this%pml_back_Hy_ii = Hz_iHy_ii + 1; this%pml_back_Hy_ij = Hz_iHy_ij + 1
       this%pml_back_Hy_ji = Hz_iHy_ji + 1; this%pml_back_Hy_jj = Hz_iHy_jj + 1
       this%pml_back_Hy_ki = Hz_iHy_ki + 1; this%pml_back_Hy_kj = Hz_iHy_kj + 1

       this%pml_back_Ezx_nx = Ez_iEz_ij - Ez_iEz_ii + 1
       this%pml_back_Ezx_ny = Ez_iEz_jj - Ez_iEz_ji + 1
       this%pml_back_Ezx_nz = Ez_iEz_kj - Ez_iEz_ki + 1
       allocate(this%pml_psi_Ezx_back(this%pml_back_Ezx_nx, this%pml_back_Ezx_ny, this%pml_back_Ezx_nz))

       this%pml_back_Eyx_nx = Ez_iEy_ij - Ez_iEy_ii + 1
       this%pml_back_Eyx_ny = Ez_iEy_jj - Ez_iEy_ji + 1
       this%pml_back_Eyx_nz = Ez_iEy_kj - Ez_iEy_ki + 1
       allocate(this%pml_psi_Eyx_back(this%pml_back_Eyx_nx, this%pml_back_Eyx_ny, this%pml_back_Eyx_nz))

       this%pml_back_Hzx_nx = Hz_iHz_ij - Hz_iHz_ii + 1
       this%pml_back_Hzx_ny = Hz_iHz_jj - Hz_iHz_ji + 1
       this%pml_back_Hzx_nz = Hz_iHz_kj - Hz_iHz_ki + 1
       allocate(this%pml_psi_Hzx_back(this%pml_back_Hzx_nx, this%pml_back_Hzx_ny, this%pml_back_Hzx_nz))

       this%pml_back_Hyx_nx = Hz_iHy_ij - Hz_iHy_ii + 1
       this%pml_back_Hyx_ny = Hz_iHy_jj - Hz_iHy_ji + 1
       this%pml_back_Hyx_nz = Hz_iHy_kj - Hz_iHy_ki + 1
       allocate(this%pml_psi_Hyx_back(this%pml_back_Hyx_nx, this%pml_back_Hyx_ny, this%pml_back_Hyx_nz))

       this%pml_psi_Ezx_back = 0.0_rkind
       this%pml_psi_Eyx_back = 0.0_rkind
       this%pml_psi_Hzx_back = 0.0_rkind
       this%pml_psi_Hyx_back = 0.0_rkind

       x_lo = lbound(P_be_x, 1)
       x_hi = ubound(P_be_x, 1)
       this%pml_coeff_x_n = x_hi - x_lo + 1
       allocate(this%pml_P_be_x_back(x_lo:x_hi))
       allocate(this%pml_P_ce_x_back(x_lo:x_hi))
       allocate(this%pml_P_bm_x_back(x_lo:x_hi))
       allocate(this%pml_P_cm_x_back(x_lo:x_hi))

       this%pml_P_be_x_back = P_be_x
       this%pml_P_ce_x_back = P_ce_x
       this%pml_P_bm_x_back = P_bm_x
       this%pml_P_cm_x_back = P_cm_x

       this%pml_back_initialized = .true.

    end subroutine gpu_init_pml_back

    !--------------------------------------------------------------------------------
    ! Initialize CPML front boundary on GPU
    !--------------------------------------------------------------------------------
    subroutine gpu_init_pml_front(this, &
                                  Ez_iEz_ii, Ez_iEz_ij, Ez_iEz_ji, Ez_iEz_jj, Ez_iEz_ki, Ez_iEz_kj, &
                                  Ez_iEy_ii, Ez_iEy_ij, Ez_iEy_ji, Ez_iEy_jj, Ez_iEy_ki, Ez_iEy_kj, &
                                  Hz_iHz_ii, Hz_iHz_ij, Hz_iHz_ji, Hz_iHz_jj, Hz_iHz_ki, Hz_iHz_kj, &
                                  Hz_iHy_ii, Hz_iHy_ij, Hz_iHy_ji, Hz_iHy_jj, Hz_iHy_ki, Hz_iHy_kj)
       class(gpu_state_t), intent(inout) :: this
       integer(kind=4), intent(in) :: Ez_iEz_ii, Ez_iEz_ij, Ez_iEz_ji, Ez_iEz_jj, Ez_iEz_ki, Ez_iEz_kj
       integer(kind=4), intent(in) :: Ez_iEy_ii, Ez_iEy_ij, Ez_iEy_ji, Ez_iEy_jj, Ez_iEy_ki, Ez_iEy_kj
       integer(kind=4), intent(in) :: Hz_iHz_ii, Hz_iHz_ij, Hz_iHz_ji, Hz_iHz_jj, Hz_iHz_ki, Hz_iHz_kj
       integer(kind=4), intent(in) :: Hz_iHy_ii, Hz_iHy_ij, Hz_iHy_ji, Hz_iHy_jj, Hz_iHy_ki, Hz_iHy_kj

       if (.not. this%initialized) return

       this%pml_front_Ez_ii = Ez_iEz_ii + 1; this%pml_front_Ez_ij = Ez_iEz_ij + 1
       this%pml_front_Ez_ji = Ez_iEz_ji + 1; this%pml_front_Ez_jj = Ez_iEz_jj + 1
       this%pml_front_Ez_ki = Ez_iEz_ki + 1; this%pml_front_Ez_kj = Ez_iEz_kj + 1

       this%pml_front_Ey_ii = Ez_iEy_ii + 1; this%pml_front_Ey_ij = Ez_iEy_ij + 1
       this%pml_front_Ey_ji = Ez_iEy_ji + 1; this%pml_front_Ey_jj = Ez_iEy_jj + 1
       this%pml_front_Ey_ki = Ez_iEy_ki + 1; this%pml_front_Ey_kj = Ez_iEy_kj + 1

       this%pml_front_Hz_ii = Hz_iHz_ii + 1; this%pml_front_Hz_ij = Hz_iHz_ij + 1
       this%pml_front_Hz_ji = Hz_iHz_ji + 1; this%pml_front_Hz_jj = Hz_iHz_jj + 1
       this%pml_front_Hz_ki = Hz_iHz_ki + 1; this%pml_front_Hz_kj = Hz_iHz_kj + 1

       this%pml_front_Hy_ii = Hz_iHy_ii + 1; this%pml_front_Hy_ij = Hz_iHy_ij + 1
       this%pml_front_Hy_ji = Hz_iHy_ji + 1; this%pml_front_Hy_jj = Hz_iHy_jj + 1
       this%pml_front_Hy_ki = Hz_iHy_ki + 1; this%pml_front_Hy_kj = Hz_iHy_kj + 1

       this%pml_front_initialized = .true.

    end subroutine gpu_init_pml_front

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
   ! Update CPML left boundary coefficients on device - called every step
   !--------------------------------------------------------------------------------
   subroutine gpu_update_pml_left_coeffs(this, P_be_y, P_ce_y, P_bm_y, P_cm_y)
      class(gpu_state_t), intent(inout) :: this
      real(kind=rkind), dimension(:), intent(in) :: P_be_y, P_ce_y, P_bm_y, P_cm_y

      if (.not. this%pml_left_initialized) return

      this%pml_P_be_y_left = P_be_y
      this%pml_P_ce_y_left = P_ce_y
      this%pml_P_bm_y_left = P_bm_y
      this%pml_P_cm_y_left = P_cm_y

   end subroutine gpu_update_pml_left_coeffs

end module gpu_core_m
