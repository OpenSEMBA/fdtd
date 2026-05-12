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

        ! MUR boundary - persistent device coefficient arrays (per media)
        integer(kind=4) :: mur_numMedia
        real(kind=rkind), pointer, device, dimension(:) :: mur_left_CAB1, mur_left_CAB3, mur_left_cab4
        real(kind=rkind), pointer, device, dimension(:) :: mur_right_CAB1, mur_right_CAB3, mur_right_cab4
        real(kind=rkind), pointer, device, dimension(:) :: mur_down_CAB1, mur_down_CAB3, mur_down_cab4
        real(kind=rkind), pointer, device, dimension(:) :: mur_up_CAB1, mur_up_CAB3, mur_up_cab4
        real(kind=rkind), pointer, device, dimension(:) :: mur_back_CAB1, mur_back_CAB3, mur_back_cab4
        real(kind=rkind), pointer, device, dimension(:) :: mur_front_CAB1, mur_front_CAB3, mur_front_cab4

        ! MUR boundary - persistent device past-field arrays
        real(kind=rkind), pointer, device, dimension(:,:,:) :: mur_past_Hx_left, mur_past_Hz_left
        real(kind=rkind), pointer, device, dimension(:,:,:) :: mur_past_Hx_right, mur_past_Hz_right
        real(kind=rkind), pointer, device, dimension(:,:,:) :: mur_past_Hy_down, mur_past_Hx_down
        real(kind=rkind), pointer, device, dimension(:,:,:) :: mur_past_Hy_up, mur_past_Hx_up
        real(kind=rkind), pointer, device, dimension(:,:,:) :: mur_past_Hz_back, mur_past_Hy_back
        real(kind=rkind), pointer, device, dimension(:,:,:) :: mur_past_Hz_front, mur_past_Hy_front

        ! MUR boundary limits
        integer(kind=4) :: mur_left_Hx_ii, mur_left_Hx_ij, mur_left_Hx_ji, mur_left_Hx_jj, mur_left_Hx_ki, mur_left_Hx_kj
        integer(kind=4) :: mur_left_Hz_ii, mur_left_Hz_ij, mur_left_Hz_ji, mur_left_Hz_jj, mur_left_Hz_ki, mur_left_Hz_kj
        integer(kind=4) :: mur_right_Hx_ii, mur_right_Hx_ij, mur_right_Hx_ji, mur_right_Hx_jj, mur_right_Hx_ki, mur_right_Hx_kj
        integer(kind=4) :: mur_right_Hz_ii, mur_right_Hz_ij, mur_right_Hz_ji, mur_right_Hz_jj, mur_right_Hz_ki, mur_right_Hz_kj
        integer(kind=4) :: mur_down_Hy_ii, mur_down_Hy_ij, mur_down_Hy_ji, mur_down_Hy_jj, mur_down_Hy_ki, mur_down_Hy_kj
        integer(kind=4) :: mur_down_Hx_ii, mur_down_Hx_ij, mur_down_Hx_ji, mur_down_Hx_jj, mur_down_Hx_ki, mur_down_Hx_kj
        integer(kind=4) :: mur_up_Hy_ii, mur_up_Hy_ij, mur_up_Hy_ji, mur_up_Hy_jj, mur_up_Hy_ki, mur_up_Hy_kj
        integer(kind=4) :: mur_up_Hx_ii, mur_up_Hx_ij, mur_up_Hx_ji, mur_up_Hx_jj, mur_up_Hx_ki, mur_up_Hx_kj
        integer(kind=4) :: mur_back_Hz_ii, mur_back_Hz_ij, mur_back_Hz_ji, mur_back_Hz_jj, mur_back_Hz_ki, mur_back_Hz_kj
        integer(kind=4) :: mur_back_Hy_ii, mur_back_Hy_ij, mur_back_Hy_ji, mur_back_Hy_jj, mur_back_Hy_ki, mur_back_Hy_kj
        integer(kind=4) :: mur_front_Hz_ii, mur_front_Hz_ij, mur_front_Hz_ji, mur_front_Hz_jj, mur_front_Hz_ki, mur_front_Hz_kj
        integer(kind=4) :: mur_front_Hy_ii, mur_front_Hy_ij, mur_front_Hy_ji, mur_front_Hy_jj, mur_front_Hy_ki, mur_front_Hy_kj

        ! MUR flags
        logical :: mur_initialized = .false.

        ! Flags
       logical :: initialized = .false.
       logical :: fields_on_device = .false.
       logical :: pml_left_initialized = .false.
       logical :: pml_right_initialized = .false.
       logical :: pml_down_initialized = .false.
       logical :: pml_up_initialized = .false.
       logical :: pml_back_initialized = .false.
       logical :: pml_front_initialized = .false.

       ! Download tracking for lazy field transfer
       integer(kind=4) :: last_download_step = 0
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
       if (.not. this%fields_on_device) return

       this%Ex = this%Ex_d
       this%Ey = this%Ey_d
       this%Ez = this%Ez_d
       this%Hx = this%Hx_d
       this%Hy = this%Hy_d
       this%Hz = this%Hz_d
       this%fields_on_device = .false.
       this%last_download_step = 0

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

        ! Deallocate MUR device memory
        if (associated(this%mur_left_CAB1)) deallocate(this%mur_left_CAB1)
        if (associated(this%mur_left_CAB3)) deallocate(this%mur_left_CAB3)
        if (associated(this%mur_left_cab4)) deallocate(this%mur_left_cab4)
        if (associated(this%mur_right_CAB1)) deallocate(this%mur_right_CAB1)
        if (associated(this%mur_right_CAB3)) deallocate(this%mur_right_CAB3)
        if (associated(this%mur_right_cab4)) deallocate(this%mur_right_cab4)
        if (associated(this%mur_down_CAB1)) deallocate(this%mur_down_CAB1)
        if (associated(this%mur_down_CAB3)) deallocate(this%mur_down_CAB3)
        if (associated(this%mur_down_cab4)) deallocate(this%mur_down_cab4)
        if (associated(this%mur_up_CAB1)) deallocate(this%mur_up_CAB1)
        if (associated(this%mur_up_CAB3)) deallocate(this%mur_up_CAB3)
        if (associated(this%mur_up_cab4)) deallocate(this%mur_up_cab4)
        if (associated(this%mur_back_CAB1)) deallocate(this%mur_back_CAB1)
        if (associated(this%mur_back_CAB3)) deallocate(this%mur_back_CAB3)
        if (associated(this%mur_back_cab4)) deallocate(this%mur_back_cab4)
        if (associated(this%mur_front_CAB1)) deallocate(this%mur_front_CAB1)
        if (associated(this%mur_front_CAB3)) deallocate(this%mur_front_CAB3)
        if (associated(this%mur_front_cab4)) deallocate(this%mur_front_cab4)
        if (associated(this%mur_past_Hx_left)) deallocate(this%mur_past_Hx_left)
        if (associated(this%mur_past_Hz_left)) deallocate(this%mur_past_Hz_left)
        if (associated(this%mur_past_Hx_right)) deallocate(this%mur_past_Hx_right)
        if (associated(this%mur_past_Hz_right)) deallocate(this%mur_past_Hz_right)
        if (associated(this%mur_past_Hy_down)) deallocate(this%mur_past_Hy_down)
        if (associated(this%mur_past_Hx_down)) deallocate(this%mur_past_Hx_down)
        if (associated(this%mur_past_Hy_up)) deallocate(this%mur_past_Hy_up)
        if (associated(this%mur_past_Hx_up)) deallocate(this%mur_past_Hx_up)
        if (associated(this%mur_past_Hz_back)) deallocate(this%mur_past_Hz_back)
        if (associated(this%mur_past_Hy_back)) deallocate(this%mur_past_Hy_back)
        if (associated(this%mur_past_Hz_front)) deallocate(this%mur_past_Hz_front)
        if (associated(this%mur_past_Hy_front)) deallocate(this%mur_past_Hy_front)

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
       this%last_download_step = -1024

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

    !--------------------------------------------------------------------------------
    ! Initialize MUR boundary coefficients on GPU - called after InitMURBorders
    !--------------------------------------------------------------------------------
    subroutine gpu_init_mur_coeffs(this, numMedia, &
                                   left_CAB1, left_CAB3, left_cab4, &
                                   right_CAB1, right_CAB3, right_cab4, &
                                   down_CAB1, down_CAB3, down_cab4, &
                                   up_CAB1, up_CAB3, up_cab4, &
                                   back_CAB1, back_CAB3, back_cab4, &
                                   front_CAB1, front_CAB3, front_cab4)
       class(gpu_state_t), intent(inout) :: this
       integer(kind=4), intent(in) :: numMedia
       real(kind=rkind), dimension(:), intent(in) :: left_CAB1, left_CAB3, left_cab4
       real(kind=rkind), dimension(:), intent(in) :: right_CAB1, right_CAB3, right_cab4
       real(kind=rkind), dimension(:), intent(in) :: down_CAB1, down_CAB3, down_cab4
       real(kind=rkind), dimension(:), intent(in) :: up_CAB1, up_CAB3, up_cab4
       real(kind=rkind), dimension(:), intent(in) :: back_CAB1, back_CAB3, back_cab4
       real(kind=rkind), dimension(:), intent(in) :: front_CAB1, front_CAB3, front_cab4

       integer(kind=4) :: lo, hi

       if (.not. this%initialized) return

       this%mur_numMedia = numMedia

       lo = lbound(left_CAB1, 1); hi = ubound(left_CAB1, 1)
       allocate(this%mur_left_CAB1(lo:hi)); this%mur_left_CAB1 = left_CAB1
       allocate(this%mur_left_CAB3(lo:hi)); this%mur_left_CAB3 = left_CAB3
       allocate(this%mur_left_cab4(lo:hi)); this%mur_left_cab4 = left_cab4
       allocate(this%mur_right_CAB1(lo:hi)); this%mur_right_CAB1 = right_CAB1
       allocate(this%mur_right_CAB3(lo:hi)); this%mur_right_CAB3 = right_CAB3
       allocate(this%mur_right_cab4(lo:hi)); this%mur_right_cab4 = right_cab4
       allocate(this%mur_down_CAB1(lo:hi)); this%mur_down_CAB1 = down_CAB1
       allocate(this%mur_down_CAB3(lo:hi)); this%mur_down_CAB3 = down_CAB3
       allocate(this%mur_down_cab4(lo:hi)); this%mur_down_cab4 = down_cab4
       allocate(this%mur_up_CAB1(lo:hi)); this%mur_up_CAB1 = up_CAB1
       allocate(this%mur_up_CAB3(lo:hi)); this%mur_up_CAB3 = up_CAB3
       allocate(this%mur_up_cab4(lo:hi)); this%mur_up_cab4 = up_cab4
       allocate(this%mur_back_CAB1(lo:hi)); this%mur_back_CAB1 = back_CAB1
       allocate(this%mur_back_CAB3(lo:hi)); this%mur_back_CAB3 = back_CAB3
       allocate(this%mur_back_cab4(lo:hi)); this%mur_back_cab4 = back_cab4
       allocate(this%mur_front_CAB1(lo:hi)); this%mur_front_CAB1 = front_CAB1
       allocate(this%mur_front_CAB3(lo:hi)); this%mur_front_CAB3 = front_CAB3
       allocate(this%mur_front_cab4(lo:hi)); this%mur_front_cab4 = front_cab4

       this%mur_initialized = .true.

    end subroutine gpu_init_mur_coeffs

    !--------------------------------------------------------------------------------
    ! Update MUR coefficients on device - called every step
    !--------------------------------------------------------------------------------
    subroutine gpu_update_mur_coeffs(this, &
                                     left_CAB1, left_CAB3, left_cab4, &
                                     right_CAB1, right_CAB3, right_cab4, &
                                     down_CAB1, down_CAB3, down_cab4, &
                                     up_CAB1, up_CAB3, up_cab4, &
                                     back_CAB1, back_CAB3, back_cab4, &
                                     front_CAB1, front_CAB3, front_cab4)
       class(gpu_state_t), intent(inout) :: this
       real(kind=rkind), dimension(:), intent(in) :: left_CAB1, left_CAB3, left_cab4
       real(kind=rkind), dimension(:), intent(in) :: right_CAB1, right_CAB3, right_cab4
       real(kind=rkind), dimension(:), intent(in) :: down_CAB1, down_CAB3, down_cab4
       real(kind=rkind), dimension(:), intent(in) :: up_CAB1, up_CAB3, up_cab4
       real(kind=rkind), dimension(:), intent(in) :: back_CAB1, back_CAB3, back_cab4
       real(kind=rkind), dimension(:), intent(in) :: front_CAB1, front_CAB3, front_cab4

       if (.not. this%mur_initialized) return

       this%mur_left_CAB1 = left_CAB1; this%mur_left_CAB3 = left_CAB3; this%mur_left_cab4 = left_cab4
       this%mur_right_CAB1 = right_CAB1; this%mur_right_CAB3 = right_CAB3; this%mur_right_cab4 = right_cab4
       this%mur_down_CAB1 = down_CAB1; this%mur_down_CAB3 = down_CAB3; this%mur_down_cab4 = down_cab4
       this%mur_up_CAB1 = up_CAB1; this%mur_up_CAB3 = up_CAB3; this%mur_up_cab4 = up_cab4
       this%mur_back_CAB1 = back_CAB1; this%mur_back_CAB3 = back_CAB3; this%mur_back_cab4 = back_cab4
       this%mur_front_CAB1 = front_CAB1; this%mur_front_CAB3 = front_CAB3; this%mur_front_cab4 = front_cab4

    end subroutine gpu_update_mur_coeffs

    !--------------------------------------------------------------------------------
    ! Initialize MUR past-field arrays on GPU - called after InitMURBorders
    !--------------------------------------------------------------------------------
    subroutine gpu_init_mur_past_fields(this, left_Hx_nx, left_Hx_ny, left_Hx_nz, left_Hz_nx, left_Hz_ny, left_Hz_nz, right_Hx_nx, right_Hx_ny, right_Hx_nz, right_Hz_nx, right_Hz_ny, right_Hz_nz, down_Hy_nx, down_Hy_ny, down_Hy_nz, down_Hx_nx, down_Hx_ny, down_Hx_nz, up_Hy_nx, up_Hy_ny, up_Hy_nz, up_Hx_nx, up_Hx_ny, up_Hx_nz, back_Hz_nx, back_Hz_ny, back_Hz_nz, back_Hy_nx, back_Hy_ny, back_Hy_nz, front_Hz_nx, front_Hz_ny, front_Hz_nz, front_Hy_nx, front_Hy_ny, front_Hy_nz, left_Hx, left_Hz, right_Hx, right_Hz, down_Hy, down_Hx, up_Hy, up_Hx, back_Hz, back_Hy, front_Hz, front_Hy)
       class(gpu_state_t), intent(inout) :: this
       integer(kind=4), intent(in) :: left_Hx_nx, left_Hx_ny, left_Hx_nz
       integer(kind=4), intent(in) :: left_Hz_nx, left_Hz_ny, left_Hz_nz
       integer(kind=4), intent(in) :: right_Hx_nx, right_Hx_ny, right_Hx_nz
       integer(kind=4), intent(in) :: right_Hz_nx, right_Hz_ny, right_Hz_nz
       integer(kind=4), intent(in) :: down_Hy_nx, down_Hy_ny, down_Hy_nz
       integer(kind=4), intent(in) :: down_Hx_nx, down_Hx_ny, down_Hx_nz
       integer(kind=4), intent(in) :: up_Hy_nx, up_Hy_ny, up_Hy_nz
       integer(kind=4), intent(in) :: up_Hx_nx, up_Hx_ny, up_Hx_nz
       integer(kind=4), intent(in) :: back_Hz_nx, back_Hz_ny, back_Hz_nz
       integer(kind=4), intent(in) :: back_Hy_nx, back_Hy_ny, back_Hy_nz
       integer(kind=4), intent(in) :: front_Hz_nx, front_Hz_ny, front_Hz_nz
       integer(kind=4), intent(in) :: front_Hy_nx, front_Hy_ny, front_Hy_nz
       real(kind=rkind), dimension(:,:,:), intent(in) :: left_Hx, left_Hz, right_Hx, right_Hz
       real(kind=rkind), dimension(:,:,:), intent(in) :: down_Hy, down_Hx, up_Hy, up_Hx
       real(kind=rkind), dimension(:,:,:), intent(in) :: back_Hz, back_Hy, front_Hz, front_Hy

       if (.not. this%initialized) return

       allocate(this%mur_past_Hx_left(left_Hx_nx, left_Hx_ny, left_Hx_nz))
       allocate(this%mur_past_Hz_left(left_Hz_nx, left_Hz_ny, left_Hz_nz))
       allocate(this%mur_past_Hx_right(right_Hx_nx, right_Hx_ny, right_Hx_nz))
       allocate(this%mur_past_Hz_right(right_Hz_nx, right_Hz_ny, right_Hz_nz))
       allocate(this%mur_past_Hy_down(down_Hy_nx, down_Hy_ny, down_Hy_nz))
       allocate(this%mur_past_Hx_down(down_Hx_nx, down_Hx_ny, down_Hx_nz))
       allocate(this%mur_past_Hy_up(up_Hy_nx, up_Hy_ny, up_Hy_nz))
       allocate(this%mur_past_Hx_up(up_Hx_nx, up_Hx_ny, up_Hx_nz))
       allocate(this%mur_past_Hz_back(back_Hz_nx, back_Hz_ny, back_Hz_nz))
       allocate(this%mur_past_Hy_back(back_Hy_nx, back_Hy_ny, back_Hy_nz))
       allocate(this%mur_past_Hz_front(front_Hz_nx, front_Hz_ny, front_Hz_nz))
       allocate(this%mur_past_Hy_front(front_Hy_nx, front_Hy_ny, front_Hy_nz))

       this%mur_past_Hx_left = left_Hx; this%mur_past_Hz_left = left_Hz
       this%mur_past_Hx_right = right_Hx; this%mur_past_Hz_right = right_Hz
       this%mur_past_Hy_down = down_Hy; this%mur_past_Hx_down = down_Hx
       this%mur_past_Hy_up = up_Hy; this%mur_past_Hx_up = up_Hx
       this%mur_past_Hz_back = back_Hz; this%mur_past_Hy_back = back_Hy
       this%mur_past_Hz_front = front_Hz; this%mur_past_Hy_front = front_Hy

       this%mur_initialized = .true.

    end subroutine gpu_init_mur_past_fields

    !--------------------------------------------------------------------------------
    ! Upload MUR past fields to device - called every step after CPU MUR update
    !--------------------------------------------------------------------------------
    subroutine gpu_upload_mur_past_fields(this, left_Hx, left_Hz, right_Hx, right_Hz, down_Hy, down_Hx, up_Hy, up_Hx, back_Hz, back_Hy, front_Hz, front_Hy)
       class(gpu_state_t), intent(inout) :: this
       real(kind=rkind), dimension(:,:,:), intent(in) :: left_Hx, left_Hz, right_Hx, right_Hz, down_Hy, down_Hx, up_Hy, up_Hx, back_Hz, back_Hy, front_Hz, front_Hy

       if (.not. this%mur_initialized) return

       this%mur_past_Hx_left = left_Hx; this%mur_past_Hz_left = left_Hz
       this%mur_past_Hx_right = right_Hx; this%mur_past_Hz_right = right_Hz
       this%mur_past_Hy_down = down_Hy; this%mur_past_Hx_down = down_Hx
       this%mur_past_Hy_up = up_Hy; this%mur_past_Hx_up = up_Hx
       this%mur_past_Hz_back = back_Hz; this%mur_past_Hy_back = back_Hy
       this%mur_past_Hz_front = front_Hz; this%mur_past_Hy_front = front_Hy

    end subroutine gpu_upload_mur_past_fields

     !--------------------------------------------------------------------------------
     ! Initialize MUR boundary limits on GPU
     !--------------------------------------------------------------------------------
     subroutine gpu_init_mur_limits(this, &
                                    left_Hx_ii, left_Hx_ij, left_Hx_ji, left_Hx_jj, left_Hx_ki, left_Hx_kj, &
                                    left_Hz_ii, left_Hz_ij, left_Hz_ji, left_Hz_jj, left_Hz_ki, left_Hz_kj, &
                                    right_Hx_ii, right_Hx_ij, right_Hx_ji, right_Hx_jj, right_Hx_ki, right_Hx_kj, &
                                    right_Hz_ii, right_Hz_ij, right_Hz_ji, right_Hz_jj, right_Hz_ki, right_Hz_kj, &
                                    down_Hy_ii, down_Hy_ij, down_Hy_ji, down_Hy_jj, down_Hy_ki, down_Hy_kj, &
                                    down_Hx_ii, down_Hx_ij, down_Hx_ji, down_Hx_jj, down_Hx_ki, down_Hx_kj, &
                                    up_Hy_ii, up_Hy_ij, up_Hy_ji, up_Hy_jj, up_Hy_ki, up_Hy_kj, &
                                    up_Hx_ii, up_Hx_ij, up_Hx_ji, up_Hx_jj, up_Hx_ki, up_Hx_kj, &
                                    back_Hz_ii, back_Hz_ij, back_Hz_ji, back_Hz_jj, back_Hz_ki, back_Hz_kj, &
                                    back_Hy_ii, back_Hy_ij, back_Hy_ji, back_Hy_jj, back_Hy_ki, back_Hy_kj, &
                                    front_Hz_ii, front_Hz_ij, front_Hz_ji, front_Hz_jj, front_Hz_ki, front_Hz_kj, &
                                    front_Hy_ii, front_Hy_ij, front_Hy_ji, front_Hy_jj, front_Hy_ki, front_Hy_kj)
        class(gpu_state_t), intent(inout) :: this
        integer(kind=4), intent(in) :: left_Hx_ii, left_Hx_ij, left_Hx_ji, left_Hx_jj, left_Hx_ki, left_Hx_kj
        integer(kind=4), intent(in) :: left_Hz_ii, left_Hz_ij, left_Hz_ji, left_Hz_jj, left_Hz_ki, left_Hz_kj
        integer(kind=4), intent(in) :: right_Hx_ii, right_Hx_ij, right_Hx_ji, right_Hx_jj, right_Hx_ki, right_Hx_kj
        integer(kind=4), intent(in) :: right_Hz_ii, right_Hz_ij, right_Hz_ji, right_Hz_jj, right_Hz_ki, right_Hz_kj
        integer(kind=4), intent(in) :: down_Hy_ii, down_Hy_ij, down_Hy_ji, down_Hy_jj, down_Hy_ki, down_Hy_kj
        integer(kind=4), intent(in) :: down_Hx_ii, down_Hx_ij, down_Hx_ji, down_Hx_jj, down_Hx_ki, down_Hx_kj
        integer(kind=4), intent(in) :: up_Hy_ii, up_Hy_ij, up_Hy_ji, up_Hy_jj, up_Hy_ki, up_Hy_kj
        integer(kind=4), intent(in) :: up_Hx_ii, up_Hx_ij, up_Hx_ji, up_Hx_jj, up_Hx_ki, up_Hx_kj
        integer(kind=4), intent(in) :: back_Hz_ii, back_Hz_ij, back_Hz_ji, back_Hz_jj, back_Hz_ki, back_Hz_kj
        integer(kind=4), intent(in) :: back_Hy_ii, back_Hy_ij, back_Hy_ji, back_Hy_jj, back_Hy_ki, back_Hy_kj
        integer(kind=4), intent(in) :: front_Hz_ii, front_Hz_ij, front_Hz_ji, front_Hz_jj, front_Hz_ki, front_Hz_kj
        integer(kind=4), intent(in) :: front_Hy_ii, front_Hy_ij, front_Hy_ji, front_Hy_jj, front_Hy_ki, front_Hy_kj

        if (.not. this%initialized) return

        this%mur_left_Hx_ii = left_Hx_ii; this%mur_left_Hx_ij = left_Hx_ij
        this%mur_left_Hx_ji = left_Hx_ji; this%mur_left_Hx_jj = left_Hx_jj
        this%mur_left_Hx_ki = left_Hx_ki; this%mur_left_Hx_kj = left_Hx_kj
        this%mur_left_Hz_ii = left_Hz_ii; this%mur_left_Hz_ij = left_Hz_ij
        this%mur_left_Hz_ji = left_Hz_ji; this%mur_left_Hz_jj = left_Hz_jj
        this%mur_left_Hz_ki = left_Hz_ki; this%mur_left_Hz_kj = left_Hz_kj
        this%mur_right_Hx_ii = right_Hx_ii; this%mur_right_Hx_ij = right_Hx_ij
        this%mur_right_Hx_ji = right_Hx_ji; this%mur_right_Hx_jj = right_Hx_jj
        this%mur_right_Hx_ki = right_Hx_ki; this%mur_right_Hx_kj = right_Hx_kj
        this%mur_right_Hz_ii = right_Hz_ii; this%mur_right_Hz_ij = right_Hz_ij
        this%mur_right_Hz_ji = right_Hz_ji; this%mur_right_Hz_jj = right_Hz_jj
        this%mur_right_Hz_ki = right_Hz_ki; this%mur_right_Hz_kj = right_Hz_kj
        this%mur_down_Hy_ii = down_Hy_ii; this%mur_down_Hy_ij = down_Hy_ij
        this%mur_down_Hy_ji = down_Hy_ji; this%mur_down_Hy_jj = down_Hy_jj
        this%mur_down_Hy_ki = down_Hy_ki; this%mur_down_Hy_kj = down_Hy_kj
        this%mur_down_Hx_ii = down_Hx_ii; this%mur_down_Hx_ij = down_Hx_ij
        this%mur_down_Hx_ji = down_Hx_ji; this%mur_down_Hx_jj = down_Hx_jj
        this%mur_down_Hx_ki = down_Hx_ki; this%mur_down_Hx_kj = down_Hx_kj
        this%mur_up_Hy_ii = up_Hy_ii; this%mur_up_Hy_ij = up_Hy_ij
        this%mur_up_Hy_ji = up_Hy_ji; this%mur_up_Hy_jj = up_Hy_jj
        this%mur_up_Hy_ki = up_Hy_ki; this%mur_up_Hy_kj = up_Hy_kj
        this%mur_up_Hx_ii = up_Hx_ii; this%mur_up_Hx_ij = up_Hx_ij
        this%mur_up_Hx_ji = up_Hx_ji; this%mur_up_Hx_jj = up_Hx_jj
        this%mur_up_Hx_ki = up_Hx_ki; this%mur_up_Hx_kj = up_Hx_kj
        this%mur_back_Hz_ii = back_Hz_ii; this%mur_back_Hz_ij = back_Hz_ij
        this%mur_back_Hz_ji = back_Hz_ji; this%mur_back_Hz_jj = back_Hz_jj
        this%mur_back_Hz_ki = back_Hz_ki; this%mur_back_Hz_kj = back_Hz_kj
        this%mur_back_Hy_ii = back_Hy_ii; this%mur_back_Hy_ij = back_Hy_ij
        this%mur_back_Hy_ji = back_Hy_ji; this%mur_back_Hy_jj = back_Hy_jj
        this%mur_back_Hy_ki = back_Hy_ki; this%mur_back_Hy_kj = back_Hy_kj
        this%mur_front_Hz_ii = front_Hz_ii; this%mur_front_Hz_ij = front_Hz_ij
        this%mur_front_Hz_ji = front_Hz_ji; this%mur_front_Hz_jj = front_Hz_jj
        this%mur_front_Hz_ki = front_Hz_ki; this%mur_front_Hz_kj = front_Hz_kj
        this%mur_front_Hy_ii = front_Hy_ii; this%mur_front_Hy_ij = front_Hy_ij
        this%mur_front_Hy_ji = front_Hy_ji; this%mur_front_Hy_jj = front_Hy_jj
        this%mur_front_Hy_ki = front_Hy_ki; this%mur_front_Hy_kj = front_Hy_kj

     end subroutine gpu_init_mur_limits

  end module gpu_core_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU CORE PROBE EXTENSION - Probe-aware selective download
!  Downloads only the cells that observation probes reference,
!  eliminating per-timestep full-field downloads.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_core_probe_m

   use FDETYPES_m
   use Report_m
   use cudafor
   use gpu_core_m

   implicit none

contains

   !--------------------------------------------------------------------------------
   ! Download only probe-relevant cells from device to host
   ! Much smaller than downloading all 6 fields
   !--------------------------------------------------------------------------------
   subroutine gpu_download_probes(this, sgg, Ex, Ey, Ez, Hx, Hy, Hz)
      class(gpu_state_t), intent(inout) :: this
      type(SGGFDTDINFO_t), intent(in) :: sgg
      real(kind=rkind), dimension(:,:,:), pointer, intent(inout) :: Ex, Ey, Ez, Hx, Hy, Hz

      integer(kind=4) :: ii, i, field
      integer(kind=4) :: I1, J1, K1, I2, J2, K2
      integer(kind=4) :: iii, jjj, kkk
      integer(kind=4) :: pointObservationCases(6), blockCurrentObservationCases(6)
      integer(kind=4) :: i1_m, i2_m, j1_m, j2_m, k1_m, k2_m
      logical :: is_point, is_block

      pointObservationCases = [iEx, iEy, iEz, iHx, iHy, iHz]
      blockCurrentObservationCases = [iBloqueJx, iBloqueJy, iBloqueJz, iBloqueMx, iBloqueMy, iBloqueMz]

      if (.not. this%initialized) return
      if (.not. this%fields_on_device) return

      ! Check if any probe has invalid bounds (0,0,0) — indicates element-based probe
      ! Element-based probes need full field download since cell locations
      ! are determined by element geometry, not by simple cell ranges
      do ii = 1, sgg%NumberRequest
         do i = 1, sgg%Observation(ii)%nP
            if (sgg%Observation(ii)%P(i)%XI == 0 .and. &
                sgg%Observation(ii)%P(i)%YI == 0 .and. &
                sgg%Observation(ii)%P(i)%ZI == 0) then
               ! Element-based probe — fall back to full download
               this%Ex = this%Ex_d
               this%Ey = this%Ey_d
               this%Ez = this%Ez_d
               this%Hx = this%Hx_d
               this%Hy = this%Hy_d
               this%Hz = this%Hz_d
               this%fields_on_device = .false.
               return
            end if
         end do
      end do

      ! Point probes: download individual cells
      do ii = 1, sgg%NumberRequest
         do i = 1, sgg%Observation(ii)%nP
            field = sgg%Observation(ii)%P(i)%what
            if (field == nothing) cycle

            I1 = sgg%Observation(ii)%P(i)%XI
            J1 = sgg%Observation(ii)%P(i)%YI
            K1 = sgg%Observation(ii)%P(i)%ZI

            is_point = any(field == pointObservationCases)
            is_block = any(field == blockCurrentObservationCases)

            if (is_point) then
               ! Download single cell from each field
               Ex(I1,J1,K1) = this%Ex_d(I1,J1,K1)
               Ey(I1,J1,K1) = this%Ey_d(I1,J1,K1)
               Ez(I1,J1,K1) = this%Ez_d(I1,J1,K1)
               Hx(I1,J1,K1) = this%Hx_d(I1,J1,K1)
               Hy(I1,J1,K1) = this%Hy_d(I1,J1,K1)
               Hz(I1,J1,K1) = this%Hz_d(I1,J1,K1)
            else if (is_block) then
               ! Block probes: download the surface/volume region
               ! Only download if bounds are valid (element-based probes may have 0 bounds)
               i1_m = I1; i2_m = sgg%Observation(ii)%P(i)%XE
               j1_m = J1; j2_m = sgg%Observation(ii)%P(i)%YE
               k1_m = K1; k2_m = sgg%Observation(ii)%P(i)%ZE
               ! Skip if bounds are invalid (element-based probes use elementIds, not cell ranges)
               if (i1_m > 0 .and. i2_m > 0 .and. j1_m > 0 .and. j2_m > 0 .and. k1_m > 0 .and. k2_m > 0) then
                  do kkk = k1_m, k2_m
                     do jjj = j1_m, j2_m
                        do iii = i1_m, i2_m
                           Ex(iii,jjj,kkk) = this%Ex_d(iii,jjj,kkk)
                           Ey(iii,jjj,kkk) = this%Ey_d(iii,jjj,kkk)
                           Ez(iii,jjj,kkk) = this%Ez_d(iii,jjj,kkk)
                           Hx(iii,jjj,kkk) = this%Hx_d(iii,jjj,kkk)
                           Hy(iii,jjj,kkk) = this%Hy_d(iii,jjj,kkk)
                           Hz(iii,jjj,kkk) = this%Hz_d(iii,jjj,kkk)
                        end do
                     end do
                  end do
               end if
            end if
         end do
      end do

   end subroutine gpu_download_probes

end module gpu_core_probe_m
