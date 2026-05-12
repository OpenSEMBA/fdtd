!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU SGBC CORE MODULE - CUDA Fortran (CUF) accelerated SGBC
!  Flattened arrays for all SGBC node data (no pointer indirection)
!  Handles init/upload/download/destroy for non-dispersive SGBC only.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_sgbc_core_m

   use FDETYPES_m
   use Report_m
   use cudafor

   implicit none

   integer, parameter :: SGBC_MAX_DEPTH = 16

   type gpu_state_sgbc_t
      integer(kind=4) :: numNodes, maxDepth

      ! Per-node constants (contiguous 1D arrays, indexed by node)
      integer(kind=4), device, allocatable, dimension(:) :: d_node
      real(kind=rkind), device, allocatable, dimension(:) :: gm2_ext
      integer(kind=4), device, allocatable, dimension(:) :: jmed_node
      real(kind=rkind), device, allocatable, dimension(:) :: transE, transH, alignedH
      real(kind=rkind), device, allocatable, dimension(:) :: g1_val, g2a_val, g2b_val
      real(kind=rkind), device, allocatable, dimension(:) :: gm2_externo
      logical, device, allocatable, dimension(:) :: correct_ha, correct_hb, crank, filo_placa
      integer(kind=4), device, allocatable, dimension(:) :: depth_node

      ! Per-node state (contiguous 2D: [2*maxDepth+1, numNodes])
      real(kind=rkind), device, allocatable, dimension(:,:) :: E_state, H_state, E_past_state
      real(kind=rkind), device, allocatable, dimension(:,:) :: D_state
      real(kind=rkind), device, allocatable, dimension(:,:) :: G1_int, G2_int, GM1_int, GM2_int
      real(kind=rkind), device, allocatable, dimension(:,:) :: a_coef, b_coef, c_coef
      real(kind=rkind), device, allocatable, dimension(:,:) :: rb_coef, rh_coef, rhm1_coef
      integer(kind=4), device, allocatable, dimension(:,:) :: capa_state
      real(kind=rkind), device, allocatable, dimension(:,:) :: delta_state

      ! Per-node tridiag boundary constants
      real(kind=rkind), device, allocatable, dimension(:) :: tridiag_a1, tridiag_b1, tridiag_c1
      real(kind=rkind), device, allocatable, dimension(:) :: tridiag_an, tridiag_bn, tridiag_cn

      ! Per-timestep field values (device, synced each step)
      real(kind=rkind), device, allocatable, dimension(:) :: Efield_val, Ha_Plus_val, Ha_Minu_val
      real(kind=rkind), device, allocatable, dimension(:) :: Hb_Plus_val, Hb_Minu_val
      real(kind=rkind), device, allocatable, dimension(:) :: Hyee_left, Hyee_right

      ! Host-side references for field sync
      real(kind=rkind), pointer, dimension(:,:,:) :: Ex, Ey, Ez, Hx, Hy, Hz
      integer(kind=4), pointer, dimension(:,:,:) :: Idxe, Idye, Idze, Idxh, Idyh, Idzh

      logical :: initialized = .false.
      logical :: fields_on_device = .false.
   end type

contains

   !--------------------------------------------------------------------------------
   ! Initialize GPU SGBC state - called once at SGBC init
   !--------------------------------------------------------------------------------
   subroutine gpu_init_sgbc(this, numNodes, maxDepth, dt, sgbcFreq, SGBCcrank, SGBCDispersive)
      class(gpu_state_sgbc_t), intent(inout) :: this
      integer(kind=4), intent(in) :: numNodes, maxDepth
      real(kind=rkind), intent(in) :: dt, sgbcFreq
      logical, intent(in) :: SGBCcrank, SGBCDispersive

      integer(kind=4) :: ndev, cuda_status, env_status
      character(len=16) :: enable_cuf
      integer(kind=4) :: i

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

      this%numNodes = numNodes
      this%maxDepth = maxDepth

      ! Allocate per-node constant arrays
      allocate(this%d_node(numNodes))
      allocate(this%gm2_ext(numNodes))
      allocate(this%jmed_node(numNodes))
      allocate(this%transE(numNodes))
      allocate(this%transH(numNodes))
      allocate(this%alignedH(numNodes))
      allocate(this%g1_val(numNodes))
      allocate(this%g2a_val(numNodes))
      allocate(this%g2b_val(numNodes))
      allocate(this%gm2_externo(numNodes))
      allocate(this%correct_ha(numNodes))
      allocate(this%correct_hb(numNodes))
      allocate(this%crank(numNodes))
      allocate(this%filo_placa(numNodes))
      allocate(this%depth_node(numNodes))

      ! Allocate per-node state arrays [2*maxDepth+1, numNodes]
      allocate(this%E_state(2*maxDepth+1, numNodes))
      allocate(this%H_state(2*maxDepth+1, numNodes))
      allocate(this%E_past_state(2*maxDepth+1, numNodes))
      allocate(this%D_state(2*maxDepth+1, numNodes))
      allocate(this%G1_int(2*maxDepth+1, numNodes))
      allocate(this%G2_int(2*maxDepth+1, numNodes))
      allocate(this%GM1_int(2*maxDepth+1, numNodes))
      allocate(this%GM2_int(2*maxDepth+1, numNodes))
      allocate(this%a_coef(2*maxDepth+1, numNodes))
      allocate(this%b_coef(2*maxDepth+1, numNodes))
      allocate(this%c_coef(2*maxDepth+1, numNodes))
      allocate(this%rb_coef(2*maxDepth+1, numNodes))
      allocate(this%rh_coef(2*maxDepth+1, numNodes))
      allocate(this%rhm1_coef(2*maxDepth+1, numNodes))
      allocate(this%capa_state(2*maxDepth+1, numNodes))
      allocate(this%delta_state(2*maxDepth+1, numNodes))

      ! Allocate tridiag boundary constants
      allocate(this%tridiag_a1(numNodes))
      allocate(this%tridiag_b1(numNodes))
      allocate(this%tridiag_c1(numNodes))
      allocate(this%tridiag_an(numNodes))
      allocate(this%tridiag_bn(numNodes))
      allocate(this%tridiag_cn(numNodes))

      ! Allocate per-timestep field arrays
      allocate(this%Efield_val(numNodes))
      allocate(this%Ha_Plus_val(numNodes))
      allocate(this%Ha_Minu_val(numNodes))
      allocate(this%Hb_Plus_val(numNodes))
      allocate(this%Hb_Minu_val(numNodes))
      allocate(this%Hyee_left(numNodes))
      allocate(this%Hyee_right(numNodes))

      ! Initialize to zero
      this%E_state = 0.0_rkind
      this%H_state = 0.0_rkind
      this%E_past_state = 0.0_rkind
      this%D_state = 0.0_rkind
      this%G1_int = 0.0_rkind
      this%G2_int = 0.0_rkind
      this%GM1_int = 0.0_rkind
      this%GM2_int = 0.0_rkind
      this%a_coef = 0.0_rkind
      this%b_coef = 0.0_rkind
      this%c_coef = 0.0_rkind
      this%rb_coef = 0.0_rkind
      this%rh_coef = 0.0_rkind
      this%rhm1_coef = 0.0_rkind
      this%capa_state = 0
      this%delta_state = 0.0_rkind

      this%initialized = .true.
      this%fields_on_device = .false.

   end subroutine gpu_init_sgbc

   !--------------------------------------------------------------------------------
   ! Upload per-node constants from host to device
   ! Called once after InitSGBCs completes
   !--------------------------------------------------------------------------------
   subroutine gpu_upload_sgbc_constants(this, depth_arr, gm2_ext_arr, jmed_arr, &
                                        transE_arr, transH_arr, alignedH_arr, &
                                        g1_arr, g2a_arr, g2b_arr, gm2_ext_val_arr, &
                                        correct_ha_arr, correct_hb_arr, crank_arr, &
                                        filo_placa_arr, depth_node_arr)
      class(gpu_state_sgbc_t), intent(inout) :: this
      integer(kind=4), intent(in) :: depth_arr(:), jmed_arr(:), depth_node_arr(:)
      real(kind=rkind), intent(in) :: gm2_ext_arr(:), transE_arr(:), transH_arr(:), alignedH_arr(:)
      real(kind=rkind), intent(in) :: g1_arr(:), g2a_arr(:), g2b_arr(:), gm2_ext_val_arr(:)
      logical, intent(in) :: correct_ha_arr(:), correct_hb_arr(:), crank_arr(:), filo_placa_arr(:)

      integer(kind=4) :: n

      if (.not. this%initialized) return

      n = this%numNodes
      this%d_node(1:n) = depth_arr(1:n)
      this%gm2_ext(1:n) = gm2_ext_arr(1:n)
      this%jmed_node(1:n) = jmed_arr(1:n)
      this%transE(1:n) = transE_arr(1:n)
      this%transH(1:n) = transH_arr(1:n)
      this%alignedH(1:n) = alignedH_arr(1:n)
      this%g1_val(1:n) = g1_arr(1:n)
      this%g2a_val(1:n) = g2a_arr(1:n)
      this%g2b_val(1:n) = g2b_arr(1:n)
      this%gm2_externo(1:n) = gm2_ext_val_arr(1:n)
      this%correct_ha(1:n) = correct_ha_arr(1:n)
      this%correct_hb(1:n) = correct_hb_arr(1:n)
      this%crank(1:n) = crank_arr(1:n)
      this%filo_placa(1:n) = filo_placa_arr(1:n)
      this%depth_node(1:n) = depth_node_arr(1:n)

   end subroutine gpu_upload_sgbc_constants

   !--------------------------------------------------------------------------------
   ! Upload per-node state arrays from host to device
   ! Called once after InitSGBCs completes
   !--------------------------------------------------------------------------------
   subroutine gpu_upload_sgbc_state(this, E_arr, H_arr, E_past_arr, D_arr, &
                                     G1_int_arr, G2_int_arr, GM1_int_arr, GM2_int_arr, &
                                     a_arr, b_arr, c_arr, rb_arr, rh_arr, rhm1_arr, &
                                     capa_arr, delta_arr, &
                                     tridiag_a1_arr, tridiag_b1_arr, tridiag_c1_arr, &
                                     tridiag_an_arr, tridiag_bn_arr, tridiag_cn_arr, &
                                     Hyee_left_arr, Hyee_right_arr, &
                                     offset)
      class(gpu_state_sgbc_t), intent(inout) :: this
      real(kind=rkind), intent(in) :: E_arr(:,:), H_arr(:,:), E_past_arr(:,:), D_arr(:,:)
      real(kind=rkind), intent(in) :: G1_int_arr(:,:), G2_int_arr(:,:), GM1_int_arr(:,:), GM2_int_arr(:,:)
      real(kind=rkind), intent(in) :: a_arr(:,:), b_arr(:,:), c_arr(:,:), rb_arr(:,:), rh_arr(:,:), rhm1_arr(:,:)
      integer(kind=4), intent(in) :: capa_arr(:,:)
      real(kind=rkind), intent(in) :: delta_arr(:,:)
      real(kind=rkind), intent(in) :: tridiag_a1_arr(:), tridiag_b1_arr(:), tridiag_c1_arr(:)
      real(kind=rkind), intent(in) :: tridiag_an_arr(:), tridiag_bn_arr(:), tridiag_cn_arr(:)
      real(kind=rkind), intent(in) :: Hyee_left_arr(:), Hyee_right_arr(:)
      integer(kind=4), intent(in) :: offset  ! array offset (-depth)

      integer(kind=4) :: n, sz, i, j

      if (.not. this%initialized) return

      n = this%numNodes
      sz = 2*this%maxDepth + 1

      do j = 1, n
         do i = 1, sz
            this%E_state(i, j) = E_arr(i + offset, j)
            this%H_state(i, j) = H_arr(i + offset, j)
            this%E_past_state(i, j) = E_past_arr(i + offset, j)
            this%D_state(i, j) = D_arr(i + offset, j)
            this%G1_int(i, j) = G1_int_arr(i + offset, j)
            this%G2_int(i, j) = G2_int_arr(i + offset, j)
            this%GM1_int(i, j) = GM1_int_arr(i + offset, j)
            this%GM2_int(i, j) = GM2_int_arr(i + offset, j)
            this%a_coef(i, j) = a_arr(i + offset, j)
            this%b_coef(i, j) = b_arr(i + offset, j)
            this%c_coef(i, j) = c_arr(i + offset, j)
            this%rb_coef(i, j) = rb_arr(i + offset, j)
            this%rh_coef(i, j) = rh_arr(i + offset, j)
            this%rhm1_coef(i, j) = rhm1_arr(i + offset, j)
            this%capa_state(i, j) = capa_arr(i + offset, j)
            this%delta_state(i, j) = delta_arr(i + offset, j)
         end do
         this%tridiag_a1(j) = tridiag_a1_arr(j)
         this%tridiag_b1(j) = tridiag_b1_arr(j)
         this%tridiag_c1(j) = tridiag_c1_arr(j)
         this%tridiag_an(j) = tridiag_an_arr(j)
         this%tridiag_bn(j) = tridiag_bn_arr(j)
         this%tridiag_cn(j) = tridiag_cn_arr(j)
         this%Hyee_left(j) = Hyee_left_arr(j)
         this%Hyee_right(j) = Hyee_right_arr(j)
      end do

      this%fields_on_device = .true.

   end subroutine gpu_upload_sgbc_state

   !--------------------------------------------------------------------------------
   ! Upload per-timestep field values from host to device
   ! Called before each SGBC E-advance
   !--------------------------------------------------------------------------------
   subroutine gpu_upload_sgbc_fields(this, Efield_arr, Ha_Plus_arr, Ha_Minu_arr, &
                                      Hb_Plus_arr, Hb_Minu_arr)
      class(gpu_state_sgbc_t), intent(inout) :: this
      real(kind=rkind), intent(in) :: Efield_arr(:), Ha_Plus_arr(:), Ha_Minu_arr(:)
      real(kind=rkind), intent(in) :: Hb_Plus_arr(:), Hb_Minu_arr(:)

      integer(kind=4) :: n

      if (.not. this%initialized) return

      n = this%numNodes
      this%Efield_val(1:n) = Efield_arr(1:n)
      this%Ha_Plus_val(1:n) = Ha_Plus_arr(1:n)
      this%Ha_Minu_val(1:n) = Ha_Minu_arr(1:n)
      this%Hb_Plus_val(1:n) = Hb_Plus_arr(1:n)
      this%Hb_Minu_val(1:n) = Hb_Minu_arr(1:n)

   end subroutine gpu_upload_sgbc_fields

   !--------------------------------------------------------------------------------
   ! Download per-timestep field values from device to host
   ! Called after each SGBC H-advance
   !--------------------------------------------------------------------------------
   subroutine gpu_download_sgbc_fields(this, Ha_Plus_arr, Ha_Minu_arr, &
                                        Hb_Plus_arr, Hb_Minu_arr, &
                                        Hyee_left_arr, Hyee_right_arr, &
                                        Efield_arr)
      class(gpu_state_sgbc_t), intent(inout) :: this
      real(kind=rkind), intent(out) :: Ha_Plus_arr(:), Ha_Minu_arr(:)
      real(kind=rkind), intent(out) :: Hb_Plus_arr(:), Hb_Minu_arr(:)
      real(kind=rkind), intent(out) :: Hyee_left_arr(:), Hyee_right_arr(:)
      real(kind=rkind), intent(out) :: Efield_arr(:)

      integer(kind=4) :: n

      if (.not. this%initialized) return

      n = this%numNodes
      Ha_Plus_arr(1:n) = this%Ha_Plus_val(1:n)
      Ha_Minu_arr(1:n) = this%Ha_Minu_val(1:n)
      Hb_Plus_arr(1:n) = this%Hb_Plus_val(1:n)
      Hb_Minu_arr(1:n) = this%Hb_Minu_val(1:n)
      Hyee_left_arr(1:n) = this%Hyee_left(1:n)
      Hyee_right_arr(1:n) = this%Hyee_right(1:n)
      Efield_arr(1:n) = this%Efield_val(1:n)

   end subroutine gpu_download_sgbc_fields

   !--------------------------------------------------------------------------------
   ! Download per-node state from device to host
   ! Called after each SGBC E-advance (for H-update) and at output times
   !--------------------------------------------------------------------------------
   subroutine gpu_download_sgbc_state(this, E_arr, H_arr, E_past_arr, D_arr, &
                                       Hyee_left_arr, Hyee_right_arr, &
                                       offset)
      class(gpu_state_sgbc_t), intent(inout) :: this
      real(kind=rkind), intent(out) :: E_arr(:,:), H_arr(:,:), E_past_arr(:,:), D_arr(:,:)
      real(kind=rkind), intent(out) :: Hyee_left_arr(:), Hyee_right_arr(:)
      integer(kind=4), intent(in) :: offset  ! array offset (-depth)

      integer(kind=4) :: n, sz, i, j

      if (.not. this%initialized) return

      n = this%numNodes
      sz = 2*this%maxDepth + 1

      do j = 1, n
         do i = 1, sz
            E_arr(i + offset, j) = this%E_state(i, j)
            H_arr(i + offset, j) = this%H_state(i, j)
            E_past_arr(i + offset, j) = this%E_past_state(i, j)
            D_arr(i + offset, j) = this%D_state(i, j)
         end do
         Hyee_left_arr(j) = this%Hyee_left(j)
         Hyee_right_arr(j) = this%Hyee_right(j)
      end do

      this%fields_on_device = .false.

   end subroutine gpu_download_sgbc_state

   !--------------------------------------------------------------------------------
   ! Upload per-node coefficients from host to device
   ! Called each timestep after calc_SGBCconstants
   !--------------------------------------------------------------------------------
   subroutine gpu_upload_sgbc_coeffs(this, G1_int_arr, G2_int_arr, GM1_int_arr, GM2_int_arr, &
                                      a_arr, b_arr, c_arr, rb_arr, rh_arr, rhm1_arr, &
                                      tridiag_a1_arr, tridiag_b1_arr, tridiag_c1_arr, &
                                      tridiag_an_arr, tridiag_bn_arr, tridiag_cn_arr, &
                                      offset)
      class(gpu_state_sgbc_t), intent(inout) :: this
      real(kind=rkind), intent(in) :: G1_int_arr(:,:), G2_int_arr(:,:), GM1_int_arr(:,:), GM2_int_arr(:,:)
      real(kind=rkind), intent(in) :: a_arr(:,:), b_arr(:,:), c_arr(:,:), rb_arr(:,:), rh_arr(:,:), rhm1_arr(:,:)
      real(kind=rkind), intent(in) :: tridiag_a1_arr(:), tridiag_b1_arr(:), tridiag_c1_arr(:)
      real(kind=rkind), intent(in) :: tridiag_an_arr(:), tridiag_bn_arr(:), tridiag_cn_arr(:)
      integer(kind=4), intent(in) :: offset

      integer(kind=4) :: n, sz, i, j

      if (.not. this%initialized) return

      n = this%numNodes
      sz = 2*this%maxDepth + 1

      do j = 1, n
         do i = 1, sz
            this%G1_int(i, j) = G1_int_arr(i + offset, j)
            this%G2_int(i, j) = G2_int_arr(i + offset, j)
            this%GM1_int(i, j) = GM1_int_arr(i + offset, j)
            this%GM2_int(i, j) = GM2_int_arr(i + offset, j)
            this%a_coef(i, j) = a_arr(i + offset, j)
            this%b_coef(i, j) = b_arr(i + offset, j)
            this%c_coef(i, j) = c_arr(i + offset, j)
            this%rb_coef(i, j) = rb_arr(i + offset, j)
            this%rh_coef(i, j) = rh_arr(i + offset, j)
            this%rhm1_coef(i, j) = rhm1_arr(i + offset, j)
         end do
         this%tridiag_a1(j) = tridiag_a1_arr(j)
         this%tridiag_b1(j) = tridiag_b1_arr(j)
         this%tridiag_c1(j) = tridiag_c1_arr(j)
         this%tridiag_an(j) = tridiag_an_arr(j)
         this%tridiag_bn(j) = tridiag_bn_arr(j)
         this%tridiag_cn(j) = tridiag_cn_arr(j)
      end do

   end subroutine gpu_upload_sgbc_coeffs

   !--------------------------------------------------------------------------------
   ! Destroy GPU SGBC state
   !--------------------------------------------------------------------------------
   subroutine gpu_destroy_sgbc(this)
      class(gpu_state_sgbc_t), intent(inout) :: this

      if (.not. this%initialized) return

      ! Deallocate device arrays
      if (this%initialized) deallocate(this%d_node)
      if (this%initialized) deallocate(this%gm2_ext)
      if (this%initialized) deallocate(this%jmed_node)
      if (this%initialized) deallocate(this%transE)
      if (this%initialized) deallocate(this%transH)
      if (this%initialized) deallocate(this%alignedH)
      if (this%initialized) deallocate(this%g1_val)
      if (this%initialized) deallocate(this%g2a_val)
      if (this%initialized) deallocate(this%g2b_val)
      if (this%initialized) deallocate(this%gm2_externo)
      if (this%initialized) deallocate(this%correct_ha)
      if (this%initialized) deallocate(this%correct_hb)
      if (this%initialized) deallocate(this%crank)
      if (this%initialized) deallocate(this%filo_placa)
      if (this%initialized) deallocate(this%depth_node)

      if (this%initialized) deallocate(this%E_state)
      if (this%initialized) deallocate(this%H_state)
      if (this%initialized) deallocate(this%E_past_state)
      if (this%initialized) deallocate(this%D_state)
      if (this%initialized) deallocate(this%G1_int)
      if (this%initialized) deallocate(this%G2_int)
      if (this%initialized) deallocate(this%GM1_int)
      if (this%initialized) deallocate(this%GM2_int)
      if (this%initialized) deallocate(this%a_coef)
      if (this%initialized) deallocate(this%b_coef)
      if (this%initialized) deallocate(this%c_coef)
      if (this%initialized) deallocate(this%rb_coef)
      if (this%initialized) deallocate(this%rh_coef)
      if (this%initialized) deallocate(this%rhm1_coef)
      if (this%initialized) deallocate(this%capa_state)
      if (this%initialized) deallocate(this%delta_state)

      if (this%initialized) deallocate(this%tridiag_a1)
      if (this%initialized) deallocate(this%tridiag_b1)
      if (this%initialized) deallocate(this%tridiag_c1)
      if (this%initialized) deallocate(this%tridiag_an)
      if (this%initialized) deallocate(this%tridiag_bn)
      if (this%initialized) deallocate(this%tridiag_cn)

      if (this%initialized) deallocate(this%Efield_val)
      if (this%initialized) deallocate(this%Ha_Plus_val)
      if (this%initialized) deallocate(this%Ha_Minu_val)
      if (this%initialized) deallocate(this%Hb_Plus_val)
      if (this%initialized) deallocate(this%Hb_Minu_val)
      if (this%initialized) deallocate(this%Hyee_left)
      if (this%initialized) deallocate(this%Hyee_right)

      nullify(this%Ex); nullify(this%Ey); nullify(this%Ez)
      nullify(this%Hx); nullify(this%Hy); nullify(this%Hz)
      nullify(this%Idxe); nullify(this%Idye); nullify(this%Idze)
      nullify(this%Idxh); nullify(this%Idyh); nullify(this%Idzh)

      this%initialized = .false.
      this%fields_on_device = .false.

   end subroutine gpu_destroy_sgbc

end module gpu_sgbc_core_m
