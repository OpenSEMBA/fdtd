!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU SGBC H-FIELD KERNELS - CUDA Fortran (CUF)
!  Non-dispersive SGBC H-field advance + boundary correction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_sgbc_h_m

   use FDETYPES_m
   use Report_m
   use cudafor
   use gpu_sgbc_core_m

   implicit none

contains

   !--------------------------------------------------------------------------------
   ! Advance SGBC H-field - GPU accelerated
   ! Each thread handles one node
   !--------------------------------------------------------------------------------
   subroutine gpu_advance_sgbc_h(this, dt)
      class(gpu_state_sgbc_t), intent(inout) :: this
      real(kind=rkind), intent(in) :: dt

      integer(kind=4) :: n

      if (.not. this%initialized) return

      n = this%numNodes

      call sgbc_h_kernel(this%H_state, this%E_state, this%E_past_state, &
                         this%GM1_int, this%GM2_int, &
                         this%crank, this%depth_node, &
                         this%delta_state, &
                         this%Ha_Plus_val, this%Ha_Minu_val, &
                         this%Hb_Plus_val, this%Hb_Minu_val, &
                         this%Efield_val, this%Hyee_left, this%Hyee_right, &
                         this%gm2_ext, this%correct_ha, this%correct_hb, &
                         n, dt)

   end subroutine gpu_advance_sgbc_h

   !--------------------------------------------------------------------------------
   ! Fused SGBC H-field advance + boundary correction kernel
   ! Each thread handles one node
   !--------------------------------------------------------------------------------
   subroutine sgbc_h_kernel(H_state, E_state, E_past_state, &
                            GM1_int, GM2_int, &
                            crank, depth_node, &
                            delta_state, &
                            Ha_Plus_val, Ha_Minu_val, &
                            Hb_Plus_val, Hb_Minu_val, &
                            Efield_val, Hyee_left, Hyee_right, &
                            gm2_ext, correct_ha, correct_hb, &
                            n, dt)
      integer(kind=4), intent(in) :: n
      real(kind=rkind), intent(in) :: dt
      real(kind=rkind), device, dimension(:,:) :: H_state, E_state, E_past_state
      real(kind=rkind), device, dimension(:,:) :: GM1_int, GM2_int
      integer(kind=4), device, dimension(:,:) :: delta_state
      logical, device, dimension(:) :: crank, correct_ha, correct_hb
      integer(kind=4), device, dimension(:) :: depth_node
      real(kind=rkind), device, dimension(:) :: Ha_Plus_val, Ha_Minu_val
      real(kind=rkind), device, dimension(:) :: Hb_Plus_val, Hb_Minu_val
      real(kind=rkind), device, dimension(:) :: Efield_val, Hyee_left, Hyee_right
      real(kind=rkind), device, dimension(:) :: gm2_ext

      integer(kind=4) :: node, i, depth, sz
      real(kind=rkind) :: gm1, gm2, delta_e
      real(kind=rkind) :: Ha_Plus, Ha_Minu, Hb_Plus, Hb_Minu, Efield
      real(kind=rkind) :: E_i1, E_i, E_past_i1, E_past_i
      logical :: cr, ha, hb

      !$cuf kernel do(1) <<<*, *>>>
      do node = 1, n
         depth = depth_node(node)
         sz = 2*depth + 1
         cr = crank(node)
         ha = correct_ha(node)
         hb = correct_hb(node)

         if (depth > 0) then
            ! Update internal H cells
            if (cr) then
               ! Crank-Nicolson: half-step advance
               do i = 1, depth
                  gm1 = GM1_int(i, node)
                  gm2 = GM2_int(i, node)
                  E_i1 = E_state(i+1, node)
                  E_i = E_state(i, node)
                  E_past_i1 = E_past_state(i+1, node)
                  E_past_i = E_past_state(i, node)
                  H_state(i, node) = gm1 * H_state(i, node) + &
                     0.5_rkind * gm2 * (E_i1 - E_i + E_past_i1 - E_past_i)
               end do
               ! Update Hyee values
               H_state(1, node) = Hyee_left(node)
               H_state(depth, node) = Hyee_right(node)
            else
               ! YEE mode
               do i = 1, depth
                  gm1 = GM1_int(i, node)
                  gm2 = GM2_int(i, node)
                  E_i1 = E_state(i+1, node)
                  E_i = E_state(i, node)
                  H_state(i, node) = gm1 * H_state(i, node) + gm2 * (E_i1 - E_i)
               end do
               Hyee_left(node) = H_state(1, node)
               Hyee_right(node) = H_state(depth, node)
            end if
         end if

         ! Boundary correction (AdvanceSGBCH)
         Ha_Plus = Ha_Plus_val(node)
         Ha_Minu = Ha_Minu_val(node)
         Hb_Plus = Hb_Plus_val(node)
         Hb_Minu = Hb_Minu_val(node)
         Efield = Efield_val(node)

         if (ha) then
            Ha_Plus_val(node) = Ha_Plus + gm2_ext(node) * (Efield - E_state(depth + 1, node))
            Ha_Minu_val(node) = Ha_Minu - gm2_ext(node) * (Efield - E_state(1, node))
         elseif (hb) then
            Hb_Plus_val(node) = Hb_Plus - gm2_ext(node) * (Efield - E_state(depth + 1, node))
            Hb_Minu_val(node) = Hb_Minu + gm2_ext(node) * (Efield - E_state(1, node))
         end if

         ! Update Efield to average
         Efield_val(node) = 0.5_rkind * (E_state(1, node) + E_state(depth + 1, node))
      end do

   end subroutine sgbc_h_kernel

end module gpu_sgbc_h_m
