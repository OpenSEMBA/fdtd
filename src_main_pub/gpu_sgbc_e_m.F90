!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU SGBC E-FIELD KERNELS - CUDA Fortran (CUF)
!  Non-dispersive SGBC E-field advance (YEE + Crank-Nicolson)
!  + Tridiagonal solver kernel (Thomas algorithm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_sgbc_e_m

   use FDETYPES_m
   use Report_m
   use cudafor
   use gpu_sgbc_core_m

   implicit none

contains

   !--------------------------------------------------------------------------------
   ! Advance SGBC E-field - GPU accelerated
   ! Each thread handles one node
   !--------------------------------------------------------------------------------
   subroutine gpu_advance_sgbc_e(this, dt)
      class(gpu_state_sgbc_t), intent(inout) :: this
      real(kind=rkind), intent(in) :: dt

      integer(kind=4) :: n

      if (.not. this%initialized) return

      n = this%numNodes

      call sgbc_e_kernel(this%E_state, this%H_state, this%E_past_state, this%D_state, &
                         this%G1_int, this%G2_int, this%GM1_int, this%GM2_int, &
                         this%a_coef, this%b_coef, this%c_coef, this%rb_coef, this%rh_coef, this%rhm1_coef, &
                         this%tridiag_a1, this%tridiag_b1, this%tridiag_c1, &
                         this%tridiag_an, this%tridiag_bn, this%tridiag_cn, &
                         this%g1_val, this%g2a_val, this%g2b_val, &
                         this%gm2_externo, this%transE, this%transH, this%alignedH, &
                         this%correct_ha, this%correct_hb, this%crank, this%filo_placa, &
                         this%depth_node, this%capa_state, this%delta_state, &
                         this%Efield_val, this%Ha_Plus_val, this%Ha_Minu_val, &
                         this%Hb_Plus_val, this%Hb_Minu_val, this%Hyee_left, this%Hyee_right, &
                         n, dt)

      this%fields_on_device = .true.

   end subroutine gpu_advance_sgbc_e

   !--------------------------------------------------------------------------------
   ! Fused SGBC E-field advance kernel
   ! Each thread handles one node
   !--------------------------------------------------------------------------------
   subroutine sgbc_e_kernel(E_state, H_state, E_past_state, D_state, &
                            G1_int, G2_int, GM1_int, GM2_int, &
                            a_coef, b_coef, c_coef, rb_coef, rh_coef, rhm1_coef, &
                            tridiag_a1, tridiag_b1, tridiag_c1, &
                            tridiag_an, tridiag_bn, tridiag_cn, &
                            g1_val, g2a_val, g2b_val, &
                            gm2_ext, transE, transH, alignedH, &
                            correct_ha, correct_hb, crank, filo_placa, &
                            depth_node, capa_state, delta_state, &
                            Efield_val, Ha_Plus_val, Ha_Minu_val, &
                            Hb_Plus_val, Hb_Minu_val, Hyee_left, Hyee_right, &
                            n, dt)
      integer(kind=4), intent(in) :: n
      real(kind=rkind), intent(in) :: dt
      real(kind=rkind), device, dimension(:,:) :: E_state, H_state, E_past_state, D_state
      real(kind=rkind), device, dimension(:,:) :: G1_int, G2_int, GM1_int, GM2_int
      real(kind=rkind), device, dimension(:,:) :: a_coef, b_coef, c_coef, rb_coef, rh_coef, rhm1_coef
      real(kind=rkind), device, dimension(:) :: tridiag_a1, tridiag_b1, tridiag_c1
      real(kind=rkind), device, dimension(:) :: tridiag_an, tridiag_bn, tridiag_cn
      real(kind=rkind), device, dimension(:) :: g1_val, g2a_val, g2b_val, gm2_ext
      real(kind=rkind), device, dimension(:) :: transE, transH, alignedH
      logical, device, dimension(:) :: correct_ha, correct_hb, crank, filo_placa
      integer(kind=4), device, dimension(:) :: depth_node
      integer(kind=4), device, dimension(:,:) :: capa_state
      real(kind=rkind), device, dimension(:,:) :: delta_state
      real(kind=rkind), device, dimension(:) :: Efield_val, Ha_Plus_val, Ha_Minu_val
      real(kind=rkind), device, dimension(:) :: Hb_Plus_val, Hb_Minu_val
      real(kind=rkind), device, dimension(:) :: Hyee_left, Hyee_right

      integer(kind=4) :: node, i, j, depth, sz
      real(kind=rkind) :: g1, g2, g1i, g2i, gm1, gm2, delta_e
      real(kind=rkind) :: Ha_Plus, Ha_Minu, Hb_Plus, Hb_Minu, Efield, Hyee_l, Hyee_r
      logical :: ha, hb, cr, fp

      !$cuf kernel do(1) <<<*, *>>>
      do node = 1, n
         depth = depth_node(node)
         sz = 2*depth + 1
         ha = correct_ha(node)
         hb = correct_hb(node)
         cr = crank(node)
         fp = filo_placa(node)
         gm2 = gm2_ext(node)
         g1 = g1_val(node)
         g2 = g2a_val(node)
         Efield = Efield_val(node)
         Ha_Plus = Ha_Plus_val(node)
         Ha_Minu = Ha_Minu_val(node)
         Hb_Plus = Hb_Plus_val(node)
         Hb_Minu = Hb_Minu_val(node)
         Hyee_l = Hyee_left(node)
         Hyee_r = Hyee_right(node)

         if (depth == 0) then
            ! depth=0 case: simple boundary update
            if (ha) then
               E_state(1, node) = g1 * E_state(1, node) + &
                  g2 * (Ha_Plus - Ha_Minu) - g2b_val(node) * (Hb_Plus - Hb_Minu)
            elseif (hb) then
               E_state(1, node) = g1 * E_state(1, node) + &
                  g2 * (Ha_Plus - Ha_Minu) - g2b_val(node) * (Hb_Plus - Hb_Minu)
            else
               E_state(1, node) = g1 * E_state(1, node) + &
                  g2 * (Ha_Plus - Ha_Minu) - g2b_val(node) * (Hb_Plus - Hb_Minu)
            end if
         else
            ! depth>0 case: boundary + interior
            ! Boundary cells: E(depth) and E(-depth)
            if (.not. cr) then
               ! YEE mode
               if (ha) then
                  E_state(depth + 1, node) = g1_val(1) * E_state(depth + 1, node) + &
                     (g2a_val(1) * (Ha_Plus - Hyee_r) - g2b_val(1) * (Hb_Plus - Hb_Minu))
                  E_state(1, node) = g1_val(0) * E_state(1, node) + &
                     (g2a_val(0) * (Hyee_l - Ha_Minu) - g2b_val(0) * (Hb_Plus - Hb_Minu))
               elseif (hb) then
                  E_state(depth + 1, node) = g1_val(1) * E_state(depth + 1, node) + &
                     (g2a_val(1) * (Ha_Plus - Ha_Minu) - g2b_val(1) * (Hb_Plus - Hyee_r))
                  E_state(1, node) = g1_val(0) * E_state(1, node) + &
                     (g2a_val(0) * (Ha_Plus - Ha_Minu) - g2b_val(0) * (Hyee_l - Hb_Minu))
               end if

               ! Interior cells
               do i = 2, depth
                  g1i = G1_int(i, node)
                  g2i = G2_int(i, node)
                  delta_e = 0.5_rkind * (delta_state(i, node) + delta_state(i-1, node))
                  E_state(i, node) = g1i * E_state(i, node) + g2i / delta_e * &
                     (H_state(i, node) - H_state(i-1, node))
               end do
            else
               ! Crank-Nicolson mode
               ! Copy E to E_past
               do i = 1, sz
                  E_past_state(i, node) = E_state(i, node)
               end do

               ! Compute D vector for boundaries
               if (ha) then
                  i = depth + 1
                  D_state(i, node) = -tridiag_an(node) * E_past_state(i-1, node) + &
                     tridiag_bn(node) * E_past_state(i, node) + &
                     g2a_val(1) * (Ha_Plus - Hyee_r) - g2b_val(1) * (Hb_Plus - Hb_Minu)
                  i = 1
                  D_state(i, node) = -tridiag_c1(node) * E_past_state(i+1, node) + &
                     tridiag_bn(node) * E_past_state(i, node) + &
                     g2a_val(0) * (Hyee_l - Ha_Minu) - g2b_val(0) * (Hb_Plus - Hb_Minu)
               elseif (hb) then
                  i = depth + 1
                  D_state(i, node) = -tridiag_an(node) * E_past_state(i-1, node) + &
                     tridiag_bn(node) * E_past_state(i, node) + &
                     g2a_val(1) * (Ha_Plus - Ha_Minu) - g2b_val(1) * (Hb_Plus - Hyee_r)
                  i = 1
                  D_state(i, node) = -tridiag_c1(node) * E_past_state(i+1, node) + &
                     tridiag_bn(node) * E_past_state(i, node) + &
                     g2a_val(0) * (Ha_Plus - Ha_Minu) - g2b_val(0) * (Hyee_l - Hb_Minu)
               end if

               ! Compute D vector for interior cells
               do i = 2, depth
                  D_state(i, node) = -a_coef(i, node) * E_past_state(i-1, node) - &
                     c_coef(i, node) * E_past_state(i+1, node) + &
                     rb_coef(i, node) * E_past_state(i, node) + &
                     rh_coef(i, node) * H_state(i, node) - &
                     rhm1_coef(i, node) * H_state(i-1, node)
               end do

               ! Solve tridiagonal system
               call sgbc_solve_tridiag(D_state, E_state, a_coef, b_coef, c_coef, &
                                       tridiag_a1(node), tridiag_b1(node), tridiag_c1(node), &
                                       tridiag_an(node), tridiag_bn(node), tridiag_cn(node), &
                                       sz, node)
            end if
         end if

         ! Update Efield = average of boundary cells
         Efield_val(node) = 0.5_rkind * (E_state(1, node) + E_state(depth + 1, node))
      end do

   end subroutine sgbc_e_kernel

   !--------------------------------------------------------------------------------
   ! Tridiagonal solver (Thomas algorithm) - GPU kernel
   ! Each thread solves one node's system
   !--------------------------------------------------------------------------------
   subroutine sgbc_solve_tridiag(d, x, a, b, c, a1, b1, c1, an, bn, cn, n, node)
      integer(kind=4), intent(in) :: n, node
      real(kind=rkind), device, dimension(:,:), intent(inout) :: d, x
      real(kind=rkind), device, dimension(:,:), intent(in) :: a, b, c
      real(kind=rkind), intent(in) :: a1, b1, c1, an, bn, cn

      real(kind=rkind), dimension(:), allocatable, device :: aa, bb, cc
      real(kind=rkind), dimension(:), allocatable, device :: cp, dp
      real(kind=rkind) :: m
      integer(kind=4) :: i

      allocate(aa(n), bb(n), cc(n))
      allocate(cp(n), dp(n))

      aa(1) = a1
      bb(1) = b1
      cc(1) = c1
      aa(n) = an
      bb(n) = bn
      cc(n) = cn
      aa(2:n-1) = a(2:n-1, node)
      bb(2:n-1) = b(2:n-1, node)
      cc(2:n-1) = c(2:n-1, node)

      ! Forward elimination
      cp(1) = cc(1) / bb(1)
      dp(1) = d(1, node) / bb(1)
      do i = 2, n
         m = bb(i) - cp(i-1) * aa(i)
         cp(i) = cc(i) / m
         dp(i) = (d(i, node) - dp(i-1) * aa(i)) / m
      end do

      ! Back substitution
      x(n, node) = dp(n)
      do i = n-1, 1, -1
         x(i, node) = dp(i) - cp(i) * x(i+1, node)
      end do

      deallocate(aa, bb, cc, cp, dp)

   end subroutine sgbc_solve_tridiag

end module gpu_sgbc_e_m
