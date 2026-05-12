!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU CPML KERNELS MODULE - CUDA Fortran (CUF) accelerated CPML kernels
!  Left + Right boundary CPML with persistent psi state on device.
!  Split from gpu_kernels_cuf.F90 to avoid NVHPC compiler file-size limit.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_cpml_m

   use FDETYPES_m
   use Report_m
   use cudafor
   use gpu_core_m

   implicit none

contains

   !--------------------------------------------------------------------------------
   ! CPML Left Boundary - GPU accelerated kernels
   ! Uses persistent device g2_d/gm2_d arrays (already on device from gpu_init)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceCPML_E_left(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_left_initialized) return

      ! Advance Ex on left boundary
      call gpu_advanceCPML_Ex_left_kernel(this%Ex_d, this%Hz_d, this%sggMiEx_d, &
                                          this%pml_psi_Exy_left, this%pml_P_be_y_left, this%pml_P_ce_y_left, &
                                          this%g2_d, &
                                          this%pml_left_Ex_ii, this%pml_left_Ex_ij, &
                                          this%pml_left_Ex_ji, this%pml_left_Ex_jj, &
                                          this%pml_left_Ex_ki, this%pml_left_Ex_kj, &
                                          b%Ex%XI-1)

      ! Advance Ez on left boundary
      call gpu_advanceCPML_Ez_left_kernel(this%Ez_d, this%Hx_d, this%sggMiEz_d, &
                                          this%pml_psi_Ezy_left, this%pml_P_be_y_left, this%pml_P_ce_y_left, &
                                          this%g2_d, &
                                          this%pml_left_Ez_ii, this%pml_left_Ez_ij, &
                                          this%pml_left_Ez_ji, this%pml_left_Ez_jj, &
                                          this%pml_left_Ez_ki, this%pml_left_Ez_kj, &
                                          b%Ez%XI-1)

   end subroutine gpu_advanceCPML_E_left

   subroutine gpu_advanceCPML_Ex_left_kernel(Ex_d, Hz_d, sggMiEx_d, psi_Exy_d, &
                                              P_be_y_d, P_ce_y_d, g2_d, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ex_d, Hz_d, psi_Exy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEx_d
      real(kind=rkind), device, dimension(:) :: P_be_y_d, P_ce_y_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEx_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Exy_d(i-ii+1,j-ji+1,k-ki+1) = P_be_y_d(j) * psi_Exy_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hz_d(i,j,k) - Hz_d(i,j-1,k)) * P_ce_y_d(j)
              Ex_d(i,j,k) = Ex_d(i,j,k) + g2_d(medio) * psi_Exy_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ex_left_kernel

   subroutine gpu_advanceCPML_Ez_left_kernel(Ez_d, Hx_d, sggMiEz_d, psi_Ezy_d, &
                                              P_be_y_d, P_ce_y_d, g2_d, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ez_d, Hx_d, psi_Ezy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEz_d
      real(kind=rkind), device, dimension(:) :: P_be_y_d, P_ce_y_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEz_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Ezy_d(i-ii+1,j-ji+1,k-ki+1) = P_be_y_d(j) * psi_Ezy_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hx_d(i,j,k) - Hx_d(i,j-1,k)) * P_ce_y_d(j)
              Ez_d(i,j,k) = Ez_d(i,j,k) - g2_d(medio) * psi_Ezy_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ez_left_kernel

   subroutine gpu_advanceCPML_H_left(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_left_initialized) return

      ! Advance Hx on left boundary
      call gpu_advanceCPML_Hx_left_kernel(this%Hx_d, this%Ez_d, this%sggMiHx_d, &
                                          this%pml_psi_Hxy_left, this%pml_P_bm_y_left, this%pml_P_cm_y_left, &
                                          this%gm2_d, &
                                          this%pml_left_Hx_ii, this%pml_left_Hx_ij, &
                                          this%pml_left_Hx_ji, this%pml_left_Hx_jj, &
                                          this%pml_left_Hx_ki, this%pml_left_Hx_kj, &
                                          b%Hx%XI-1)

      ! Advance Hz on left boundary
      call gpu_advanceCPML_Hz_left_kernel(this%Hz_d, this%Ex_d, this%sggMiHz_d, &
                                          this%pml_psi_Hzy_left, this%pml_P_bm_y_left, this%pml_P_cm_y_left, &
                                          this%gm2_d, &
                                          this%pml_left_Hz_ii, this%pml_left_Hz_ij, &
                                          this%pml_left_Hz_ji, this%pml_left_Hz_jj, &
                                          this%pml_left_Hz_ki, this%pml_left_Hz_kj, &
                                          b%Hz%XI-1)

   end subroutine gpu_advanceCPML_H_left

   subroutine gpu_advanceCPML_Hx_left_kernel(Hx_d, Ez_d, sggMiHx_d, psi_Hxy_d, &
                                             P_bm_y_d, P_cm_y_d, gm2_d, &
                                             ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hx_d, Ez_d, psi_Hxy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), device, dimension(:) :: P_bm_y_d, P_cm_y_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHx_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Hxy_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_y_d(j) * psi_Hxy_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ez_d(i,j+1,k) - Ez_d(i,j,k)) * P_cm_y_d(j)
              Hx_d(i,j,k) = Hx_d(i,j,k) - gm2_d(medio) * psi_Hxy_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hx_left_kernel

   subroutine gpu_advanceCPML_Hz_left_kernel(Hz_d, Ex_d, sggMiHz_d, psi_Hzy_d, &
                                             P_bm_y_d, P_cm_y_d, gm2_d, &
                                             ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hz_d, Ex_d, psi_Hzy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), device, dimension(:) :: P_bm_y_d, P_cm_y_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHz_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Hzy_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_y_d(j) * psi_Hzy_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ex_d(i,j+1,k) - Ex_d(i,j,k)) * P_cm_y_d(j)
              Hz_d(i,j,k) = Hz_d(i,j,k) + gm2_d(medio) * psi_Hzy_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hz_left_kernel

   !--------------------------------------------------------------------------------
   ! CPML Right boundary - wrapper subroutines
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceCPML_E_right(this, b, numMedia)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b
      integer, intent(in) :: numMedia
      if (.not. this%pml_right_initialized) return
      call gpu_advanceCPML_Ex_right_kernel(this%Ex_d, this%Hz_d, this%sggMiEx_d, &
                                          this%pml_psi_Exy_left, this%pml_P_be_y_left, this%pml_P_ce_y_left, &
                                          this%g2_d, numMedia, &
                                          this%pml_right_Ex_ii, this%pml_right_Ex_ij, &
                                          this%pml_right_Ex_ji, this%pml_right_Ex_jj, &
                                          this%pml_right_Ex_ki, this%pml_right_Ex_kj, &
                                          b%Ex%XI-1)
      call gpu_advanceCPML_Ez_right_kernel(this%Ez_d, this%Hx_d, this%sggMiEz_d, &
                                          this%pml_psi_Ezy_left, this%pml_P_be_y_left, this%pml_P_ce_y_left, &
                                          this%g2_d, numMedia, &
                                          this%pml_right_Ez_ii, this%pml_right_Ez_ij, &
                                          this%pml_right_Ez_ji, this%pml_right_Ez_jj, &
                                          this%pml_right_Ez_ki, this%pml_right_Ez_kj, &
                                          b%Ez%XI-1)
   end subroutine gpu_advanceCPML_E_right

   subroutine gpu_advanceCPML_H_right(this, b, numMedia)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b
      integer, intent(in) :: numMedia
      if (.not. this%pml_right_initialized) return
      call gpu_advanceCPML_Hx_right_kernel(this%Hx_d, this%Ez_d, this%sggMiHx_d, &
                                          this%pml_psi_Hxy_left, this%pml_P_bm_y_left, this%pml_P_cm_y_left, &
                                          this%gm2_d, numMedia, &
                                          this%pml_right_Hx_ii, this%pml_right_Hx_ij, &
                                          this%pml_right_Hx_ji, this%pml_right_Hx_jj, &
                                          this%pml_right_Hx_ki, this%pml_right_Hx_kj, &
                                          b%Hx%XI-1)
      call gpu_advanceCPML_Hz_right_kernel(this%Hz_d, this%Ex_d, this%sggMiHz_d, &
                                          this%pml_psi_Hzy_left, this%pml_P_bm_y_left, this%pml_P_cm_y_left, &
                                          this%gm2_d, numMedia, &
                                          this%pml_right_Hz_ii, this%pml_right_Hz_ij, &
                                          this%pml_right_Hz_ji, this%pml_right_Hz_jj, &
                                          this%pml_right_Hz_ki, this%pml_right_Hz_kj, &
                                          b%Hz%XI-1)
   end subroutine gpu_advanceCPML_H_right

   subroutine gpu_advanceCPML_Ex_right_kernel(Ex_d, Hz_d, sggMiEx_d, psi_Exy_d, &
                                              P_be_y_d, P_ce_y_d, g2_d, numMedia, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      integer, intent(in) :: numMedia
      real(kind=rkind), device, dimension(:,:,:) :: Ex_d, Hz_d, psi_Exy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEx_d
      real(kind=rkind), device, dimension(:) :: P_be_y_d, P_ce_y_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEx_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Exy_d(i-ii+1,j-ji+1,k-ki+1) = P_be_y_d(j) * psi_Exy_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hz_d(i,j,k) - Hz_d(i,j-1,k)) * P_ce_y_d(j)
              Ex_d(i,j,k) = Ex_d(i,j,k) + g2_d(medio) * psi_Exy_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ex_right_kernel

   subroutine gpu_advanceCPML_Ez_right_kernel(Ez_d, Hx_d, sggMiEz_d, psi_Ezy_d, &
                                              P_be_y_d, P_ce_y_d, g2_d, numMedia, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      integer, intent(in) :: numMedia
      real(kind=rkind), device, dimension(:,:,:) :: Ez_d, Hx_d, psi_Ezy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEz_d
      real(kind=rkind), device, dimension(:) :: P_be_y_d, P_ce_y_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEz_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Ezy_d(i-ii+1,j-ji+1,k-ki+1) = P_be_y_d(j) * psi_Ezy_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hx_d(i,j,k) - Hx_d(i,j-1,k)) * P_ce_y_d(j)
              Ez_d(i,j,k) = Ez_d(i,j,k) - g2_d(medio) * psi_Ezy_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ez_right_kernel

   subroutine gpu_advanceCPML_Hx_right_kernel(Hx_d, Ez_d, sggMiHx_d, psi_Hxy_d, &
                                              P_bm_y_d, P_cm_y_d, gm2_d, numMedia, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      integer, intent(in) :: numMedia
      real(kind=rkind), device, dimension(:,:,:) :: Hx_d, Ez_d, psi_Hxy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), device, dimension(:) :: P_bm_y_d, P_cm_y_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHx_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Hxy_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_y_d(j) * psi_Hxy_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ez_d(i,j+1,k) - Ez_d(i,j,k)) * P_cm_y_d(j)
              Hx_d(i,j,k) = Hx_d(i,j,k) - gm2_d(medio) * psi_Hxy_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hx_right_kernel

   subroutine gpu_advanceCPML_Hz_right_kernel(Hz_d, Ex_d, sggMiHz_d, psi_Hzy_d, &
                                              P_bm_y_d, P_cm_y_d, gm2_d, numMedia, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      integer, intent(in) :: numMedia
      real(kind=rkind), device, dimension(:,:,:) :: Hz_d, Ex_d, psi_Hzy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), device, dimension(:) :: P_bm_y_d, P_cm_y_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHz_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Hzy_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_y_d(j) * psi_Hzy_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ex_d(i,j+1,k) - Ex_d(i,j,k)) * P_cm_y_d(j)
              Hz_d(i,j,k) = Hz_d(i,j,k) + gm2_d(medio) * psi_Hzy_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hz_right_kernel

   !--------------------------------------------------------------------------------
   ! CPML Down Boundary - GPU accelerated kernels (z-dependent coefficients)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceCPML_E_down(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_down_initialized) return

      call gpu_advanceCPML_Ey_down_kernel(this%Ey_d, this%Hx_d, this%sggMiEy_d, &
                                          this%pml_psi_Eyz_down, this%pml_P_be_z_down, this%pml_P_ce_z_down, &
                                          this%g2_d, &
                                          this%pml_down_Ey_ii, this%pml_down_Ey_ij, &
                                          this%pml_down_Ey_ji, this%pml_down_Ey_jj, &
                                          this%pml_down_Ey_ki, this%pml_down_Ey_kj, &
                                          b%Ey%ZI-1)

      call gpu_advanceCPML_Ex_down_kernel(this%Ex_d, this%Hy_d, this%sggMiEx_d, &
                                          this%pml_psi_Exz_down, this%pml_P_be_z_down, this%pml_P_ce_z_down, &
                                          this%g2_d, &
                                          this%pml_down_Ex_ii, this%pml_down_Ex_ij, &
                                          this%pml_down_Ex_ji, this%pml_down_Ex_jj, &
                                          this%pml_down_Ex_ki, this%pml_down_Ex_kj, &
                                          b%Ex%ZI-1)

   end subroutine gpu_advanceCPML_E_down

   subroutine gpu_advanceCPML_Ey_down_kernel(Ey_d, Hx_d, sggMiEy_d, psi_Eyz_d, &
                                             P_be_z_d, P_ce_z_d, g2_d, &
                                             ii, ij, ji, jj, ki, kj, zi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, zi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ey_d, Hx_d, psi_Eyz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEy_d
      real(kind=rkind), device, dimension(:) :: P_be_z_d, P_ce_z_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEy_d(i-zi_offset,j-zi_offset,k-zi_offset)
              psi_Eyz_d(i-ii+1,j-ji+1,k-ki+1) = P_be_z_d(k) * psi_Eyz_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hx_d(i,j,k) - Hx_d(i,j,k-1)) * P_ce_z_d(k)
              Ey_d(i,j,k) = Ey_d(i,j,k) + g2_d(medio) * psi_Eyz_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ey_down_kernel

   subroutine gpu_advanceCPML_Ex_down_kernel(Ex_d, Hy_d, sggMiEx_d, psi_Exz_d, &
                                             P_be_z_d, P_ce_z_d, g2_d, &
                                             ii, ij, ji, jj, ki, kj, zi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, zi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ex_d, Hy_d, psi_Exz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEx_d
      real(kind=rkind), device, dimension(:) :: P_be_z_d, P_ce_z_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEx_d(i-zi_offset,j-zi_offset,k-zi_offset)
              psi_Exz_d(i-ii+1,j-ji+1,k-ki+1) = P_be_z_d(k) * psi_Exz_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hy_d(i,j,k) - Hy_d(i,j,k-1)) * P_ce_z_d(k)
              Ex_d(i,j,k) = Ex_d(i,j,k) - g2_d(medio) * psi_Exz_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ex_down_kernel

   subroutine gpu_advanceCPML_H_down(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_down_initialized) return

      call gpu_advanceCPML_Hy_down_kernel(this%Hy_d, this%Ex_d, this%sggMiHy_d, &
                                          this%pml_psi_Hyz_down, this%pml_P_bm_z_down, this%pml_P_cm_z_down, &
                                          this%gm2_d, &
                                          this%pml_down_Hy_ii, this%pml_down_Hy_ij, &
                                          this%pml_down_Hy_ji, this%pml_down_Hy_jj, &
                                          this%pml_down_Hy_ki, this%pml_down_Hy_kj, &
                                          b%Hy%ZI-1)

      call gpu_advanceCPML_Hx_down_kernel(this%Hx_d, this%Ey_d, this%sggMiHx_d, &
                                          this%pml_psi_Hxz_down, this%pml_P_bm_z_down, this%pml_P_cm_z_down, &
                                          this%gm2_d, &
                                          this%pml_down_Hx_ii, this%pml_down_Hx_ij, &
                                          this%pml_down_Hx_ji, this%pml_down_Hx_jj, &
                                          this%pml_down_Hx_ki, this%pml_down_Hx_kj, &
                                          b%Hx%ZI-1)

   end subroutine gpu_advanceCPML_H_down

   subroutine gpu_advanceCPML_Hy_down_kernel(Hy_d, Ex_d, sggMiHy_d, psi_Hyz_d, &
                                             P_bm_z_d, P_cm_z_d, gm2_d, &
                                             ii, ij, ji, jj, ki, kj, zi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, zi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hy_d, Ex_d, psi_Hyz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), device, dimension(:) :: P_bm_z_d, P_cm_z_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHy_d(i-zi_offset,j-zi_offset,k-zi_offset)
              psi_Hyz_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_z_d(k) * psi_Hyz_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ex_d(i,j,k+1) - Ex_d(i,j,k)) * P_cm_z_d(k)
              Hy_d(i,j,k) = Hy_d(i,j,k) - gm2_d(medio) * psi_Hyz_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hy_down_kernel

   subroutine gpu_advanceCPML_Hx_down_kernel(Hx_d, Ey_d, sggMiHx_d, psi_Hxz_d, &
                                             P_bm_z_d, P_cm_z_d, gm2_d, &
                                             ii, ij, ji, jj, ki, kj, zi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, zi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hx_d, Ey_d, psi_Hxz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), device, dimension(:) :: P_bm_z_d, P_cm_z_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHx_d(i-zi_offset,j-zi_offset,k-zi_offset)
              psi_Hxz_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_z_d(k) * psi_Hxz_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ey_d(i,j,k+1) - Ey_d(i,j,k)) * P_cm_z_d(k)
              Hx_d(i,j,k) = Hx_d(i,j,k) + gm2_d(medio) * psi_Hxz_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hx_down_kernel

   !--------------------------------------------------------------------------------
   ! CPML Up Boundary - GPU accelerated kernels (z-dependent, same as down)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceCPML_E_up(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_up_initialized) return

      call gpu_advanceCPML_Ey_up_kernel(this%Ey_d, this%Hx_d, this%sggMiEy_d, &
                                        this%pml_psi_Eyz_down, this%pml_P_be_z_down, this%pml_P_ce_z_down, &
                                        this%g2_d, &
                                        this%pml_up_Ey_ii, this%pml_up_Ey_ij, &
                                        this%pml_up_Ey_ji, this%pml_up_Ey_jj, &
                                        this%pml_up_Ey_ki, this%pml_up_Ey_kj, &
                                        b%Ey%ZI-1)

      call gpu_advanceCPML_Ex_up_kernel(this%Ex_d, this%Hy_d, this%sggMiEx_d, &
                                        this%pml_psi_Exz_down, this%pml_P_be_z_down, this%pml_P_ce_z_down, &
                                        this%g2_d, &
                                        this%pml_up_Ex_ii, this%pml_up_Ex_ij, &
                                        this%pml_up_Ex_ji, this%pml_up_Ex_jj, &
                                        this%pml_up_Ex_ki, this%pml_up_Ex_kj, &
                                        b%Ex%ZI-1)

   end subroutine gpu_advanceCPML_E_up

   subroutine gpu_advanceCPML_Ey_up_kernel(Ey_d, Hx_d, sggMiEy_d, psi_Eyz_d, &
                                           P_be_z_d, P_ce_z_d, g2_d, &
                                           ii, ij, ji, jj, ki, kj, zi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, zi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ey_d, Hx_d, psi_Eyz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEy_d
      real(kind=rkind), device, dimension(:) :: P_be_z_d, P_ce_z_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEy_d(i-zi_offset,j-zi_offset,k-zi_offset)
              psi_Eyz_d(i-ii+1,j-ji+1,k-ki+1) = P_be_z_d(k) * psi_Eyz_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hx_d(i,j,k) - Hx_d(i,j,k-1)) * P_ce_z_d(k)
              Ey_d(i,j,k) = Ey_d(i,j,k) + g2_d(medio) * psi_Eyz_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ey_up_kernel

   subroutine gpu_advanceCPML_Ex_up_kernel(Ex_d, Hy_d, sggMiEx_d, psi_Exz_d, &
                                           P_be_z_d, P_ce_z_d, g2_d, &
                                           ii, ij, ji, jj, ki, kj, zi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, zi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ex_d, Hy_d, psi_Exz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEx_d
      real(kind=rkind), device, dimension(:) :: P_be_z_d, P_ce_z_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEx_d(i-zi_offset,j-zi_offset,k-zi_offset)
              psi_Exz_d(i-ii+1,j-ji+1,k-ki+1) = P_be_z_d(k) * psi_Exz_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hy_d(i,j,k) - Hy_d(i,j,k-1)) * P_ce_z_d(k)
              Ex_d(i,j,k) = Ex_d(i,j,k) - g2_d(medio) * psi_Exz_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ex_up_kernel

   subroutine gpu_advanceCPML_H_up(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_up_initialized) return

      call gpu_advanceCPML_Hy_up_kernel(this%Hy_d, this%Ex_d, this%sggMiHy_d, &
                                        this%pml_psi_Hyz_down, this%pml_P_bm_z_down, this%pml_P_cm_z_down, &
                                        this%gm2_d, &
                                        this%pml_up_Hy_ii, this%pml_up_Hy_ij, &
                                        this%pml_up_Hy_ji, this%pml_up_Hy_jj, &
                                        this%pml_up_Hy_ki, this%pml_up_Hy_kj, &
                                        b%Hy%ZI-1)

      call gpu_advanceCPML_Hx_up_kernel(this%Hx_d, this%Ey_d, this%sggMiHx_d, &
                                        this%pml_psi_Hxz_down, this%pml_P_bm_z_down, this%pml_P_cm_z_down, &
                                        this%gm2_d, &
                                        this%pml_up_Hx_ii, this%pml_up_Hx_ij, &
                                        this%pml_up_Hx_ji, this%pml_up_Hx_jj, &
                                        this%pml_up_Hx_ki, this%pml_up_Hx_kj, &
                                        b%Hx%ZI-1)

   end subroutine gpu_advanceCPML_H_up

   subroutine gpu_advanceCPML_Hy_up_kernel(Hy_d, Ex_d, sggMiHy_d, psi_Hyz_d, &
                                           P_bm_z_d, P_cm_z_d, gm2_d, &
                                           ii, ij, ji, jj, ki, kj, zi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, zi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hy_d, Ex_d, psi_Hyz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), device, dimension(:) :: P_bm_z_d, P_cm_z_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHy_d(i-zi_offset,j-zi_offset,k-zi_offset)
              psi_Hyz_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_z_d(k) * psi_Hyz_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ex_d(i,j,k+1) - Ex_d(i,j,k)) * P_cm_z_d(k)
              Hy_d(i,j,k) = Hy_d(i,j,k) - gm2_d(medio) * psi_Hyz_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hy_up_kernel

   subroutine gpu_advanceCPML_Hx_up_kernel(Hx_d, Ey_d, sggMiHx_d, psi_Hxz_d, &
                                           P_bm_z_d, P_cm_z_d, gm2_d, &
                                           ii, ij, ji, jj, ki, kj, zi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, zi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hx_d, Ey_d, psi_Hxz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), device, dimension(:) :: P_bm_z_d, P_cm_z_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHx_d(i-zi_offset,j-zi_offset,k-zi_offset)
              psi_Hxz_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_z_d(k) * psi_Hxz_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ey_d(i,j,k+1) - Ey_d(i,j,k)) * P_cm_z_d(k)
              Hx_d(i,j,k) = Hx_d(i,j,k) + gm2_d(medio) * psi_Hxz_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hx_up_kernel

   !--------------------------------------------------------------------------------
   ! CPML Back Boundary - GPU accelerated kernels (x-dependent coefficients)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceCPML_E_back(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_back_initialized) return

      call gpu_advanceCPML_Ez_back_kernel(this%Ez_d, this%Hy_d, this%sggMiEz_d, &
                                          this%pml_psi_Ezx_back, this%pml_P_be_x_back, this%pml_P_ce_x_back, &
                                          this%g2_d, &
                                          this%pml_back_Ez_ii, this%pml_back_Ez_ij, &
                                          this%pml_back_Ez_ji, this%pml_back_Ez_jj, &
                                          this%pml_back_Ez_ki, this%pml_back_Ez_kj, &
                                          b%Ez%XI-1)

      call gpu_advanceCPML_Ey_back_kernel(this%Ey_d, this%Hz_d, this%sggMiEy_d, &
                                          this%pml_psi_Eyx_back, this%pml_P_be_x_back, this%pml_P_ce_x_back, &
                                          this%g2_d, &
                                          this%pml_back_Ey_ii, this%pml_back_Ey_ij, &
                                          this%pml_back_Ey_ji, this%pml_back_Ey_jj, &
                                          this%pml_back_Ey_ki, this%pml_back_Ey_kj, &
                                          b%Ey%XI-1)

   end subroutine gpu_advanceCPML_E_back

   subroutine gpu_advanceCPML_Ez_back_kernel(Ez_d, Hy_d, sggMiEz_d, psi_Ezx_d, &
                                             P_be_x_d, P_ce_x_d, g2_d, &
                                             ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ez_d, Hy_d, psi_Ezx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEz_d
      real(kind=rkind), device, dimension(:) :: P_be_x_d, P_ce_x_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEz_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Ezx_d(i-ii+1,j-ji+1,k-ki+1) = P_be_x_d(i) * psi_Ezx_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hy_d(i,j,k) - Hy_d(i-1,j,k)) * P_ce_x_d(i)
              Ez_d(i,j,k) = Ez_d(i,j,k) + g2_d(medio) * psi_Ezx_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ez_back_kernel

   subroutine gpu_advanceCPML_Ey_back_kernel(Ey_d, Hz_d, sggMiEy_d, psi_Eyx_d, &
                                             P_be_x_d, P_ce_x_d, g2_d, &
                                             ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ey_d, Hz_d, psi_Eyx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEy_d
      real(kind=rkind), device, dimension(:) :: P_be_x_d, P_ce_x_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEy_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Eyx_d(i-ii+1,j-ji+1,k-ki+1) = P_be_x_d(i) * psi_Eyx_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hz_d(i,j,k) - Hz_d(i-1,j,k)) * P_ce_x_d(i)
              Ey_d(i,j,k) = Ey_d(i,j,k) - g2_d(medio) * psi_Eyx_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ey_back_kernel

   subroutine gpu_advanceCPML_H_back(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_back_initialized) return

      call gpu_advanceCPML_Hz_back_kernel(this%Hz_d, this%Hy_d, this%sggMiHz_d, &
                                          this%pml_psi_Hzx_back, this%pml_P_bm_x_back, this%pml_P_cm_x_back, &
                                          this%gm2_d, &
                                          this%pml_back_Hz_ii, this%pml_back_Hz_ij, &
                                          this%pml_back_Hz_ji, this%pml_back_Hz_jj, &
                                          this%pml_back_Hz_ki, this%pml_back_Hz_kj, &
                                          b%Hz%XI-1)

      call gpu_advanceCPML_Hy_back_kernel(this%Hy_d, this%Ez_d, this%sggMiHy_d, &
                                          this%pml_psi_Hyx_back, this%pml_P_bm_x_back, this%pml_P_cm_x_back, &
                                          this%gm2_d, &
                                          this%pml_back_Hy_ii, this%pml_back_Hy_ij, &
                                          this%pml_back_Hy_ji, this%pml_back_Hy_jj, &
                                          this%pml_back_Hy_ki, this%pml_back_Hy_kj, &
                                          b%Hy%XI-1)

   end subroutine gpu_advanceCPML_H_back

   subroutine gpu_advanceCPML_Hz_back_kernel(Hz_d, Hy_d, sggMiHz_d, psi_Hzx_d, &
                                             P_bm_x_d, P_cm_x_d, gm2_d, &
                                             ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hz_d, Hy_d, psi_Hzx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), device, dimension(:) :: P_bm_x_d, P_cm_x_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHz_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Hzx_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_x_d(i) * psi_Hzx_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hy_d(i,j,k) - Hy_d(i-1,j,k)) * P_cm_x_d(i)
              Hz_d(i,j,k) = Hz_d(i,j,k) - gm2_d(medio) * psi_Hzx_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hz_back_kernel

   subroutine gpu_advanceCPML_Hy_back_kernel(Hy_d, Ez_d, sggMiHy_d, psi_Hyx_d, &
                                             P_bm_x_d, P_cm_x_d, gm2_d, &
                                             ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hy_d, Ez_d, psi_Hyx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), device, dimension(:) :: P_bm_x_d, P_cm_x_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHy_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Hyx_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_x_d(i) * psi_Hyx_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ez_d(i,j,k) - Ez_d(i-1,j,k)) * P_cm_x_d(i)
              Hy_d(i,j,k) = Hy_d(i,j,k) - gm2_d(medio) * psi_Hyx_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hy_back_kernel

   !--------------------------------------------------------------------------------
   ! CPML Front Boundary - GPU accelerated kernels (x-dependent, same as back)
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceCPML_E_front(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_front_initialized) return

      call gpu_advanceCPML_Ez_front_kernel(this%Ez_d, this%Hy_d, this%sggMiEz_d, &
                                           this%pml_psi_Ezx_back, this%pml_P_be_x_back, this%pml_P_ce_x_back, &
                                           this%g2_d, &
                                           this%pml_front_Ez_ii, this%pml_front_Ez_ij, &
                                           this%pml_front_Ez_ji, this%pml_front_Ez_jj, &
                                           this%pml_front_Ez_ki, this%pml_front_Ez_kj, &
                                           b%Ez%XI-1)

      call gpu_advanceCPML_Ey_front_kernel(this%Ey_d, this%Hz_d, this%sggMiEy_d, &
                                           this%pml_psi_Eyx_back, this%pml_P_be_x_back, this%pml_P_ce_x_back, &
                                           this%g2_d, &
                                           this%pml_front_Ey_ii, this%pml_front_Ey_ij, &
                                           this%pml_front_Ey_ji, this%pml_front_Ey_jj, &
                                           this%pml_front_Ey_ki, this%pml_front_Ey_kj, &
                                           b%Ey%XI-1)

   end subroutine gpu_advanceCPML_E_front

   subroutine gpu_advanceCPML_Ez_front_kernel(Ez_d, Hy_d, sggMiEz_d, psi_Ezx_d, &
                                              P_be_x_d, P_ce_x_d, g2_d, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ez_d, Hy_d, psi_Ezx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEz_d
      real(kind=rkind), device, dimension(:) :: P_be_x_d, P_ce_x_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEz_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Ezx_d(i-ii+1,j-ji+1,k-ki+1) = P_be_x_d(i) * psi_Ezx_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hy_d(i,j,k) - Hy_d(i-1,j,k)) * P_ce_x_d(i)
              Ez_d(i,j,k) = Ez_d(i,j,k) + g2_d(medio) * psi_Ezx_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ez_front_kernel

   subroutine gpu_advanceCPML_Ey_front_kernel(Ey_d, Hz_d, sggMiEy_d, psi_Eyx_d, &
                                              P_be_x_d, P_ce_x_d, g2_d, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Ey_d, Hz_d, psi_Eyx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiEy_d
      real(kind=rkind), device, dimension(:) :: P_be_x_d, P_ce_x_d, g2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiEy_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Eyx_d(i-ii+1,j-ji+1,k-ki+1) = P_be_x_d(i) * psi_Eyx_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hz_d(i,j,k) - Hz_d(i-1,j,k)) * P_ce_x_d(i)
              Ey_d(i,j,k) = Ey_d(i,j,k) - g2_d(medio) * psi_Eyx_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Ey_front_kernel

   subroutine gpu_advanceCPML_H_front(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%pml_front_initialized) return

      call gpu_advanceCPML_Hz_front_kernel(this%Hz_d, this%Hy_d, this%sggMiHz_d, &
                                           this%pml_psi_Hzx_back, this%pml_P_bm_x_back, this%pml_P_cm_x_back, &
                                           this%gm2_d, &
                                           this%pml_front_Hz_ii, this%pml_front_Hz_ij, &
                                           this%pml_front_Hz_ji, this%pml_front_Hz_jj, &
                                           this%pml_front_Hz_ki, this%pml_front_Hz_kj, &
                                           b%Hz%XI-1)

      call gpu_advanceCPML_Hy_front_kernel(this%Hy_d, this%Ez_d, this%sggMiHy_d, &
                                           this%pml_psi_Hyx_back, this%pml_P_bm_x_back, this%pml_P_cm_x_back, &
                                           this%gm2_d, &
                                           this%pml_front_Hy_ii, this%pml_front_Hy_ij, &
                                           this%pml_front_Hy_ji, this%pml_front_Hy_jj, &
                                           this%pml_front_Hy_ki, this%pml_front_Hy_kj, &
                                           b%Hy%XI-1)

   end subroutine gpu_advanceCPML_H_front

   subroutine gpu_advanceCPML_Hz_front_kernel(Hz_d, Hy_d, sggMiHz_d, psi_Hzx_d, &
                                              P_bm_x_d, P_cm_x_d, gm2_d, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hz_d, Hy_d, psi_Hzx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), device, dimension(:) :: P_bm_x_d, P_cm_x_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHz_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Hzx_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_x_d(i) * psi_Hzx_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Hy_d(i,j,k) - Hy_d(i-1,j,k)) * P_cm_x_d(i)
              Hz_d(i,j,k) = Hz_d(i,j,k) - gm2_d(medio) * psi_Hzx_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hz_front_kernel

   subroutine gpu_advanceCPML_Hy_front_kernel(Hy_d, Ez_d, sggMiHy_d, psi_Hyx_d, &
                                              P_bm_x_d, P_cm_x_d, gm2_d, &
                                              ii, ij, ji, jj, ki, kj, xi_offset)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj, xi_offset
      real(kind=rkind), device, dimension(:,:,:) :: Hy_d, Ez_d, psi_Hyx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), device, dimension(:) :: P_bm_x_d, P_cm_x_d, gm2_d

      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(3) <<<*, *>>>
      do k=ki,kj
         do j=ji,jj
            do i=ii,ij
              medio = sggMiHy_d(i-xi_offset,j-xi_offset,k-xi_offset)
              psi_Hyx_d(i-ii+1,j-ji+1,k-ki+1) = P_bm_x_d(i) * psi_Hyx_d(i-ii+1,j-ji+1,k-ki+1) + &
                 (Ez_d(i,j,k) - Ez_d(i-1,j,k)) * P_cm_x_d(i)
              Hy_d(i,j,k) = Hy_d(i,j,k) - gm2_d(medio) * psi_Hyx_d(i-ii+1,j-ji+1,k-ki+1)
            end do
         end do
      end do
   end subroutine gpu_advanceCPML_Hy_front_kernel

end module gpu_cpml_m
