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

end module gpu_cpml_m
