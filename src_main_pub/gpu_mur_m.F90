!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  GPU MUR KERNELS MODULE - CUDA Fortran (CUF) accelerated MUR boundaries
!  First-order Mur ABC for all 6 boundaries with persistent state on device.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module gpu_mur_m

   use FDETYPES_m
   use Report_m
   use cudafor
   use gpu_core_m

   implicit none

contains

   !--------------------------------------------------------------------------------
   ! MUR Left boundary - GPU accelerated kernels
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceMUR_H_left(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%mur_initialized) return

      call gpu_advanceMUR_Hx_left_kernel(this%Hx_d, this%sggMiHx_d, &
                                         this%mur_past_Hx_left, this%mur_left_CAB1, &
                                         this%mur_left_Hx_ii, this%mur_left_Hx_ij, &
                                         this%mur_left_Hx_ji, this%mur_left_Hx_jj, &
                                         this%mur_left_Hx_ki, this%mur_left_Hx_kj)

      call gpu_advanceMUR_Hz_left_kernel(this%Hz_d, this%sggMiHz_d, &
                                         this%mur_past_Hz_left, this%mur_left_CAB1, &
                                         this%mur_left_Hz_ii, this%mur_left_Hz_ij, &
                                         this%mur_left_Hz_ji, this%mur_left_Hz_jj, &
                                         this%mur_left_Hz_ki, this%mur_left_Hz_kj)

   end subroutine gpu_advanceMUR_H_left

   subroutine gpu_advanceMUR_Hx_left_kernel(Hx_d, sggMiHx_d, past_Hx_d, CAB1_d, &
                                            ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hx_d, past_Hx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do k=ki,kj
         do i=ii,ij
            medio = sggMiHx_d(i, ji+1, k)
            Hx_d(i, ji, k) = past_Hx_d(i, ji+1, k) + &
                             CAB1_d(medio) * (Hx_d(i, ji+1, k) - past_Hx_d(i, ji, k))
         end do
      end do
   end subroutine gpu_advanceMUR_Hx_left_kernel

   subroutine gpu_advanceMUR_Hz_left_kernel(Hz_d, sggMiHz_d, past_Hz_d, CAB1_d, &
                                            ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hz_d, past_Hz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do k=ki,kj
         do i=ii,ij
            medio = sggMiHz_d(i, ji+1, k)
            Hz_d(i, ji, k) = past_Hz_d(i, ji+1, k) + &
                             CAB1_d(medio) * (Hz_d(i, ji+1, k) - past_Hz_d(i, ji, k))
         end do
      end do
   end subroutine gpu_advanceMUR_Hz_left_kernel

   !--------------------------------------------------------------------------------
   ! MUR Right boundary - GPU accelerated kernels
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceMUR_H_right(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%mur_initialized) return

      call gpu_advanceMUR_Hx_right_kernel(this%Hx_d, this%sggMiHx_d, &
                                          this%mur_past_Hx_right, this%mur_right_CAB1, &
                                          this%mur_right_Hx_ii, this%mur_right_Hx_ij, &
                                          this%mur_right_Hx_ji, this%mur_right_Hx_jj, &
                                          this%mur_right_Hx_ki, this%mur_right_Hx_kj)

      call gpu_advanceMUR_Hz_right_kernel(this%Hz_d, this%sggMiHz_d, &
                                          this%mur_past_Hz_right, this%mur_right_CAB1, &
                                          this%mur_right_Hz_ii, this%mur_right_Hz_ij, &
                                          this%mur_right_Hz_ji, this%mur_right_Hz_jj, &
                                          this%mur_right_Hz_ki, this%mur_right_Hz_kj)

   end subroutine gpu_advanceMUR_H_right

   subroutine gpu_advanceMUR_Hx_right_kernel(Hx_d, sggMiHx_d, past_Hx_d, CAB1_d, &
                                             ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hx_d, past_Hx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do k=ki,kj
         do i=ii,ij
            medio = sggMiHx_d(i, ji-1, k)
            Hx_d(i, ji, k) = past_Hx_d(i, ji-1, k) + &
                             CAB1_d(medio) * (Hx_d(i, ji-1, k) - past_Hx_d(i, ji, k))
         end do
      end do
   end subroutine gpu_advanceMUR_Hx_right_kernel

   subroutine gpu_advanceMUR_Hz_right_kernel(Hz_d, sggMiHz_d, past_Hz_d, CAB1_d, &
                                             ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hz_d, past_Hz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, k
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do k=ki,kj
         do i=ii,ij
            medio = sggMiHz_d(i, ji-1, k)
            Hz_d(i, ji, k) = past_Hz_d(i, ji-1, k) + &
                             CAB1_d(medio) * (Hz_d(i, ji-1, k) - past_Hz_d(i, ji, k))
         end do
      end do
   end subroutine gpu_advanceMUR_Hz_right_kernel

   !--------------------------------------------------------------------------------
   ! MUR Down boundary - GPU accelerated kernels
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceMUR_H_down(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%mur_initialized) return

      call gpu_advanceMUR_Hy_down_kernel(this%Hy_d, this%sggMiHy_d, &
                                         this%mur_past_Hy_down, this%mur_down_CAB1, &
                                         this%mur_down_Hy_ii, this%mur_down_Hy_ij, &
                                         this%mur_down_Hy_ji, this%mur_down_Hy_jj, &
                                         this%mur_down_Hy_ki, this%mur_down_Hy_kj)

      call gpu_advanceMUR_Hx_down_kernel(this%Hx_d, this%sggMiHx_d, &
                                         this%mur_past_Hx_down, this%mur_down_CAB1, &
                                         this%mur_down_Hx_ii, this%mur_down_Hx_ij, &
                                         this%mur_down_Hx_ji, this%mur_down_Hx_jj, &
                                         this%mur_down_Hx_ki, this%mur_down_Hx_kj)

   end subroutine gpu_advanceMUR_H_down

   subroutine gpu_advanceMUR_Hy_down_kernel(Hy_d, sggMiHy_d, past_Hy_d, CAB1_d, &
                                            ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hy_d, past_Hy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, j
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do j=ji,jj
         do i=ii,ij
            medio = sggMiHy_d(i, j, ki+1)
            Hy_d(i, j, ki) = past_Hy_d(i, j, ki+1) + &
                             CAB1_d(medio) * (Hy_d(i, j, ki+1) - past_Hy_d(i, j, ki))
         end do
      end do
   end subroutine gpu_advanceMUR_Hy_down_kernel

   subroutine gpu_advanceMUR_Hx_down_kernel(Hx_d, sggMiHx_d, past_Hx_d, CAB1_d, &
                                            ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hx_d, past_Hx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, j
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do j=ji,jj
         do i=ii,ij
            medio = sggMiHx_d(i, j, ki+1)
            Hx_d(i, j, ki) = past_Hx_d(i, j, ki+1) + &
                             CAB1_d(medio) * (Hx_d(i, j, ki+1) - past_Hx_d(i, j, ki))
         end do
      end do
   end subroutine gpu_advanceMUR_Hx_down_kernel

   !--------------------------------------------------------------------------------
   ! MUR Up boundary - GPU accelerated kernels
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceMUR_H_up(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%mur_initialized) return

      call gpu_advanceMUR_Hy_up_kernel(this%Hy_d, this%sggMiHy_d, &
                                       this%mur_past_Hy_up, this%mur_up_CAB1, &
                                       this%mur_up_Hy_ii, this%mur_up_Hy_ij, &
                                       this%mur_up_Hy_ji, this%mur_up_Hy_jj, &
                                       this%mur_up_Hy_ki, this%mur_up_Hy_kj)

      call gpu_advanceMUR_Hx_up_kernel(this%Hx_d, this%sggMiHx_d, &
                                       this%mur_past_Hx_up, this%mur_up_CAB1, &
                                       this%mur_up_Hx_ii, this%mur_up_Hx_ij, &
                                       this%mur_up_Hx_ji, this%mur_up_Hx_jj, &
                                       this%mur_up_Hx_ki, this%mur_up_Hx_kj)

   end subroutine gpu_advanceMUR_H_up

   subroutine gpu_advanceMUR_Hy_up_kernel(Hy_d, sggMiHy_d, past_Hy_d, CAB1_d, &
                                          ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hy_d, past_Hy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, j
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do j=ji,jj
         do i=ii,ij
            medio = sggMiHy_d(i, j, kj-1)
            Hy_d(i, j, kj) = past_Hy_d(i, j, kj-1) + &
                             CAB1_d(medio) * (Hy_d(i, j, kj-1) - past_Hy_d(i, j, kj))
         end do
      end do
   end subroutine gpu_advanceMUR_Hy_up_kernel

   subroutine gpu_advanceMUR_Hx_up_kernel(Hx_d, sggMiHx_d, past_Hx_d, CAB1_d, &
                                          ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hx_d, past_Hx_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, j
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do j=ji,jj
         do i=ii,ij
            medio = sggMiHx_d(i, j, kj-1)
            Hx_d(i, j, kj) = past_Hx_d(i, j, kj-1) + &
                             CAB1_d(medio) * (Hx_d(i, j, kj-1) - past_Hx_d(i, j, kj))
         end do
      end do
   end subroutine gpu_advanceMUR_Hx_up_kernel

   !--------------------------------------------------------------------------------
   ! MUR Back boundary - GPU accelerated kernels
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceMUR_H_back(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%mur_initialized) return

      call gpu_advanceMUR_Hz_back_kernel(this%Hz_d, this%sggMiHz_d, &
                                         this%mur_past_Hz_back, this%mur_back_CAB1, &
                                         this%mur_back_Hz_ii, this%mur_back_Hz_ij, &
                                         this%mur_back_Hz_ji, this%mur_back_Hz_jj, &
                                         this%mur_back_Hz_ki, this%mur_back_Hz_kj)

      call gpu_advanceMUR_Hy_back_kernel(this%Hy_d, this%sggMiHy_d, &
                                         this%mur_past_Hy_back, this%mur_back_CAB1, &
                                         this%mur_back_Hy_ii, this%mur_back_Hy_ij, &
                                         this%mur_back_Hy_ji, this%mur_back_Hy_jj, &
                                         this%mur_back_Hy_ki, this%mur_back_Hy_kj)

   end subroutine gpu_advanceMUR_H_back

   subroutine gpu_advanceMUR_Hz_back_kernel(Hz_d, sggMiHz_d, past_Hz_d, CAB1_d, &
                                            ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hz_d, past_Hz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, j
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do j=ji,jj
         do i=ii,ij
            medio = sggMiHz_d(i, j, ki-1)
            Hz_d(i, j, ki) = past_Hz_d(i, j, ki-1) + &
                             CAB1_d(medio) * (Hz_d(i, j, ki-1) - past_Hz_d(i, j, ki))
         end do
      end do
   end subroutine gpu_advanceMUR_Hz_back_kernel

   subroutine gpu_advanceMUR_Hy_back_kernel(Hy_d, sggMiHy_d, past_Hy_d, CAB1_d, &
                                            ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hy_d, past_Hy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, j
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do j=ji,jj
         do i=ii,ij
            medio = sggMiHy_d(i, j, ki-1)
            Hy_d(i, j, ki) = past_Hy_d(i, j, ki-1) + &
                             CAB1_d(medio) * (Hy_d(i, j, ki-1) - past_Hy_d(i, j, ki))
         end do
      end do
   end subroutine gpu_advanceMUR_Hy_back_kernel

   !--------------------------------------------------------------------------------
   ! MUR Front boundary - GPU accelerated kernels
   !--------------------------------------------------------------------------------
   subroutine gpu_advanceMUR_H_front(this, b)
      class(gpu_state_t), intent(inout) :: this
      type(bounds_t), intent(in) :: b

      if (.not. this%mur_initialized) return

      call gpu_advanceMUR_Hz_front_kernel(this%Hz_d, this%sggMiHz_d, &
                                          this%mur_past_Hz_front, this%mur_front_CAB1, &
                                          this%mur_front_Hz_ii, this%mur_front_Hz_ij, &
                                          this%mur_front_Hz_ji, this%mur_front_Hz_jj, &
                                          this%mur_front_Hz_ki, this%mur_front_Hz_kj)

      call gpu_advanceMUR_Hy_front_kernel(this%Hy_d, this%sggMiHy_d, &
                                          this%mur_past_Hy_front, this%mur_front_CAB1, &
                                          this%mur_front_Hy_ii, this%mur_front_Hy_ij, &
                                          this%mur_front_Hy_ji, this%mur_front_Hy_jj, &
                                          this%mur_front_Hy_ki, this%mur_front_Hy_kj)

   end subroutine gpu_advanceMUR_H_front

   subroutine gpu_advanceMUR_Hz_front_kernel(Hz_d, sggMiHz_d, past_Hz_d, CAB1_d, &
                                             ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hz_d, past_Hz_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, j
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do j=ji,jj
         do i=ii,ij
            medio = sggMiHz_d(i, j, kj+1)
            Hz_d(i, j, kj) = past_Hz_d(i, j, kj+1) + &
                             CAB1_d(medio) * (Hz_d(i, j, kj+1) - past_Hz_d(i, j, kj))
         end do
      end do
   end subroutine gpu_advanceMUR_Hz_front_kernel

   subroutine gpu_advanceMUR_Hy_front_kernel(Hy_d, sggMiHy_d, past_Hy_d, CAB1_d, &
                                             ii, ij, ji, jj, ki, kj)
      integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
      real(kind=rkind), device, dimension(:,:,:) :: Hy_d, past_Hy_d
      integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
      real(kind=rkind), device, dimension(:) :: CAB1_d

      integer(kind=4) :: i, j
      integer(kind=integersizeofmediamatrices) :: medio

      !$cuf kernel do(2) <<<*, *>>>
      do j=ji,jj
         do i=ii,ij
            medio = sggMiHy_d(i, j, kj+1)
            Hy_d(i, j, kj) = past_Hy_d(i, j, kj+1) + &
                             CAB1_d(medio) * (Hy_d(i, j, kj+1) - past_Hy_d(i, j, kj))
         end do
      end do
   end subroutine gpu_advanceMUR_Hy_front_kernel

end module gpu_mur_m
