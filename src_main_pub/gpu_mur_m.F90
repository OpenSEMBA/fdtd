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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
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

      !$cuf kernel do(2) 
      do j=ji,jj
         do i=ii,ij
            medio = sggMiHy_d(i, j, kj+1)
            Hy_d(i, j, kj) = past_Hy_d(i, j, kj+1) + &
                             CAB1_d(medio) * (Hy_d(i, j, kj+1) - past_Hy_d(i, j, kj))
         end do
      end do
   end subroutine gpu_advanceMUR_Hy_front_kernel

    !--------------------------------------------------------------------------------
    ! Update MUR past fields on device - copy current Hx/Hy/Hz to past arrays
    ! Called after each MUR step to prepare for next timestep
    !--------------------------------------------------------------------------------
    subroutine gpu_update_mur_past_left(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call gpu_update_mur_past_Hx_left_kernel(this%Hx_d, this%mur_past_Hx_left, &
                                               this%mur_left_Hx_ii, this%mur_left_Hx_ij, &
                                               this%mur_left_Hx_ji, this%mur_left_Hx_jj, &
                                               this%mur_left_Hx_ki, this%mur_left_Hx_kj)

       call gpu_update_mur_past_Hz_left_kernel(this%Hz_d, this%mur_past_Hz_left, &
                                               this%mur_left_Hz_ii, this%mur_left_Hz_ij, &
                                               this%mur_left_Hz_ji, this%mur_left_Hz_jj, &
                                               this%mur_left_Hz_ki, this%mur_left_Hz_kj)
    end subroutine gpu_update_mur_past_left

    subroutine gpu_update_mur_past_Hx_left_kernel(Hx_d, past_Hx_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hx_d, past_Hx_d

       integer(kind=4) :: i, k

       !$cuf kernel do(2) 
       do k=ki,kj
          do i=ii,ij
             past_Hx_d(i, ji, k) = Hx_d(i, ji, k)
          end do
       end do
    end subroutine gpu_update_mur_past_Hx_left_kernel

    subroutine gpu_update_mur_past_Hz_left_kernel(Hz_d, past_Hz_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hz_d, past_Hz_d

       integer(kind=4) :: i, k

       !$cuf kernel do(2) 
       do k=ki,kj
          do i=ii,ij
             past_Hz_d(i, ji, k) = Hz_d(i, ji, k)
          end do
       end do
    end subroutine gpu_update_mur_past_Hz_left_kernel

    subroutine gpu_update_mur_past_right(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call gpu_update_mur_past_Hx_right_kernel(this%Hx_d, this%mur_past_Hx_right, &
                                                this%mur_right_Hx_ii, this%mur_right_Hx_ij, &
                                                this%mur_right_Hx_ji, this%mur_right_Hx_jj, &
                                                this%mur_right_Hx_ki, this%mur_right_Hx_kj)

       call gpu_update_mur_past_Hz_right_kernel(this%Hz_d, this%mur_past_Hz_right, &
                                                this%mur_right_Hz_ii, this%mur_right_Hz_ij, &
                                                this%mur_right_Hz_ji, this%mur_right_Hz_jj, &
                                                this%mur_right_Hz_ki, this%mur_right_Hz_kj)
    end subroutine gpu_update_mur_past_right

    subroutine gpu_update_mur_past_Hx_right_kernel(Hx_d, past_Hx_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hx_d, past_Hx_d

       integer(kind=4) :: i, k

       !$cuf kernel do(2) 
       do k=ki,kj
          do i=ii,ij
             past_Hx_d(i, ji, k) = Hx_d(i, ji, k)
          end do
       end do
    end subroutine gpu_update_mur_past_Hx_right_kernel

    subroutine gpu_update_mur_past_Hz_right_kernel(Hz_d, past_Hz_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hz_d, past_Hz_d

       integer(kind=4) :: i, k

       !$cuf kernel do(2) 
       do k=ki,kj
          do i=ii,ij
             past_Hz_d(i, ji, k) = Hz_d(i, ji, k)
          end do
       end do
    end subroutine gpu_update_mur_past_Hz_right_kernel

    subroutine gpu_update_mur_past_down(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call gpu_update_mur_past_Hy_down_kernel(this%Hy_d, this%mur_past_Hy_down, &
                                               this%mur_down_Hy_ii, this%mur_down_Hy_ij, &
                                               this%mur_down_Hy_ji, this%mur_down_Hy_jj, &
                                               this%mur_down_Hy_ki, this%mur_down_Hy_kj)

       call gpu_update_mur_past_Hx_down_kernel(this%Hx_d, this%mur_past_Hx_down, &
                                               this%mur_down_Hx_ii, this%mur_down_Hx_ij, &
                                               this%mur_down_Hx_ji, this%mur_down_Hx_jj, &
                                               this%mur_down_Hx_ki, this%mur_down_Hx_kj)
    end subroutine gpu_update_mur_past_down

    subroutine gpu_update_mur_past_Hy_down_kernel(Hy_d, past_Hy_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hy_d, past_Hy_d

       integer(kind=4) :: i, k

       !$cuf kernel do(2) 
       do k=ki,kj
          do i=ii,ij
             past_Hy_d(i, ji, k) = Hy_d(i, ji, k)
          end do
       end do
    end subroutine gpu_update_mur_past_Hy_down_kernel

    subroutine gpu_update_mur_past_Hx_down_kernel(Hx_d, past_Hx_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hx_d, past_Hx_d

       integer(kind=4) :: i, k

       !$cuf kernel do(2) 
       do k=ki,kj
          do i=ii,ij
             past_Hx_d(i, ji, k) = Hx_d(i, ji, k)
          end do
       end do
    end subroutine gpu_update_mur_past_Hx_down_kernel

    subroutine gpu_update_mur_past_up(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call gpu_update_mur_past_Hy_up_kernel(this%Hy_d, this%mur_past_Hy_up, &
                                             this%mur_up_Hy_ii, this%mur_up_Hy_ij, &
                                             this%mur_up_Hy_ji, this%mur_up_Hy_jj, &
                                             this%mur_up_Hy_ki, this%mur_up_Hy_kj)

       call gpu_update_mur_past_Hx_up_kernel(this%Hx_d, this%mur_past_Hx_up, &
                                             this%mur_up_Hx_ii, this%mur_up_Hx_ij, &
                                             this%mur_up_Hx_ji, this%mur_up_Hx_jj, &
                                             this%mur_up_Hx_ki, this%mur_up_Hx_kj)
    end subroutine gpu_update_mur_past_up

    subroutine gpu_update_mur_past_Hy_up_kernel(Hy_d, past_Hy_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hy_d, past_Hy_d

       integer(kind=4) :: i, k

       !$cuf kernel do(2) 
       do k=ki,kj
          do i=ii,ij
             past_Hy_d(i, ji, k) = Hy_d(i, ji, k)
          end do
       end do
    end subroutine gpu_update_mur_past_Hy_up_kernel

    subroutine gpu_update_mur_past_Hx_up_kernel(Hx_d, past_Hx_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hx_d, past_Hx_d

       integer(kind=4) :: i, k

       !$cuf kernel do(2) 
       do k=ki,kj
          do i=ii,ij
             past_Hx_d(i, ji, k) = Hx_d(i, ji, k)
          end do
       end do
    end subroutine gpu_update_mur_past_Hx_up_kernel

    subroutine gpu_update_mur_past_back(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call gpu_update_mur_past_Hz_back_kernel(this%Hz_d, this%mur_past_Hz_back, &
                                               this%mur_back_Hz_ii, this%mur_back_Hz_ij, &
                                               this%mur_back_Hz_ji, this%mur_back_Hz_jj, &
                                               this%mur_back_Hz_ki, this%mur_back_Hz_kj)

       call gpu_update_mur_past_Hy_back_kernel(this%Hy_d, this%mur_past_Hy_back, &
                                               this%mur_back_Hy_ii, this%mur_back_Hy_ij, &
                                               this%mur_back_Hy_ji, this%mur_back_Hy_jj, &
                                               this%mur_back_Hy_ki, this%mur_back_Hy_kj)
    end subroutine gpu_update_mur_past_back

    subroutine gpu_update_mur_past_Hz_back_kernel(Hz_d, past_Hz_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hz_d, past_Hz_d

       integer(kind=4) :: i, j

       !$cuf kernel do(2) 
       do j=ji,jj
          do i=ii,ij
             past_Hz_d(i, j, kj) = Hz_d(i, j, kj)
          end do
       end do
    end subroutine gpu_update_mur_past_Hz_back_kernel

    subroutine gpu_update_mur_past_Hy_back_kernel(Hy_d, past_Hy_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hy_d, past_Hy_d

       integer(kind=4) :: i, j

       !$cuf kernel do(2) 
       do j=ji,jj
          do i=ii,ij
             past_Hy_d(i, j, kj) = Hy_d(i, j, kj)
          end do
       end do
    end subroutine gpu_update_mur_past_Hy_back_kernel

    subroutine gpu_update_mur_past_front(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call gpu_update_mur_past_Hz_front_kernel(this%Hz_d, this%mur_past_Hz_front, &
                                                this%mur_front_Hz_ii, this%mur_front_Hz_ij, &
                                                this%mur_front_Hz_ji, this%mur_front_Hz_jj, &
                                                this%mur_front_Hz_ki, this%mur_front_Hz_kj)

       call gpu_update_mur_past_Hy_front_kernel(this%Hy_d, this%mur_past_Hy_front, &
                                                this%mur_front_Hy_ii, this%mur_front_Hy_ij, &
                                                this%mur_front_Hy_ji, this%mur_front_Hy_jj, &
                                                this%mur_front_Hy_ki, this%mur_front_Hy_kj)
    end subroutine gpu_update_mur_past_front

    subroutine gpu_update_mur_past_Hz_front_kernel(Hz_d, past_Hz_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hz_d, past_Hz_d

       integer(kind=4) :: i, j

       !$cuf kernel do(2) 
       do j=ji,jj
          do i=ii,ij
             past_Hz_d(i, j, kj) = Hz_d(i, j, kj)
          end do
       end do
    end subroutine gpu_update_mur_past_Hz_front_kernel

    subroutine gpu_update_mur_past_Hy_front_kernel(Hy_d, past_Hy_d, ii, ij, ji, jj, ki, kj)
       integer(kind=4), intent(in) :: ii, ij, ji, jj, ki, kj
       real(kind=rkind), device, dimension(:,:,:) :: Hy_d, past_Hy_d

       integer(kind=4) :: i, j

       !$cuf kernel do(2) 
       do j=ji,jj
          do i=ii,ij
             past_Hy_d(i, j, kj) = Hy_d(i, j, kj)
          end do
       end do
   end subroutine gpu_update_mur_past_Hy_front_kernel

    !--------------------------------------------------------------------------------
    ! FUSED MUR KERNELS — 3 kernels instead of 24 (advance + update_past per field)
    ! Reduces kernel launches from 300K to ~30K (3 per timestep)
    !--------------------------------------------------------------------------------

    ! Fused advance for Hx: handles left, right, down, up boundaries
    subroutine gpu_fused_mur_advance_hx(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call fused_mur_advance_hx_kernel( &
          this%Hx_d, this%sggMiHx_d, &
          this%mur_past_Hx_left, this%mur_past_Hx_right, &
          this%mur_past_Hx_down, this%mur_past_Hx_up, &
          this%mur_left_Hx_ii, this%mur_left_Hx_ij, this%mur_left_Hx_ji, this%mur_left_Hx_jj, this%mur_left_Hx_ki, this%mur_left_Hx_kj, &
          this%mur_right_Hx_ii, this%mur_right_Hx_ij, this%mur_right_Hx_ji, this%mur_right_Hx_jj, this%mur_right_Hx_ki, this%mur_right_Hx_kj, &
          this%mur_down_Hx_ii, this%mur_down_Hx_ij, this%mur_down_Hx_ji, this%mur_down_Hx_jj, this%mur_down_Hx_ki, this%mur_down_Hx_kj, &
          this%mur_up_Hx_ii, this%mur_up_Hx_ij, this%mur_up_Hx_ji, this%mur_up_Hx_jj, this%mur_up_Hx_ki, this%mur_up_Hx_kj, &
          this%mur_left_CAB1, this%mur_right_CAB1, this%mur_down_CAB1, this%mur_up_CAB1)

    end subroutine gpu_fused_mur_advance_hx

    ! Fused advance for Hy: handles down, up, back, front boundaries
    subroutine gpu_fused_mur_advance_hy(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call fused_mur_advance_hy_kernel( &
          this%Hy_d, this%sggMiHy_d, &
          this%mur_past_Hy_down, this%mur_past_Hy_up, &
          this%mur_past_Hy_back, this%mur_past_Hy_front, &
          this%mur_down_Hy_ii, this%mur_down_Hy_ij, this%mur_down_Hy_ji, this%mur_down_Hy_jj, this%mur_down_Hy_ki, this%mur_down_Hy_kj, &
          this%mur_up_Hy_ii, this%mur_up_Hy_ij, this%mur_up_Hy_ji, this%mur_up_Hy_jj, this%mur_up_Hy_ki, this%mur_up_Hy_kj, &
          this%mur_back_Hy_ii, this%mur_back_Hy_ij, this%mur_back_Hy_ji, this%mur_back_Hy_jj, this%mur_back_Hy_ki, this%mur_back_Hy_kj, &
          this%mur_front_Hy_ii, this%mur_front_Hy_ij, this%mur_front_Hy_ji, this%mur_front_Hy_jj, this%mur_front_Hy_ki, this%mur_front_Hy_kj, &
          this%mur_down_CAB1, this%mur_up_CAB1, this%mur_back_CAB1, this%mur_front_CAB1)

    end subroutine gpu_fused_mur_advance_hy

    ! Fused advance for Hz: handles left, right, back, front boundaries
    subroutine gpu_fused_mur_advance_hz(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call fused_mur_advance_hz_kernel( &
          this%Hz_d, this%sggMiHz_d, &
          this%mur_past_Hz_left, this%mur_past_Hz_right, &
          this%mur_past_Hz_back, this%mur_past_Hz_front, &
          this%mur_left_Hz_ii, this%mur_left_Hz_ij, this%mur_left_Hz_ji, this%mur_left_Hz_jj, this%mur_left_Hz_ki, this%mur_left_Hz_kj, &
          this%mur_right_Hz_ii, this%mur_right_Hz_ij, this%mur_right_Hz_ji, this%mur_right_Hz_jj, this%mur_right_Hz_ki, this%mur_right_Hz_kj, &
          this%mur_back_Hz_ii, this%mur_back_Hz_ij, this%mur_back_Hz_ji, this%mur_back_Hz_jj, this%mur_back_Hz_ki, this%mur_back_Hz_kj, &
          this%mur_front_Hz_ii, this%mur_front_Hz_ij, this%mur_front_Hz_ji, this%mur_front_Hz_jj, this%mur_front_Hz_ki, this%mur_front_Hz_kj, &
          this%mur_left_CAB1, this%mur_right_CAB1, this%mur_back_CAB1, this%mur_front_CAB1)

    end subroutine gpu_fused_mur_advance_hz

    ! Fused update_past for Hx: left, right, down, up
    subroutine gpu_fused_mur_update_past_hx(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call fused_mur_update_past_hx_kernel( &
          this%Hx_d, &
          this%mur_past_Hx_left, this%mur_past_Hx_right, &
          this%mur_past_Hx_down, this%mur_past_Hx_up, &
          this%mur_left_Hx_ii, this%mur_left_Hx_ij, this%mur_left_Hx_ji, this%mur_left_Hx_jj, this%mur_left_Hx_ki, this%mur_left_Hx_kj, &
          this%mur_right_Hx_ii, this%mur_right_Hx_ij, this%mur_right_Hx_ji, this%mur_right_Hx_jj, this%mur_right_Hx_ki, this%mur_right_Hx_kj, &
          this%mur_down_Hx_ii, this%mur_down_Hx_ij, this%mur_down_Hx_ji, this%mur_down_Hx_jj, this%mur_down_Hx_ki, this%mur_down_Hx_kj, &
          this%mur_up_Hx_ii, this%mur_up_Hx_ij, this%mur_up_Hx_ji, this%mur_up_Hx_jj, this%mur_up_Hx_ki, this%mur_up_Hx_kj)

    end subroutine gpu_fused_mur_update_past_hx

    ! Fused update_past for Hy: down, up, back, front
    subroutine gpu_fused_mur_update_past_hy(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call fused_mur_update_past_hy_kernel( &
          this%Hy_d, &
          this%mur_past_Hy_down, this%mur_past_Hy_up, &
          this%mur_past_Hy_back, this%mur_past_Hy_front, &
          this%mur_down_Hy_ii, this%mur_down_Hy_ij, this%mur_down_Hy_ji, this%mur_down_Hy_jj, this%mur_down_Hy_ki, this%mur_down_Hy_kj, &
          this%mur_up_Hy_ii, this%mur_up_Hy_ij, this%mur_up_Hy_ji, this%mur_up_Hy_jj, this%mur_up_Hy_ki, this%mur_up_Hy_kj, &
          this%mur_back_Hy_ii, this%mur_back_Hy_ij, this%mur_back_Hy_ji, this%mur_back_Hy_jj, this%mur_back_Hy_ki, this%mur_back_Hy_kj, &
          this%mur_front_Hy_ii, this%mur_front_Hy_ij, this%mur_front_Hy_ji, this%mur_front_Hy_jj, this%mur_front_Hy_ki, this%mur_front_Hy_kj)

    end subroutine gpu_fused_mur_update_past_hy

    ! Fused update_past for Hz: left, right, back, front
    subroutine gpu_fused_mur_update_past_hz(this, b)
       class(gpu_state_t), intent(inout) :: this
       type(bounds_t), intent(in) :: b

       if (.not. this%mur_initialized) return

       call fused_mur_update_past_hz_kernel( &
          this%Hz_d, &
          this%mur_past_Hz_left, this%mur_past_Hz_right, &
          this%mur_past_Hz_back, this%mur_past_Hz_front, &
          this%mur_left_Hz_ii, this%mur_left_Hz_ij, this%mur_left_Hz_ji, this%mur_left_Hz_jj, this%mur_left_Hz_ki, this%mur_left_Hz_kj, &
          this%mur_right_Hz_ii, this%mur_right_Hz_ij, this%mur_right_Hz_ji, this%mur_right_Hz_jj, this%mur_right_Hz_ki, this%mur_right_Hz_kj, &
          this%mur_back_Hz_ii, this%mur_back_Hz_ij, this%mur_back_Hz_ji, this%mur_back_Hz_jj, this%mur_back_Hz_ki, this%mur_back_Hz_kj, &
          this%mur_front_Hz_ii, this%mur_front_Hz_ij, this%mur_front_Hz_ji, this%mur_front_Hz_jj, this%mur_front_Hz_ki, this%mur_front_Hz_kj)

    end subroutine gpu_fused_mur_update_past_hz

    !===============================================================================
    ! Fused advance Hx kernel — 4 boundaries: left, right, down, up
    !===============================================================================
    subroutine fused_mur_advance_hx_kernel(Hx_d, sggMiHx_d, &
                                           past_Hx_left, past_Hx_right, past_Hx_down, past_Hx_up, &
                                           l_ii, l_ij, l_ji, l_jj, l_ki, l_kj, &
                                           r_ii, r_ij, r_ji, r_jj, r_ki, r_kj, &
                                           d_ii, d_ij, d_ji, d_jj, d_ki, d_kj, &
                                           u_ii, u_ij, u_ji, u_jj, u_ki, u_kj, &
                                           left_CAB1, right_CAB1, down_CAB1, up_CAB1)

       real(kind=rkind), device, dimension(:,:,:) :: Hx_d
       integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHx_d
       real(kind=rkind), device, dimension(:,:,:) :: past_Hx_left, past_Hx_right, past_Hx_down, past_Hx_up
       real(kind=rkind), device, dimension(:) :: left_CAB1, right_CAB1, down_CAB1, up_CAB1
       integer(kind=4), intent(in) :: l_ii, l_ij, l_ji, l_jj, l_ki, l_kj
       integer(kind=4), intent(in) :: r_ii, r_ij, r_ji, r_jj, r_ki, r_kj
       integer(kind=4), intent(in) :: d_ii, d_ij, d_ji, d_jj, d_ki, d_kj
       integer(kind=4), intent(in) :: u_ii, u_ij, u_ji, u_jj, u_ki, u_kj

       integer(kind=4) :: i, j, k
       integer(kind=integersizeofmediamatrices) :: medio

       !$cuf kernel do(2) 
       do k=l_ki,l_kj
          do i=l_ii,l_ij
             medio = sggMiHx_d(i, l_ji+1, k)
             Hx_d(i, l_ji, k) = past_Hx_left(i, l_ji+1, k) + &
                  left_CAB1(medio) * (Hx_d(i, l_ji+1, k) - past_Hx_left(i, l_ji, k))
          end do
       end do

       !$cuf kernel do(2) 
       do k=r_ki,r_kj
          do i=r_ii,r_ij
             medio = sggMiHx_d(i, r_ji-1, k)
             Hx_d(i, r_ji, k) = past_Hx_right(i, r_ji-1, k) + &
                  right_CAB1(medio) * (Hx_d(i, r_ji-1, k) - past_Hx_right(i, r_ji, k))
          end do
       end do

       !$cuf kernel do(2) 
       do j=d_ji,d_jj
          do i=d_ii,d_ij
             medio = sggMiHx_d(i, j, d_ki+1)
             Hx_d(i, j, d_ki) = past_Hx_down(i, j, d_ki+1) + &
                  down_CAB1(medio) * (Hx_d(i, j, d_ki+1) - past_Hx_down(i, j, d_ki))
          end do
       end do

       !$cuf kernel do(2) 
       do j=u_ji,u_jj
          do i=u_ii,u_ij
             medio = sggMiHx_d(i, j, u_kj-1)
             Hx_d(i, j, u_kj) = past_Hx_up(i, j, u_kj-1) + &
                  up_CAB1(medio) * (Hx_d(i, j, u_kj-1) - past_Hx_up(i, j, u_kj))
          end do
       end do

    end subroutine fused_mur_advance_hx_kernel

    !===============================================================================
    ! Fused advance Hy kernel — 4 boundaries: down, up, back, front
    !===============================================================================
    subroutine fused_mur_advance_hy_kernel(Hy_d, sggMiHy_d, &
                                           past_Hy_down, past_Hy_up, past_Hy_back, past_Hy_front, &
                                           d_ii, d_ij, d_ji, d_jj, d_ki, d_kj, &
                                           u_ii, u_ij, u_ji, u_jj, u_ki, u_kj, &
                                           b_ii, b_ij, b_ji, b_jj, b_ki, b_kj, &
                                           f_ii, f_ij, f_ji, f_jj, f_ki, f_kj, &
                                           down_CAB1, up_CAB1, back_CAB1, front_CAB1)

       real(kind=rkind), device, dimension(:,:,:) :: Hy_d
       integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHy_d
       real(kind=rkind), device, dimension(:,:,:) :: past_Hy_down, past_Hy_up, past_Hy_back, past_Hy_front
       real(kind=rkind), device, dimension(:) :: down_CAB1, up_CAB1, back_CAB1, front_CAB1
       integer(kind=4), intent(in) :: d_ii, d_ij, d_ji, d_jj, d_ki, d_kj
       integer(kind=4), intent(in) :: u_ii, u_ij, u_ji, u_jj, u_ki, u_kj
       integer(kind=4), intent(in) :: b_ii, b_ij, b_ji, b_jj, b_ki, b_kj
       integer(kind=4), intent(in) :: f_ii, f_ij, f_ji, f_jj, f_ki, f_kj

       integer(kind=4) :: i, j, k
       integer(kind=integersizeofmediamatrices) :: medio

       !$cuf kernel do(2) 
       do j=d_ji,d_jj
          do i=d_ii,d_ij
             medio = sggMiHy_d(i, j, d_ki+1)
             Hy_d(i, j, d_ki) = past_Hy_down(i, j, d_ki+1) + &
                  down_CAB1(medio) * (Hy_d(i, j, d_ki+1) - past_Hy_down(i, j, d_ki))
          end do
       end do

       !$cuf kernel do(2) 
       do j=u_ji,u_jj
          do i=u_ii,u_ij
             medio = sggMiHy_d(i, j, u_kj-1)
             Hy_d(i, j, u_kj) = past_Hy_up(i, j, u_kj-1) + &
                  up_CAB1(medio) * (Hy_d(i, j, u_kj-1) - past_Hy_up(i, j, u_kj))
          end do
       end do

       !$cuf kernel do(2) 
       do j=b_ji,b_jj
          do i=b_ii,b_ij
             medio = sggMiHy_d(i, j, b_ki-1)
             Hy_d(i, j, b_ki) = past_Hy_back(i, j, b_ki-1) + &
                  back_CAB1(medio) * (Hy_d(i, j, b_ki-1) - past_Hy_back(i, j, b_ki))
          end do
       end do

       !$cuf kernel do(2) 
       do j=f_ji,f_jj
          do i=f_ii,f_ij
             medio = sggMiHy_d(i, j, f_kj+1)
             Hy_d(i, j, f_kj) = past_Hy_front(i, j, f_kj+1) + &
                  front_CAB1(medio) * (Hy_d(i, j, f_kj+1) - past_Hy_front(i, j, f_kj))
          end do
       end do

    end subroutine fused_mur_advance_hy_kernel

    !===============================================================================
    ! Fused advance Hz kernel — 4 boundaries: left, right, back, front
    !===============================================================================
    subroutine fused_mur_advance_hz_kernel(Hz_d, sggMiHz_d, &
                                           past_Hz_left, past_Hz_right, past_Hz_back, past_Hz_front, &
                                           l_ii, l_ij, l_ji, l_jj, l_ki, l_kj, &
                                           r_ii, r_ij, r_ji, r_jj, r_ki, r_kj, &
                                           b_ii, b_ij, b_ji, b_jj, b_ki, b_kj, &
                                           f_ii, f_ij, f_ji, f_jj, f_ki, f_kj, &
                                           left_CAB1, right_CAB1, back_CAB1, front_CAB1)

       real(kind=rkind), device, dimension(:,:,:) :: Hz_d
       integer(kind=integersizeofmediamatrices), device, dimension(:,:,:) :: sggMiHz_d
       real(kind=rkind), device, dimension(:,:,:) :: past_Hz_left, past_Hz_right, past_Hz_back, past_Hz_front
       real(kind=rkind), device, dimension(:) :: left_CAB1, right_CAB1, back_CAB1, front_CAB1
       integer(kind=4), intent(in) :: l_ii, l_ij, l_ji, l_jj, l_ki, l_kj
       integer(kind=4), intent(in) :: r_ii, r_ij, r_ji, r_jj, r_ki, r_kj
       integer(kind=4), intent(in) :: b_ii, b_ij, b_ji, b_jj, b_ki, b_kj
       integer(kind=4), intent(in) :: f_ii, f_ij, f_ji, f_jj, f_ki, f_kj

       integer(kind=4) :: i, j, k
       integer(kind=integersizeofmediamatrices) :: medio

       !$cuf kernel do(2) 
       do k=l_ki,l_kj
          do i=l_ii,l_ij
             medio = sggMiHz_d(i, l_ji+1, k)
             Hz_d(i, l_ji, k) = past_Hz_left(i, l_ji+1, k) + &
                  left_CAB1(medio) * (Hz_d(i, l_ji+1, k) - past_Hz_left(i, l_ji, k))
          end do
       end do

       !$cuf kernel do(2) 
       do k=r_ki,r_kj
          do i=r_ii,r_ij
             medio = sggMiHz_d(i, r_ji-1, k)
             Hz_d(i, r_ji, k) = past_Hz_right(i, r_ji-1, k) + &
                  right_CAB1(medio) * (Hz_d(i, r_ji-1, k) - past_Hz_right(i, r_ji, k))
          end do
       end do

       !$cuf kernel do(2) 
       do j=b_ji,b_jj
          do i=b_ii,b_ij
             medio = sggMiHz_d(i, j, b_ki-1)
             Hz_d(i, j, b_ki) = past_Hz_back(i, j, b_ki-1) + &
                  back_CAB1(medio) * (Hz_d(i, j, b_ki-1) - past_Hz_back(i, j, b_ki))
          end do
       end do

       !$cuf kernel do(2) 
       do j=f_ji,f_jj
          do i=f_ii,f_ij
             medio = sggMiHz_d(i, j, f_kj+1)
             Hz_d(i, j, f_kj) = past_Hz_front(i, j, f_kj+1) + &
                  front_CAB1(medio) * (Hz_d(i, j, f_kj+1) - past_Hz_front(i, j, f_kj))
          end do
       end do

    end subroutine fused_mur_advance_hz_kernel

    !===============================================================================
    ! Fused update_past Hx kernel — 4 boundaries: left, right, down, up
    !===============================================================================
    subroutine fused_mur_update_past_hx_kernel(Hx_d, &
                                               past_Hx_left, past_Hx_right, past_Hx_down, past_Hx_up, &
                                               l_ii, l_ij, l_ji, l_jj, l_ki, l_kj, &
                                               r_ii, r_ij, r_ji, r_jj, r_ki, r_kj, &
                                               d_ii, d_ij, d_ji, d_jj, d_ki, d_kj, &
                                               u_ii, u_ij, u_ji, u_jj, u_ki, u_kj)

       real(kind=rkind), device, dimension(:,:,:) :: Hx_d
       real(kind=rkind), device, dimension(:,:,:) :: past_Hx_left, past_Hx_right, past_Hx_down, past_Hx_up
       integer(kind=4), intent(in) :: l_ii, l_ij, l_ji, l_jj, l_ki, l_kj
       integer(kind=4), intent(in) :: r_ii, r_ij, r_ji, r_jj, r_ki, r_kj
       integer(kind=4), intent(in) :: d_ii, d_ij, d_ji, d_jj, d_ki, d_kj
       integer(kind=4), intent(in) :: u_ii, u_ij, u_ji, u_jj, u_ki, u_kj

       integer(kind=4) :: i, j, k

       !$cuf kernel do(2) 
       do k=l_ki,l_kj
          do i=l_ii,l_ij
             past_Hx_left(i, l_ji, k) = Hx_d(i, l_ji, k)
          end do
       end do

       !$cuf kernel do(2) 
       do k=r_ki,r_kj
          do i=r_ii,r_ij
             past_Hx_right(i, r_ji, k) = Hx_d(i, r_ji, k)
          end do
       end do

       !$cuf kernel do(2) 
       do k=d_ki,d_kj
          do i=d_ii,d_ij
             past_Hx_down(i, d_ji, k) = Hx_d(i, d_ji, k)
          end do
       end do

       !$cuf kernel do(2) 
       do k=u_ki,u_kj
          do i=u_ii,u_ij
             past_Hx_up(i, u_ji, k) = Hx_d(i, u_ji, k)
          end do
       end do

    end subroutine fused_mur_update_past_hx_kernel

    !===============================================================================
    ! Fused update_past Hy kernel — 4 boundaries: down, up, back, front
    !===============================================================================
    subroutine fused_mur_update_past_hy_kernel(Hy_d, &
                                               past_Hy_down, past_Hy_up, past_Hy_back, past_Hy_front, &
                                               d_ii, d_ij, d_ji, d_jj, d_ki, d_kj, &
                                               u_ii, u_ij, u_ji, u_jj, u_ki, u_kj, &
                                               b_ii, b_ij, b_ji, b_jj, b_ki, b_kj, &
                                               f_ii, f_ij, f_ji, f_jj, f_ki, f_kj)

       real(kind=rkind), device, dimension(:,:,:) :: Hy_d
       real(kind=rkind), device, dimension(:,:,:) :: past_Hy_down, past_Hy_up, past_Hy_back, past_Hy_front
       integer(kind=4), intent(in) :: d_ii, d_ij, d_ji, d_jj, d_ki, d_kj
       integer(kind=4), intent(in) :: u_ii, u_ij, u_ji, u_jj, u_ki, u_kj
       integer(kind=4), intent(in) :: b_ii, b_ij, b_ji, b_jj, b_ki, b_kj
       integer(kind=4), intent(in) :: f_ii, f_ij, f_ji, f_jj, f_ki, f_kj

       integer(kind=4) :: i, j, k

       !$cuf kernel do(2) 
       do k=d_ki,d_kj
          do i=d_ii,d_ij
             past_Hy_down(i, d_ji, k) = Hy_d(i, d_ji, k)
          end do
       end do

       !$cuf kernel do(2) 
       do k=u_ki,u_kj
          do i=u_ii,u_ij
             past_Hy_up(i, u_ji, k) = Hy_d(i, u_ji, k)
          end do
       end do

       !$cuf kernel do(2) 
       do j=b_ji,b_jj
          do i=b_ii,b_ij
             past_Hy_back(i, j, b_kj) = Hy_d(i, j, b_kj)
          end do
       end do

       !$cuf kernel do(2) 
       do j=f_ji,f_jj
          do i=f_ii,f_ij
             past_Hy_front(i, j, f_kj) = Hy_d(i, j, f_kj)
          end do
       end do

    end subroutine fused_mur_update_past_hy_kernel

    !===============================================================================
    ! Fused update_past Hz kernel — 4 boundaries: left, right, back, front
    !===============================================================================
    subroutine fused_mur_update_past_hz_kernel(Hz_d, &
                                               past_Hz_left, past_Hz_right, past_Hz_back, past_Hz_front, &
                                               l_ii, l_ij, l_ji, l_jj, l_ki, l_kj, &
                                               r_ii, r_ij, r_ji, r_jj, r_ki, r_kj, &
                                               b_ii, b_ij, b_ji, b_jj, b_ki, b_kj, &
                                               f_ii, f_ij, f_ji, f_jj, f_ki, f_kj)

       real(kind=rkind), device, dimension(:,:,:) :: Hz_d
       real(kind=rkind), device, dimension(:,:,:) :: past_Hz_left, past_Hz_right, past_Hz_back, past_Hz_front
       integer(kind=4), intent(in) :: l_ii, l_ij, l_ji, l_jj, l_ki, l_kj
       integer(kind=4), intent(in) :: r_ii, r_ij, r_ji, r_jj, r_ki, r_kj
       integer(kind=4), intent(in) :: b_ii, b_ij, b_ji, b_jj, b_ki, b_kj
       integer(kind=4), intent(in) :: f_ii, f_ij, f_ji, f_jj, f_ki, f_kj

       integer(kind=4) :: i, j, k

       !$cuf kernel do(2) 
       do k=l_ki,l_kj
          do i=l_ii,l_ij
             past_Hz_left(i, l_ji, k) = Hz_d(i, l_ji, k)
          end do
       end do

       !$cuf kernel do(2) 
       do k=r_ki,r_kj
          do i=r_ii,r_ij
             past_Hz_right(i, r_ji, k) = Hz_d(i, r_ji, k)
          end do
       end do

       !$cuf kernel do(2) 
       do j=b_ji,b_jj
          do i=b_ii,b_ij
             past_Hz_back(i, j, b_kj) = Hz_d(i, j, b_kj)
          end do
       end do

       !$cuf kernel do(2) 
       do j=f_ji,f_jj
          do i=f_ii,f_ij
             past_Hz_front(i, j, f_kj) = Hz_d(i, j, f_kj)
          end do
       end do

    end subroutine fused_mur_update_past_hz_kernel

  end module gpu_mur_m
