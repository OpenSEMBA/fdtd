! Device first-order magnetic Mur ABC (CompileWithCUDA).

#ifdef CompileWithCUDA
logical function cuda_mur_is_ready()
   cuda_mur_is_ready = cuda_mur_ready
end function

subroutine InitMURBorders_cuda(b)
   use fdtd_cuda_m
   type(bounds_t), intent(in) :: b
   type(fdtd_mur_job_c) :: jobs(12)
   integer :: n, irc, slot, cabn

   cuda_mur_ready = .false.
   if (.not. fdtd_cuda_ok_f()) return
   if (.not. allocated(left_CAB1)) return

   cabn = size(left_CAB1)
   irc = fdtd_cuda_upload_mur_cab_f(0, cabn, left_CAB1);  if (irc == 0) return
   irc = fdtd_cuda_upload_mur_cab_f(1, cabn, right_CAB1); if (irc == 0) return
   irc = fdtd_cuda_upload_mur_cab_f(2, cabn, down_CAB1);  if (irc == 0) return
   irc = fdtd_cuda_upload_mur_cab_f(3, cabn, up_CAB1);    if (irc == 0) return
   irc = fdtd_cuda_upload_mur_cab_f(4, cabn, back_CAB1);  if (irc == 0) return
   irc = fdtd_cuda_upload_mur_cab_f(5, cabn, front_CAB1); if (irc == 0) return

   n = 0
   ! left: ghost at YI, interior +Y
   call try_add(n, jobs, 3, 0, 1, +1, MURc(iHx)%YI(left), &
      MURc(iHx)%XI(left), MURc(iHx)%XE(left), MURc(iHx)%ZI(left), MURc(iHx)%ZE(left), &
      b%Hx%XI, b%Hx%YI, b%Hx%ZI, regLR(left)%Past_Hx, &
      MURc(iHx)%XI(left), MURc(iHx)%XE(left), MURc(iHx)%YI(left), MURc(iHx)%YE(left), &
      MURc(iHx)%ZI(left), MURc(iHx)%ZE(left))
   call try_add(n, jobs, 5, 0, 1, +1, MURc(iHz)%YI(left), &
      MURc(iHz)%XI(left), MURc(iHz)%XE(left), MURc(iHz)%ZI(left), MURc(iHz)%ZE(left), &
      b%Hz%XI, b%Hz%YI, b%Hz%ZI, regLR(left)%Past_Hz, &
      MURc(iHz)%XI(left), MURc(iHz)%XE(left), MURc(iHz)%YI(left), MURc(iHz)%YE(left), &
      MURc(iHz)%ZI(left), MURc(iHz)%ZE(left))
   ! right: ghost at YE, interior -Y
   call try_add(n, jobs, 3, 1, 1, -1, MURc(iHx)%YE(right), &
      MURc(iHx)%XI(right), MURc(iHx)%XE(right), MURc(iHx)%ZI(right), MURc(iHx)%ZE(right), &
      b%Hx%XI, b%Hx%YI, b%Hx%ZI, regLR(right)%Past_Hx, &
      MURc(iHx)%XI(right), MURc(iHx)%XE(right), MURc(iHx)%YI(right), MURc(iHx)%YE(right), &
      MURc(iHx)%ZI(right), MURc(iHx)%ZE(right))
   call try_add(n, jobs, 5, 1, 1, -1, MURc(iHz)%YE(right), &
      MURc(iHz)%XI(right), MURc(iHz)%XE(right), MURc(iHz)%ZI(right), MURc(iHz)%ZE(right), &
      b%Hz%XI, b%Hz%YI, b%Hz%ZI, regLR(right)%Past_Hz, &
      MURc(iHz)%XI(right), MURc(iHz)%XE(right), MURc(iHz)%YI(right), MURc(iHz)%YE(right), &
      MURc(iHz)%ZI(right), MURc(iHz)%ZE(right))
   ! down: ghost at ZI, interior +Z
   call try_add(n, jobs, 4, 2, 2, +1, MURc(iHy)%ZI(down), &
      MURc(iHy)%XI(down), MURc(iHy)%XE(down), MURc(iHy)%YI(down), MURc(iHy)%YE(down), &
      b%Hy%XI, b%Hy%YI, b%Hy%ZI, regDU(down)%Past_Hy, &
      MURc(iHy)%XI(down), MURc(iHy)%XE(down), MURc(iHy)%YI(down), MURc(iHy)%YE(down), &
      MURc(iHy)%ZI(down), MURc(iHy)%ZE(down))
   call try_add(n, jobs, 3, 2, 2, +1, MURc(iHx)%ZI(down), &
      MURc(iHx)%XI(down), MURc(iHx)%XE(down), MURc(iHx)%YI(down), MURc(iHx)%YE(down), &
      b%Hx%XI, b%Hx%YI, b%Hx%ZI, regDU(down)%Past_Hx, &
      MURc(iHx)%XI(down), MURc(iHx)%XE(down), MURc(iHx)%YI(down), MURc(iHx)%YE(down), &
      MURc(iHx)%ZI(down), MURc(iHx)%ZE(down))
   ! up: ghost at ZE, interior -Z
   call try_add(n, jobs, 4, 3, 2, -1, MURc(iHy)%ZE(up), &
      MURc(iHy)%XI(up), MURc(iHy)%XE(up), MURc(iHy)%YI(up), MURc(iHy)%YE(up), &
      b%Hy%XI, b%Hy%YI, b%Hy%ZI, regDU(up)%Past_Hy, &
      MURc(iHy)%XI(up), MURc(iHy)%XE(up), MURc(iHy)%YI(up), MURc(iHy)%YE(up), &
      MURc(iHy)%ZI(up), MURc(iHy)%ZE(up))
   call try_add(n, jobs, 3, 3, 2, -1, MURc(iHx)%ZE(up), &
      MURc(iHx)%XI(up), MURc(iHx)%XE(up), MURc(iHx)%YI(up), MURc(iHx)%YE(up), &
      b%Hx%XI, b%Hx%YI, b%Hx%ZI, regDU(up)%Past_Hx, &
      MURc(iHx)%XI(up), MURc(iHx)%XE(up), MURc(iHx)%YI(up), MURc(iHx)%YE(up), &
      MURc(iHx)%ZI(up), MURc(iHx)%ZE(up))
   ! back: ghost at XI, interior +X
   call try_add(n, jobs, 5, 4, 0, +1, MURc(iHz)%XI(back), &
      MURc(iHz)%YI(back), MURc(iHz)%YE(back), MURc(iHz)%ZI(back), MURc(iHz)%ZE(back), &
      b%Hz%XI, b%Hz%YI, b%Hz%ZI, regBF(back)%Past_Hz, &
      MURc(iHz)%XI(back), MURc(iHz)%XE(back), MURc(iHz)%YI(back), MURc(iHz)%YE(back), &
      MURc(iHz)%ZI(back), MURc(iHz)%ZE(back))
   call try_add(n, jobs, 4, 4, 0, +1, MURc(iHy)%XI(back), &
      MURc(iHy)%YI(back), MURc(iHy)%YE(back), MURc(iHy)%ZI(back), MURc(iHy)%ZE(back), &
      b%Hy%XI, b%Hy%YI, b%Hy%ZI, regBF(back)%Past_Hy, &
      MURc(iHy)%XI(back), MURc(iHy)%XE(back), MURc(iHy)%YI(back), MURc(iHy)%YE(back), &
      MURc(iHy)%ZI(back), MURc(iHy)%ZE(back))
   ! front: ghost at XE, interior -X
   call try_add(n, jobs, 5, 5, 0, -1, MURc(iHz)%XE(front), &
      MURc(iHz)%YI(front), MURc(iHz)%YE(front), MURc(iHz)%ZI(front), MURc(iHz)%ZE(front), &
      b%Hz%XI, b%Hz%YI, b%Hz%ZI, regBF(front)%Past_Hz, &
      MURc(iHz)%XI(front), MURc(iHz)%XE(front), MURc(iHz)%YI(front), MURc(iHz)%YE(front), &
      MURc(iHz)%ZI(front), MURc(iHz)%ZE(front))
   call try_add(n, jobs, 4, 5, 0, -1, MURc(iHy)%XE(front), &
      MURc(iHy)%YI(front), MURc(iHy)%YE(front), MURc(iHy)%ZI(front), MURc(iHy)%ZE(front), &
      b%Hy%XI, b%Hy%YI, b%Hy%ZI, regBF(front)%Past_Hy, &
      MURc(iHy)%XI(front), MURc(iHy)%XE(front), MURc(iHy)%YI(front), MURc(iHy)%YE(front), &
      MURc(iHy)%ZI(front), MURc(iHy)%ZE(front))

   if (n < 1) return
   irc = fdtd_cuda_set_mur_jobs_f(jobs, n)
   if (irc == 0) return

   slot = 0
   call try_upload(slot, irc, regLR(left)%Past_Hx);  if (irc == 0) return
   call try_upload(slot, irc, regLR(left)%Past_Hz);  if (irc == 0) return
   call try_upload(slot, irc, regLR(right)%Past_Hx); if (irc == 0) return
   call try_upload(slot, irc, regLR(right)%Past_Hz); if (irc == 0) return
   call try_upload(slot, irc, regDU(down)%Past_Hy);  if (irc == 0) return
   call try_upload(slot, irc, regDU(down)%Past_Hx);  if (irc == 0) return
   call try_upload(slot, irc, regDU(up)%Past_Hy);    if (irc == 0) return
   call try_upload(slot, irc, regDU(up)%Past_Hx);    if (irc == 0) return
   call try_upload(slot, irc, regBF(back)%Past_Hz);  if (irc == 0) return
   call try_upload(slot, irc, regBF(back)%Past_Hy);  if (irc == 0) return
   call try_upload(slot, irc, regBF(front)%Past_Hz); if (irc == 0) return
   call try_upload(slot, irc, regBF(front)%Past_Hy); if (irc == 0) return
   if (slot /= n) return

   cuda_mur_ready = .true.
contains
   subroutine try_add(njobs, jobs, field_comp, cab_which, wall_axis, neigh_sign, plane_abs, &
                      a0, a1, b0, b1, e_xi, e_yi, e_zi, past, &
                      sxi, sxe, syi, sye, szi, sze)
      integer, intent(inout) :: njobs
      type(fdtd_mur_job_c), intent(inout) :: jobs(:)
      integer, intent(in) :: field_comp, cab_which, wall_axis, neigh_sign, plane_abs
      integer, intent(in) :: a0, a1, b0, b1, e_xi, e_yi, e_zi
      real(kind=RKIND), pointer, intent(in) :: past(:,:,:)
      integer, intent(in) :: sxi, sxe, syi, sye, szi, sze
      integer :: lb(3), ub(3)
      if (.not. associated(past)) return
      if (size(past, 1) < 1 .or. size(past, 2) < 1 .or. size(past, 3) < 1) return
      if (szi > sze) return
      if (a1 < a0 .or. b1 < b0) return
      njobs = njobs + 1
      lb = lbound(past)
      ub = ubound(past)
      jobs(njobs)%field_comp = field_comp
      jobs(njobs)%media_comp = field_comp
      jobs(njobs)%cab_which = cab_which
      jobs(njobs)%wall_axis = wall_axis
      jobs(njobs)%neigh_sign = neigh_sign
      jobs(njobs)%plane_abs = plane_abs
      jobs(njobs)%a0 = a0
      jobs(njobs)%a1 = a1
      jobs(njobs)%b0 = b0
      jobs(njobs)%b1 = b1
      jobs(njobs)%e_xi = e_xi
      jobs(njobs)%e_yi = e_yi
      jobs(njobs)%e_zi = e_zi
      jobs(njobs)%p_xi = lb(1)
      jobs(njobs)%p_yi = lb(2)
      jobs(njobs)%p_zi = lb(3)
      jobs(njobs)%nx_p = ub(1) - lb(1) + 1
      jobs(njobs)%ny_p = ub(2) - lb(2) + 1
      jobs(njobs)%nz_p = ub(3) - lb(3) + 1
      jobs(njobs)%store_xi = sxi
      jobs(njobs)%store_xe = sxe
      jobs(njobs)%store_yi = syi
      jobs(njobs)%store_ye = sye
      jobs(njobs)%store_zi = szi
      jobs(njobs)%store_ze = sze
   end subroutine

   subroutine try_upload(slot, irc, past)
      integer, intent(inout) :: slot, irc
      real(kind=RKIND), pointer, intent(in) :: past(:,:,:)
      irc = 1
      if (.not. associated(past)) return
      if (size(past, 1) < 1 .or. size(past, 2) < 1 .or. size(past, 3) < 1) return
      irc = fdtd_cuda_upload_mur_past_f(slot, past, size(past))
      if (irc /= 0) slot = slot + 1
   end subroutine
end subroutine InitMURBorders_cuda

subroutine AdvanceMagneticMUR_cuda()
   use fdtd_cuda_m
   integer :: irc
   if (.not. cuda_mur_ready) return
   irc = fdtd_cuda_advance_mur_f()
   if (irc == 0) then
      cuda_mur_ready = .false.
      call stoponerror(0, 0, 'CUDA Mur ABC advance failed')
   end if
end subroutine AdvanceMagneticMUR_cuda
#endif
