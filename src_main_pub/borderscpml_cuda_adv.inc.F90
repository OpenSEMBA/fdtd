! Appended CUDA CPML advance routines for BORDERS_CPML_m

#ifdef CompileWithCUDA
logical function cuda_cpml_is_ready()
   cuda_cpml_is_ready = cuda_cpml_ready
end function

subroutine InitCPMLBorders_cuda()
   use fdtd_cuda_m
   integer :: region, n, irc, slot
   integer :: lb(3), ub(3)

   cuda_cpml_ready = .false.
   if (.not. fdtd_cuda_ok_f()) return

   irc = fdtd_cuda_alloc_cpml_1d_f(0, size(P_be_x), P_be_x); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(1, size(P_be_y), P_be_y); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(2, size(P_be_z), P_be_z); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(6, size(P_ce_x), P_ce_x); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(7, size(P_ce_y), P_ce_y); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(8, size(P_ce_z), P_ce_z); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(12, size(P_bm_x), P_bm_x); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(13, size(P_bm_y), P_bm_y); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(14, size(P_bm_z), P_bm_z); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(18, size(P_cm_x), P_cm_x); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(19, size(P_cm_y), P_cm_y); if (irc==0) return
   irc = fdtd_cuda_alloc_cpml_1d_f(20, size(P_cm_z), P_cm_z); if (irc==0) return

   do region = left, right
      slot = 0 + (region - left) * 4
      call cuda_upload_psi_assoc(slot+0, regLR(region)%Psi_Exy)
      call cuda_upload_psi_assoc(slot+1, regLR(region)%Psi_Ezy)
      call cuda_upload_psi_assoc(slot+2, regLR(region)%Psi_Hxy)
      call cuda_upload_psi_assoc(slot+3, regLR(region)%Psi_Hzy)
   end do
   do region = down, up
      slot = 8 + (region - down) * 4
      call cuda_upload_psi_assoc(slot+0, regDU(region)%Psi_Eyz)
      call cuda_upload_psi_assoc(slot+1, regDU(region)%Psi_Exz)
      call cuda_upload_psi_assoc(slot+2, regDU(region)%Psi_Hyz)
      call cuda_upload_psi_assoc(slot+3, regDU(region)%Psi_Hxz)
   end do
   do region = back, front
      slot = 16 + (region - back) * 4
      call cuda_upload_psi_assoc(slot+0, regBF(region)%Psi_Ezx)
      call cuda_upload_psi_assoc(slot+1, regBF(region)%Psi_Eyx)
      call cuda_upload_psi_assoc(slot+2, regBF(region)%Psi_Hzx)
      call cuda_upload_psi_assoc(slot+3, regBF(region)%Psi_Hyx)
   end do

   cuda_cpml_ready = .true.
contains
   subroutine cuda_upload_psi_assoc(slot, psi)
      integer, intent(in) :: slot
      real(kind=RKIND), pointer, intent(in) :: psi(:,:,:)
      integer :: n, irc, lb(3), ub(3)
      if (.not. associated(psi)) return
      lb = lbound(psi); ub = ubound(psi)
      n = size(psi)
      irc = fdtd_cuda_alloc_psi_f(slot, n)
      if (irc == 0) return
      irc = fdtd_cuda_upload_psi_f(slot, psi, n)
      cuda_psi_origin(slot,:) = lb
      cuda_psi_shape(slot,:) = ub - lb + 1
   end subroutine
end subroutine InitCPMLBorders_cuda

subroutine fill_cpml_job(job, field_comp, psi_slot, media_comp, h_comp, h_diff_axis, &
                         p_b_which, p_c_which, p_base, free_axis, &
                         xi, xe, yi, ye, zi, ze, e_xi, e_yi, e_zi, h_sign)
   use fdtd_cuda_m
   type(fdtd_cpml_job_c), intent(out) :: job
   integer, intent(in) :: field_comp, psi_slot, media_comp, h_comp, h_diff_axis
   integer, intent(in) :: p_b_which, p_c_which, p_base, free_axis
   integer, intent(in) :: xi, xe, yi, ye, zi, ze, e_xi, e_yi, e_zi, h_sign
   job%field_comp = field_comp
   job%psi_slot = psi_slot
   job%media_comp = media_comp
   job%h_comp_a = h_comp
   job%h_comp_b = h_comp
   job%h_diff_axis = h_diff_axis
   job%p_b_which = p_b_which
   job%p_c_which = p_c_which
   job%p_base = p_base
   job%free_axis = free_axis
   job%xi = xi; job%xe = xe; job%yi = yi; job%ye = ye; job%zi = zi; job%ze = ze
   job%e_xi = e_xi; job%e_yi = e_yi; job%e_zi = e_zi
   job%psi_xi = cuda_psi_origin(psi_slot, 1)
   job%psi_yi = cuda_psi_origin(psi_slot, 2)
   job%psi_zi = cuda_psi_origin(psi_slot, 3)
   job%nx_psi = cuda_psi_shape(psi_slot, 1)
   job%ny_psi = cuda_psi_shape(psi_slot, 2)
   job%nz_psi = cuda_psi_shape(psi_slot, 3)
   job%h_sign = h_sign
   job%use_fixed_medio = 0
   job%medio_fixed = 1
end subroutine

subroutine AdvanceelectricCPML_cuda(NumMedia, b)
   use fdtd_cuda_m
   integer, intent(in) :: NumMedia
   type(bounds_t), intent(in) :: b
   type(fdtd_cpml_job_c) :: jobs(12)
   integer :: region, slot, n, irc
   integer :: pby, pcy, pbz, pcz, pbx, pcx
   if (NumMedia < 0) return
   pby = lbound(P_be_y, 1); pcy = lbound(P_ce_y, 1)
   pbz = lbound(P_be_z, 1); pcz = lbound(P_ce_z, 1)
   pbx = lbound(P_be_x, 1); pcx = lbound(P_ce_x, 1)
   n = 0

   do region = left, right
      slot = 0 + (region - left) * 4
      ! Ex += G2 * Psi_Exy ; dH = Hz - Hz(j-1) ; free axis y
      n = n + 1
      call fill_cpml_job(jobs(n), 0, slot+0, 0, 5, 1, 1, 7, pby, 1, &
         PMLc(iEx)%XI(region), PMLc(iEx)%XE(region), PMLc(iEx)%YI(region), PMLc(iEx)%YE(region), &
         PMLc(iEx)%ZI(region), PMLc(iEx)%ZE(region), b%Ex%XI, b%Ex%YI, b%Ex%ZI, +1)
      ! Ez -= G2 * Psi_Ezy ; dH = Hx - Hx(j-1)
      n = n + 1
      call fill_cpml_job(jobs(n), 2, slot+1, 2, 3, 1, 1, 7, pby, 1, &
         PMLc(iEz)%XI(region), PMLc(iEz)%XE(region), PMLc(iEz)%YI(region), PMLc(iEz)%YE(region), &
         PMLc(iEz)%ZI(region), PMLc(iEz)%ZE(region), b%Ez%XI, b%Ez%YI, b%Ez%ZI, -1)
   end do

   do region = down, up
      slot = 8 + (region - down) * 4
      ! Ey += ; dH = Hx-Hx(k-1) free z
      n = n + 1
      call fill_cpml_job(jobs(n), 1, slot+0, 1, 3, 2, 2, 8, pbz, 2, &
         PMLc(iEy)%XI(region), PMLc(iEy)%XE(region), PMLc(iEy)%YI(region), PMLc(iEy)%YE(region), &
         PMLc(iEy)%ZI(region), PMLc(iEy)%ZE(region), b%Ey%XI, b%Ey%YI, b%Ey%ZI, +1)
      ! Ex -= ; dH = Hy-Hy(k-1)
      n = n + 1
      call fill_cpml_job(jobs(n), 0, slot+1, 0, 4, 2, 2, 8, pbz, 2, &
         PMLc(iEx)%XI(region), PMLc(iEx)%XE(region), PMLc(iEx)%YI(region), PMLc(iEx)%YE(region), &
         PMLc(iEx)%ZI(region), PMLc(iEx)%ZE(region), b%Ex%XI, b%Ex%YI, b%Ex%ZI, -1)
   end do

   do region = back, front
      slot = 16 + (region - back) * 4
      ! Ez += ; dH = Hy-Hy(i-1) free x
      n = n + 1
      call fill_cpml_job(jobs(n), 2, slot+0, 2, 4, 0, 0, 6, pbx, 0, &
         PMLc(iEz)%XI(region), PMLc(iEz)%XE(region), PMLc(iEz)%YI(region), PMLc(iEz)%YE(region), &
         PMLc(iEz)%ZI(region), PMLc(iEz)%ZE(region), b%Ez%XI, b%Ez%YI, b%Ez%ZI, +1)
      ! Ey -= ; dH = Hz-Hz(i-1)
      n = n + 1
      call fill_cpml_job(jobs(n), 1, slot+1, 1, 5, 0, 0, 6, pbx, 0, &
         PMLc(iEy)%XI(region), PMLc(iEy)%XE(region), PMLc(iEy)%YI(region), PMLc(iEy)%YE(region), &
         PMLc(iEy)%ZI(region), PMLc(iEy)%ZE(region), b%Ey%XI, b%Ey%YI, b%Ey%ZI, -1)
   end do
   irc = fdtd_cuda_cpml_apply_n_f(jobs, n)
   if (irc == 0) call stoponerror(0, 0, 'CUDA electric CPML advance failed')
end subroutine AdvanceelectricCPML_cuda

subroutine AdvanceMagneticCPML_cuda(NumMedia, b)
   use fdtd_cuda_m
   integer, intent(in) :: NumMedia
   type(bounds_t), intent(in) :: b
   type(fdtd_cpml_job_c) :: jobs(12)
   integer :: region, slot, n, irc
   integer :: pby, pcy, pbz, pcz, pbx, pcx
   if (NumMedia < 0) return
   pby = lbound(P_bm_y, 1); pcy = lbound(P_cm_y, 1)
   pbz = lbound(P_bm_z, 1); pcz = lbound(P_cm_z, 1)
   pbx = lbound(P_bm_x, 1); pcx = lbound(P_cm_x, 1)
   n = 0

   do region = left, right
      slot = 0 + (region - left) * 4
      ! Hx += Gm2*Psi_Hxy ; dE = Ez - Ez(j-1) — match host signs from AdvanceMagneticCPML
      n = n + 1
      call fill_cpml_job(jobs(n), 3, slot+2, 3, 2, 1, 13, 19, pby, 1, &
         PMLc(iHx)%XI(region), PMLc(iHx)%XE(region), PMLc(iHx)%YI(region), PMLc(iHx)%YE(region), &
         PMLc(iHx)%ZI(region), PMLc(iHx)%ZE(region), b%Hx%XI, b%Hx%YI, b%Hx%ZI, +1)
      n = n + 1
      call fill_cpml_job(jobs(n), 5, slot+3, 5, 0, 1, 13, 19, pby, 1, &
         PMLc(iHz)%XI(region), PMLc(iHz)%XE(region), PMLc(iHz)%YI(region), PMLc(iHz)%YE(region), &
         PMLc(iHz)%ZI(region), PMLc(iHz)%ZE(region), b%Hz%XI, b%Hz%YI, b%Hz%ZI, -1)
   end do
   do region = down, up
      slot = 8 + (region - down) * 4
      n = n + 1
      call fill_cpml_job(jobs(n), 4, slot+2, 4, 0, 2, 14, 20, pbz, 2, &
         PMLc(iHy)%XI(region), PMLc(iHy)%XE(region), PMLc(iHy)%YI(region), PMLc(iHy)%YE(region), &
         PMLc(iHy)%ZI(region), PMLc(iHy)%ZE(region), b%Hy%XI, b%Hy%YI, b%Hy%ZI, +1)
      n = n + 1
      call fill_cpml_job(jobs(n), 3, slot+3, 3, 1, 2, 14, 20, pbz, 2, &
         PMLc(iHx)%XI(region), PMLc(iHx)%XE(region), PMLc(iHx)%YI(region), PMLc(iHx)%YE(region), &
         PMLc(iHx)%ZI(region), PMLc(iHx)%ZE(region), b%Hx%XI, b%Hx%YI, b%Hx%ZI, -1)
   end do
   do region = back, front
      slot = 16 + (region - back) * 4
      n = n + 1
      call fill_cpml_job(jobs(n), 5, slot+2, 5, 1, 0, 12, 18, pbx, 0, &
         PMLc(iHz)%XI(region), PMLc(iHz)%XE(region), PMLc(iHz)%YI(region), PMLc(iHz)%YE(region), &
         PMLc(iHz)%ZI(region), PMLc(iHz)%ZE(region), b%Hz%XI, b%Hz%YI, b%Hz%ZI, +1)
      n = n + 1
      call fill_cpml_job(jobs(n), 4, slot+3, 4, 2, 0, 12, 18, pbx, 0, &
         PMLc(iHy)%XI(region), PMLc(iHy)%XE(region), PMLc(iHy)%YI(region), PMLc(iHy)%YE(region), &
         PMLc(iHy)%ZI(region), PMLc(iHy)%ZE(region), b%Hy%XI, b%Hy%YI, b%Hy%ZI, -1)
   end do
   irc = fdtd_cuda_cpml_apply_n_f(jobs, n)
   if (irc == 0) call stoponerror(0, 0, 'CUDA magnetic CPML advance failed')
end subroutine AdvanceMagneticCPML_cuda
#endif
