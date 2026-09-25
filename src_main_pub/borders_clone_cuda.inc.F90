! Device PMC (sign -1) then periodic (sign +1) magnetic ghost clones.

logical function cuda_clone_is_ready()
   cuda_clone_is_ready = cuda_clone_ready
end function

subroutine InitMagneticClone_cuda(sgg, b)
   type(SGGFDTDINFO_t), intent(in) :: sgg
   type(bounds_t), intent(in) :: b
   type(fdtd_clone_job_c) :: jobs(24)
   integer :: n, irc

   cuda_clone_ready = .false.
   if (.not. fdtd_cuda_ok_f()) return
   n = 0
   call add_pmc(n, jobs, sgg, b)
   call add_periodic(n, jobs, sgg, b)
   if (n < 1) return
   irc = fdtd_cuda_set_clone_jobs_f(jobs, n)
   if (irc == 0) return
   cuda_clone_ready = .true.
end subroutine InitMagneticClone_cuda

subroutine AdvanceMagneticClone_cuda()
   integer :: irc
   if (.not. cuda_clone_ready) return
   irc = fdtd_cuda_advance_clones_f()
   if (irc == 0) call stoponerror(0, 0, 'CUDA magnetic clone failed')
end subroutine AdvanceMagneticClone_cuda

subroutine add_one(n, jobs, comp, wall, ghost, source, sgn, a0, a1, b0, b1, ex, ey, ez)
   integer, intent(inout) :: n
   type(fdtd_clone_job_c), intent(inout) :: jobs(24)
   integer, intent(in) :: comp, wall, ghost, source, sgn, a0, a1, b0, b1, ex, ey, ez
   if (n >= 24) return
   n = n + 1
   jobs(n)%field_comp = comp
   jobs(n)%wall_axis = wall
   jobs(n)%ghost = ghost
   jobs(n)%source = source
   jobs(n)%sign = sgn
   jobs(n)%a0 = a0
   jobs(n)%a1 = a1
   jobs(n)%b0 = b0
   jobs(n)%b1 = b1
   jobs(n)%e_xi = ex
   jobs(n)%e_yi = ey
   jobs(n)%e_zi = ez
end subroutine

subroutine add_pmc(n, jobs, sgg, b)
   integer, intent(inout) :: n
   type(fdtd_clone_job_c), intent(inout) :: jobs(24)
   type(SGGFDTDINFO_t), intent(in) :: sgg
   type(bounds_t), intent(in) :: b
   ! Same face order as MinusCloneMagneticPMC.
   if (sgg%Border%IsDownPMC) &
      call add_one(n, jobs, 3, 2, sgg%Sweep(iHx)%ZI - 1, sgg%Sweep(iHx)%ZI, -1, &
         b%Hx%XI, b%Hx%XE, b%Hx%YI, b%Hx%YE, b%Hx%XI, b%Hx%YI, b%Hx%ZI)
   if (sgg%Border%IsUpPMC) &
      call add_one(n, jobs, 3, 2, sgg%Sweep(iHx)%ZE + 1, sgg%Sweep(iHx)%ZE, -1, &
         b%Hx%XI, b%Hx%XE, b%Hx%YI, b%Hx%YE, b%Hx%XI, b%Hx%YI, b%Hx%ZI)
   if (sgg%Border%IsLeftPMC) &
      call add_one(n, jobs, 3, 1, sgg%Sweep(iHx)%YI - 1, sgg%Sweep(iHx)%YI, -1, &
         b%Hx%XI, b%Hx%XE, b%Hx%ZI, b%Hx%ZE, b%Hx%XI, b%Hx%YI, b%Hx%ZI)
   if (sgg%Border%IsRightPMC) &
      call add_one(n, jobs, 3, 1, sgg%Sweep(iHx)%YE + 1, sgg%Sweep(iHx)%YE, -1, &
         b%Hx%XI, b%Hx%XE, b%Hx%ZI, b%Hx%ZE, b%Hx%XI, b%Hx%YI, b%Hx%ZI)
   if (sgg%Border%IsBackPMC) &
      call add_one(n, jobs, 4, 0, sgg%Sweep(iHy)%XI - 1, sgg%Sweep(iHy)%XI, -1, &
         b%Hy%YI, b%Hy%YE, b%Hy%ZI, b%Hy%ZE, b%Hy%XI, b%Hy%YI, b%Hy%ZI)
   if (sgg%Border%IsFrontPMC) &
      call add_one(n, jobs, 4, 0, sgg%Sweep(iHy)%XE + 1, sgg%Sweep(iHy)%XE, -1, &
         b%Hy%YI, b%Hy%YE, b%Hy%ZI, b%Hy%ZE, b%Hy%XI, b%Hy%YI, b%Hy%ZI)
   if (sgg%Border%IsDownPMC) &
      call add_one(n, jobs, 4, 2, sgg%Sweep(iHy)%ZI - 1, sgg%Sweep(iHy)%ZI, -1, &
         b%Hy%XI, b%Hy%XE, b%Hy%YI, b%Hy%YE, b%Hy%XI, b%Hy%YI, b%Hy%ZI)
   if (sgg%Border%IsUpPMC) &
      call add_one(n, jobs, 4, 2, sgg%Sweep(iHy)%ZE + 1, sgg%Sweep(iHy)%ZE, -1, &
         b%Hy%XI, b%Hy%XE, b%Hy%YI, b%Hy%YE, b%Hy%XI, b%Hy%YI, b%Hy%ZI)
   if (sgg%Border%IsBackPMC) &
      call add_one(n, jobs, 5, 0, sgg%Sweep(iHz)%XI - 1, sgg%Sweep(iHz)%XI, -1, &
         b%Hz%YI, b%Hz%YE, b%Hz%ZI, b%Hz%ZE, b%Hz%XI, b%Hz%YI, b%Hz%ZI)
   if (sgg%Border%IsFrontPMC) &
      call add_one(n, jobs, 5, 0, sgg%Sweep(iHz)%XE + 1, sgg%Sweep(iHz)%XE, -1, &
         b%Hz%YI, b%Hz%YE, b%Hz%ZI, b%Hz%ZE, b%Hz%XI, b%Hz%YI, b%Hz%ZI)
   if (sgg%Border%IsLeftPMC) &
      call add_one(n, jobs, 5, 1, sgg%Sweep(iHz)%YI - 1, sgg%Sweep(iHz)%YI, -1, &
         b%Hz%XI, b%Hz%XE, b%Hz%ZI, b%Hz%ZE, b%Hz%XI, b%Hz%YI, b%Hz%ZI)
   if (sgg%Border%IsRightPMC) &
      call add_one(n, jobs, 5, 1, sgg%Sweep(iHz)%YE + 1, sgg%Sweep(iHz)%YE, -1, &
         b%Hz%XI, b%Hz%XE, b%Hz%ZI, b%Hz%ZE, b%Hz%XI, b%Hz%YI, b%Hz%ZI)
end subroutine

subroutine add_periodic(n, jobs, sgg, b)
   integer, intent(inout) :: n
   type(fdtd_clone_job_c), intent(inout) :: jobs(24)
   type(SGGFDTDINFO_t), intent(in) :: sgg
   type(bounds_t), intent(in) :: b
   ! Same face order as CloneMagneticPeriodic. Source is the opposite interior plane.
   if (sgg%Border%IsDownPeriodic) &
      call add_one(n, jobs, 3, 2, sgg%Sweep(iHx)%ZI - 1, sgg%Sweep(iHx)%ZE, +1, &
         b%Hx%XI, b%Hx%XE, b%Hx%YI, b%Hx%YE, b%Hx%XI, b%Hx%YI, b%Hx%ZI)
   if (sgg%Border%IsUpPeriodic) &
      call add_one(n, jobs, 3, 2, sgg%Sweep(iHx)%ZE + 1, sgg%Sweep(iHx)%ZI, +1, &
         b%Hx%XI, b%Hx%XE, b%Hx%YI, b%Hx%YE, b%Hx%XI, b%Hx%YI, b%Hx%ZI)
   if (sgg%Border%IsLeftPeriodic) &
      call add_one(n, jobs, 3, 1, sgg%Sweep(iHx)%YI - 1, sgg%Sweep(iHx)%YE, +1, &
         b%Hx%XI, b%Hx%XE, b%Hx%ZI, b%Hx%ZE, b%Hx%XI, b%Hx%YI, b%Hx%ZI)
   if (sgg%Border%IsRightPeriodic) &
      call add_one(n, jobs, 3, 1, sgg%Sweep(iHx)%YE + 1, sgg%Sweep(iHx)%YI, +1, &
         b%Hx%XI, b%Hx%XE, b%Hx%ZI, b%Hx%ZE, b%Hx%XI, b%Hx%YI, b%Hx%ZI)
   if (sgg%Border%IsBackPeriodic) &
      call add_one(n, jobs, 4, 0, sgg%Sweep(iHy)%XI - 1, sgg%Sweep(iHy)%XE, +1, &
         b%Hy%YI, b%Hy%YE, b%Hy%ZI, b%Hy%ZE, b%Hy%XI, b%Hy%YI, b%Hy%ZI)
   if (sgg%Border%IsFrontPeriodic) &
      call add_one(n, jobs, 4, 0, sgg%Sweep(iHy)%XE + 1, sgg%Sweep(iHy)%XI, +1, &
         b%Hy%YI, b%Hy%YE, b%Hy%ZI, b%Hy%ZE, b%Hy%XI, b%Hy%YI, b%Hy%ZI)
   if (sgg%Border%IsDownPeriodic) &
      call add_one(n, jobs, 4, 2, sgg%Sweep(iHy)%ZI - 1, sgg%Sweep(iHy)%ZE, +1, &
         b%Hy%XI, b%Hy%XE, b%Hy%YI, b%Hy%YE, b%Hy%XI, b%Hy%YI, b%Hy%ZI)
   if (sgg%Border%IsUpPeriodic) &
      call add_one(n, jobs, 4, 2, sgg%Sweep(iHy)%ZE + 1, sgg%Sweep(iHy)%ZI, +1, &
         b%Hy%XI, b%Hy%XE, b%Hy%YI, b%Hy%YE, b%Hy%XI, b%Hy%YI, b%Hy%ZI)
   if (sgg%Border%IsBackPeriodic) &
      call add_one(n, jobs, 5, 0, sgg%Sweep(iHz)%XI - 1, sgg%Sweep(iHz)%XE, +1, &
         b%Hz%YI, b%Hz%YE, b%Hz%ZI, b%Hz%ZE, b%Hz%XI, b%Hz%YI, b%Hz%ZI)
   if (sgg%Border%IsFrontPeriodic) &
      call add_one(n, jobs, 5, 0, sgg%Sweep(iHz)%XE + 1, sgg%Sweep(iHz)%XI, +1, &
         b%Hz%YI, b%Hz%YE, b%Hz%ZI, b%Hz%ZE, b%Hz%XI, b%Hz%YI, b%Hz%ZI)
   if (sgg%Border%IsLeftPeriodic) &
      call add_one(n, jobs, 5, 1, sgg%Sweep(iHz)%YI - 1, sgg%Sweep(iHz)%YE, +1, &
         b%Hz%XI, b%Hz%XE, b%Hz%ZI, b%Hz%ZE, b%Hz%XI, b%Hz%YI, b%Hz%ZI)
   if (sgg%Border%IsRightPeriodic) &
      call add_one(n, jobs, 5, 1, sgg%Sweep(iHz)%YE + 1, sgg%Sweep(iHz)%YI, +1, &
         b%Hz%XI, b%Hz%XE, b%Hz%ZI, b%Hz%ZE, b%Hz%XI, b%Hz%YI, b%Hz%ZI)
end subroutine
