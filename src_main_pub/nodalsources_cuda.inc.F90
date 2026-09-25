! Device nodal sources (CompileWithCUDA). Host AdvanceNodalE/H stays the CPU path.

#ifdef CompileWithCUDA
logical function cuda_nodal_is_ready()
   cuda_nodal_is_ready = cuda_nodal_ready
end function

subroutine InitNodalSources_cuda(sgg, b)
   type(SGGFDTDINFO_t), intent(in) :: sgg
   type(bounds_t), intent(in) :: b
   type(fdtd_nodal_job_c), allocatable :: jobs(:)
   real(kind=RKIND), allocatable :: samples(:)
   integer(c_int), allocatable :: skip_e(:), skip_h(:)
   integer :: n_jobs, n_samples, m, irc, cap

   cuda_nodal_ready = .false.
   if (.not. fdtd_cuda_ok_f()) return

   cap = Nodal_Ex%numHard + Nodal_Ex%numSoft + Nodal_Ey%numHard + Nodal_Ey%numSoft + &
         Nodal_Ez%numHard + Nodal_Ez%numSoft + Nodal_Hx%numHard + Nodal_Hx%numSoft + &
         Nodal_Hy%numHard + Nodal_Hy%numSoft + Nodal_Hz%numHard + Nodal_Hz%numSoft
   if (cap < 1) return

   allocate(jobs(cap))
   n_samples = count_samples()
   if (n_samples < 1) return
   allocate(samples(0:n_samples-1))
   n_jobs = 0
   n_samples = 0
   call push_comp(Nodal_Ex, 0, b%Ex%XI, b%Ex%YI, b%Ex%ZI)
   call push_comp(Nodal_Ey, 1, b%Ey%XI, b%Ey%YI, b%Ey%ZI)
   call push_comp(Nodal_Ez, 2, b%Ez%XI, b%Ez%YI, b%Ez%ZI)
   call push_comp(Nodal_Hx, 3, b%Hx%XI, b%Hx%YI, b%Hx%ZI)
   call push_comp(Nodal_Hy, 4, b%Hy%XI, b%Hy%YI, b%Hy%ZI)
   call push_comp(Nodal_Hz, 5, b%Hz%XI, b%Hz%YI, b%Hz%ZI)
   if (n_jobs < 1) return

   allocate(skip_e(0:sgg%NumMedia), skip_h(0:sgg%NumMedia))
   do m = 0, sgg%NumMedia
      skip_e(m) = 0
      skip_h(m) = 0
      if (sgg%Med(m)%Is%PEC) skip_e(m) = 1
      if (sgg%Med(m)%Is%PMC) skip_h(m) = 1
   end do

   irc = fdtd_cuda_upload_nodal_f(jobs, n_jobs, samples, n_samples, skip_e, skip_h, sgg%NumMedia + 1)
   if (irc == 0) return
   cuda_nodal_ready = .true.

contains

   integer function count_samples()
      count_samples = span(Nodal_Ex) + span(Nodal_Ey) + span(Nodal_Ez) + &
                      span(Nodal_Hx) + span(Nodal_Hy) + span(Nodal_Hz)
   end function

   integer function span(nod)
      type(nodsou_t), intent(in) :: nod
      integer :: ii
      span = 0
      if (associated(nod%nodHard)) then
         do ii = 1, nod%numHard
            span = span + nod%nodHard(ii)%numus + 1
         end do
      end if
      if (associated(nod%nodSoft)) then
         do ii = 1, nod%numSoft
            span = span + nod%nodSoft(ii)%numus + 1
         end do
      end if
   end function

   subroutine push_comp(nod, comp, ox, oy, oz)
      type(nodsou_t), intent(in) :: nod
      integer, intent(in) :: comp, ox, oy, oz
      integer :: ii
      if (associated(nod%nodHard)) then
         do ii = 1, nod%numHard
            call push_one(nod%nodHard(ii), comp, 1, ox, oy, oz)
         end do
      end if
      if (associated(nod%nodSoft)) then
         do ii = 1, nod%numSoft
            call push_one(nod%nodSoft(ii), comp, 0, ox, oy, oz)
         end do
      end if
   end subroutine

   subroutine push_one(loc, comp, is_hard, ox, oy, oz)
      type(NodalLocal_t), intent(in) :: loc
      integer, intent(in) :: comp, is_hard, ox, oy, oz
      integer :: k
      if (loc%punto%XI > loc%punto%XE) return
      if (loc%punto%YI > loc%punto%YE) return
      if (loc%punto%ZI > loc%punto%ZE) return
      n_jobs = n_jobs + 1
      jobs(n_jobs)%field_comp = comp
      jobs(n_jobs)%hard = is_hard
      jobs(n_jobs)%initial_only = 0
      if (loc%IsInitialValue) jobs(n_jobs)%initial_only = 1
      jobs(n_jobs)%xi = loc%punto%XI
      jobs(n_jobs)%xe = loc%punto%XE
      jobs(n_jobs)%yi = loc%punto%YI
      jobs(n_jobs)%ye = loc%punto%YE
      jobs(n_jobs)%zi = loc%punto%ZI
      jobs(n_jobs)%ze = loc%punto%ZE
      jobs(n_jobs)%e_xi = ox
      jobs(n_jobs)%e_yi = oy
      jobs(n_jobs)%e_zi = oz
      jobs(n_jobs)%evol_off = n_samples
      jobs(n_jobs)%numus = loc%numus
      jobs(n_jobs)%amplitude = real(loc%punto%amplitude, c_float)
      jobs(n_jobs)%deltaevol = real(loc%deltaevol, c_float)
      do k = 0, loc%numus
         samples(n_samples + k) = loc%evol(k)
      end do
      n_samples = n_samples + loc%numus + 1
   end subroutine
end subroutine InitNodalSources_cuda
#endif
