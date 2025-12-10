module mod_domain
   use FDETYPES
   use outputTypes
   implicit none



   interface domain_t
     module procedure new_domain_time, new_domain_freq, new_domain_both
   end interface domain_t

contains
   function new_domain_time(tstart, tstop, tstep) result(new_domain)
      real(kind=RKIND_tiempo), intent(in) :: tstart, tstop, tstep
      type(domain_t) :: new_domain

      new_domain%tstart = tstart
      new_domain%tstop = tstop
      new_domain%tstep = tstep
      new_domain%domainType = TIME_DOMAIN
   end function new_domain_time

   function new_domain_freq(fstart, fstop, fnum, logarithmicSpacing) result(new_domain)
      real(kind=RKIND), intent(in)     :: fstart, fstop
      integer(kind=SINGLE), intent(in) :: fnum
      logical, intent(in)    :: logarithmicSpacing
      type(domain_t) :: new_domain

      new_domain%fstart = fstart
      new_domain%fstop = fstop
      new_domain%fnum = fnum
      new_domain%fstep = (fstop - fstart) / fnum
      new_domain%logarithmicSpacing = logarithmicSpacing
      
      new_domain%domainType = FREQUENCY_DOMAIN

      
   end function new_domain_freq

   function new_domain_both(tstart, tstop, tstep, fstart, fstop, fnum, logarithmicSpacing) result(new_domain)
      real(kind=RKIND_tiempo), intent(in) :: tstart, tstop, tstep
      real(kind=RKIND), intent(in)     :: fstart, fstop
      integer(kind=SINGLE), intent(in) :: fnum
      logical, intent(in)   :: logarithmicSpacing
      type(domain_t) :: new_domain

      new_domain%tstart = tstart
      new_domain%tstop = tstop
      new_domain%tstep = tstep

      new_domain%fstart = fstart
      new_domain%fstop = fstop
      new_domain%fnum = fnum
      new_domain%fstep = (fstop - fstart) / fnum
      new_domain%logarithmicSpacing = logarithmicSpacing

      new_domain%domainType = BOTH_DOMAIN

      
      
   end function new_domain_both

end module mod_domain
