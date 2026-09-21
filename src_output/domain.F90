module domain_m
   use FDETYPES_m
   use outputTypes_m
   implicit none

   private
   public :: domain_t

   interface domain_t
      module procedure new_domain_time, new_domain_freq, new_domain_both, null_domain
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
      if (fnum > 1_SINGLE) then
         new_domain%fstep = (fstop - fstart)/real(fnum - 1_SINGLE, RKIND)
      else
         new_domain%fstep = 0.0_RKIND
      end if
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
      if (fnum > 1_SINGLE) then
         new_domain%fstep = (fstop - fstart)/real(fnum - 1_SINGLE, RKIND)
      else
         new_domain%fstep = 0.0_RKIND
      end if
      new_domain%logarithmicSpacing = logarithmicSpacing

      new_domain%domainType = BOTH_DOMAIN

   end function new_domain_both

   function null_domain() result(new_domain)
      type(domain_t) :: new_domain
      new_domain%domainType = UNDEFINED_DOMAIN
   end function

end module domain_m
