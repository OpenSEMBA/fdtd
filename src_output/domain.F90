module mod_domain
   use FDETYPES
   implicit none

   integer, parameter :: UNDEFINED_DOMAIN = -1
   integer, parameter :: TIME_DOMAIN = 0
   integer, parameter :: FREQUENCY_DOMAIN = 1
   integer, parameter :: BOTH_DOMAIN = 2

   interface domain_t
     module procedure new_domain_time, new_domain_freq, new_domain_both
   end interface domain_t
   type :: domain_t
      real(kind=RKIND_tiempo) :: tstart = 0.0_RKIND_tiempo, tstop = 0.0_RKIND_tiempo, tstep = 0.0_RKIND_tiempo
      real(kind=RKIND)        :: fstart = 0.0_RKIND, fstop = 0.0_RKIND, fstep
      integer(kind=SINGLE)    :: fnum = 0
      integer(kind=SINGLE)    :: domainType = UNDEFINED_DOMAIN
      logical                 :: logarithmicSpacing = .false.
   end type domain_t

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
      logical, intent(in), optional    :: logarithmicSpacing
      type(domain_t) :: new_domain

      new_domain%fstart = fstart
      new_domain%fstop = fstop
      new_domain%fnum = fnum
      new_domain%fstep = (fstop - fstart) / fnum
      
      new_domain%domainType = FREQUENCY_DOMAIN

      if (present(logarithmicSpacing)) then
         new_domain%logarithmicSpacing = logarithmicSpacing
      end if
   end function new_domain_freq

   function new_domain_both(tstart, tstop, tstep, fstart, fstop, fnum, logarithmicSpacing) result(new_domain)
      real(kind=RKIND_tiempo), intent(in) :: tstart, tstop, tstep
      real(kind=RKIND), intent(in)     :: fstart, fstop
      integer(kind=SINGLE), intent(in) :: fnum
      logical, intent(in), optional    :: logarithmicSpacing
      type(domain_t) :: new_domain

      new_domain%tstart = tstart
      new_domain%tstop = tstop
      new_domain%tstep = tstep

      new_domain%fstart = fstart
      new_domain%fstop = fstop
      new_domain%fnum = fnum
      new_domain%fstep = (fstop - fstart) / fnum

      new_domain%domainType = BOTH_DOMAIN

      if (present(logarithmicSpacing)) then
         new_domain%logarithmicSpacing = logarithmicSpacing
      end if
   end function new_domain_both

end module mod_domain
