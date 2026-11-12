module lorentz_testingTools_m
   use FDETYPES_m
   use SEMBA_FDTD_m
   implicit none
   private
   public :: test_lorentz_res2_s11, test_lorentz_plasma_s11

   ! Yee node (Ex) observed in the committed 1D-equivalent Lorentz slab cases.
   integer, parameter :: PROBE_I = 1, PROBE_J = 1, PROBE_K = 200

   real(kind=RKIND), dimension(:), allocatable, save :: incident_samples
   integer(kind=4), save :: number_of_steps = 0
   real(kind=RKIND), save :: time_step = 0.0_RKIND
   logical, save :: incident_cached = .false.

   real(kind=RKIND), parameter :: MAX_ERROR_DB = 0.5_RKIND
   real(kind=RKIND), parameter :: MAX_ERROR_DEG = 12.0_RKIND

contains

   ! Run one case step by step and return the Ex samples at the observation node.
   subroutine run_case(case_file, samples)
      character(len=*), intent(in) :: case_file
      real(kind=RKIND), dimension(:), allocatable, intent(out) :: samples
      type(semba_fdtd_t) :: semba
      type(solver_t) :: solver

      call semba%init("-i "//trim(case_file))
      solver = semba%create_solver()
#ifdef CompileWithMTLN
      solver%mtln_parsed = semba%mtln_parsed
#endif
      call solver%init()

      number_of_steps = solver%control%finaltimestep
      time_step = solver%sgg%dt
      allocate(samples(0:number_of_steps))
      samples = 0.0_RKIND

      do while (solver%n <= solver%control%finaltimestep)
         call solver%step()
         samples(solver%n) = solver%get_field_value(IEX, PROBE_I, PROBE_J, PROBE_K)
         solver%n = solver%n + 1
      end do

      call solver%destroy_and_deallocate()
   end subroutine

   ! The incident field comes from an all-vacuum run with the same mesh, source and node.
   subroutine cache_incident_field()
      if (incident_cached) return
      call run_case('vacuum_reference.fdtd.json', incident_samples)
      incident_cached = .true.
   end subroutine

   ! S11(f) = DFT(E_total - E_incident) / DFT(E_incident) at one frequency.
   subroutine s11_at_frequency(total, incident, frequency, value)
      real(kind=RKIND), intent(in) :: total(0:), incident(0:), frequency
      complex(kind=CKIND), intent(out) :: value
      real(kind=8), parameter :: two_pi = 6.283185307179586476925286766559d0
      real(kind=8) :: phase
      complex(kind=CKIND) :: weight, numerator, denominator
      integer :: n

      numerator = (0.0_RKIND, 0.0_RKIND)
      denominator = (0.0_RKIND, 0.0_RKIND)
      do n = 0, ubound(total, 1)
         phase = two_pi*real(frequency, 8)*real(n, 8)*real(time_step, 8)
         weight = cmplx(cos(phase), -sin(phase), CKIND)
         numerator = numerator + cmplx(total(n) - incident(n), 0.0_RKIND, CKIND)*weight
         denominator = denominator + cmplx(incident(n), 0.0_RKIND, CKIND)*weight
      end do
      value = numerator/denominator
   end subroutine

   real(kind=RKIND) function wrapped_degrees(angle) result(res)
      real(kind=RKIND), intent(in) :: angle
      res = modulo(angle + 180.0_RKIND, 360.0_RKIND) - 180.0_RKIND
   end function

   ! Compare the simulated S11 against the committed dgtd and scikit-rf vectors.
   integer function compare_case(case_name) result(errors)
      character(len=*), intent(in) :: case_name
      real(kind=RKIND), dimension(:), allocatable :: total
      real(kind=RKIND) :: frequency, skrf_db, skrf_phase, dgtd_db, dgtd_phase
      real(kind=RKIND) :: db, phase, phase_error
      complex(kind=CKIND) :: value
      character(len=256) :: header
      character(len=256) :: reference_file
      integer :: unit, ios

      errors = 0
      call cache_incident_field()
      call run_case(trim(case_name)//'.fdtd.json', total)
      if (size(total) /= size(incident_samples)) then
         write (*, *) 'Lorentz test: case and reference time windows differ'
         errors = errors + 1
         return
      end if

      reference_file = 'references/'//trim(case_name)//'_s11.csv'
      open (newunit=unit, file=trim(reference_file), status='old', action='read', iostat=ios)
      if (ios /= 0) then
         write (*, *) 'Lorentz test: cannot open '//trim(reference_file)
         errors = errors + 1
         return
      end if
      read (unit, '(a)') header

      do
         read (unit, *, iostat=ios) frequency, skrf_db, skrf_phase, dgtd_db, dgtd_phase
         if (ios /= 0) exit
         call s11_at_frequency(total, incident_samples, frequency, value)
         db = 20.0_RKIND*log10(abs(value))
         phase = atan2(aimag(value), real(value))*180.0_RKIND/3.14159265358979323846_RKIND

         if (abs(db - skrf_db) > MAX_ERROR_DB) then
            write (*, '(a,f9.3,a,f8.3,a,f8.3)') 'Lorentz scikit-rf |S11| mismatch at ', &
               frequency/1.0e6_RKIND, ' MHz: fdtd ', db, ' dB, reference ', skrf_db
            errors = errors + 1
         end if
         if (abs(db - dgtd_db) > MAX_ERROR_DB) then
            write (*, '(a,f9.3,a,f8.3,a,f8.3)') 'Lorentz dgtd |S11| mismatch at ', &
               frequency/1.0e6_RKIND, ' MHz: fdtd ', db, ' dB, reference ', dgtd_db
            errors = errors + 1
         end if
         phase_error = wrapped_degrees(phase - skrf_phase)
         if (abs(phase_error) > MAX_ERROR_DEG) then
            write (*, '(a,f9.3,a,f8.2,a,f8.2)') 'Lorentz scikit-rf S11 phase mismatch at ', &
               frequency/1.0e6_RKIND, ' MHz: fdtd ', phase, ' deg, reference ', skrf_phase
            errors = errors + 1
         end if
         phase_error = wrapped_degrees(phase - dgtd_phase)
         if (abs(phase_error) > MAX_ERROR_DEG) then
            write (*, '(a,f9.3,a,f8.2,a,f8.2)') 'Lorentz dgtd S11 phase mismatch at ', &
               frequency/1.0e6_RKIND, ' MHz: fdtd ', phase, ' deg, reference ', dgtd_phase
            errors = errors + 1
         end if
      end do
      close (unit)
   end function

   ! The slab cases are PEC-backed, so S11 fully characterizes them (S21 = 0).
   integer function test_lorentz_res2_s11() bind(C) result(errors)
      errors = 0
      call chdir("./testData/cases/lorentz/")
      errors = compare_case('lorentz_res2')
      call chdir("../../..")
   end function

   integer function test_lorentz_plasma_s11() bind(C) result(errors)
      errors = 0
      call chdir("./testData/cases/lorentz/")
      errors = compare_case('lorentz_plasma')
      call chdir("../../..")
   end function

end module
