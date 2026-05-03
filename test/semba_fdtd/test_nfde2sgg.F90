module test_nfde2sgg_m
   use SEMBA_FDTD_m
   use smbjson_m, only: fdtdjson_parser_t => parser_t
   use NFDETypes_m
   use nfde_rotate_m
   use Report_m, only: OnPrint, resetFatalError
   use FDETYPES_m, only: RKIND, set_priorities, setglobal, &
                          EPSILON_VACUUM, MU_VACUUM
   use interpreta_switches_m, only: default_flags

#ifndef __GFORTRAN__
   use IFPORT
#endif

   implicit none

contains

   subroutine change_to_dir(dir, original_dir, ierr)
      character(len=*), intent(in) :: dir
      character(len=1024), intent(out) :: original_dir
      integer, intent(out) :: ierr
#ifdef __GFORTRAN__
      call getcwd(original_dir)
      call chdir(dir, ierr)
#else
      integer :: s
      s = getcwd(original_dir)
      ierr = chdir(dir)
#endif
   end subroutine

   subroutine restore_dir(original_dir)
      character(len=1024), intent(in) :: original_dir
      integer :: s
#ifdef __GFORTRAN__
      call chdir(trim(original_dir), s)
#else
      s = chdir(trim(original_dir))
#endif
   end subroutine

   subroutine setup_fdtd(fdtd)
      type(semba_fdtd_t), intent(inout) :: fdtd
      call initEntrada(fdtd%l)
      call default_flags(fdtd%l)
      fdtd%eps0 = EPSILON_VACUUM
      fdtd%mu0  = MU_VACUUM
      fdtd%cluz = 1.0_RKIND / sqrt(fdtd%eps0 * fdtd%mu0)
      call OnPrint
      fdtd%l%num_procs = 1
      fdtd%l%layoutnumber = 0
      call setglobal(fdtd%l%layoutnumber, fdtd%l%num_procs)
      call set_priorities(.false., .false., .false.)
   end subroutine

   subroutine run_nfde2sgg_on_case(case_dir, json_file, basename, err)
      character(len=*), intent(in) :: case_dir, json_file, basename
      integer, intent(inout) :: err

      type(semba_fdtd_t) :: fdtd
      type(fdtdjson_parser_t) :: json_parser
      type(Parseador_t), pointer :: parser => null()
      real(kind=RKIND) :: dt_before
      character(len=1024) :: original_dir
      integer :: chdir_status

      call change_to_dir(case_dir, original_dir, chdir_status)
      if (chdir_status /= 0) then
         write(*,*) 'FAIL: could not chdir to ', case_dir
         err = err + 1
         return
      end if

      call setup_fdtd(fdtd)
      call resetFatalError()

      json_parser = fdtdjson_parser_t(json_file)
      allocate(parser)
      parser = json_parser%readProblemDescription()
      call nfde_rotate(parser, fdtd%l%mpidir)

      fdtd%l%fichin = basename
      fdtd%l%nEntradaRoot = basename
      fdtd%l%extension = '.json'
      fdtd%sgg%nEntradaRoot = basename

      dt_before = nfde2sgg(fdtd, parser)

      if (fdtd%l%fatalerror .or. fdtd%l%fatalerrornfde2sgg) then
         err = err + 1
         write(*,*) 'FAIL: fatal error in nfde2sgg for ', json_file
      end if
      if (fdtd%sgg%dt <= 0.0_RKIND) then
         err = err + 1
         write(*,*) 'FAIL: sgg%dt is not positive:', fdtd%sgg%dt
      end if
      if (dt_before <= 0.0_RKIND) then
         err = err + 1
         write(*,*) 'FAIL: dt_before is not positive:', dt_before
      end if
      if (fdtd%sgg%dt > dt_before) then
         err = err + 1
         write(*,*) 'FAIL: dt was increased (CFL violation):', fdtd%sgg%dt, '>', dt_before
      end if

      if (associated(parser)) deallocate(parser)
      call restore_dir(original_dir)
   end subroutine

end module test_nfde2sgg_m


integer function test_nfde2sgg_planewave() bind(C) result(err)
   use test_nfde2sgg_m
   err = 0
   call run_nfde2sgg_on_case( &
      'testData/cases/planewave/', &
      'pw-in-box.fdtd.json', &
      'pw-in-box', err)
end function

integer function test_nfde2sgg_planewave_periodic() bind(C) result(err)
   use test_nfde2sgg_m
   err = 0
   call run_nfde2sgg_on_case( &
      'testData/cases/planewave/', &
      'pw-with-periodic.fdtd.json', &
      'pw-with-periodic', err)
end function

integer function test_nfde2sgg_holland() bind(C) result(err)
   use test_nfde2sgg_m
   err = 0
   call run_nfde2sgg_on_case( &
      'testData/cases/holland/', &
      'holland1981.fdtd.json', &
      'holland1981', err)
end function

integer function test_nfde2sgg_sources_voltage() bind(C) result(err)
   use test_nfde2sgg_m
   err = 0
   call run_nfde2sgg_on_case( &
      'testData/cases/sources/', &
      'sources_voltage.fdtd.json', &
      'sources_voltage', err)
end function
