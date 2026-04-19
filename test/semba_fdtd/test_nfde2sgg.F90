integer function test_nfde2sgg_planewave() bind(C) result(err)
   use SEMBA_FDTD_m
   use smbjson_m, only: fdtdjson_parser_t => parser_t
   use NFDETypes_m
   use nfde_rotate_m
   use Report_m, only: OnPrint
   use FDETYPES_m, only: RKIND, set_priorities, setglobal
   use interpreta_switches_m, only: default_flags

   implicit none

   character(len=*), parameter :: TEST_CASE_DIR = 'testData/cases/planewave/'
   character(len=*), parameter :: testfile = 'pw-in-box.fdtd.json'
   character(len=*), parameter :: testbase = 'pw-in-box'

   type(semba_fdtd_t) :: fdtd
   type(fdtdjson_parser_t) :: json_parser
   type(Parseador_t), pointer :: parser => null()
   real(kind=RKIND) :: dtantesdecorregir
   character(len=1024) :: original_dir
   integer :: chdir_status

   err = 0

   call getcwd(original_dir)
   call chdir(TEST_CASE_DIR, chdir_status)
   if (chdir_status /= 0) then
      write(*,*) 'FAIL: could not chdir to ', TEST_CASE_DIR
      err = 1
      return
   end if

   call initEntrada(fdtd%l)
   call default_flags(fdtd%l)

   fdtd%eps0 = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
   fdtd%mu0  = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
   fdtd%cluz = 1.0_RKIND / sqrt(fdtd%eps0 * fdtd%mu0)

   call OnPrint
   fdtd%l%num_procs = 1
   fdtd%l%layoutnumber = 0
   call setglobal(fdtd%l%layoutnumber, fdtd%l%num_procs)
   call set_priorities(.false., .false., .false.)

   json_parser = fdtdjson_parser_t(testfile)
   allocate(parser)
   parser = json_parser%readProblemDescription()

   call nfde_rotate(parser, fdtd%l%mpidir)

   fdtd%l%fichin = testbase
   fdtd%l%nEntradaRoot = testbase
   fdtd%l%extension = '.json'
   fdtd%sgg%nEntradaRoot = testbase

   call fdtd%nfde2sgg(parser, dtantesdecorregir)

   if (fdtd%l%fatalerror .or. fdtd%l%fatalerrornfde2sgg) then
      err = err + 1
      write(*,*) 'FAIL: fatal error in nfde2sgg'
   end if

   if (fdtd%sgg%dt <= 0.0_RKIND) then
      err = err + 1
      write(*,*) 'FAIL: sgg%dt is not positive:', fdtd%sgg%dt
   end if

   if (dtantesdecorregir <= 0.0_RKIND) then
      err = err + 1
      write(*,*) 'FAIL: dtantesdecorregir is not positive:', dtantesdecorregir
   end if

   if (associated(parser)) deallocate(parser)

   call chdir(trim(original_dir), chdir_status)

end function
