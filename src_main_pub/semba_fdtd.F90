module SEMBA_FDTD_m

   use version_m
   use Report_m
   use Getargs_m
   !
   use FDETYPES_m
   use Solver_m         
   use resuming_m
   !nfde parser stuff
   use NFDETypes_m                
   use nfde_rotate_m           


#ifdef CompilePrivateVersion  
   use ParseadorClass
#endif

#ifdef CompileWithSMBJSON
   use smbjson_m, only: fdtdjson_parser_t => parser_t
#endif

   use Preprocess_m
   use storeData_m
   use xdmf_h5_m
   !
#ifdef CompileWithMPI
   use MPIcomm_m
   use build_t_linea_mpi_m
#ifdef CompileWithStochastic
   use MPI_stochastic
#endif
#endif

   use EpsMuTimeScale_m

   use interpreta_switches_m
   use, intrinsic:: iso_fortran_env, only: stdin=>input_unit

   implicit none

   ! should eps0 and mu0 be global variables?

   type, public :: semba_fdtd_t 
      type(entrada_t) :: l
      type(tiempo_t) :: time_comienzo
      real(kind=8) time_desdelanzamiento
      type(media_matrices_t) :: media
      type(SGGFDTDINFO_t) :: sgg
      type(limit_t), dimension(1:6) :: fullsize, SINPML_fullsize
      real(kind=RKIND) :: eps0,mu0,cluz
      real(kind=RKIND) :: maxSourceValue
      character(len=BUFSIZE) :: whoami, whoamishort
#ifdef CompileWithMTLN
      type(mtln_t) :: mtln_parsed
#endif
      type(taglist_t) :: tag_numbers
      type(tagtype_t) :: tagtype
      logical :: finishedwithsuccess

   contains
      procedure :: init => semba_init
      procedure :: launch => semba_launch
      procedure :: end => semba_end
      procedure :: create_solver => semba_create_solver
      procedure :: update_after_simulation => semba_update_after_simulation
   end type semba_fdtd_t 

   
contains

   subroutine semba_init(this, input_flags)
      class(semba_fdtd_t) :: this
      character(len=*), optional :: input_flags

      real(kind=RKIND) :: dtantesdecorregir
      real(kind=RKIND) :: dxmin,dymin,dzmin,dtlay
      
      logical :: dummylog,l_auxinput, l_auxoutput, ThereArethinslots
      logical :: hayinput
      logical :: lexis

      character(len=BUFSIZE) :: f= ' ', chain = ' ', chain3 = ' ',chain4 = ' ', chaindummy= ' '
      character(len=BUFSIZE_LONG) :: slices = ' '
      character(len=BUFSIZE) :: dubuf
      character(len=BUFSIZE) :: buff
      character(len=BUFSIZE) :: filename_h5bin ! File name

      integer(kind=4) :: myunit,jmed
      integer(kind=4) :: finaltimestepantesdecorregir,NEWfinaltimestep,thefileno
      integer(kind=4) :: statuse
      integer(kind=4) :: status, i, field
      integer(kind=4) :: my_iostat


      type(Parseador_t), pointer :: parser
      type(t_NFDE_FILE_t), pointer :: NFDE_FILE
      type(solver_t) :: solver 
         
#ifdef CompileWithMPI
      LOGICAL :: fatalerror_aux
      type(XYZlimit_t), dimension(1:6) :: tempalloc
#endif

      integer(kind=4) :: conf_err

      call initEntrada(this%l) 

      this%eps0= 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
      this%mu0 = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
      this%cluz=1.0_RKIND/sqrt(this%eps0*this%mu0)
      
      call OnPrint

#ifdef CompileWithMPI
      call InitGeneralMPI (this%l%layoutnumber, this%l%num_procs)
      SUBCOMM_MPI=MPI_COMM_WORLD !default el this%l%stochastic es el global a menos que luego se divida
#else
      this%l%num_procs = 1
      this%l%layoutnumber = 0
#endif
      call setglobal(this%l%layoutnumber,this%l%num_procs) !para crear variables globales con info MPI
         
      write(this%whoamishort, '(i5)') this%l%layoutnumber + 1
      write(this%whoami, '(a,i5,a,i5,a)') '(', this%l%layoutnumber + 1, '/', this%l%num_procs, ') '
         
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      call get_secnds(this%l%time_out2)
      this%time_desdelanzamiento= this%l%time_out2%segundos
#ifndef keeppause
      if (this%l%layoutnumber==0) then
         open(38, file='running')
         write (38,*) '!END'
         CLOSE (38,status='delete')
         open(38, file='pause')
         write (38,*) '!END'
         CLOSE (38,status='delete')
         open(38, file='relaunch')
         write (38,*) '!END'
         CLOSE (38,status='delete')
         open(38, file='forcestop')
         write (38,*) '!END'
         CLOSE (38,status='delete')
      end if
#endif

   if (this%l%layoutnumber==0) then
         my_iostat=0
   3443  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' 
         open(11, file='SEMBA_FDTD_temp.log',err=3443,iostat=my_iostat,action='write')
         write (11,*) '!END'
         CLOSE (11,status='delete')
         my_iostat=0
   3447  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',3447,'.',this%l%layoutnumber,'SEMBA_FDTD_temp.log' 
         open(11, file='SEMBA_FDTD_temp.log',err=3447,iostat=my_iostat,status='new',action='write')
         call print_credits(this%l)
         CLOSE (11)
   end if

#ifdef CompileWithMPI
      call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif


   652 continue

      call CLOSEWARNINGFILE(this%l%layoutnumber,this%l%num_procs,dummylog,.false.,.false.) !aqui ya no se tiene en cuenta el this%l%fatalerror

      write(this%l%opcionespararesumeo, '(a,i4,a)') 'mpirun -n ', this%l%num_procs,' '
      call default_flags(this%l)    !set all default flags

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      call get_secnds(this%time_comienzo)
      !temporarily until later
      if (this%l%layoutnumber == 0) then
         open(11, file='SEMBA_FDTD_temp.log',position='append')
         this%l%file11isopen=.true.
      end if
      !

#ifdef CompileWithMPI
      !wait until everything comes out
      call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif

      !see if there is semaphore to pause continuing
      inquire(file='pause', EXIST=this%l%pausar)
#ifdef CompileWithMPI
      this%l%l_aux = this%l%pausar
      call MPI_AllReduce (this%l%l_aux, this%l%pausar, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, this%l%ierr)
#endif
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      call get_secnds (this%l%time_out2)
      this%l%time_begin = this%l%time_out2%segundos
      write(dubuf,*) 'Paused at              ', this%l%time_out2%fecha(7:8), '/', this%l%time_out2%fecha(5:6), '/', &
      &                this%l%time_out2%fecha(1:4), '  ', this%l%time_out2%hora(1:2), ':', this%l%time_out2%hora(3:4)
      if (this%l%pausar) call print11 (this%l%layoutnumber, dubuf)
      do while (this%l%pausar)
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
         call get_secnds (this%l%time_out2)
         this%l%time_end = this%l%time_out2%segundos
         if (this%l%time_end-this%l%time_begin > 10.0_RKIND) then
            inquire(file='pause', EXIST=this%l%pausar)
#ifdef CompileWithMPI
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
            this%l%l_aux = this%l%pausar
            call MPI_AllReduce (this%l%l_aux, this%l%pausar, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, this%l%ierr)
            call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
            call get_secnds (this%l%time_out2)
            this%l%time_begin = this%l%time_out2%segundos
            write(dubuf,*) 'Paused at              ', this%l%time_out2%fecha(7:8), '/', this%l%time_out2%fecha(5:6), '/', &
            &                this%l%time_out2%fecha(1:4), ' ', this%l%time_out2%hora(1:2), ':', this%l%time_out2%hora(3:4)
            if (this%l%pausar) call print11 (this%l%layoutnumber, dubuf)
         end if
      end do
      !fin del semaphoro

#ifdef keeppause   
      inquire(file='forcestop', EXIST=this%l%forcestop)
      if (this%l%forcestop) then
         if (this%l%layoutnumber==0) then
            open(38, file='running')
            write (38,*) '!END'
            CLOSE (38,status='delete')
            open(38, file='pause')
            write (38,*) '!END'
            CLOSE (38,status='delete')
            open(38, file='relaunch')
            write (38,*) '!END'
            CLOSE (38,status='delete')
            open(38, file='forcestop')
            write (38,*) '!END'
            CLOSE (38,status='delete')
         end if
#ifdef CompileWithMPI
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
         call MPI_FINALIZE (this%l%ierr)
#endif
         STOP
      end if
#endif

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      call get_secnds (this%l%time_out2)
      
   
      if (present(input_flags)) then 
         this%l%read_command_line = .false.
         this%l%chain2 = input_flags
         this%l%length = len(input_flags)
      else
      ! mira el command_line y el fichero launch 251022
         call get_command (this%l%chain2, this%l%length, status)
         if (status /= 0) then
            call stoponerror (this%l%layoutnumber, this%l%num_procs, 'General error',.true.); goto 652
         end if
      end if

      this%l%chain2=trim(adjustl(this%l%chain2))
      !concatena con lo que haya en launch
      inquire(file='launch', EXIST=hayinput)
      if (hayinput) then
         open(9, file='launch', FORM='formatted',action='read')
         READ (9, '(a)') chain3
         chain3=trim(adjustl(chain3))
         CLOSE (9)               
         print *,'----> launch input file '//trim(adjustl(chain3))
      end if
#ifdef CompileWithMPI
      call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif


      this%l%chain2=trim(adjustl(this%l%chain2))//' '//trim(adjustl(chain3))

      call buscaswitchficheroinput(this%l)
      

   if (status /= 0) then
       call stoponerror (this%l%layoutnumber, this%l%num_procs, 'Error in searching input file. Correct and remove pause file',.true.); goto 652
   end if
!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!
   call print_credits(this%l)

#ifdef CompileWithMPI
   call initialize_MPI_process(this%l%filefde,this%l%extension) 
#else
#ifdef CompilePrivateVersion
   if (trim(adjustl(this%l%extension))=='.nfde') then 
#ifdef CompileWithMTLN
      NFDE_FILE => cargar_NFDE_FILE (this%l%filefde)
#else
      call WarnErrReport(".nfde files are not supported when compiling with MTLN.", .true.)
#endif
   else
      allocate (NFDE_FILE)
   end if
#else
   allocate (NFDE_FILE)
#endif
#endif

   call data_loader(this%l%filefde, parser)

   this%sgg%extraswitches=parser%switches
!!!da preferencia a los switches por linea de comando
   call getcommandargument (this%l%chain2, 1, chaindummy, this%l%length, statuse, getBinaryPath())

   this%l%chain2=trim(adjustl(this%l%chain2))
   chaindummy=trim(adjustl(chaindummy))
   this%l%length=len(trim(adjustl(chaindummy)))
   this%l%chain2=trim(adjustl(chaindummy))//' '//trim(adjustl(this%sgg%extraswitches))//' '//trim(adjustl(this%l%chain2(this%l%length+1:)))               
   this%l%chaininput=trim(adjustl(this%l%chain2))
!!!!
   call interpreta(this%l,status )      
   this%sgg%nEntradaRoot=trim (adjustl(this%l%nEntradaRoot))

#ifdef CompileWithMPI            
   call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif

   call nfde_rotate (parser,this%l%mpidir)

#ifdef CompileWithMPI            
   call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif

#ifdef CompileWithMTLN   
   if (parser%general%mtlnProblem) then 
      call solver%launch_mtln_simulation(parser%mtln, this%l%nEntradaRoot, this%l%layoutnumber) 
      STOP
   end if
#endif

#ifdef CompileWithHDF
   !!!!tunel a lo bestia para crear el .h5 a 021219
      if (this%l%createh5filefromsinglebin) then
      if (this%l%layoutnumber==0) then
         inquire(file=trim(adjustl(this%sgg%nEntradaRoot))//'_h5bin.txt',exist=lexis)
         if (.not.lexis) goto 9083
         open(newunit=myunit,file=trim(adjustl(this%sgg%nEntradaRoot))//'_h5bin.txt',form='formatted',err=9083) !lista de todos los .h5bin
         do 
            read (myunit,'(a)',end=84552) filename_h5bin
            call createh5filefromsinglebin(filename_h5bin,this%l%vtkindex) 
            print *, 'Processed '//trim(adjustl(filename_h5bin))
         end do
   84552  close(myunit)
         print *, 'END: SUCCESS creating '//trim(adjustl(this%sgg%nEntradaRoot))//'_h5bin.txt'
         stop
   9083   call stoponerror (0, this%l%num_procs, 'Invalid _h5bin.txt file',.true.); statuse=-1; !return
      end if
#ifdef CompileWithMPI
         !wait until everything comes out
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         stop
      end if
#endif

      if (status /= 0) then
         call print11(this%l%layoutnumber,'Remove running and pause files. If error persists check switches for error.  '//this%l%chain2,.true.)
         call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' ');  goto 652
      end if

      call set_priorities(this%l%prioritizeCOMPOoverPEC,this%l%prioritizeISOTROPICBODYoverall,this%l%prioritizeTHINWIRE) !!! asigna las prioridades
      if (this%l%finaltimestep /= -2) then
         ! nfde part
         call print11 (this%l%layoutnumber, 'INIT conversion internal ASCII => Binary')
         call print11 (this%l%layoutnumber, SEPARADOR//SEPARADOR//SEPARADOR)

         call print11 (this%l%layoutnumber, SEPARADOR//SEPARADOR//SEPARADOR)
         !!!!!!!!!!!!!!!!!!!!!!
         call NFDE2sgg
         this%l%fatalerror=this%l%fatalerror.or.this%l%fatalerrornfde2sgg
         !!!!!!!!!!!!!!!!!!!!!
#ifdef CompileWithMPI
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         call print11 (this%l%layoutnumber, '[OK] Ended conversion internal ASCII => Binary')
         !release memory created by newPARSER
         if (this%l%fatalerror) then
            if (allocated(this%media%sggMiEx)) deallocate(this%media%sggMiEx, this%media%sggMiEy, this%media%sggMiEz,this%media%sggMiHx, this%media%sggMiHy, this%media%sggMiHz,this%media%sggMiNo,this%media%sggMtag)
            call stoponerror (this%l%layoutnumber, this%l%num_procs, 'Error in .nfde file syntax. Check all *Warnings* and *tmpWarnings* files, correct and remove pause file if any',.true.); goto 652
         end if

         if (allocated(this%media%sggMiEx)) then !para el this%l%skindepthpre no se allocatea nada
         call AssigLossyOrPECtoNodes(this%sgg,this%media)

         if (this%l%createmap) call store_geomData (this%sgg,this%media, this%l%geomfile)
         end if
         !
#ifdef CompileWithMPI
         !wait until everything comes out
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
      end if
      write(dubuf,*) '[OK] Ended Conformal Mesh';  call print11(this%l%layoutnumber,dubuf)
      if (this%l%finaltimestep==0) this%l%finaltimestep=this%sgg%TimeSteps !no quitar
      if (this%l%forcesteps) then
         this%sgg%TimeSteps = this%l%finaltimestep
#ifdef CompileWithMTLN
         this%mtln_parsed%number_of_steps = this%l%finaltimestep 
#endif
      else
         this%l%finaltimestep = this%sgg%TimeSteps
      end if
      if (.not.this%l%forcesteps) then
         finaltimestepantesdecorregir=this%l%finaltimestep
         if (dtantesdecorregir /= 0.0) then
            this%l%finaltimestep=int(dtantesdecorregir/this%sgg%dt*finaltimestepantesdecorregir)
         end if
#ifdef CompileWithMPI
         call MPI_AllReduce( this%l%finaltimestep, NEWfinaltimestep, 1_4, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, this%l%ierr)
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
         this%l%finaltimestep=NEWfinaltimestep
#endif
#ifdef CompileWithMTLN
            this%mtln_parsed%number_of_steps = this%l%finaltimestep 
#endif
            if (finaltimestepantesdecorregir/=this%l%finaltimestep) then
               write(dubuf,*) SEPARADOR//separador//separador
               call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) 'Original Final Time Step= ',finaltimestepantesdecorregir
               if (this%l%layoutnumber==0) call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) 'Corrected Final Time Step= ',this%l%finaltimestep
               if (this%l%layoutnumber==0) call print11(this%l%layoutnumber,dubuf)
            end if
      end if
      !check that simulation can actually be done for the kind of media requested
      do i = 1, this%sgg%nummedia
         if (this%sgg%Med(i)%Is%ThinWire) then
#ifndef CompileWithBerengerWires
      if  ((this%l%wiresflavor=='berenger')) then
            call stoponerror (this%l%layoutnumber, this%l%num_procs, 'Berenger Wires without support. Recompile!')
      end if
#endif
#ifndef CompileWithSlantedWires
      if  ((this%l%wiresflavor=='slanted').or.(this%l%wiresflavor=='semistructured')) then
            call stoponerror (this%l%layoutnumber, this%l%num_procs, 'slanted Wires without support. Recompile!')
      end if
#endif
            continue
         end if
         !
         if ((this%sgg%Med(i)%Is%AnisMultiport) .OR. (this%sgg%Med(i)%Is%multiport).OR. (this%sgg%Med(i)%Is%SGBC)) then
#ifndef CompileWithNIBC
            if (this%l%mibc) call stoponerror (this%l%layoutnumber, this%l%num_procs, 'this%l%mibc Multiports without support. Recompile!')
#endif
            continue
         end if
   !altair no conformal sgbc 201119
#ifdef NoConformalSGBC
         if (this%sgg%Med(i)%Is%sgbc .and. this%l%input_conformal_flag) then
            call stoponerror (this%l%layoutnumber, this%l%num_procs, 'Conformal sgbc not allowed. ')
         end if
#endif
   !    
      end do
      
      
      if (this%l%thereare_stoch.and.(.not.this%l%chosenyesornostochastic)) then
         call stoponerror (this%l%layoutnumber, this%l%num_procs, '!STOCH found in .nfde. Specify either -stoch or -nostoch')
      end if
#ifndef CompileWithSlantedWires
      if (this%l%hay_slanted_wires) then
         call stoponerror (this%l%layoutnumber, this%l%num_procs, 'slanted wires without slanted support. Recompile ()')
      end if
#endif   
      if (this%l%hay_slanted_wires .AND. ((trim(adjustl(this%l%wiresflavor))/='slanted').AND.(trim(adjustl(this%l%wiresflavor))/='semistructured'))) then
         call stoponerror (this%l%layoutnumber, this%l%num_procs, 'slanted wires require -this%l%wiresflavor Slanted/semistructured')
      end if

      
      !Error abrezanjas y no this%l%resume conformal
      ThereArethinslots=.FALSE.
      do jmed=1,this%sgg%NumMedia
         if (this%sgg%Med(jmed)%Is%ThinSlot) ThereArethinslots=.true.
      end do
      if (this%l%resume.and.this%l%run_with_abrezanjas.and.ThereArethinslots) then   
            call stoponerror (this%l%layoutnumber, this%l%num_procs, 'this%l%resume -r currently unsupported by conformal solver',.true.); statuse=-1; !return
      end if
      !
   !!!SOME FINAL REPORTING

      if (this%l%layoutnumber==0) then
         write(dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         call print11 (this%l%layoutnumber, dubuf)
         call print11 (this%l%layoutnumber, 'Solver launched with options:')
         write(dubuf,*) this%l%mibc          
         call print11 (this%l%layoutnumber, '---> this%l%mibc    solver for NIBC multilayer: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%ade         
         call print11 (this%l%layoutnumber, '---> this%l%ade     solver for ADC multilayer: '//trim(adjustl(dubuf)))
         Write(dubuf,*) this%l%sgbc    
         call print11 (this%l%layoutnumber, '---> sgbc    solver for multilayer: '//trim(adjustl(dubuf)))
         if (this%l%sgbc) then
               write(dubuf,*) this%l%sgbcDispersive      
               call print11 (this%l%layoutnumber, '---> sgbc DISPERSIVE solver for multilayer: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%sgbccrank     
               call print11 (this%l%layoutnumber, '---> sgbc Crank-Nicolson solver for multilayer: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%sgbcdepth
               call print11 (this%l%layoutnumber, '---> sgbc Depth: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%sgbcfreq
               call print11 (this%l%layoutnumber, '---> sgbc Freq: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%sgbcresol
               call print11 (this%l%layoutnumber, '---> sgbc Resol: '//trim(adjustl(dubuf)))
         end if
         write(dubuf,*) this%l%skindepthpre
         call print11 (this%l%layoutnumber, '---> this%l%skindepthpre preprocessing for multilayer: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%flag_conf_sgg
         call print11 (this%l%layoutnumber, '---> Conformal file external: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%input_conformal_flag      
         call print11 (this%l%layoutnumber, '---> Conformal solver: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%run_with_abrezanjas
         call print11 (this%l%layoutnumber, '---> Conformal thin-gap solver: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%run_with_dmma
         call print11 (this%l%layoutnumber, '---> DMMA thin-gap solver: '//trim(adjustl(dubuf)))
#ifdef CompileWithMTLN
         write(dubuf,'(a)') 'MTLN wires'
         call print11 (this%l%layoutnumber, '---> Wire model: '//trim(adjustl(dubuf)))
#else
         write(dubuf,'(a)') this%l%wiresflavor
         call print11 (this%l%layoutnumber, '---> Wire model: '//trim(adjustl(dubuf)))
         write(dubuf,'(a)') this%l%inductance_model
         call print11 (this%l%layoutnumber, '---> Inductance model: '//trim(adjustl(dubuf)))
         if (trim(adjustl(this%l%wiresflavor))=='berenger') then
               write(dubuf,*) this%l%mindistwires
               call print11 (this%l%layoutnumber, '---> Berenger minimum distance between wires: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%mtlnberenger
               call print11 (this%l%layoutnumber, '---> Berenger -this%l%mtlnberenger MTLN switch: '//trim(adjustl(dubuf)))
         end if
         if (trim(adjustl(this%l%wiresflavor))=='holland') then
               write(dubuf,*) this%l%stableradholland                 
               call print11 (this%l%layoutnumber, '---> Holland -this%l%stableradholland automatic correction switch: '//trim(adjustl(dubuf)))
         end if
         write(dubuf,*) this%l%TAPARRABOS                
         call print11 (this%l%layoutnumber, '---> Thin-wire double-tails removed: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%fieldtotl                
         call print11 (this%l%layoutnumber, '---> Thin-wire -this%l%fieldtotl experimental switch: '//trim(adjustl(dubuf)))

         write(dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         call print11 (this%l%layoutnumber, dubuf)
#endif
      end if
      if (this%l%layoutnumber == 0) then
         call erasesignalingfiles(this%l%simu_devia)
      end if
      
      if (this%l%layoutnumber==0) then
         
         open(newunit=thefileno,FILE = trim(adjustl(this%l%nEntradaRoot))//'_tag_paraviewfilters.txt')
               write(thefileno,'(a)') trim(adjustl('### FOR SLICE CURRENT VTK PROBES select the "current_t" or "current_f"                           '))   
               write(thefileno,'(a)') trim(adjustl('### FOR MAP VTK PROBES select the "mediatype" layer                                               '))             
               write(thefileno,'(a)') trim(adjustl('### For Paraview versions over 5.10 just use the Threshold exisiting filter to select the interval'))           
               write(thefileno,'(a)') trim(adjustl('### ######################'))
               write(thefileno,'(a)') trim(adjustl('### For Paraview versions under 5.10 Copy and paste the next as a programmable filter to select only one interval of tags'))
               write(thefileno,'(a)') trim(adjustl('import vtk                                                                                        '))
               write(thefileno,'(a)') trim(adjustl('inp = self.GetInputDataObject(0, 0)                                                               '))
               write(thefileno,'(a)') trim(adjustl('outp = self.GetOutputDataObject(0)                                                                '))
               write(thefileno,'(a)') trim(adjustl('thresh = vtk.vtkThreshold()                                                                       '))
               write(thefileno,'(a)') trim(adjustl('thresh.SetInputData(inp)                                                                          '))
               write(thefileno,'(a)') trim(adjustl('thresh.SetInputArrayToProcess(0, 0, 0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "tagnumber")     '))
               write(thefileno,'(a)') trim(adjustl('thresh.ThresholdBetween(64,127)                                                              '))
               write(thefileno,'(a)') trim(adjustl('thresh.Update()                                                              '))
               write(thefileno,'(a)') trim(adjustl('outp.ShallowCopy(thresh.GetOutput())    '))        
               write(thefileno,'(a)') trim(adjustl( '# Replace the thresh.ThresholdBetween numbers by tag intervals below to filter by tags           '))
               write(thefileno,'(a)')               '# ( -1e21    , -1e-3    ) '//trim(adjustl('Candidates for undesired free-space slots'))
               write(thefileno,'(a,i9,a,i9,a)')     '# (  0       ,  63      ) '//trim(adjustl('Nodal sources, etc.'))
               do i=1,this%tagtype%numertags
                  write(thefileno,'(a,i9,a,i9,a)') '# (',i*64,' , ',i*64+63,') '//trim(adjustl(this%tagtype%tag(i))) !los shifteo 6 bits y les sumo 2**campo ! idea de los 3 bits de 151020
               end do
               !!
               write(thefileno,'(a)') trim(adjustl( '###    '))   
               write(thefileno,'(a)') trim(adjustl( '###    '))   
               write(thefileno,'(a)') trim(adjustl( '### FOR MAP VTK PROBES select the "mediatype" layer                                               '))                
               write(thefileno,'(a)') trim(adjustl( '### For Paraview versions over 5.10 just use the Threshold exisiting filter to select the interval'))           
               write(thefileno,'(a)') trim(adjustl( '### ######################'))
               write(thefileno,'(a)') trim(adjustl( '### For Paraview versions under 5.10Copy and paste the next as a programmable filter to select only one types of media'))
               write(thefileno,'(a)') trim(adjustl( 'import vtk                                                                                        '))
               write(thefileno,'(a)') trim(adjustl( 'inp = self.GetInputDataObject(0, 0)                                                               '))
               write(thefileno,'(a)') trim(adjustl( 'outp = self.GetOutputDataObject(0)                                                                '))
               write(thefileno,'(a)') trim(adjustl( 'thresh = vtk.vtkThreshold()                                                                       '))
               write(thefileno,'(a)') trim(adjustl( 'thresh.SetInputData(inp)                                                                          '))
               write(thefileno,'(a)') trim(adjustl( 'thresh.SetInputArrayToProcess(0, 0, 0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, "mediatype")     '))
               write(thefileno,'(a)') trim(adjustl( 'thresh.ThresholdBetween(0.0,0.5)                                                              '))
               write(thefileno,'(a)') trim(adjustl( 'thresh.Update()                                                              '))
               write(thefileno,'(a)') trim(adjustl( 'outp.ShallowCopy(thresh.GetOutput())  '))
               write(thefileno,'(a)') trim(adjustl( '# Replace the thresh.ThresholdBetween numbers by media types below to filter by media types           '))
               write(thefileno,'(a)') '# ( -100 , -100 ) '//trim(adjustl('Candidates for undesired free-space slots                               (Surface)'))
               write(thefileno,'(a)') '# (  0.0 ,  0.0 ) '//trim(adjustl('PEC                                                                     (Surface)'))
               write(thefileno,'(a)') '# (  0.5 ,  0.5 ) '//trim(adjustl('PEC                                                                     (Line)'))
               write(thefileno,'(a)') '# (  1.5 ,  1.5 ) '//trim(adjustl('Dispersive electric or magnetic isotropic or anisotropic                (Line)'))
               write(thefileno,'(a)') '# (  100 ,  199 ) '//trim(adjustl('Dispersive electric/magnetic isotropic/anisotropic (+indexmedium)       (Surface) '))
               write(thefileno,'(a)') '# (  2.5 ,  2.5 ) '//trim(adjustl('Dielectric isotropic or anisotropic                                     (Line)'))
               write(thefileno,'(a)') '# (  200 ,  299 ) '//trim(adjustl('Dielectric isotropic or anisotropic (+indexmedium)                      (Surface)'))
               write(thefileno,'(a)') '# (  3.5 ,  3.5 ) '//trim(adjustl('sgbc/this%l%mibc Isotropic/anisotropic Multiport                        (Line)'))
               write(thefileno,'(a)') '# (  300 ,  399 ) '//trim(adjustl('sgbc/this%l%mibc Isotropic/anisotropic Multiport (+indexmedium)         (Surface)'))
               write(thefileno,'(a)') '# (  4.5 ,  4.5 ) '//trim(adjustl('Thin slot                                                               (Line)'))
               write(thefileno,'(a)') '# (  5.0 ,  5.0 ) '//trim(adjustl('Already_YEEadvanced_byconformal                                         (Surface)'))
               write(thefileno,'(a)') '# (  5.5 ,  5.5 ) '//trim(adjustl('Already_YEEadvanced_byconformal                                         (Line)'))
               write(thefileno,'(a)') '# (  6.0 ,  6.0 ) '//trim(adjustl('Split_and_useless                                                       (Surface)'))
               write(thefileno,'(a)') '# (  6.5 ,  6.5 ) '//trim(adjustl('Split_and_useless                                                       (Line)'))
               write(thefileno,'(a)') '# (  7.0 ,  7.0 ) '//trim(adjustl('Edge Not colliding thin wires                                           (Line)'))
               write(thefileno,'(a)') '# (  8.0 ,  8.0 ) '//trim(adjustl('Thin wire segments colliding with structure                             (Line)'))
               write(thefileno,'(a)') '# (  8.5 ,  8.5 ) '//trim(adjustl('Soft/Hard Nodal CURRENT/FIELD ELECTRIC DENSITY SOURCE                   (Line)'))
               write(thefileno,'(a)') '# (  9.0 ,  9.0 ) '//trim(adjustl('Soft/Hard Nodal CURRENT/FIELD MAGNETIC DENSITY SOURCE                   (Line)'))
               write(thefileno,'(a)') '# (   10 ,   11 ) '//trim(adjustl('LeftEnd/RightEnd/Ending wire segment                                    (Wire)'))
               write(thefileno,'(a)') '# (   20 ,   20 ) '//trim(adjustl('Intermediate wire segment +number_holland_parallel or +number_berenger  (Wire) '))
               write(thefileno,'(a)') '# (   12 ,   12 ) '//trim(adjustl('Edge Not colliding multiwires                                           (Multiwire)'))
               write(thefileno,'(a)') '# (   13 ,   13 ) '//trim(adjustl('Multiwire segments colliding with structure                             (Multiwire)'))
               write(thefileno,'(a)') '# (   14 ,   15 ) '//trim(adjustl('LeftEnd/RightEnd/Ending multiwire segment                               (Multiwire)'))
               write(thefileno,'(a)') '# (   60 ,   60 ) '//trim(adjustl('Intermediate multiwire segment + number parallel segments               (Multiwire) '))
               write(thefileno,'(a)') '# (  400 ,  499 ) '//trim(adjustl('Thin slot (+indexmedium)                                                (Surface)'))
               write(thefileno,'(a)') '# ( 1000 , 1999 ) '//trim(adjustl('Conformal Volume PEC (+indexmedium)                                     (Surface)'))
               write(thefileno,'(a)') '# ( 2000 , 2999 ) '//trim(adjustl('Conformal Volume PEC (+indexmedium)                                     (Line)'))
               write(thefileno,'(a)') '# ( -0.5 , -0.5 ) '//trim(adjustl('Other types of media                                                    (Line)'))
               write(thefileno,'(a)') '# ( -1.0 , -1.0 ) '//trim(adjustl('Other types of media                                                    (Surface)'))
         close(thefileno)
      end if

contains 
   subroutine NFDE2sgg     
   !!!!!!!!!      
         real(kind=rkind) :: dt,finaldt
         logical fatalerror
         ! parser now holds all the .nfde info
         !first read the limits
#ifdef CompileWithMPI
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         call read_limits_nogeom (this%l%layoutnumber,this%l%num_procs, this%sgg, this%fullsize, this%SINPML_fullsize, parser,this%l%MurAfterPML,this%l%mur_exist)
      
         dtantesdecorregir=this%sgg%dt

         dxmin=minval(this%sgg%DX)
         dymin=minval(this%sgg%DY)
         dzmin=minval(this%sgg%DZ)
         !!!
         dtlay=(1.0_RKIND/(this%cluz*sqrt(((1.0_RKIND / dxmin)**2.0_RKIND )+((1.0_RKIND / dymin)**2.0_RKIND )+((1.0_RKIND / dzmin)**2.0_RKIND ))))
         dt=dtlay
#ifdef CompileWithMPI
         call MPIupdateMin(dtlay,dt)
#endif

         if (this%l%forcecfl) then
            this%sgg%dt=dt*this%l%cfl
            write(dubuf,*) SEPARADOR//separador//separador
            call print11(this%l%layoutnumber,dubuf)
            write(dubuf,*) 'Correcting sgg%dt with -this%l%cfl switch. New time step: ',this%sgg%dt
            call print11(this%l%layoutnumber,dubuf)
            write(dubuf,*) SEPARADOR//separador//separador
            call print11(this%l%layoutnumber,dubuf)
         else
            if (dtantesdecorregir == 0.0 .or. this%sgg%dt > dt*heurCFL) then
               write(dubuf,*) SEPARADOR//separador//separador
               call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) 'Automatically correcting dt for stability reasons: '
               call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) 'Original dt: ', this%sgg%dt
               call print11(this%l%layoutnumber,dubuf)
               this%sgg%dt=dt*heurCFL
               write(dubuf,*) 'New dt: ', this%sgg%dt
               call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) SEPARADOR//separador//separador
               call print11(this%l%layoutnumber,dubuf)
            end if
         end if
         !!!!!!!!!!!!No es preciso re-sincronizar pero lo hago !!!!!!!!!!!!!!!!!!!!!!!!!!
         finaldt=this%sgg%dt
#ifdef CompileWithMPI
         call MPIupdateMin(real(this%sgg%dt,RKIND),finaldt)
#endif
         !!!!!!!!!!!!!!
         this%l%cfl=this%sgg%dt/dtlay
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(this%l%layoutnumber,dubuf)
         write(dubuf,*) 'CFLN= ',this%l%cfl
         call print11(this%l%layoutnumber,dubuf)
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(this%l%layoutnumber,dubuf)

         write(dubuf,*) SEPARADOR//separador//separador
         call print11(this%l%layoutnumber,dubuf)
         write(dubuf,*) 'Deltat= ',this%sgg%dt
         if (this%l%layoutnumber==0) call print11(this%l%layoutnumber,dubuf)
#ifdef CompileWithMPI
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(this%l%layoutnumber,dubuf)
         if (this%l%mur_exist.and.this%l%mur_first) then
            this%l%mur_second=.false.
         else
            this%l%mur_second=.false. !arreglar cuando se arregle el bug de las mur second
            this%l%mur_first=.true. !arreglar cuando se arregle el bug de las mur second
         end if
#ifdef CompileWithMPI
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         !LATER OVERRRIDEN BY MPI
         !ALLOCATED ONE MORE TO KEEP PMC INFO FOR THE HX,HY,HZ FIELDS
         this%sgg%Alloc(1:6)%XI = this%fullsize(1:6)%XI - 1
         this%sgg%Alloc(1:6)%XE = this%fullsize(1:6)%XE + 1
         this%sgg%Alloc(1:6)%YI = this%fullsize(1:6)%YI - 1
         this%sgg%Alloc(1:6)%YE = this%fullsize(1:6)%YE + 1
         !REDUCE THE SWEEP AREA BY 1
         this%sgg%Sweep(1:6)%XI = this%fullsize(1:6)%XI
         this%sgg%Sweep(1:6)%XE = this%fullsize(1:6)%XE
         this%sgg%Sweep(1:6)%YI = this%fullsize(1:6)%YI
         this%sgg%Sweep(1:6)%YE = this%fullsize(1:6)%YE
         !
         if (this%l%num_procs == 1) then
            this%sgg%Alloc(1:6)%ZI = this%fullsize(1:6)%ZI - 1
            this%sgg%Alloc(1:6)%ZE = this%fullsize(1:6)%ZE + 1
            !REDUCE THE SWEEP AREA BY 1
            this%sgg%Sweep(1:6)%ZI = this%fullsize(1:6)%ZI
            this%sgg%Sweep(1:6)%ZE = this%fullsize(1:6)%ZE
            !!incluido aqui pq se precisa para clip 16/07/15
            do field = iEx, iHz
               this%sgg%SINPMLSweep(field)%XI = Max (this%SINPML_fullsize(field)%XI, this%sgg%Sweep(field)%XI)
               this%sgg%SINPMLSweep(field)%XE = Min (this%SINPML_fullsize(field)%XE, this%sgg%Sweep(field)%XE)
               this%sgg%SINPMLSweep(field)%YI = Max (this%SINPML_fullsize(field)%YI, this%sgg%Sweep(field)%YI)
               this%sgg%SINPMLSweep(field)%YE = Min (this%SINPML_fullsize(field)%YE, this%sgg%Sweep(field)%YE)
               this%sgg%SINPMLSweep(field)%ZI = Max (this%SINPML_fullsize(field)%ZI, this%sgg%Sweep(field)%ZI)
               this%sgg%SINPMLSweep(field)%ZE = Min (this%SINPML_fullsize(field)%ZE, this%sgg%Sweep(field)%ZE)
            end do
            !!fin 16/07/15
            write(dubuf,*) 'INIT NFDE --------> GEOM'
            call print11 (this%l%layoutnumber, dubuf)
            call read_geomData (this%sgg,this%media,this%tag_numbers, this%l%fichin, this%l%layoutnumber, this%l%num_procs, this%SINPML_fullsize, this%fullsize, parser, &
            this%l%groundwires,this%l%attfactorc,this%l%mibc,this%l%sgbc,this%l%sgbcDispersive,this%l%MEDIOEXTRA,this%maxSourceValue,this%l%skindepthpre,this%l%createmapvtk,this%l%input_conformal_flag,this%l%CLIPREGION,this%l%boundwireradius,this%l%maxwireradius,this%l%updateshared,this%l%run_with_dmma, this%eps0, &
            this%mu0,.false.,this%l%hay_slanted_wires,this%l%verbose,this%l%ignoresamplingerrors,this%tagtype,this%l%wiresflavor)            
            ! call read_geomData (this%sgg,this%sggMtag,this%tag_numbers, this%sggMiNo,this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz, this%l%fichin, this%l%layoutnumber, this%l%num_procs, this%SINPML_fullsize, this%fullsize, parser, &
            ! this%l%groundwires,this%l%attfactorc,this%l%mibc,this%l%sgbc,this%l%sgbcDispersive,this%l%MEDIOEXTRA,this%maxSourceValue,this%l%skindepthpre,this%l%createmapvtk,this%l%input_conformal_flag,this%l%CLIPREGION,this%l%boundwireradius,this%l%maxwireradius,this%l%updateshared,this%l%run_with_dmma, this%eps0, &
            ! this%mu0,.false.,this%l%hay_slanted_wires,this%l%verbose,this%l%ignoresamplingerrors,this%tagtype,this%l%wiresflavor)
#ifdef CompileWithMTLN
            if (trim(adjustl(this%l%extension))=='.json')  then 
               this%mtln_parsed = parser%mtln
               this%mtln_parsed%time_step = this%sgg%dt
            end if
            ! if (trim(adjustl(this%l%extension))=='.json')  mtln_solver = mtlnCtor(parser%mtln)   
#endif
            write(dubuf,*) '[OK] ENDED NFDE --------> GEOM'
            call print11 (this%l%layoutnumber, dubuf)
            !writing
            slices = '!SLICES'
            write(buff, '(i7)') this%sgg%Sweep(iHz)%ZE - this%sgg%Sweep(iHz)%ZI
            slices = trim (adjustl(slices)) // '_' // trim (adjustl(buff))
            if (this%l%resume .AND. (slices /= this%l%slicesoriginales)) then
               buff='Different resumed/original MPI slices: '//trim(adjustl(slices))//' '//&
               & trim(adjustl(this%l%slicesoriginales))
               call stoponerror (this%l%layoutnumber, this%l%num_procs, buff)
            end if
            call print11 (this%l%layoutnumber, trim(adjustl(slices)))
            !end writing
            write(buff, '(a,i7,a,i7)') '_________Spanning from z=', this%sgg%Sweep(iHz)%ZI, ' to z=', this%sgg%Sweep(iHz)%ZE
            call print11 (this%l%layoutnumber, trim(adjustl(buff)))
#ifdef CompileWithMPI
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#ifdef CompileWithStochastic
            if (this%l%stochastic) then
               buff='this%l%stochastic uncompatible with MPI this%l%num_procs smaller than 2'
               call stoponerror (this%l%layoutnumber, this%l%num_procs, buff)
            end if
#endif
#endif
         ELSE !del this%l%num_procs==1       
#ifdef CompileWithMPI
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#ifdef CompileWithStochastic
            if (this%l%stochastic) then
               call HalvesStochasticMPI(this%l%layoutnumber,this%l%num_procs,this%l%simu_devia)
            end if
#endif
                     
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)   
   !!!ahora divide el espacio computacional
            call MPIdivide (this%sgg, this%fullsize, this%SINPML_fullsize, this%l%layoutnumber, this%l%num_procs, this%l%forcing, this%l%forced, this%l%slicesoriginales, this%l%resume,this%l%fatalerror)
            !
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)   
            if (this%l%fatalerror) then
   !intenta recuperarte
               return
            end if
      
            ! if the layout is pure PML then take at least a line of non PML to build the PML data insider read_geomDAta
            ! Uses extra memory but later matrix sggm is deallocated in favor of smaller sggMIEX, etc
            do field = iEx, iHz
               tempalloc(field)%ZE = this%sgg%Alloc(field)%ZE
               tempalloc(field)%ZI = this%sgg%Alloc(field)%ZI
               this%sgg%Alloc(field)%ZE = Max (this%sgg%Alloc(field)%ZE, this%SINPML_fullsize(field)%ZI+1)
               this%sgg%Alloc(field)%ZI = Min (this%sgg%Alloc(field)%ZI, this%SINPML_fullsize(field)%ZE-1)
            end do
            !   
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)  
            !!incluido aqui pq se precisa para clip 16/07/15
            do field = iEx, iHz
               this%sgg%SINPMLSweep(field)%XI = Max (this%SINPML_fullsize(field)%XI, this%sgg%Sweep(field)%XI)
               this%sgg%SINPMLSweep(field)%XE = Min (this%SINPML_fullsize(field)%XE, this%sgg%Sweep(field)%XE)
               this%sgg%SINPMLSweep(field)%YI = Max (this%SINPML_fullsize(field)%YI, this%sgg%Sweep(field)%YI)
               this%sgg%SINPMLSweep(field)%YE = Min (this%SINPML_fullsize(field)%YE, this%sgg%Sweep(field)%YE)
               this%sgg%SINPMLSweep(field)%ZI = Max (this%SINPML_fullsize(field)%ZI, this%sgg%Sweep(field)%ZI)
               this%sgg%SINPMLSweep(field)%ZE = Min (this%SINPML_fullsize(field)%ZE, this%sgg%Sweep(field)%ZE)
            end do
            !!fin 16/07/15
            write(dubuf,*) 'INIT NFDE --------> GEOM'
            call print11 (this%l%layoutnumber, dubuf)           

            call read_geomData (this%sgg,this%media,this%tag_numbers, this%l%fichin, this%l%layoutnumber, this%l%num_procs, this%SINPML_fullsize, this%fullsize, parser, &
            this%l%groundwires,this%l%attfactorc,this%l%mibc,this%l%sgbc,this%l%sgbcDispersive,this%l%MEDIOEXTRA,this%maxSourceValue,this%l%skindepthpre,this%l%createmapvtk,this%l%input_conformal_flag,this%l%CLIPREGION,this%l%boundwireradius,this%l%maxwireradius,this%l%updateshared,this%l%run_with_dmma, &
            this%eps0,this%mu0,this%l%simu_devia,this%l%hay_slanted_wires,this%l%verbose,this%l%ignoresamplingerrors,this%tagtype,this%l%wiresflavor)
            ! call read_geomData (this%sgg,this%sggMtag,this%tag_numbers, this%sggMiNo,this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz, this%l%fichin, this%l%layoutnumber, this%l%num_procs, this%SINPML_fullsize, this%fullsize, parser, &
            ! this%l%groundwires,this%l%attfactorc,this%l%mibc,this%l%sgbc,this%l%sgbcDispersive,this%l%MEDIOEXTRA,this%maxSourceValue,this%l%skindepthpre,this%l%createmapvtk,this%l%input_conformal_flag,this%l%CLIPREGION,this%l%boundwireradius,this%l%maxwireradius,this%l%updateshared,this%l%run_with_dmma, &
            ! this%eps0,this%mu0,this%l%simu_devia,this%l%hay_slanted_wires,this%l%verbose,this%l%ignoresamplingerrors,this%tagtype,this%l%wiresflavor)


#ifdef CompileWithMPI
            !wait until everything comes out
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
#ifdef CompileWithMTLN
            if (trim(adjustl(this%l%extension))=='.json')  then 
               this%mtln_parsed = parser%mtln
               this%mtln_parsed%time_step = this%sgg%dt
            end if
#endif
            write(dubuf,*) '[OK] ENDED NFDE --------> GEOM'
            call print11 (this%l%layoutnumber, dubuf)
            !restore back the indexes
            do field = iEx, iHz
               this%sgg%Alloc(field)%ZE = tempalloc(field)%ZE
               this%sgg%Alloc(field)%ZI = tempalloc(field)%ZI
            end do
#endif
            continue
         end if !del this%l%num_procs==1
         !
#ifdef CompileWithMPI
         !wait until everything comes out
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         !!!!!!!!!!!!!lo dejo aqui debajo tambien aunque ya se ha calculado antes para lo del clipping
         do field = iEx, iHz
            this%sgg%SINPMLSweep(field)%XI = Max (this%SINPML_fullsize(field)%XI, this%sgg%Sweep(field)%XI)
            this%sgg%SINPMLSweep(field)%XE = Min (this%SINPML_fullsize(field)%XE, this%sgg%Sweep(field)%XE)
            this%sgg%SINPMLSweep(field)%YI = Max (this%SINPML_fullsize(field)%YI, this%sgg%Sweep(field)%YI)
            this%sgg%SINPMLSweep(field)%YE = Min (this%SINPML_fullsize(field)%YE, this%sgg%Sweep(field)%YE)
            this%sgg%SINPMLSweep(field)%ZI = Max (this%SINPML_fullsize(field)%ZI, this%sgg%Sweep(field)%ZI)
            this%sgg%SINPMLSweep(field)%ZE = Min (this%SINPML_fullsize(field)%ZE, this%sgg%Sweep(field)%ZE)
         end do
         return
      end subroutine

#ifdef CompileWithMPI
   subroutine initialize_MPI_process(filename, extension)
      character(len=BUFSIZE), intent(in) :: filename, extension
      integer(kind=4) :: mpi_t_linea_t,longitud4
      integer(kind=8) :: rawInfoBuffer, numeroLineasFichero, i8, longitud8
      type(t_NFDE_FILE_t), pointer :: rawFileInfo

      write (dubuf,*) 'INIT Reading file '//trim (adjustl(this%whoami))//' ', trim (adjustl(filename))

      call print11 (this%l%layoutnumber, dubuf)

      if (this%l%layoutnumber==0) then
#ifdef CompilePrivateVersion
         if (trim(adjustl(extension))=='.nfde') then 
#ifdef CompileWithMTLN
            call stoponerror(this%l%layoutnumber, this%l%num_procs, &
               'NFDE files are not supported when compiling with MTLN', .true.)
#endif
            NFDE_FILE => cargar_NFDE_FILE (filename)
         else
            call carga_raw_info(rawFileInfo, filename, extension)
            NFDE_FILE => rawFileInfo
         end if
#else
         call carga_raw_info(rawFileInfo, filename, extension)
         NFDE_FILE => rawFileInfo
#endif
      else
        allocate(NFDE_FILE)
      end if

      write(dubuf,*) '[OK]';  call print11(this%l%layoutnumber,dubuf)

      write(dubuf,*) 'INIT Sharing file through MPI'; call print11 (this%l%layoutnumber, dubuf)
      !
      call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
      !
      numeroLineasFichero=NFDE_FILE%numero
      call MPI_BCAST(numeroLineasFichero, 1_4, MPI_INTEGER8, 0_4, SUBCOMM_MPI, this%l%ierr)      
      if (this%l%layoutnumber/=0) then
         NFDE_FILE%targ = 1
         NFDE_FILE%numero=numeroLineasFichero
        allocate(NFDE_FILE%lineas(NFDE_FILE%numero))
      end if
      call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)

      call build_derived_t_linea(mpi_t_linea_t)

      rawInfoBuffer=ceiling(maxmpibytes*1.0_8/(BUFSIZE*1.0_8+8.0_8),8)

      do i8=1, numeroLineasFichero, rawInfoBuffer
                  longitud8=min(rawInfoBuffer, numeroLineasFichero - i8 + 1)
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
            if ((longitud8>huge(1_4)).or.(longitud8>maxmpibytes)) then
               print *,'Stop. Buggy error: MPI longitud greater that greatest integer*4'
               stop
            else
               longitud4=int(longitud8,4)
            end if
            call MPI_BCAST(NFDE_FILE%lineas(i8),longitud4,mpi_t_linea_t,0_4,SUBCOMM_MPI,this%l%ierr)    
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
      end do
   end subroutine initialize_MPI_process

#endif
   subroutine data_loader(filename, parsedProblem)
      character(len=1024), intent(in) :: filename
      type(Parseador_t), pointer :: parsedProblem
      type(fdtdjson_parser_t) :: parsed_t

      write (dubuf,*) 'INIT interpreting geometrical data from ', trim (adjustl(filename))
      call print11 (this%l%layoutnumber, dubuf)

   
      if (trim(adjustl(this%l%extension))=='.nfde') then 
#ifdef CompilePrivateVersion   
         parsedProblem => newparser (NFDE_FILE)
         ! this%l%mpidir = NFDE_FILE%mpidir
         this%l%thereare_stoch=NFDE_FILE%thereare_stoch
#else
         print *,'Not compiled with cargaNFDEINDEX'
         stop
#endif
      
#ifdef CompileWithSMBJSON
      elseif (trim(adjustl(this%l%extension))=='.json') then
         parsed_t = fdtdjson_parser_t(filename)   
         allocate(parsedProblem)
         parsedProblem = parsed_t%readProblemDescription()
#endif

      else
         print *, 'Neither .nfde nor .json files used as input after -i'
         stop
      end if

      write(dubuf,*) '[OK] '//trim(adjustl(this%whoami))//' Parser still working ';  call print11(this%l%layoutnumber,dubuf)       
#ifdef CompileWithMPI            
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
      return
   end subroutine data_loader

   function countLinesInJSONOneLiner(filename, unit) result(res)
      character(len=*), intent(in) :: filename
      integer(kind=4), intent(in) :: unit
      integer(kind=4) :: res
      character(len=BUFSIZE) :: l_aux
      integer :: size_read, pos, d, io
      res = 0
      open(UNIT=unit, FILE=trim(adjustl(filename)), STATUS='old',form='formatted')
      DO
         READ (unit, '(A)', advance='no', iostat = io, size = size_read) l_aux
         if (size_read == 0) exit
         res = res + 1
      end do
      CLOSE (unit)

   end function

   subroutine readLines(rInfo, filename, unit)
      type(t_NFDE_FILE_t), pointer :: rInfo
      character(len=*), intent(in) :: filename
      integer(kind=4), intent(in) :: unit

      type(t_linea_t), pointer :: linea
      character(len=BUFSIZE) :: l_aux
      character(len=BUFSIZE) :: buffer

     allocate(rInfo%lineas(rInfo%numero))
      rInfo%numero = 0
      open(UNIT=unit, FILE=trim(adjustl(filename)), STATUS='old',form='formatted')
      DO
         READ (unit, '(A)', end=2010) l_aux
         if (len_trim (adjustl(l_aux))>=BUFSIZE) then
            write(buffer,*) 'Line in .nfde larger than ',BUFSIZE,'Recompile '
            call warnerrreport(buffer,.TRUE.) !ABORTA
         end if
         rInfo%numero = rInfo%numero + 1
         linea => rInfo%lineas (rInfo%numero)
         linea%dato = adjustl(l_aux)
         linea%LEN=len_trim (linea%dato)
      end do
   2010   CLOSE (unit)

   end subroutine

   subroutine readLinesFromJSONOneLiner(rInfo, filename, unit)
      type(t_NFDE_FILE_t), pointer :: rInfo
      character(len=*), intent(in) :: filename
      integer(kind=4), intent(in) :: unit

      integer(kind=4) :: io, size_read, pos, d
      type(t_linea_t), pointer :: linea
      character(len=BUFSIZE) :: l_aux
      character(len=BUFSIZE) :: buffer

     allocate(rInfo%lineas(rInfo%numero))
      rInfo%numero = 0
      open(UNIT=unit, FILE=trim(adjustl(filename)), STATUS='old',form='formatted')
      DO
         READ (unit, '(A)', advance='no', iostat = io, size = size_read) l_aux
         if (size_read == 0) exit
         rInfo%numero = rInfo%numero + 1
         linea => rInfo%lineas (rInfo%numero)
         linea%dato = adjustl(l_aux)
         linea%LEN=len_trim (linea%dato)
      end do
      CLOSE (unit)

   end subroutine

   subroutine carga_raw_info (rawFileInfo, filename, extension)
      character(len=*), intent(in) :: filename, extension
      type(t_NFDE_FILE_t), pointer :: rawFileInfo
      
      type(t_linea_t), pointer :: linea
      LOGICAL :: ok
      character(len=BUFSIZE) :: l_aux
      character(len=BUFSIZE) :: buffer
      integer(kind=4) :: i,tamanio,i0,ascii,offset,ascii_menos1,j,k
      Character (Len=:), Allocatable :: fichero
      integer(kind=4), parameter :: UNIT_EF = 10

      integer(kind=4) :: prelines = 0, io
     allocate(rawFileInfo)
      rawFileInfo%numero = 0
      rawFileInfo%targ = 1

      !precount
      open(UNIT=UNIT_EF, FILE=trim(adjustl(filename)), STATUS='old',form='formatted')
      DO
         READ (UNIT_EF, '(A)', iostat=io) l_aux
         if (io/=0) exit
         prelines = prelines + 1
      end do
      CLOSE (UNIT_EF)

      if (prelines == 1 .and. trim(adjustl(extension))=='.json') then
         rawFileInfo%numero = countLinesInJSONOneLiner(filename, UNIT_EF)      
         call readLinesFromJSONOneLiner(rawFileInfo, filename, UNIT_EF)
      else 
         rawFileInfo%numero = prelines
         call readLines(rawFileInfo, filename, UNIT_EF)
      end if

      do k=1,rawFileInfo%numero
          linea => rawFileInfo%lineas (k)
          do j=1,linea%len
              i=j
              buscaespa: do while ((ichar(linea%dato(i:i))==32).or.(ichar(linea%dato(i:i))==9))
                 if ((ichar(linea%dato(i+1:i+1))==32).or.(ichar(linea%dato(i+1:i+1))==9)) then
                     linea%dato = trim (adjustl(linea%dato(1:i)))//' '//trim (adjustl(linea%dato(i+2:linea%len)))
                 end if
                 i=i+1
                 if (i>linea%len) exit buscaespa
              end do buscaespa
          end do
          !update
          linea%dato =  trim (adjustl(linea%dato))
          linea%LEN=len_trim (adjustl(linea%dato))   
     end do


      return
   end subroutine carga_raw_info


   end subroutine semba_init  


   function semba_create_solver(this) result (res)
      class(semba_fdtd_t) :: this
      type(solver_t) :: res
      res = solver_ctor(this%sgg,this%media,this%tag_numbers,& 
                        this%SINPML_fullsize,this%fullsize, & 
                        this%finishedwithsuccess, this%eps0,this%mu0, & 
                        this%tagtype,this%l, this%maxSourceValue, & 
                        this%time_desdelanzamiento)
   end function

   subroutine semba_update_after_simulation(this, success, sgg, eps, mu, media)
      class(semba_fdtd_t) :: this
      logical :: success
      type(SGGFDTDINFO_t) :: sgg
      type(media_matrices_t) :: media
      real(kind=rkind) :: eps ,mu
      this%finishedwithsuccess = success
      this%sgg = sgg
      this%eps0 = eps
      this%mu0 = mu
      this%media = media
   end subroutine
   
   subroutine semba_launch(this)
      class(semba_fdtd_t) :: this
      type(solver_t) :: solver
      character(len=BUFSIZE) :: dubuf
      logical :: dummylog

      ! call each simulation   !ojo que los layoutnumbers empiezan en 0
      if (this%l%finaltimestep /= 0) then
#ifdef CompileWithMPI
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         this%finishedwithsuccess=.false.
         solver = this%create_solver()
#ifdef CompileWithMTLN
         solver%mtln_parsed =  this%mtln_parsed
#endif
         if ((this%l%finaltimestep >= 0).and.(.not.this%l%skindepthpre)) then
            call solver%launch_simulation()
            call this%update_after_simulation(solver%finishedwithsuccess, solver%sgg, solver%eps0,solver%mu0,solver%media)

            deallocate(this%media%sggMiEx, this%media%sggMiEy, this%media%sggMiEz,this%media%sggMiHx, this%media%sggMiHy, this%media%sggMiHz,this%media%sggMiNo,this%media%sggMtag)
         else
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
            call get_secnds (this%l%time_out2)
            if (this%l%layoutnumber == 0) then
               call print_credits(this%l)
               write(dubuf,*) 'BEGUN '//trim (adjustl(this%l%nEntradaRoot)),' at ', this%time_comienzo%fecha(7:8), &
               & '/', this%time_comienzo%fecha(5:6), '/', this%time_comienzo%fecha(1:4),' , ',  &
               & this%time_comienzo%hora(1:2), ':', this%time_comienzo%hora(3:4)
               call print11 (this%l%layoutnumber, dubuf)
               write(dubuf,*) 'ENDED '//trim (adjustl(this%l%nEntradaRoot)),' at ', this%l%time_out2%fecha(7:8), &
               & '/', this%l%time_out2%fecha(5:6), '/', this%l%time_out2%fecha(1:4),' , ',  &
               & this%l%time_out2%hora(1:2), ':', this%l%time_out2%hora(3:4)
               call print11 (this%l%layoutnumber, dubuf)
               write(dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
               call print11 (this%l%layoutnumber, dubuf)
               call print11 (this%l%layoutnumber, dubuf)
            end if
            call CLOSEWARNINGFILE(this%l%layoutnumber,this%l%num_procs,dummylog,this%l%stochastic,this%l%simu_devia) !aqui ya no se tiene en cuenta el this%l%fatalerror
#ifdef CompileWithMPI
            call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
#ifdef CompileWithMPI
            call MPI_FINALIZE (this%l%ierr)
#endif
            stop
         end if
      end if
      !
#ifdef CompileWithMPI
      call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
   end subroutine semba_launch

   subroutine semba_end(this)
      class(semba_fdtd_t) :: this
      character(len=BUFSIZE) :: dubuf
      logical :: existe  
      character(len=BUFSIZE) :: filenombre= ' '

      if (this%l%layoutnumber == 0) then
         if (this%l%run) then
            open(38, file='running')
            write(38, '(a)') '!END'
            CLOSE (38,status='delete')
         end if
         write(dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         call print11 (this%l%layoutnumber, dubuf)
         write(dubuf,*) 'DONE :  ', trim (adjustl(this%l%nEntradaRoot)), ' UNTIL n=', this%l%finaltimestep
         call print11 (this%l%layoutnumber, dubuf)
         write(dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         call print11 (this%l%layoutnumber, dubuf)
         call erasesignalingfiles(this%l%simu_devia)

      end if

#ifdef CompileWithMPI
      call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
      !
      if (this%l%deleteintermediates) then
         write(dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         call print11 (this%l%layoutnumber, dubuf)
         write(dubuf,*) 'Attempting to delete all intermediate data files'
         call print11 (this%l%layoutnumber, dubuf)
         write(dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         call print11 (this%l%layoutnumber, dubuf)
         inquire(file=trim(adjustl(this%l%nEntradaRoot))//'_Outputrequests_'//trim(adjustl(this%whoamishort))//'.txt', EXIST=existe)
         if (existe) then
            open(19, file=trim(adjustl(this%l%nEntradaRoot))//'_Outputrequests_'//trim(adjustl(this%whoamishort))//'.txt')
            buscafile: DO
               READ (19, '(a)', end=76) filenombre
               if (trim(adjustl(filenombre)) == '!END') then
                  EXIT buscafile
               ELSE
                  open(34, file=trim(adjustl(filenombre)))
                  write(34,*) '!END'
                  CLOSE (34, STATUS='delete')
               end if
            end do buscafile
   76       continue
            CLOSE (19, STATUS='delete')
            if (this%l%layoutnumber == 0) then
               open(33, file=trim(adjustl(this%l%nEntradaRoot))//'_Outputlists.dat')
               write(33,*) '!END'
               CLOSE (33, STATUS='delete')
            end if
         end if
      end if
      !
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      call get_secnds (this%l%time_out2)
      if (this%l%layoutnumber == 0) then
         call print_credits(this%l)
         write(dubuf,*) 'BEGUN '//trim (adjustl(this%l%nEntradaRoot)),' at ', this%time_comienzo%fecha(7:8), &
         & '/', this%time_comienzo%fecha(5:6), '/', this%time_comienzo%fecha(1:4),' , ',  &
         & this%time_comienzo%hora(1:2), ':', this%time_comienzo%hora(3:4)
         call print11 (this%l%layoutnumber, dubuf)
         write(dubuf,*) 'ENDED '//trim (adjustl(this%l%nEntradaRoot)),' at ', this%l%time_out2%fecha(7:8), &
         & '/', this%l%time_out2%fecha(5:6), '/', this%l%time_out2%fecha(1:4),' , ',  &
         & this%l%time_out2%hora(1:2), ':', this%l%time_out2%hora(3:4)
         call print11 (this%l%layoutnumber, dubuf)
         write(dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         call print11 (this%l%layoutnumber, dubuf)
         call print11 (this%l%layoutnumber, dubuf)
      end if
      inquire(file='relaunch', EXIST=this%l%relaunching)
#ifdef CompileWithMPI
      call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif

#ifdef keeppause
      if (this%l%fatalerror) then
         fatalerror_aux=.true.
#ifdef CompileWithMPI
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
         call MPI_AllReduce(fatalerror_aux, this%l%fatalerror, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, this%l%ierr)
#else
         this%l%fatalerror = fatalerror_aux
#endif
      if (this%l%fatalerror) this%l%relaunching=.true.
#ifdef CompileWithMPI
      call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
   end if
#endif

      if (this%l%relaunching.and.(.not.this%finishedwithsuccess)) then
         if (this%l%layoutnumber == 0) then
            call print11 (this%l%layoutnumber, SEPARADOR//SEPARADOR)
            call print11 (this%l%layoutnumber, 'Not finishing solicited either manually or by an error condition. Edit of create launch file and remove pause file ')
            call print11 (this%l%layoutnumber, SEPARADOR//SEPARADOR)
            open(9, file='pause', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9)
            open(9, file='relaunch', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9,status='delete')
         end if
         !!!!!
#ifdef CompileWithMPI
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         if (this%l%layoutnumber == 0) then
            call CloseReportingFiles
         end if
         ! GO TO 652
      end if
   !si ha acabado con exito sal borrando signal files
      if (this%finishedwithsuccess) then
         if (this%l%layoutnumber == 0) then
            open(9, file='pause', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9,status='delete')
            open(9, file='relaunch', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9,status='delete')
            open(9, file='running', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9,status='delete')
      end if
      end if

#ifdef CompileWithMPI
         call MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif

      if (this%l%layoutnumber == 0) then
         call CloseReportingFiles
      end if
      !**************************************************************************************************

#ifdef CompileWithMPI
      call MPI_FINALIZE (this%l%ierr)
#endif
   end subroutine semba_end

   subroutine initEntrada(input)
      type(entrada_t), intent(inout) :: input
      input%geomfile = ' ';
      input%prefix = ' ';
      input%fichin = ' ';
      input%chain2 = ' ';
      input%opcionestotales = ' ' 
      input%nEntradaRoot = ' ';
      input%fileFDE = ' ';
      input%fileH5 = ' '
      input%prefixopci = ' ';
      input%prefixopci1 = ' ';
      input%opcionespararesumeo = ' ';
      input%opcionesoriginales = ' '
      input%slicesoriginales = ' ';
      input%chdummy = ' ';
      input%flushsecondsFields=0.;
      input%flushsecondsData=0.;
      input%time_end=0. 
      input%existeNFDE=.false.;
      input%existeh5=.false.
      input%creditosyaprinteados=.false.
      call input%EpsMuTimeScale_input_parameters%init0()

   end subroutine

end module SEMBA_FDTD_m