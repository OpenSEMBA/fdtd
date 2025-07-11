module SEMBA_FDTD_mod

   USE version
   USE Report
   USE Getargs
   !
   USE fdetypes
   USE Solver_mod         
   USE Resuming
   !nfde parser stuff
   USE NFDETypes                
   use nfde_rotate_m           


#ifdef CompilePrivateVersion  
   USE ParseadorClass
#endif

#ifdef CompileWithSMBJSON
   USE smbjson, only: fdtdjson_parser_t => parser_t
#endif

   USE Preprocess_m
   USE storeData
   USE xdmf_h5
   !
#ifdef CompileWithMPI
   USE MPIcomm
   USE build_t_linea_mpi
#ifdef CompileWithStochastic
   use MPI_stochastic
#endif
#endif

#ifdef CompileWithConformal
   USE CONFORMAL_INI_CLASS
   USE CONFORMAL_TOOLS
   USE CONFORMAL_MAPPED
   USE CONFORMAL_TYPES
   USE Conformal_TimeSteps_m
#endif
   use EpsMuTimeScale_m

   use interpreta_switches_m
   use, intrinsic:: iso_fortran_env, only: stdin=>input_unit

   IMPLICIT NONE

   ! should eps0 and mu0 be global variables?

   type, public :: semba_fdtd_t 
      type (entrada_t) :: l
      TYPE (tiempo_t) :: time_comienzo
      real (KIND=8) time_desdelanzamiento
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES) , allocatable , dimension(:,:,:) ::  sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz
      integer (KIND=IKINDMTAG) , allocatable , dimension(:,:,:) :: sggMtag
      type (SGGFDTDINFO)   :: sgg
      type (limit_t), DIMENSION (1:6) :: fullsize, SINPML_fullsize
      real (KIND=RKIND) ::  eps0,mu0,cluz
      real (KIND=RKIND) :: maxSourceValue
      character (LEN=BUFSIZE) :: whoami, whoamishort
#ifdef CompileWithMTLN
      type(mtln_t) :: mtln_parsed
#endif
      type (taglist_t) :: tag_numbers
      type (tagtype_t) :: tagtype
      logical :: finishedwithsuccess

   contains
      procedure :: init => semba_init   
      procedure :: launch => semba_launch   
      procedure :: end => semba_end   
   end type semba_fdtd_t 

   
contains

   subroutine semba_init(this, input_flags)
      class(semba_fdtd_t) :: this
      character (len=*), optional :: input_flags

      real (KIND=RKIND) :: dtantesdecorregir
      real (KIND=RKIND)   ::  dxmin,dymin,dzmin,dtlay
      
      logical :: dummylog,l_auxinput, l_auxoutput, ThereArethinslots
      logical :: hayinput
      logical :: lexis
      logical :: newrotate !300124 tiramos con el rotador antiguo

      character (LEN=BUFSIZE) ::  f= ' ', chain = ' ', chain3 = ' ',chain4 = ' ', chaindummy= ' '
      character (LEN=BUFSIZE_LONG) :: slices = ' '
      character (LEN=BUFSIZE) :: dubuf
      character (LEN=BUFSIZE) :: buff
      character (LEN=BUFSIZE) :: filename_h5bin ! File name

      integer (KIND=4) :: myunit,jmed
      integer (kind=4) :: finaltimestepantesdecorregir,NEWfinaltimestep,thefileno
      integer (kind=4) :: statuse
      integer (KIND=4) ::  status, i, field
      INTEGER (KIND=4) ::  verdadero_mpidir
      integer (kind=4) :: my_iostat


      type (Parseador), POINTER :: parser
      type (t_NFDE_FILE), POINTER :: NFDE_FILE
      type(solver_t) :: solver 
         
#ifdef CompileWithMPI
      LOGICAL :: fatalerror_aux
      TYPE (XYZlimit_t), DIMENSION (1:6) :: tempalloc
#endif

   integer (kind=4) :: conf_err
#ifdef CompileWithConformal
      type (conf_conflicts_t), pointer  :: conf_conflicts
#endif
      ! call sleep(5)
      call initEntrada(this%l) 
      newrotate=.false.       !!ojo tocar luego                     
   !!200918 !!!si se lanza con -pscal se overridea esto
      this%eps0= 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
      this%mu0 = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
      this%cluz=1.0_RKIND/sqrt(this%eps0*this%mu0)
      
      CALL OnPrint

#ifdef CompileWithMPI
      CALL InitGeneralMPI (this%l%layoutnumber, this%l%size)
      SUBCOMM_MPI=MPI_COMM_WORLD !default el this%l%stochastic es el global a menos que luego se divida
#else
      this%l%size = 1
      this%l%layoutnumber = 0
#endif
      call setglobal(this%l%layoutnumber,this%l%size) !para crear variables globales con info MPI
         
      WRITE (this%whoamishort, '(i5)') this%l%layoutnumber + 1
      WRITE (this%whoami, '(a,i5,a,i5,a)') '(', this%l%layoutnumber + 1, '/', this%l%size, ') '
         
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      call get_secnds(this%l%time_out2)
      this%time_desdelanzamiento= this%l%time_out2%segundos
#ifndef keeppause
      if (this%l%layoutnumber==0) then
         OPEN (38, file='running')
         write (38,*) '!END'
         CLOSE (38,status='delete')
         OPEN (38, file='pause')
         write (38,*) '!END'
         CLOSE (38,status='delete')
         OPEN (38, file='relaunch')
         write (38,*) '!END'
         CLOSE (38,status='delete')
         OPEN (38, file='forcestop')
         write (38,*) '!END'
         CLOSE (38,status='delete')
      endif
#endif

   if (this%l%layoutnumber==0) then
         my_iostat=0
   3443  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' 
         OPEN (11, file='SEMBA_FDTD_temp.log',err=3443,iostat=my_iostat,action='write')
         write (11,*) '!END'
         CLOSE (11,status='delete')
         my_iostat=0
   3447  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',3447,'.',this%l%layoutnumber,'SEMBA_FDTD_temp.log' 
         OPEN (11, file='SEMBA_FDTD_temp.log',err=3447,iostat=my_iostat,status='new',action='write')
         call print_credits(this%l)
         CLOSE (11)
   endif

#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif


   652 continue

      CALL CLOSEWARNINGFILE(this%l%layoutnumber,this%l%size,dummylog,.false.,.false.) !aqui ya no se tiene en cuenta el this%l%fatalerror

      WRITE (this%l%opcionespararesumeo, '(a,i4,a)') 'mpirun -n ', this%l%size,' '
      call default_flags(this%l)    !set all default flags

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      call get_secnds(this%time_comienzo)
      !temporarily until later
      IF (this%l%layoutnumber == 0) THEN
         OPEN (11, file='SEMBA_FDTD_temp.log',position='append')
         this%l%file11isopen=.true.
      END IF
      !

#ifdef CompileWithMPI
      !wait until everything comes out
      CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif

      !see if there is semaphore to pause continuing
      INQUIRE (file='pause', EXIST=this%l%pausar)
#ifdef CompileWithMPI
      this%l%l_aux = this%l%pausar
      CALL MPI_AllReduce (this%l%l_aux, this%l%pausar, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, this%l%ierr)
#endif
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      CALL get_secnds (this%l%time_out2)
      this%l%time_begin = this%l%time_out2%segundos
      WRITE (dubuf,*) 'Paused at              ', this%l%time_out2%fecha(7:8), '/', this%l%time_out2%fecha(5:6), '/', &
      &                this%l%time_out2%fecha(1:4), '  ', this%l%time_out2%hora(1:2), ':', this%l%time_out2%hora(3:4)
      IF (this%l%pausar) CALL print11 (this%l%layoutnumber, dubuf)
      DO while (this%l%pausar)
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
         CALL get_secnds (this%l%time_out2)
         this%l%time_end = this%l%time_out2%segundos
         IF (this%l%time_end-this%l%time_begin > 10.0_RKIND) THEN
            INQUIRE (file='pause', EXIST=this%l%pausar)
#ifdef CompileWithMPI
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
            this%l%l_aux = this%l%pausar
            CALL MPI_AllReduce (this%l%l_aux, this%l%pausar, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, this%l%ierr)
            call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
            CALL get_secnds (this%l%time_out2)
            this%l%time_begin = this%l%time_out2%segundos
            WRITE (dubuf,*) 'Paused at              ', this%l%time_out2%fecha(7:8), '/', this%l%time_out2%fecha(5:6), '/', &
            &                this%l%time_out2%fecha(1:4), ' ', this%l%time_out2%hora(1:2), ':', this%l%time_out2%hora(3:4)
            IF (this%l%pausar) CALL print11 (this%l%layoutnumber, dubuf)
         END IF
      END DO
      !fin del semaphoro

#ifdef keeppause   
      INQUIRE (file='forcestop', EXIST=this%l%forcestop)
      if (this%l%forcestop) then
         if (this%l%layoutnumber==0) then
            OPEN (38, file='running')
            write (38,*) '!END'
            CLOSE (38,status='delete')
            OPEN (38, file='pause')
            write (38,*) '!END'
            CLOSE (38,status='delete')
            OPEN (38, file='relaunch')
            write (38,*) '!END'
            CLOSE (38,status='delete')
            OPEN (38, file='forcestop')
            write (38,*) '!END'
            CLOSE (38,status='delete')
         endif
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
         CALL MPI_FINALIZE (this%l%ierr)
#endif
         STOP
      endif
#endif

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      CALL get_secnds (this%l%time_out2)
      
   
      if (present(input_flags)) then 
         this%l%read_command_line = .false.
         this%l%chain2 = input_flags
         this%l%length = len(input_flags)
      else
      ! mira el command_line y el fichero launch 251022
         CALL get_command (this%l%chain2, this%l%length, status)
         IF (status /= 0) then
            CALL stoponerror (this%l%layoutnumber, this%l%size, 'General error',.true.); goto 652
         endif
      end if

      this%l%chain2=trim(adjustl(this%l%chain2))
      !concatena con lo que haya en launch
      INQUIRE (file='launch', EXIST=hayinput)
      if (hayinput) then
         OPEN (9, file='launch', FORM='formatted',action='read')
         READ (9, '(a)') chain3
         chain3=trim(adjustl(chain3))
         CLOSE (9)               
         print *,'----> launch input file '//trim(adjustl(chain3))
      endif
#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif


      this%l%chain2=trim(adjustl(this%l%chain2))//' '//trim(adjustl(chain3))

      call buscaswitchficheroinput(this%l)
      

      IF (status /= 0) then
         CALL stoponerror (this%l%layoutnumber, this%l%size, 'Error in searching input file. Correct and remove pause file',.true.); goto 652
      endif
   !!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!
      call print_credits(this%l)
      if (trim(adjustl(this%l%extension))=='.nfde') then   
#ifdef CompilePrivateVersion   
         call cargaNFDE(this%l%filefde,parser)
#else
         print *,'Not compiled with cargaNFDEINDEX'
         stop
#endif
#ifdef CompileWithSMBJSON
      elseif (trim(adjustl(this%l%extension))=='.json') then
         call cargaFDTDJSON(this%l%fichin, parser)
#endif
      else
         print *, 'Neither .nfde nor .json files used as input after -i'
         stop
      endif
      


      this%sgg%extraswitches=parser%switches
   !!!da preferencia a los switches por linea de comando
      CALL getcommandargument (this%l%chain2, 1, chaindummy, this%l%length, statuse, getBinaryPath())

      this%l%chain2=trim(adjustl(this%l%chain2))
      chaindummy=trim(adjustl(chaindummy))
      this%l%length=len(trim(adjustl(chaindummy)))
      this%l%chain2=trim(adjustl(chaindummy))//' '//trim(adjustl(this%sgg%extraswitches))//' '//trim(adjustl(this%l%chain2(this%l%length+1:)))               
      this%l%chaininput=trim(adjustl(this%l%chain2))
   !!!!
      

      call interpreta(this%l,status )      
      this%sgg%nEntradaRoot=trim (adjustl(this%l%nEntradaRoot))

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
   9083   CALL stoponerror (0, this%l%size, 'Invalid _h5bin.txt file',.true.); statuse=-1; !return
      endif
#ifdef CompileWithMPI
         !wait until everything comes out
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         stop
      endif
#endif

      IF (status /= 0) then
         call print11(this%l%layoutnumber,'Remove running and pause files. If error persists check switches for error.  '//this%l%chain2,.true.)
         call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' '); call print11(this%l%layoutnumber,' ');  goto 652
      endif

      call set_priorities(this%l%prioritizeCOMPOoverPEC,this%l%prioritizeISOTROPICBODYoverall,this%l%prioritizeTHINWIRE) !!! asigna las prioridades
      if (this%l%finaltimestep /= -2) then
         ! nfde part
         CALL print11 (this%l%layoutnumber, 'INIT conversion internal ASCII => Binary')
         CALL print11 (this%l%layoutnumber, SEPARADOR//SEPARADOR//SEPARADOR)

         CALL print11 (this%l%layoutnumber, SEPARADOR//SEPARADOR//SEPARADOR)
         !!!!!!!!!!!!!!!!!!!!!!
         call NFDE2sgg
         this%l%fatalerror=this%l%fatalerror.or.this%l%fatalerrornfde2sgg
         !!!!!!!!!!!!!!!!!!!!!
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         CALL print11 (this%l%layoutnumber, '[OK] Ended conversion internal ASCII => Binary')
         !release memory created by newPARSER
         if (this%l%fatalerror) then
            if (allocated(this%sggMiEx)) deallocate (this%sggMiEx, this%sggMiEy, this%sggMiEz,this%sggMiHx, this%sggMiHy, this%sggMiHz,this%sggMiNo,this%sggMtag)
            CALL stoponerror (this%l%layoutnumber, this%l%size, 'Error in .nfde file syntax. Check all *Warnings* and *tmpWarnings* files, correct and remove pause file if any',.true.); goto 652
         endif

         !*************************************************************************
         !***[conformal] ******************************************
         !*************************************************************************
         !conformal conformal ini          ref: ##Confini##
#ifdef CompileWithConformal
      if (this%l%input_conformal_flag) then

            !md notes:
            ![1]      Todos los procesos parsean el archivo -conf completo.
            ![2]      El parseador es INDEPENDIENTE de del resto del problema (dimensiones,
            !         particion MPI, ... )
            ![3]      Posteriormente conf_mesh obtenido por el parseador sera tratado por cada
            !         proceso atendiedo al resto del porblema y la particion MPI

            conf_parameter%output_file_report_id = 47;
            !......................................................................
         write(dubuf,*) 'Init Searching for Conformal Mesh ...';  call print11(this%l%layoutnumber,dubuf)
#ifdef CompileWithMPI
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
            CALL conformal_ini (TRIM(this%l%conformal_file_input_name),trim(this%l%fileFDE),parser,&
               &this%sgg, this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz,this%l%run_with_abrezanjas,&
               &this%fullsize,this%l%layoutnumber,this%l%mpidir, this%l%input_conformal_flag,conf_err,this%l%verbose)
#endif
            !......................................................................
#ifndef CompileWithMPI
            !CALL conformal_ini (TRIM(this%l%conformal_file_input_name),trim(this%l%fileFDE),sgg,fullsize,0,conf_err,this%l%verbose)
         CALL conformal_ini (TRIM(this%l%conformal_file_input_name),trim(this%l%fileFDE),parser,&
               &this%sgg, this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz,&
               &this%l%run_with_abrezanjas,this%fullsize,0,this%l%mpidir,this%l%input_conformal_flag,conf_err,this%l%verbose)
#endif
            if(conf_err/=0)then
               call WarnErrReport(Trim(buff),.true.)
            end if

#ifdef CompilePrivateVersion  
         if (trim(adjustl(this%l%extension))=='.nfde') then
         CALL Destroy_Parser (parser)  
         DEALLOCATE (NFDE_FILE%lineas)
         DEALLOCATE (NFDE_FILE)
         nullify (NFDE_FILE)
         endif
#endif      
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CompileWithMPI
         !wait until everything comes out
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
            l_auxinput = this%l%input_conformal_flag
            call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
            call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, this%l%ierr)
            this%l%input_conformal_flag = l_auxoutput
#endif
            !......................................................................
#ifdef CompileWithMPI
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif       
            if (this%l%resume.and.this%l%flag_conf_sgg) then
                  CALL stoponerror (this%l%layoutnumber, this%l%size, 'this%l%resume -r currently unsupported by conformal solver',.true.); statuse=-1; !return
            end if
            if (this%l%input_conformal_flag.and.this%l%flag_conf_sgg) then
               write(dubuf,*) '----> Conformal Mesh found';  call print11(this%l%layoutnumber,dubuf)
            else   
               write(dubuf,*) '----> No Conformal Mesh found';  call print11(this%l%layoutnumber,dubuf)
            endif
      end if !FIN DEL: if (this%l%input_conformal_flag) then
      
#endif

         !*************************************************************************
         !*************************************************************************
         !*************************************************************************

#ifdef CompileWithConformal
         !*************************************************************************
         !***[conformal] ******************************************
         !*************************************************************************
         !conformal mapped reff: ##Confmapped##

         !call creamatricesdedibujoencadaslabmpi(sgg%alloc(iEx)%XI,....,sgg%Sweep(iEx)%...)

         if (this%l%input_conformal_flag) then
               write(dubuf,*) '----> this%l%input_conformal_flag True and init';  call print11(this%l%layoutnumber,dubuf)
            call conf_geometry_mapped_for_UGRDTD (&
            &conf_conflicts, &
            &this%sgg,this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz, &
            &this%fullsize, this%SINPML_fullsize,this%l%layoutnumber,conf_err,this%l%verbose);
            !call conf_geometry_mapped_for_UGRDTD (sgg, fullsize, this%SINPML_fullsize,this%l%layoutnumber,conf_err,this%l%verbose); //refactor JUL15
            if(conf_err==0)then
            else
               buff=''; buff = 'Program aborted.';
               call WarnErrReport(Trim(buff),.true.)
            end if
               write(dubuf,*) '----> this%l%input_conformal_flag True and exit';  call print11(this%l%layoutnumber,dubuf)
         end if

#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         !*************************************************************************
         !*************************************************************************
         !*************************************************************************
#endif

         if (allocated(this%sggMiEx)) then !para el this%l%skindepthpre no se allocatea nada
#ifdef CompileWithConformal
         call AssigLossyOrPECtoNodes(this%sgg,this%sggMiNo,this%sggMiEx,this%sggMiEy,this%sggMiEz,&
                                       &conf_conflicts,this%l%input_conformal_flag)
#else
         call AssigLossyOrPECtoNodes(this%sgg,this%sggMiNo,this%sggMiEx,this%sggMiEy,this%sggMiEz)
#endif
         IF (this%l%createmap) CALL store_geomData (this%sgg,this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz, this%l%geomfile)
         endif
         !
#ifdef CompileWithMPI
         !wait until everything comes out
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
      endif
      write(dubuf,*) '[OK] Ended Conformal Mesh';  call print11(this%l%layoutnumber,dubuf)
      if (this%l%finaltimestep==0) this%l%finaltimestep=this%sgg%TimeSteps !no quitar
      IF (this%l%forcesteps) then
         this%sgg%TimeSteps = this%l%finaltimestep
      else
         this%l%finaltimestep = this%sgg%TimeSteps
      endif
      IF (.not.this%l%forcesteps) then
            finaltimestepantesdecorregir=this%l%finaltimestep
            this%l%finaltimestep=int(dtantesdecorregir/this%sgg%dt*finaltimestepantesdecorregir)
#ifdef CompileWithMPI
            call MPI_AllReduce( this%l%finaltimestep, NEWfinaltimestep, 1_4, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, this%l%ierr)
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
            this%l%finaltimestep=NEWfinaltimestep
#endif
            if (finaltimestepantesdecorregir/=this%l%finaltimestep) then
               write(dubuf,*) SEPARADOR//separador//separador
               call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) 'Original Final Time Step= ',finaltimestepantesdecorregir
               if (this%l%layoutnumber==0) call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) 'Corrected Final Time Step= ',this%l%finaltimestep
               if (this%l%layoutnumber==0) call print11(this%l%layoutnumber,dubuf)
            endif
      endif
      !check that simulation can actually be done for the kind of media requested
      DO i = 1, this%sgg%nummedia
         IF (this%sgg%Med(i)%Is%ThinWire) THEN
#ifndef CompileWithBerengerWires
      if  ((this%l%wiresflavor=='berenger')) then
            CALL stoponerror (this%l%layoutnumber, this%l%size, 'Berenger Wires without support. Recompile!')
      endif
#endif
#ifndef CompileWithSlantedWires
      if  ((this%l%wiresflavor=='slanted').or.(this%l%wiresflavor=='semistructured')) then
            CALL stoponerror (this%l%layoutnumber, this%l%size, 'slanted Wires without support. Recompile!')
      endif
#endif
            CONTINUE
         END IF
         !
         IF ((this%sgg%Med(i)%Is%AnisMultiport) .OR. (this%sgg%Med(i)%Is%multiport).OR. (this%sgg%Med(i)%Is%SGBC)) THEN
#ifndef CompileWithNIBC
            if (this%l%mibc) CALL stoponerror (this%l%layoutnumber, this%l%size, 'this%l%mibc Multiports without support. Recompile!')
#endif
            CONTINUE
         END IF
   !altair no conformal sgbc 201119
#ifdef NoConformalSGBC
         IF (this%sgg%Med(i)%Is%sgbc .and. this%l%input_conformal_flag) THEN
            CALL stoponerror (this%l%layoutnumber, this%l%size, 'Conformal sgbc not allowed. ')
         END IF
#endif
   !    
      END DO
      
      
      IF (this%l%thereare_stoch.and.(.not.this%l%chosenyesornostochastic)) THEN
         CALL stoponerror (this%l%layoutnumber, this%l%size, '!STOCH found in .nfde. Specify either -stoch or -nostoch')
      END IF
#ifndef CompileWithSlantedWires
      IF (this%l%hay_slanted_wires) THEN
         CALL stoponerror (this%l%layoutnumber, this%l%size, 'slanted wires without slanted support. Recompile ()')
      END IF
#endif   
      IF (this%l%hay_slanted_wires .AND. ((trim(adjustl(this%l%wiresflavor))/='slanted').AND.(trim(adjustl(this%l%wiresflavor))/='semistructured'))) THEN
         CALL stoponerror (this%l%layoutnumber, this%l%size, 'slanted wires require -this%l%wiresflavor Slanted/semistructured')
      endif

      
      !Error abrezanjas y no this%l%resume conformal
      ThereArethinslots=.FALSE.
      do jmed=1,this%sgg%NumMedia
         if (this%sgg%Med(jmed)%Is%ThinSlot) ThereArethinslots=.true.
      end do
      if (this%l%resume.and.this%l%run_with_abrezanjas.and.ThereArethinslots) then   
            CALL stoponerror (this%l%layoutnumber, this%l%size, 'this%l%resume -r currently unsupported by conformal solver',.true.); statuse=-1; !return
      end if
      !
   !!!SOME FINAL REPORTING

      if (this%l%layoutnumber==0) then
         WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         CALL print11 (this%l%layoutnumber, dubuf)
         CALL print11 (this%l%layoutnumber, 'Solver launched with options:')
         write(dubuf,*) this%l%mibc          
         CALL print11 (this%l%layoutnumber, '---> this%l%mibc    solver for NIBC multilayer: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%ade         
         CALL print11 (this%l%layoutnumber, '---> this%l%ade     solver for ADC multilayer: '//trim(adjustl(dubuf)))
         Write(dubuf,*) this%l%sgbc    
         CALL print11 (this%l%layoutnumber, '---> sgbc    solver for multilayer: '//trim(adjustl(dubuf)))
         if (this%l%sgbc) then
               write(dubuf,*) this%l%sgbcDispersive      
               CALL print11 (this%l%layoutnumber, '---> sgbc DISPERSIVE solver for multilayer: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%sgbccrank     
               CALL print11 (this%l%layoutnumber, '---> sgbc Crank-Nicolson solver for multilayer: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%sgbcdepth
               CALL print11 (this%l%layoutnumber, '---> sgbc Depth: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%sgbcfreq
               CALL print11 (this%l%layoutnumber, '---> sgbc Freq: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%sgbcresol
               CALL print11 (this%l%layoutnumber, '---> sgbc Resol: '//trim(adjustl(dubuf)))
         endif
         write(dubuf,*) this%l%skindepthpre
         CALL print11 (this%l%layoutnumber, '---> this%l%skindepthpre preprocessing for multilayer: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%flag_conf_sgg
         CALL print11 (this%l%layoutnumber, '---> Conformal file external: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%input_conformal_flag      
         CALL print11 (this%l%layoutnumber, '---> Conformal solver: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%run_with_abrezanjas
         CALL print11 (this%l%layoutnumber, '---> Conformal thin-gap solver: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%run_with_dmma
         CALL print11 (this%l%layoutnumber, '---> DMMA thin-gap solver: '//trim(adjustl(dubuf)))
         write(dubuf,'(a)') this%l%wiresflavor
         CALL print11 (this%l%layoutnumber, '---> Wire model: '//trim(adjustl(dubuf)))
         write(dubuf,'(a)') this%l%inductance_model
         CALL print11 (this%l%layoutnumber, '---> Inductance model: '//trim(adjustl(dubuf)))
         if (trim(adjustl(this%l%wiresflavor))=='berenger') then
               write(dubuf,*) this%l%mindistwires
               CALL print11 (this%l%layoutnumber, '---> Berenger minimum distance between wires: '//trim(adjustl(dubuf)))
               write(dubuf,*) this%l%mtlnberenger
               CALL print11 (this%l%layoutnumber, '---> Berenger -this%l%mtlnberenger MTLN switch: '//trim(adjustl(dubuf)))
         endif
         if (trim(adjustl(this%l%wiresflavor))=='holland') then
               write(dubuf,*) this%l%stableradholland                 
               CALL print11 (this%l%layoutnumber, '---> Holland -this%l%stableradholland automatic correction switch: '//trim(adjustl(dubuf)))
         endif
         write(dubuf,*) this%l%TAPARRABOS                
         CALL print11 (this%l%layoutnumber, '---> Thin-wire double-tails removed: '//trim(adjustl(dubuf)))
         write(dubuf,*) this%l%fieldtotl                
         CALL print11 (this%l%layoutnumber, '---> Thin-wire -this%l%fieldtotl experimental switch: '//trim(adjustl(dubuf)))
         WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         CALL print11 (this%l%layoutnumber, dubuf)
      endif

      IF (this%l%layoutnumber == 0) THEN
         call erasesignalingfiles(this%l%simu_devia)
      endif
      
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
               write(thefileno,'(a)') '# (  3.5 ,  3.5 ) '//trim(adjustl('sgbc/this%l%mibc Isotropic/anisotropic Multiport                               (Line)'))
               write(thefileno,'(a)') '# (  300 ,  399 ) '//trim(adjustl('sgbc/this%l%mibc Isotropic/anisotropic Multiport (+indexmedium)                (Surface)'))
               write(thefileno,'(a)') '# (  4.5 ,  4.5 ) '//trim(adjustl('Thin slot                                                               (Line)'))
               write(thefileno,'(a)') '# (  5.0 ,  5.0 ) '//trim(adjustl('Already_YEEadvanced_byconformal                                         (Surface)'))
               write(thefileno,'(a)') '# (  5.5 ,  5.5 ) '//trim(adjustl('Already_YEEadvanced_byconformal                                         (Line)'))
               write(thefileno,'(a)') '# (  6.0 ,  6.0 ) '//trim(adjustl('Split_and_useless                                                       (Surface)'))
               write(thefileno,'(a)') '# (  6.5 ,  6.5 ) '//trim(adjustl('Split_and_useless                                                       (Line)'))
               write(thefileno,'(a)') '# (  7.0 ,  7.0 ) '//trim(adjustl('Edge Not colliding thin wires                                           (Line)'))
               write(thefileno,'(a)') '# (  8.0 ,  8.0 ) '//trim(adjustl('Thin wire segments colliding with structure                             (Line)'))
               write(thefileno,'(a)') '# (  8.5 ,  8.5 ) '//trim(adjustl('Soft/Hard Nodal CURRENT/FIELD ELECTRIC DENSITY SOURCE                   (Line)'))
               write(thefileno,'(a)') '# (  9.0 ,  9.0 ) '//trim(adjustl('Soft/Hard Nodal CURRENT/FIELD MAGNETIC DENSITY SOURCE                   (Line)'))
               write(thefileno,'(a)') '# (   10 ,   11 ) '//trim(adjustl('LeftEnd/RightEnd/Ending wire segment                                                 (Wire)'))
               write(thefileno,'(a)') '# (   20 ,   20 ) '//trim(adjustl('Intermediate wire segment +number_holland_parallel or +number_berenger       (Wire) '))
               write(thefileno,'(a)') '# (  400 ,  499 ) '//trim(adjustl('Thin slot (+indexmedium)                                                (Surface)'))
               write(thefileno,'(a)') '# ( -0.5 , -0.5 ) '//trim(adjustl('Other types of media                                                    (Line)'))
               write(thefileno,'(a)') '# ( -1.0 , -1.0 ) '//trim(adjustl('Other types of media                                                    (Surface)'))
         close(thefileno)
      endif

contains 
   subroutine NFDE2sgg     
   !!!!!!!!!      
         real (kind=rkind) :: dt,finaldt
         logical fatalerror
         ! parser now holds all the .nfde info
         !first read the limits
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         CALL read_limits_nogeom (this%l%layoutnumber,this%l%size, this%sgg, this%fullsize, this%SINPML_fullsize, parser,this%l%MurAfterPML,this%l%mur_exist)
      
         dtantesdecorregir=this%sgg%dt
         !!!!!corrige el delta de t si es necesario !sgg15 310715 bug distintos sgg%dt !!!!!!!!!!

         dxmin=minval(this%sgg%DX)
         dymin=minval(this%sgg%DY)
         dzmin=minval(this%sgg%DZ)
         !!!
         dtlay=(1.0_RKIND/(this%cluz*sqrt(((1.0_RKIND / dxmin)**2.0_RKIND )+((1.0_RKIND / dymin)**2.0_RKIND )+((1.0_RKIND / dzmin)**2.0_RKIND ))))
         dt=dtlay
#ifdef CompileWithMPI
         call MPIupdateMin(dtlay,dt)
#endif

         !!!write(dubuf,*) SEPARADOR//separador//separador
         !!!call print11(this%l%layoutnumber,dubuf)
         !!!write(dubuf,*) '--->dt,dxmin,dymin,dzmin,sgg%dt  ',dt,dxmin,dymin,dzmin,sgg%dt
         !!!call print11(this%l%layoutnumber,dubuf)
         !!!write(dubuf,*) SEPARADOR//separador//separador
         !!!call print11(this%l%layoutnumber,dubuf)

         if (this%l%forcecfl) then
            this%sgg%dt=dt*this%l%cfl
            write(dubuf,*) SEPARADOR//separador//separador
            call print11(this%l%layoutnumber,dubuf)
            write(dubuf,*) 'Correcting sgg%dt with -this%l%cfl switch. New time step: ',this%sgg%dt
            call print11(this%l%layoutnumber,dubuf)
            write(dubuf,*) SEPARADOR//separador//separador
            call print11(this%l%layoutnumber,dubuf)
         else
            if (this%sgg%dt > dt*heurCFL) then
               write(dubuf,*) SEPARADOR//separador//separador
               call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) 'Automatically correcting dt for stability reasons: '
               call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) 'Original dt: ',this%sgg%dt
               call print11(this%l%layoutnumber,dubuf)
               this%sgg%dt=dt*heurCFL
               write(dubuf,*) 'New dt: ',this%sgg%dt
               call print11(this%l%layoutnumber,dubuf)
               write(dubuf,*) SEPARADOR//separador//separador
               call print11(this%l%layoutnumber,dubuf)
            endif
         endif
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
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(this%l%layoutnumber,dubuf)
         if (this%l%mur_exist.and.this%l%mur_first) then
            this%l%mur_second=.false.
         else
            this%l%mur_second=.false. !arreglar cuando se arregle el bug de las mur second
            this%l%mur_first=.true. !arreglar cuando se arregle el bug de las mur second
         endif
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
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
         IF (this%l%size == 1) THEN
            this%sgg%Alloc(1:6)%ZI = this%fullsize(1:6)%ZI - 1
            this%sgg%Alloc(1:6)%ZE = this%fullsize(1:6)%ZE + 1
            !REDUCE THE SWEEP AREA BY 1
            this%sgg%Sweep(1:6)%ZI = this%fullsize(1:6)%ZI
            this%sgg%Sweep(1:6)%ZE = this%fullsize(1:6)%ZE
            !!incluido aqui pq se precisa para clip 16/07/15
            DO field = iEx, iHz
               this%sgg%SINPMLSweep(field)%XI = Max (this%SINPML_fullsize(field)%XI, this%sgg%Sweep(field)%XI)
               this%sgg%SINPMLSweep(field)%XE = Min (this%SINPML_fullsize(field)%XE, this%sgg%Sweep(field)%XE)
               this%sgg%SINPMLSweep(field)%YI = Max (this%SINPML_fullsize(field)%YI, this%sgg%Sweep(field)%YI)
               this%sgg%SINPMLSweep(field)%YE = Min (this%SINPML_fullsize(field)%YE, this%sgg%Sweep(field)%YE)
               this%sgg%SINPMLSweep(field)%ZI = Max (this%SINPML_fullsize(field)%ZI, this%sgg%Sweep(field)%ZI)
               this%sgg%SINPMLSweep(field)%ZE = Min (this%SINPML_fullsize(field)%ZE, this%sgg%Sweep(field)%ZE)
            END DO
            !!fin 16/07/15
            WRITE (dubuf,*) 'INIT NFDE --------> GEOM'
            CALL print11 (this%l%layoutnumber, dubuf)
            CALL read_geomData (this%sgg,this%sggMtag,this%tag_numbers, this%sggMiNo,this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz, this%l%fichin, this%l%layoutnumber, this%l%size, this%SINPML_fullsize, this%fullsize, parser, &
            this%l%groundwires,this%l%attfactorc,this%l%mibc,this%l%sgbc,this%l%sgbcDispersive,this%l%MEDIOEXTRA,this%maxSourceValue,this%l%skindepthpre,this%l%createmapvtk,this%l%input_conformal_flag,this%l%CLIPREGION,this%l%boundwireradius,this%l%maxwireradius,this%l%updateshared,this%l%run_with_dmma, this%eps0, &
            this%mu0,.false.,this%l%hay_slanted_wires,this%l%verbose,this%l%ignoresamplingerrors,this%tagtype,this%l%wiresflavor)
#ifdef CompileWithMTLN
            if (trim(adjustl(this%l%extension))=='.json')  then 
               this%mtln_parsed = parser%mtln
               this%mtln_parsed%time_step = this%sgg%dt
            end if
            ! if (trim(adjustl(this%l%extension))=='.json')  mtln_solver = mtlnCtor(parser%mtln)   
#endif
            WRITE (dubuf,*) '[OK] ENDED NFDE --------> GEOM'
            CALL print11 (this%l%layoutnumber, dubuf)
            !writing
            slices = '!SLICES'
            WRITE (buff, '(i7)') this%sgg%Sweep(iHz)%ZE - this%sgg%Sweep(iHz)%ZI
            slices = trim (adjustl(slices)) // '_' // trim (adjustl(buff))
            IF (this%l%resume .AND. (slices /= this%l%slicesoriginales)) THEN
               buff='Different resumed/original MPI slices: '//trim(adjustl(slices))//' '//&
               & trim(adjustl(this%l%slicesoriginales))
               CALL stoponerror (this%l%layoutnumber, this%l%size, buff)
            END IF
            CALL print11 (this%l%layoutnumber, trim(adjustl(slices)))
            !end writing
            WRITE (buff, '(a,i7,a,i7)') '_________Spanning from z=', this%sgg%Sweep(iHz)%ZI, ' to z=', this%sgg%Sweep(iHz)%ZE
            CALL print11 (this%l%layoutnumber, trim(adjustl(buff)))
#ifdef CompileWithMPI
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#ifdef CompileWithStochastic
            if (this%l%stochastic) then
               buff='this%l%stochastic uncompatible with MPI this%l%size smaller than 2'
               CALL stoponerror (this%l%layoutnumber, this%l%size, buff)
            endif
#endif
#endif
         ELSE !del this%l%size==1       
#ifdef CompileWithMPI
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#ifdef CompileWithStochastic
            if (this%l%stochastic) then
               call HalvesStochasticMPI(this%l%layoutnumber,this%l%size,this%l%simu_devia)
            endif
#endif
                     
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)   
   !!!ahora divide el espacio computacional
            CALL MPIdivide (this%sgg, this%fullsize, this%SINPML_fullsize, this%l%layoutnumber, this%l%size, this%l%forcing, this%l%forced, this%l%slicesoriginales, this%l%resume,this%l%fatalerror)
            !
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)   
            if (this%l%fatalerror) then
   !intenta recuperarte
               return
            endif
      
            ! if the layout is pure PML then take at least a line of non PML to build the PML data insider read_geomDAta
            ! Uses extra memory but later matrix sggm is deallocated in favor of smaller sggMIEX, etc
            DO field = iEx, iHz
               tempalloc(field)%ZE = this%sgg%Alloc(field)%ZE
               tempalloc(field)%ZI = this%sgg%Alloc(field)%ZI
               this%sgg%Alloc(field)%ZE = Max (this%sgg%Alloc(field)%ZE, this%SINPML_fullsize(field)%ZI+1)
               this%sgg%Alloc(field)%ZI = Min (this%sgg%Alloc(field)%ZI, this%SINPML_fullsize(field)%ZE-1)
            END DO
            !   
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)  
            !!incluido aqui pq se precisa para clip 16/07/15
            DO field = iEx, iHz
               this%sgg%SINPMLSweep(field)%XI = Max (this%SINPML_fullsize(field)%XI, this%sgg%Sweep(field)%XI)
               this%sgg%SINPMLSweep(field)%XE = Min (this%SINPML_fullsize(field)%XE, this%sgg%Sweep(field)%XE)
               this%sgg%SINPMLSweep(field)%YI = Max (this%SINPML_fullsize(field)%YI, this%sgg%Sweep(field)%YI)
               this%sgg%SINPMLSweep(field)%YE = Min (this%SINPML_fullsize(field)%YE, this%sgg%Sweep(field)%YE)
               this%sgg%SINPMLSweep(field)%ZI = Max (this%SINPML_fullsize(field)%ZI, this%sgg%Sweep(field)%ZI)
               this%sgg%SINPMLSweep(field)%ZE = Min (this%SINPML_fullsize(field)%ZE, this%sgg%Sweep(field)%ZE)
            END DO
            !!fin 16/07/15
            WRITE (dubuf,*) 'INIT NFDE --------> GEOM'
            CALL print11 (this%l%layoutnumber, dubuf)           

            CALL read_geomData (this%sgg,this%sggMtag,this%tag_numbers, this%sggMiNo,this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz, this%l%fichin, this%l%layoutnumber, this%l%size, this%SINPML_fullsize, this%fullsize, parser, &
            this%l%groundwires,this%l%attfactorc,this%l%mibc,this%l%sgbc,this%l%sgbcDispersive,this%l%MEDIOEXTRA,this%maxSourceValue,this%l%skindepthpre,this%l%createmapvtk,this%l%input_conformal_flag,this%l%CLIPREGION,this%l%boundwireradius,this%l%maxwireradius,this%l%updateshared,this%l%run_with_dmma, &
            this%eps0,this%mu0,this%l%simu_devia,this%l%hay_slanted_wires,this%l%verbose,this%l%ignoresamplingerrors,this%tagtype,this%l%wiresflavor)


#ifdef CompileWithMPI
            !wait until everything comes out
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
#ifdef CompileWithMTLN
            if (trim(adjustl(this%l%extension))=='.json')  then 
               this%mtln_parsed = parser%mtln
               this%mtln_parsed%time_step = this%sgg%dt
            end if
#endif
            WRITE (dubuf,*) '[OK] ENDED NFDE --------> GEOM'
            CALL print11 (this%l%layoutnumber, dubuf)
            !restore back the indexes
            DO field = iEx, iHz
               this%sgg%Alloc(field)%ZE = tempalloc(field)%ZE
               this%sgg%Alloc(field)%ZI = tempalloc(field)%ZI
            END DO
#endif
            CONTINUE
         END IF !del this%l%size==1
         !
#ifdef CompileWithMPI
         !wait until everything comes out
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         !!!!!!!!!!!!!lo dejo aqui debajo tambien aunque ya se ha calculado antes para lo del clipping
         DO field = iEx, iHz
            this%sgg%SINPMLSweep(field)%XI = Max (this%SINPML_fullsize(field)%XI, this%sgg%Sweep(field)%XI)
            this%sgg%SINPMLSweep(field)%XE = Min (this%SINPML_fullsize(field)%XE, this%sgg%Sweep(field)%XE)
            this%sgg%SINPMLSweep(field)%YI = Max (this%SINPML_fullsize(field)%YI, this%sgg%Sweep(field)%YI)
            this%sgg%SINPMLSweep(field)%YE = Min (this%SINPML_fullsize(field)%YE, this%sgg%Sweep(field)%YE)
            this%sgg%SINPMLSweep(field)%ZI = Max (this%SINPML_fullsize(field)%ZI, this%sgg%Sweep(field)%ZI)
            this%sgg%SINPMLSweep(field)%ZE = Min (this%SINPML_fullsize(field)%ZE, this%sgg%Sweep(field)%ZE)
         END DO
         return
      end subroutine

   end subroutine semba_init  


   subroutine semba_launch(this)
      class(semba_fdtd_t) :: this
      type(solver_t) :: solver
      character (LEN=BUFSIZE) :: dubuf
      logical :: dummylog

      ! call each simulation   !ojo que los layoutnumbers empiezan en 0
      IF (this%l%finaltimestep /= 0) THEN
#ifdef CompileWithMPI
         !wait until everything comes out
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         this%finishedwithsuccess=.false.


#ifdef CompileWithMTLN
         solver%mtln_parsed =  this%mtln_parsed
#endif

         if ((this%l%finaltimestep >= 0).and.(.not.this%l%skindepthpre)) then
            CALL solver%launch_simulation (this%sgg,this%sggMtag,this%tag_numbers,this%sggMiNo, this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz,&
                                           this%SINPML_fullsize,this%fullsize,this%finishedwithsuccess,this%eps0,this%mu0,this%tagtype,&
                                           this%l, this%maxSourceValue, this%time_desdelanzamiento)

            deallocate (this%sggMiEx, this%sggMiEy, this%sggMiEz,this%sggMiHx, this%sggMiHy, this%sggMiHz,this%sggMiNo,this%sggMtag)
         else
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
            CALL get_secnds (this%l%time_out2)
            IF (this%l%layoutnumber == 0) THEN
               call print_credits(this%l)
               WRITE (dubuf,*) 'BEGUN '//trim (adjustl(this%l%nEntradaRoot)),' at ', this%time_comienzo%fecha(7:8), &
               & '/', this%time_comienzo%fecha(5:6), '/', this%time_comienzo%fecha(1:4),' , ',  &
               & this%time_comienzo%hora(1:2), ':', this%time_comienzo%hora(3:4)
               CALL print11 (this%l%layoutnumber, dubuf)
               WRITE (dubuf,*) 'ENDED '//trim (adjustl(this%l%nEntradaRoot)),' at ', this%l%time_out2%fecha(7:8), &
               & '/', this%l%time_out2%fecha(5:6), '/', this%l%time_out2%fecha(1:4),' , ',  &
               & this%l%time_out2%hora(1:2), ':', this%l%time_out2%hora(3:4)
               CALL print11 (this%l%layoutnumber, dubuf)
               WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
               CALL print11 (this%l%layoutnumber, dubuf)
               CALL print11 (this%l%layoutnumber, dubuf)
            ENDIF
            !!!!!!!        CALL CLOSEdxfFILE(this%l%layoutnumber,this%l%size)
            CALL CLOSEWARNINGFILE(this%l%layoutnumber,this%l%size,dummylog,this%l%stochastic,this%l%simu_devia) !aqui ya no se tiene en cuenta el this%l%fatalerror
#ifdef CompileWithMPI
            !wait until everything comes out
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
#ifdef CompileWithMPI
            CALL MPI_FINALIZE (this%l%ierr)
#endif
            stop
         endif
      END IF
      !
#ifdef CompileWithMPI
      !wait until everything comes out
      CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif

   end subroutine semba_launch

   subroutine semba_end(this)
      class(semba_fdtd_t) :: this
      character (LEN=BUFSIZE) :: dubuf
      logical :: existe  
      character (LEN=BUFSIZE) :: filenombre= ' '

      IF (this%l%layoutnumber == 0) THEN
         if (this%l%run) then
            OPEN (38, file='running')
            WRITE (38, '(a)') '!END'
            CLOSE (38,status='delete')
         endif
         WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         CALL print11 (this%l%layoutnumber, dubuf)
         WRITE (dubuf,*) 'DONE :  ', trim (adjustl(this%l%nEntradaRoot)), ' UNTIL n=', this%l%finaltimestep
         CALL print11 (this%l%layoutnumber, dubuf)
         WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         CALL print11 (this%l%layoutnumber, dubuf)
         call erasesignalingfiles(this%l%simu_devia)

      END IF

#ifdef CompileWithMPI
      !wait until everything comes out
      CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
      !
      IF (this%l%deleteintermediates) THEN
         WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         CALL print11 (this%l%layoutnumber, dubuf)
         WRITE (dubuf,*) 'Attempting to delete all intermediate data files'
         CALL print11 (this%l%layoutnumber, dubuf)
         WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         CALL print11 (this%l%layoutnumber, dubuf)
         INQUIRE (file=trim(adjustl(this%l%nEntradaRoot))//'_Outputrequests_'//trim(adjustl(this%whoamishort))//'.txt', EXIST=existe)
         IF (existe) THEN
            OPEN (19, file=trim(adjustl(this%l%nEntradaRoot))//'_Outputrequests_'//trim(adjustl(this%whoamishort))//'.txt')
            buscafile: DO
               READ (19, '(a)', end=76) filenombre
               IF (trim(adjustl(filenombre)) == '!END') THEN
                  EXIT buscafile
               ELSE
                  OPEN (34, file=trim(adjustl(filenombre)))
                  WRITE (34,*) '!END'
                  CLOSE (34, STATUS='delete')
               END IF
            END DO buscafile
   76       CONTINUE
            CLOSE (19, STATUS='delete')
            IF (this%l%layoutnumber == 0) THEN
               OPEN (33, file=trim(adjustl(this%l%nEntradaRoot))//'_Outputlists.dat')
               WRITE (33,*) '!END'
               CLOSE (33, STATUS='delete')
            END IF
         END IF
      END IF
      !

      !**************************************************************************************************
      !***[conformal] *******************************************************************
      !**************************************************************************************************
      !delete conformal memory   reff: ##Conf_end##
#ifdef CompileWithConformal
      if(this%l%input_conformal_flag)then
         call conf_sMesh%delete
         call conf_timeSteps%delete;
         call delete_conf_tools();
      end if
#endif
      !**************************************************************************************************
      !**************************************************************************************************
      !**************************************************************************************************

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,this%l%ierr)
#endif
      CALL get_secnds (this%l%time_out2)
      IF (this%l%layoutnumber == 0) THEN
         call print_credits(this%l)
         WRITE (dubuf,*) 'BEGUN '//trim (adjustl(this%l%nEntradaRoot)),' at ', this%time_comienzo%fecha(7:8), &
         & '/', this%time_comienzo%fecha(5:6), '/', this%time_comienzo%fecha(1:4),' , ',  &
         & this%time_comienzo%hora(1:2), ':', this%time_comienzo%hora(3:4)
         CALL print11 (this%l%layoutnumber, dubuf)
         WRITE (dubuf,*) 'ENDED '//trim (adjustl(this%l%nEntradaRoot)),' at ', this%l%time_out2%fecha(7:8), &
         & '/', this%l%time_out2%fecha(5:6), '/', this%l%time_out2%fecha(1:4),' , ',  &
         & this%l%time_out2%hora(1:2), ':', this%l%time_out2%hora(3:4)
         CALL print11 (this%l%layoutnumber, dubuf)
         WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
         CALL print11 (this%l%layoutnumber, dubuf)
         CALL print11 (this%l%layoutnumber, dubuf)
      ENDIF
      INQUIRE (file='relaunch', EXIST=this%l%relaunching)
#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
      ! Error reading check

#ifdef keeppause
      if (this%l%fatalerror) then
         fatalerror_aux=.true.
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
         call MPI_AllReduce(fatalerror_aux, this%l%fatalerror, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, this%l%ierr)
#else
         this%l%fatalerror = fatalerror_aux
#endif
      if (this%l%fatalerror) this%l%relaunching=.true.
#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
   endif
#endif

      IF (this%l%relaunching.and.(.not.this%finishedwithsuccess)) THEN
         IF (this%l%layoutnumber == 0) THEN
            CALL print11 (this%l%layoutnumber, SEPARADOR//SEPARADOR)
            CALL print11 (this%l%layoutnumber, 'Not finishing solicited either manually or by an error condition. Edit of create launch file and remove pause file ')
            CALL print11 (this%l%layoutnumber, SEPARADOR//SEPARADOR)
            OPEN (9, file='pause', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9)
            OPEN (9, file='relaunch', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9,status='delete')
         endif
         !!!!!
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
         IF (this%l%layoutnumber == 0) THEN
            CALL CloseReportingFiles
         endif
         ! GO TO 652
      END IF
   !si ha acabado con exito sal borrando signal files
      IF (this%finishedwithsuccess) THEN
         IF (this%l%layoutnumber == 0) THEN
            OPEN (9, file='pause', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9,status='delete')
            OPEN (9, file='relaunch', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9,status='delete')
            OPEN (9, file='running', FORM='formatted')
            write (9, '(a)') ' '
            CLOSE (9,status='delete')
      endif
      endif

#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif

      IF (this%l%layoutnumber == 0) THEN
         CALL CloseReportingFiles
      endif
      !**************************************************************************************************

#ifdef CompileWithMPI
      CALL MPI_FINALIZE (this%l%ierr)
#endif
      ! STOP
      !

   end subroutine semba_end

   subroutine initEntrada(input)
      type(entrada_t), intent(inout) :: input
#ifdef CompileWithConformal
      input%conformal_file_input_name=char(0);  
#endif
      input%geomfile = ' ';
      input%prefix = ' ';input%fichin = ' '; input%chain2 = ' '; input%opcionestotales = ' ' 
      input%nEntradaRoot = ' '; input%fileFDE = ' '; input%fileH5 = ' '
      input%prefixopci = ' '; input%prefixopci1 = ' ';input%opcionespararesumeo = ' '; input%opcionesoriginales = ' '
      input%slicesoriginales = ' '; ; input%chdummy = ' '
      input%flushsecondsFields=0.; input%flushsecondsData=0.; input%time_end=0. 
      input%existeNFDE=.false.; input%existeconf=.false.; input%existecmsh=.false.; input%existeh5=.false.
      input%creditosyaprinteados=.false.
      call input%EpsMuTimeScale_input_parameters%init0()

   end subroutine

#ifdef CompileWithSMBJSON
   subroutine cargaFDTDJSON(filename, parsed)
      character(len=1024), intent(in) :: filename
      type(Parseador), pointer :: parsed
      
      character(len=:), allocatable :: usedFilename    
      type(fdtdjson_parser_t) :: parser
      
      usedFilename = adjustl(trim(filename)) // ".json"
      parser = fdtdjson_parser_t(usedFilename)
      
      allocate(parsed)
      parsed = parser%readProblemDescription()
   end subroutine cargaFDTDJSON
#endif

#ifdef CompilePrivateVersion 
   subroutine cargaNFDE(local_nfde,local_parser)
      CHARACTER (LEN=BUFSIZE) :: local_nfde
      TYPE (Parseador), POINTER :: local_parser
      INTEGER (KIND=8) :: numero,i8,troncho,longitud
      integer (kind=4) :: mpi_t_linea_t,longitud4
      IF (this%l%existeNFDE) THEN
         WRITE (dubuf,*) 'INIT Reading file '//trim (adjustl(this%whoami))//' ', trim (adjustl(local_nfde))
         CALL print11 (this%l%layoutnumber, dubuf)
   !!!!!!!!!!!!!!!!!!!!!!!
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
         if (this%l%layoutnumber==0) then
            NFDE_FILE => cargar_NFDE_FILE (local_nfde)
         !!!ya se allocatea dentro
         else
            ALLOCATE (NFDE_FILE)
         endif
         !
         write(dubuf,*) '[OK]';  call print11(this%l%layoutnumber,dubuf)
         !--->
         WRITE (dubuf,*) 'INIT Sharing file through MPI'; CALL print11 (this%l%layoutnumber, dubuf)
         !
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
         !
         numero=NFDE_FILE%numero
         call MPI_BCAST(numero, 1_4, MPI_INTEGER8, 0_4, SUBCOMM_MPI, this%l%ierr)      
         if (this%l%layoutnumber/=0) then
            NFDE_FILE%targ = 1
            NFDE_FILE%numero=numero
            ALLOCATE (NFDE_FILE%lineas(NFDE_FILE%numero))
         endif
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
         !CREAMOS EL DERIVED TYPE y lo enviamos !para evitar el error de Marconi asociado a PSM2_MQ_RECVREQS_MAX 100617

         CALL build_derived_t_linea(mpi_t_linea_t)

         !problema del limite de mandar mas de 2^29 bytes con MPI!!!  Los soluciono partiendo en maxmpibytes (2^27) (algo menos por prudencia)! 040716
         troncho=ceiling(maxmpibytes*1.0_8/(BUFSIZE*1.0_8+8.0_8),8)
         do i8=1,numero,troncho
            longitud=min(troncho,numero-i8+1)
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
            if ((longitud>huge(1_4)).or.(longitud>maxmpibytes)) then
               print *,'Stop. Buggy error: MPI longitud greater that greatest integer*4'
               stop
            else
               longitud4=int(longitud,4)
            endif
            call MPI_BCAST(NFDE_FILE%lineas(i8),longitud4,mpi_t_linea_t,0_4,SUBCOMM_MPI,this%l%ierr)    
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
            CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
         end do
#else
      NFDE_FILE => cargar_NFDE_FILE (local_nfde)
#endif
         write(dubuf,*) '[OK]';  call print11(this%l%layoutnumber,dubuf)
         !--->
      END IF    
      NFDE_FILE%mpidir=this%l%mpidir
      WRITE (dubuf,*) 'INIT interpreting geometrical data from ', trim (adjustl(local_nfde))
      CALL print11 (this%l%layoutnumber, dubuf)
      if(newrotate) then
         verdadero_mpidir=NFDE_FILE%mpidir
         NFDE_FILE%mpidir=3
      endif
      local_parser => newparser (NFDE_FILE)         
#ifdef CompileWithMPI            
      CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
      if(newrotate) then      
         NFDE_FILE%mpidir=verdadero_mpidir
         call nfde_rotate (local_parser,NFDE_FILE%mpidir)
#ifdef CompileWithMPI            
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
      endif
      this%l%thereare_stoch=NFDE_FILE%thereare_stoch
      this%l%mpidir=NFDE_FILE%mpidir !bug 100419
      write(dubuf,*) '[OK] '//trim(adjustl(this%whoami))//' newparser (NFDE_FILE)';  call print11(this%l%layoutnumber,dubuf)       
#ifdef CompileWithMPI            
         CALL MPI_Barrier (SUBCOMM_MPI, this%l%ierr)
#endif
      return

   end subroutine cargaNFDE
#endif


end module SEMBA_FDTD_mod
!
