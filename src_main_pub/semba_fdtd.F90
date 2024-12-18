!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MIT License
! 
! Copyright (c) 2023 University of Granada
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SEMBA_FDTD LAUNCHER MODULE
!  Creation date Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! in Unix Systems checkt that ulimit -n gives high numbers (more than 2e6 files may be open at the same time!!!)
! in Ubuntu this is in /etc/security/limits.conf
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

PROGRAM SEMBA_FDTD_launcher

   USE version
   USE Report
   USE Getargs
   !
   USE fdetypes
   USE Solver         !!!el timestepping.F90
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

   !*************************************************
   !***[conformal] ******************
   !*************************************************
#ifdef CompileWithConformal
   USE CONFORMAL_INI_CLASS
   USE CONFORMAL_TOOLS
   USE CONFORMAL_MAPPED
   USE CONFORMAL_TYPES
   USE Conformal_TimeSteps_m
#endif
   use EpsMuTimeScale_m
   !*************************************************
   !*************************************************
   !*************************************************
   !
   use interpreta_switches_m
   IMPLICIT NONE
   !
!!!24118 pscaling
   REAL (KIND=RKIND)              ::  eps0,mu0,cluz
!!!241018 fin pscaling
   integer (KIND=IKINDMTAG) , allocatable , dimension(:,:,:) ::  sggMtag
   type(taglist_t) :: tag_numbers
   integer (KIND=INTEGERSIZEOFMEDIAMATRICES) , allocatable , dimension(:,:,:) ::  sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz
   
   LOGICAL :: dummylog,finishedwithsuccess,l_auxinput, l_auxoutput, ThereArethinslots
   integer (KIND=4) :: myunit,jmed
   REAL (KIND=RKIND) :: maxSourceValue
   !
   REAL (KIND=RKIND) :: dtantesdecorregir
   integer (kind=4) :: finaltimestepantesdecorregir,NEWfinaltimestep,thefileno
   TYPE (Parseador), POINTER :: parser
   type (SGGFDTDINFO)   :: sgg
   TYPE (limit_t), DIMENSION (1:6) :: fullsize, SINPML_fullsize
   !
   LOGICAL :: existe  
   INTEGER (KIND=4) ::  status, i, field
   CHARACTER (LEN=BUFSIZE) ::  f, chain, chain3,chain4, chaindummy, filenombre

   CHARACTER (LEN=BUFSIZE_LONG) :: slices
   CHARACTER (LEN=BUFSIZE) :: whoami, whoamishort
   CHARACTER (LEN=BUFSIZE) :: dubuf
   integer (kind=4) :: statuse

   LOGICAL :: hayinput
   !
   TYPE (t_NFDE_FILE), POINTER :: NFDE_FILE

   type (tagtype_t) :: tagtype
   REAL (KIND=RKIND)   ::  dxmin,dymin,dzmin,dtlay
      
#ifdef CompileWithMPI
   LOGICAL :: fatalerror_aux
   TYPE (XYZlimit_t), DIMENSION (1:6) :: tempalloc
#endif
   TYPE (tiempo_t) :: time_comienzo
   CHARACTER (LEN=BUFSIZE) :: buff
   REAL (KIND=8) time_desdelanzamiento
   CHARACTER (LEN=BUFSIZE) :: filename_h5bin ! File name

   !****************************************************************************
   !****************************************************************************
   !conformal existence flags   ref: ##Confflag##
   integer (kind=4) :: conf_err
#ifdef CompileWithConformal
   type (conf_conflicts_t), pointer  :: conf_conflicts
#endif

   !****************************************************************************
   !****************************************************************************
   !****************************************************************************

   type (entrada_t) :: l    
#ifdef CompileWithMTLN
   type(mtln_t) :: mtln_parsed
#endif
   logical :: lexis
   integer (kind=4) :: my_iostat
   
   INTEGER (KIND=4) ::  verdadero_mpidir
   logical :: newrotate !300124 tiramos con el rotador antiguo

   newrotate=.false.       !!ojo tocar luego                     
!!200918 !!!si se lanza con -pscal se overridea esto
   Eps0= 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
   Mu0 = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
   cluz=1.0_RKIND/sqrt(eps0*mu0)
!!!   
#ifdef CompileWithConformal
   l%conformal_file_input_name=char(0);  
#endif
   slices = ' '; chain3 = ' ';chain4 = ' ' ;chaindummy=' '
    l%geomfile = ' '; filenombre = ' '
!!!   
   l%prefix = ' ';l%fichin = ' '; f = ' '; chain = ' '; l%chain2 = ' '; l%opcionestotales = ' ' 
   l%nEntradaRoot = ' '; l%fileFDE = ' '; l%fileH5 = ' '
   l%prefixopci = ' '; l%prefixopci1 = ' ';l%opcionespararesumeo = ' '; l%opcionesoriginales = ' '
   l%slicesoriginales = ' '; ; l%chdummy = ' '
   l%flushsecondsFields=0.; l%flushsecondsData=0.; l%time_end=0. 
   l%existeNFDE=.false.; l%existeconf=.false.; l%existecmsh=.false.; l%existeh5=.false.
   l%creditosyaprinteados=.false.
   !activate printing through screen
   CALL OnPrint
   !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!
   call l%EpsMuTimeScale_input_parameters%init0()

#ifdef CompileWithMPI
   CALL InitGeneralMPI (l%layoutnumber, l%size)
   SUBCOMM_MPI=MPI_COMM_WORLD !default el l%stochastic es el global a menos que luego se divida
#else
   l%size = 1
   l%layoutnumber = 0
#endif
    call setglobal(l%layoutnumber,l%size) !para crear variables globales con info MPI
    
   if (l%size.gt.maxcores) then
       print *,'Maximum cores ',maxcores,' reached.  to recompile'
       stop
   endif
       
   WRITE (whoamishort, '(i5)') l%layoutnumber + 1
   WRITE (whoami, '(a,i5,a,i5,a)') '(', l%layoutnumber + 1, '/', l%size, ') '

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
   call get_secnds(l%time_out2)
   time_desdelanzamiento= l%time_out2%segundos
#ifndef keeppause
   if (l%layoutnumber==0) then
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

  if (l%layoutnumber==0) then
      my_iostat=0
3443  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' 
      OPEN (11, file='SEMBA_FDTD_temp.log',err=3443,iostat=my_iostat,action='write')
      write (11,*) '!END'
      CLOSE (11,status='delete')
      my_iostat=0
3447  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',3447,'.',l%layoutnumber,'SEMBA_FDTD_temp.log' 
      OPEN (11, file='SEMBA_FDTD_temp.log',err=3447,iostat=my_iostat,status='new',action='write')
      call print_credits(l)
      CLOSE (11)
  endif

#ifdef CompileWithMPI
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif


652 continue

   CALL CLOSEWARNINGFILE(l%layoutnumber,l%size,dummylog,.false.,.false.) !aqui ya no se tiene en cuenta el l%fatalerror

   WRITE (l%opcionespararesumeo, '(a,i4,a)') 'mpirun -n ', l%size,' '
   call default_flags(l)    !set all default flags

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
   call get_secnds(time_comienzo)
   !temporarily until later
   IF (l%layoutnumber == 0) THEN
      OPEN (11, file='SEMBA_FDTD_temp.log',position='append')
      l%file11isopen=.true.
   END IF
   !

#ifdef CompileWithMPI
   !wait until everything comes out
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif

   !see if there is semaphore to pause continuing
   INQUIRE (file='pause', EXIST=l%pausar)
#ifdef CompileWithMPI
   l%l_aux = l%pausar
   CALL MPI_AllReduce (l%l_aux, l%pausar, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, l%ierr)
#endif
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
   CALL get_secnds (l%time_out2)
   l%time_begin = l%time_out2%segundos
   WRITE (dubuf,*) 'Paused at              ', l%time_out2%fecha(7:8), '/', l%time_out2%fecha(5:6), '/', &
   &                l%time_out2%fecha(1:4), '  ', l%time_out2%hora(1:2), ':', l%time_out2%hora(3:4)
   IF (l%pausar) CALL print11 (l%layoutnumber, dubuf)
   DO while (l%pausar)
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
      CALL get_secnds (l%time_out2)
      l%time_end = l%time_out2%segundos
      IF (l%time_end-l%time_begin > 10.0_RKIND) THEN
         INQUIRE (file='pause', EXIST=l%pausar)
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
         l%l_aux = l%pausar
         CALL MPI_AllReduce (l%l_aux, l%pausar, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, l%ierr)
         call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
         CALL get_secnds (l%time_out2)
         l%time_begin = l%time_out2%segundos
         WRITE (dubuf,*) 'Paused at              ', l%time_out2%fecha(7:8), '/', l%time_out2%fecha(5:6), '/', &
         &                l%time_out2%fecha(1:4), ' ', l%time_out2%hora(1:2), ':', l%time_out2%hora(3:4)
         IF (l%pausar) CALL print11 (l%layoutnumber, dubuf)
      END IF
   END DO
   !fin del semaphoro

#ifdef keeppause   
   INQUIRE (file='forcestop', EXIST=l%forcestop)
   if (l%forcestop) then
       if (l%layoutnumber==0) then
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
        CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
        CALL MPI_FINALIZE (l%ierr)
#endif
        STOP
   endif
#endif

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
   CALL get_secnds (l%time_out2)
   !
   ! mira el command_line y el fichero launch 251022
   CALL get_command (l%chain2, l%length, status)
   IF (status /= 0) then
      CALL stoponerror (l%layoutnumber, l%size, 'General error',.true.); goto 652
   endif

   l%chain2=trim(adjustl(l%chain2))
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
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif


   l%chain2=trim(adjustl(l%chain2))//' '//trim(adjustl(chain3))

   call buscaswitchficheroinput(l)
   

   IF (status /= 0) then
       CALL stoponerror (l%layoutnumber, l%size, 'Error in searching input file. Correct and remove pause file',.true.); goto 652
   endif
!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!
   call print_credits(l)
   if (trim(adjustl(l%extension))=='.nfde') then   
#ifdef CompilePrivateVersion   
       call cargaNFDE(l%filefde,parser)
#else
       print *,'Not compiled with cargaNFDEINDEX'
       stop
#endif
#ifdef CompileWithSMBJSON
   elseif (trim(adjustl(l%extension))=='.json') then
        call cargaFDTDJSON(l%fichin, parser)
#endif
   else
       print *, 'Neither .nfde nor .json files used as input after -i'
       stop
   endif
   


!!!!!!!!!!!!!!!!!!!!!!!
   sgg%extraswitches=parser%switches
!!!da preferencia a los switches por linea de comando
   CALL getcommandargument (l%chain2, 1, chaindummy, l%length, statuse)

   l%chain2=trim(adjustl(l%chain2))
   chaindummy=trim(adjustl(chaindummy))
   l%length=len(trim(adjustl(chaindummy)))
   l%chain2=trim(adjustl(chaindummy))//' '//trim(adjustl(sgg%extraswitches))//' '//trim(adjustl(l%chain2(l%length+1:)))               
   l%chaininput=trim(adjustl(l%chain2))
!!!!
   

   call interpreta(l,status )      
   sgg%nEntradaRoot=trim (adjustl(l%nEntradaRoot))

#ifdef CompileWithMTLN   
   if (parser%general%mtlnProblem) then 
      call launch_mtln_simulation(parser%mtln, l%nEntradaRoot, l%layoutnumber) 
      STOP
   end if
#endif

#ifdef CompileWithHDF
!!!!tunel a lo bestia para crear el .h5 a 021219
   if (l%createh5filefromsinglebin) then
     if (l%layoutnumber==0) then
       inquire(file=trim(adjustl(sgg%nEntradaRoot))//'_h5bin.txt',exist=lexis)
       if (.not.lexis) goto 9083
       open(newunit=myunit,file=trim(adjustl(sgg%nEntradaRoot))//'_h5bin.txt',form='formatted',err=9083) !lista de todos los .h5bin
       do 
           read (myunit,'(a)',end=84552) filename_h5bin
           call createh5filefromsinglebin(filename_h5bin,l%vtkindex) 
           print *, 'Processed '//trim(adjustl(filename_h5bin))
       end do
84552  close(myunit)
       print *, 'END: SUCCESS creating '//trim(adjustl(sgg%nEntradaRoot))//'_h5bin.txt'
       stop
9083   CALL stoponerror (0, l%size, 'Invalid _h5bin.txt file',.true.); statuse=-1; !return
     endif
#ifdef CompileWithMPI
      !wait until everything comes out
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      stop
   endif
#endif

   IF (status /= 0) then
      call print11(l%layoutnumber,'Remove running and pause files. If error persists check switches for error.  '//l%chain2,.true.)
      call print11(l%layoutnumber,' '); call print11(l%layoutnumber,' '); call print11(l%layoutnumber,' '); call print11(l%layoutnumber,' '); call print11(l%layoutnumber,' '); call print11(l%layoutnumber,' ');  goto 652
   endif
!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!

   call set_priorities(l%prioritizeCOMPOoverPEC,l%prioritizeISOTROPICBODYoverall,l%prioritizeTHINWIRE) !!! asigna las prioridades
   if (l%finaltimestep /= -2) then
      ! nfde part
      CALL print11 (l%layoutnumber, 'INIT conversion internal ASCII => Binary')
      CALL print11 (l%layoutnumber, SEPARADOR//SEPARADOR//SEPARADOR)

      CALL print11 (l%layoutnumber, SEPARADOR//SEPARADOR//SEPARADOR)
      !!!!!!!!!!!!!!!!!!!!!!
      call NFDE2sgg
      l%fatalerror=l%fatalerror.or.l%fatalerrornfde2sgg
      !!!!!!!!!!!!!!!!!!!!!
      !NOTE: md: necesito parser vivo hata el conformal ini, lo paso abajo
      ! CALL Destroy_Parser (parser)
      ! DEALLOCATE (NFDE_FILE%lineas)
      ! DEALLOCATE (NFDE_FILE)
      ! nullify (NFDE_FILE)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CompileWithMPI
      !wait until everything comes out
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      CALL print11 (l%layoutnumber, '[OK] Ended conversion internal ASCII => Binary')
      !release memory created by newPARSER
      if (l%fatalerror) then
!!intenta recuperarte
         if (allocated(sggMiEx)) deallocate (sggMiEx, sggMiEy, sggMiEz,sggMiHx, sggMiHy, sggMiHz,sggMiNo,sggMtag)
         CALL stoponerror (l%layoutnumber, l%size, 'Error in .nfde file syntax. Check all *Warnings* and *tmpWarnings* files, correct and remove pause file if any',.true.); goto 652
      endif
      !**********************************************
      !INIT DXF OUTPUT
      !
      !!!!    CALL INITdxfFILE (l%layoutnumber, l%size, l%nEntradaRoot)
      !************************************************
      !
      ! IF (l%createmap) CALL store_geomData (sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, l%geomfile)
      ! Se hace otra vez luego (sgg 220817)

      !*************************************************************************
      !***[conformal] ******************************************
      !*************************************************************************
      !conformal conformal ini          ref: ##Confini##
#ifdef CompileWithConformal
    if (l%input_conformal_flag) then

         !md notes:
         ![1]      Todos los procesos parsean el archivo -conf completo.
         ![2]      El parseador es INDEPENDIENTE de del resto del problema (dimensiones,
         !         particion MPI, ... )
         ![3]      Posteriormente conf_mesh obtenido por el parseador sera tratado por cada
         !         proceso atendiedo al resto del porblema y la particion MPI

         conf_parameter%output_file_report_id = 47;
         !......................................................................
        write(dubuf,*) 'Init Searching for Conformal Mesh ...';  call print11(l%layoutnumber,dubuf)
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
         CALL conformal_ini (TRIM(l%conformal_file_input_name),trim(l%fileFDE),parser,&
            &sgg, sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,l%run_with_abrezanjas,&
            &fullsize,l%layoutnumber,l%mpidir, l%input_conformal_flag,conf_err,l%verbose)
#endif
         !......................................................................
#ifndef CompileWithMPI
         !CALL conformal_ini (TRIM(l%conformal_file_input_name),trim(l%fileFDE),sgg,fullsize,0,conf_err,l%verbose)
        CALL conformal_ini (TRIM(l%conformal_file_input_name),trim(l%fileFDE),parser,&
            &sgg, sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,&
            &l%run_with_abrezanjas,fullsize,0,l%mpidir,l%input_conformal_flag,conf_err,l%verbose)
#endif
         if(conf_err/=0)then
            call WarnErrReport(Trim(buff),.true.)
         end if
      !NOTE: md: lo necesito despues del conformal init (antes se borraba mas arriba)
      !REVIEW: sgg

#ifdef CompilePrivateVersion  
      if (trim(adjustl(l%extension))=='.nfde') then
        CALL Destroy_Parser (parser)  
        DEALLOCATE (NFDE_FILE%lineas)
        DEALLOCATE (NFDE_FILE)
        nullify (NFDE_FILE)
      endif
#endif      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CompileWithMPI
       !wait until everything comes out
       CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
         !tocado para reducir en MPI 1119 con sgg pa no se que demonios quiere
         !l%input_conformal_flag = .True. ! lo fuerzo para evitar deadlocks tengo que revisarlo
         !lo que sigue debe resolve los deadlocks 1119 sgg
#ifdef CompileWithMPI
         l_auxinput = l%input_conformal_flag
         call MPI_Barrier(SUBCOMM_MPI,l%ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, l%ierr)
         l%input_conformal_flag = l_auxoutput
#endif
         !......................................................................
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif       
         if (l%resume.and.l%flag_conf_sgg) then
               CALL stoponerror (l%layoutnumber, l%size, 'l%resume -r currently unsupported by conformal solver',.true.); statuse=-1; !return
         end if
         if (l%input_conformal_flag.and.l%flag_conf_sgg) then
             write(dubuf,*) '----> Conformal Mesh found';  call print11(l%layoutnumber,dubuf)
         else   
             write(dubuf,*) '----> No Conformal Mesh found';  call print11(l%layoutnumber,dubuf)
         endif
    end if !FIN DEL: if (l%input_conformal_flag) then
    
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

      if (l%input_conformal_flag) then
             write(dubuf,*) '----> l%input_conformal_flag True and init';  call print11(l%layoutnumber,dubuf)
         call conf_geometry_mapped_for_UGRDTD (&
         &conf_conflicts, &
         &sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, &
         &fullsize, SINPML_fullsize,l%layoutnumber,conf_err,l%verbose);
         !call conf_geometry_mapped_for_UGRDTD (sgg, fullsize, SINPML_fullsize,l%layoutnumber,conf_err,l%verbose); //refactor JUL15
         if(conf_err==0)then
         else
            buff=''; buff = 'Program aborted.';
            call WarnErrReport(Trim(buff),.true.)
         end if
             write(dubuf,*) '----> l%input_conformal_flag True and exit';  call print11(l%layoutnumber,dubuf)
      end if

#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      !*************************************************************************
      !*************************************************************************
      !*************************************************************************
#endif

      !310715

      if (allocated(sggMiEx)) then !para el l%skindepthpre no se allocatea nada
#ifdef CompileWithConformal
        call AssigLossyOrPECtoNodes(sgg,sggMiNo,sggMiEx,sggMiEy,sggMiEz,&
                                    &conf_conflicts,l%input_conformal_flag)
#else
        call AssigLossyOrPECtoNodes(sgg,sggMiNo,sggMiEx,sggMiEy,sggMiEz)
#endif
        IF (l%createmap) CALL store_geomData (sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, l%geomfile)
      endif
      !
#ifdef CompileWithMPI
      !wait until everything comes out
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   endif
   write(dubuf,*) '[OK] Ended Conformal Mesh';  call print11(l%layoutnumber,dubuf)
   if (l%finaltimestep==0) l%finaltimestep=sgg%TimeSteps !no quitar
   IF (l%forcesteps) then
      sgg%TimeSteps = l%finaltimestep
   else
      l%finaltimestep = sgg%TimeSteps
   endif
!aniadido correcion timesteps finales si no hay l%forcesteps 250417
   IF (.not.l%forcesteps) then
         finaltimestepantesdecorregir=l%finaltimestep
         l%finaltimestep=int(dtantesdecorregir/sgg%dt*finaltimestepantesdecorregir)
#ifdef CompileWithMPI
         call MPI_AllReduce( l%finaltimestep, NEWfinaltimestep, 1_4, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, l%ierr)
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
         l%finaltimestep=NEWfinaltimestep
#endif
         if (finaltimestepantesdecorregir/=l%finaltimestep) then
             write(dubuf,*) SEPARADOR//separador//separador
             call print11(l%layoutnumber,dubuf)
             write(dubuf,*) 'Original Final Time Step= ',finaltimestepantesdecorregir
             if (l%layoutnumber==0) call print11(l%layoutnumber,dubuf)
             write(dubuf,*) 'Corrected Final Time Step= ',l%finaltimestep
             if (l%layoutnumber==0) call print11(l%layoutnumber,dubuf)
         endif
   endif
!fin aniadido 250417
   !check that simulation can actually be done for the kind of media requested
   DO i = 1, sgg%nummedia
      IF (sgg%Med(i)%Is%ThinWire) THEN
#ifndef CompileWithBerengerWires
    if  ((l%wiresflavor=='berenger')) then
         CALL stoponerror (l%layoutnumber, l%size, 'Berenger Wires without support. Recompile!')
    endif
#endif
#ifndef CompileWithSlantedWires
    if  ((l%wiresflavor=='slanted').or.(l%wiresflavor=='semistructured')) then
         CALL stoponerror (l%layoutnumber, l%size, 'slanted Wires without support. Recompile!')
    endif
#endif
         CONTINUE
      END IF
      !
      IF ((sgg%Med(i)%Is%AnisMultiport) .OR. (sgg%Med(i)%Is%multiport).OR. (sgg%Med(i)%Is%SGBC)) THEN
#ifndef CompileWithNIBC
         if (l%mibc) CALL stoponerror (l%layoutnumber, l%size, 'l%mibc Multiports without support. Recompile!')
#endif
         CONTINUE
      END IF
!altair no conformal sgbc 201119
#ifdef NoConformalSGBC
      IF (sgg%Med(i)%Is%sgbc .and. l%input_conformal_flag) THEN
         CALL stoponerror (l%layoutnumber, l%size, 'Conformal sgbc not allowed. ')
      END IF
#endif
!    
   END DO
   
   
   IF (l%thereare_stoch.and.(.not.l%chosenyesornostochastic)) THEN
      CALL stoponerror (l%layoutnumber, l%size, '!STOCH found in .nfde. Specify either -stoch or -nostoch')
   END IF
#ifndef CompileWithSlantedWires
   IF (l%hay_slanted_wires) THEN
      CALL stoponerror (l%layoutnumber, l%size, 'slanted wires without slanted support. Recompile ()')
   END IF
#endif   
   IF (l%hay_slanted_wires .AND. ((trim(adjustl(l%wiresflavor))/='slanted').AND.(trim(adjustl(l%wiresflavor))/='semistructured'))) THEN
      CALL stoponerror (l%layoutnumber, l%size, 'slanted wires require -l%wiresflavor Slanted/semistructured')
   endif

   
!punietero error abrezanjas y no l%resume conformal !niapa 121020
          ThereArethinslots=.FALSE.
          do jmed=1,sgg%NumMedia
             if (sgg%Med(jmed)%Is%ThinSlot) ThereArethinslots=.true.
          end do
         if (l%resume.and.l%run_with_abrezanjas.and.ThereArethinslots) then   
             CALL stoponerror (l%layoutnumber, l%size, 'l%resume -r currently unsupported by conformal solver',.true.); statuse=-1; !return
         end if
!fin niapa  
   !
!!!SOME FINAL REPORTING

   if (l%layoutnumber==0) then
        WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
        CALL print11 (l%layoutnumber, dubuf)
        CALL print11 (l%layoutnumber, 'Solver launched with options:')
        write(dubuf,*) l%mibc          
        CALL print11 (l%layoutnumber, '---> l%mibc    solver for NIBC multilayer: '//trim(adjustl(dubuf)))
        write(dubuf,*) l%ade         
        CALL print11 (l%layoutnumber, '---> l%ade     solver for ADC multilayer: '//trim(adjustl(dubuf)))
        Write(dubuf,*) l%sgbc    
        CALL print11 (l%layoutnumber, '---> sgbc    solver for multilayer: '//trim(adjustl(dubuf)))
        if (l%sgbc) then
            write(dubuf,*) l%sgbcDispersive      
            CALL print11 (l%layoutnumber, '---> sgbc DISPERSIVE solver for multilayer: '//trim(adjustl(dubuf)))
            write(dubuf,*) l%sgbccrank     
            CALL print11 (l%layoutnumber, '---> sgbc Crank-Nicolson solver for multilayer: '//trim(adjustl(dubuf)))
            write(dubuf,*) l%sgbcdepth
            CALL print11 (l%layoutnumber, '---> sgbc Depth: '//trim(adjustl(dubuf)))
            write(dubuf,*) l%sgbcfreq
            CALL print11 (l%layoutnumber, '---> sgbc Freq: '//trim(adjustl(dubuf)))
            write(dubuf,*) l%sgbcresol
            CALL print11 (l%layoutnumber, '---> sgbc Resol: '//trim(adjustl(dubuf)))
        endif
        write(dubuf,*) l%skindepthpre
        CALL print11 (l%layoutnumber, '---> l%skindepthpre preprocessing for multilayer: '//trim(adjustl(dubuf)))
        write(dubuf,*) l%flag_conf_sgg
        CALL print11 (l%layoutnumber, '---> Conformal file external: '//trim(adjustl(dubuf)))
        write(dubuf,*) l%input_conformal_flag      
        CALL print11 (l%layoutnumber, '---> Conformal solver: '//trim(adjustl(dubuf)))
        write(dubuf,*) l%run_with_abrezanjas
        CALL print11 (l%layoutnumber, '---> Conformal thin-gap solver: '//trim(adjustl(dubuf)))
        write(dubuf,*) l%run_with_dmma
        CALL print11 (l%layoutnumber, '---> DMMA thin-gap solver: '//trim(adjustl(dubuf)))
        write(dubuf,'(a)') l%wiresflavor
        CALL print11 (l%layoutnumber, '---> Wire model: '//trim(adjustl(dubuf)))
        write(dubuf,'(a)') l%inductance_model
        CALL print11 (l%layoutnumber, '---> Inductance model: '//trim(adjustl(dubuf)))
        if (trim(adjustl(l%wiresflavor))=='berenger') then
            write(dubuf,*) l%mindistwires
            CALL print11 (l%layoutnumber, '---> Berenger minimum distance between wires: '//trim(adjustl(dubuf)))
            write(dubuf,*) l%mtlnberenger
            CALL print11 (l%layoutnumber, '---> Berenger -l%mtlnberenger MTLN switch: '//trim(adjustl(dubuf)))
        endif
        if (trim(adjustl(l%wiresflavor))=='holland') then
            write(dubuf,*) l%stableradholland                 
            CALL print11 (l%layoutnumber, '---> Holland -l%stableradholland automatic correction switch: '//trim(adjustl(dubuf)))
        endif
        write(dubuf,*) l%TAPARRABOS                
        CALL print11 (l%layoutnumber, '---> Thin-wire double-tails removed: '//trim(adjustl(dubuf)))
        write(dubuf,*) l%fieldtotl                
        CALL print11 (l%layoutnumber, '---> Thin-wire -l%fieldtotl experimental switch: '//trim(adjustl(dubuf)))
        WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
        CALL print11 (l%layoutnumber, dubuf)
    endif

   IF (l%layoutnumber == 0) THEN
      call erasesignalingfiles(l%simu_devia)
   endif
   
    if (l%layoutnumber==0) then
        
        open(newunit=thefileno,FILE = trim(adjustl(l%nEntradaRoot))//'_tag_paraviewfilters.txt')
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
            do i=1,tagtype%numertags
                write(thefileno,'(a,i9,a,i9,a)') '# (',i*64,' , ',i*64+63,') '//trim(adjustl(tagtype%tag(i))) !los shifteo 6 bits y les sumo 2**campo ! idea de los 3 bits de 151020
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
            write(thefileno,'(a)') '# (  3.5 ,  3.5 ) '//trim(adjustl('sgbc/l%mibc Isotropic/anisotropic Multiport                               (Line)'))
            write(thefileno,'(a)') '# (  300 ,  399 ) '//trim(adjustl('sgbc/l%mibc Isotropic/anisotropic Multiport (+indexmedium)                (Surface)'))
            write(thefileno,'(a)') '# (  4.5 ,  4.5 ) '//trim(adjustl('Thin slot                                                               (Line)'))
            write(thefileno,'(a)') '# (  5.0 ,  5.0 ) '//trim(adjustl('Already_YEEadvanced_byconformal                                         (Surface)'))
            write(thefileno,'(a)') '# (  5.5 ,  5.5 ) '//trim(adjustl('Already_YEEadvanced_byconformal                                         (Line)'))
            write(thefileno,'(a)') '# (  6.0 ,  6.0 ) '//trim(adjustl('Split_and_useless                                                       (Surface)'))
            write(thefileno,'(a)') '# (  6.5 ,  6.5 ) '//trim(adjustl('Split_and_useless                                                       (Line)'))
            write(thefileno,'(a)') '# (  7.0 ,  7.0 ) '//trim(adjustl('Edge Not colliding thin wires                                           (Line)'))
            write(thefileno,'(a)') '# (  8.0 ,  8.0 ) '//trim(adjustl('Thin wire segments colliding with structure                             (Line)'))
            write(thefileno,'(a)') '# (  8.5 ,  8.5 ) '//trim(adjustl('Soft/Hard Nodal CURRENT/FIELD ELECTRIC DENSITY SOURCE                   (Line)'))
            write(thefileno,'(a)') '# (  9.0 ,  9.0 ) '//trim(adjustl('Soft/Hard Nodal CURRENT/FIELD MAGNETIC DENSITY SOURCE                   (Line)'))
            write(thefileno,'(a)') '# (   10 ,   11 ) '//trim(adjustl('LeftEnd/RightEnd/Ending  segment                                                 (Wire)'))
            write(thefileno,'(a)') '# (   20 ,   20 ) '//trim(adjustl('Intermediate segment +number_holland_parallel or +number_berenger       (Wire) '))
            write(thefileno,'(a)') '# (  400 ,  499 ) '//trim(adjustl('Thin slot (+indexmedium)                                                (Surface)'))
            write(thefileno,'(a)') '# ( -0.5 , -0.5 ) '//trim(adjustl('Other types of media                                                    (Line)'))
            write(thefileno,'(a)') '# ( -1.0 , -1.0 ) '//trim(adjustl('Other types of media                                                    (Surface)'))
        close(thefileno)
    endif

   ! call each simulation   !ojo que los layoutnumbers empiezan en 0
   IF (l%finaltimestep /= 0) THEN
#ifdef CompileWithMPI
      !wait until everything comes out
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      finishedwithsuccess=.false.
      if ((l%finaltimestep >= 0).and.(.not.l%skindepthpre)) then
#ifdef CompileWithMTLN
         CALL launch_simulation (sgg,sggMtag,tag_numbers, sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,&
           SINPML_fullsize,fullsize,finishedwithsuccess,Eps0,Mu0,tagtype, &
!los del tipo l%             
           l%simu_devia,l%cfl, l%nEntradaRoot, l%finaltimestep, l%resume, l%saveall,l%makeholes, &
           l%connectendings, l%isolategroupgroups,l%stableradholland, l%flushsecondsFields,l%mtlnberenger, &
           l%flushsecondsData, l%layoutnumber, l%size, l%createmap, &
           l%inductance_model, l%inductance_order, l%wirethickness, l%maxCPUtime,time_desdelanzamiento, &
           l%nresumeable2, l%resume_fromold,l%groundwires,l%noSlantedcrecepelo,l%sgbc,l%sgbcDispersive,l%mibc,l%attfactorc,l%attfactorw,&
           l%alphamaxpar,l%alphaOrden,l%kappamaxpar,l%mur_second,l%MurAfterPML,l%MEDIOEXTRA,&
           l%singlefilewrite,maxSourceValue,l%NOcompomur,l%ade, &
           l%conformalskin,l%strictOLD,l%TAPARRABOS,l%wiresflavor,l%mindistwires,l%facesNF2FF,l%NF2FFDecim,l%vtkindex,&
           l%createh5bin,l%wirecrank, &
           l%opcionestotales,l%sgbcfreq,l%sgbcresol,l%sgbccrank,l%sgbcdepth,l%fatalerror,l%fieldtotl,l%permitscaling, &
           l%EpsMuTimeScale_input_parameters, &
           l%stochastic,l%mpidir,l%verbose,l%precision,l%hopf,l%ficherohopf,l%niapapostprocess,l%planewavecorr, &
           l%dontwritevtk,l%experimentalVideal,l%forceresampled,l%factorradius,l%factordelta,l%noconformalmapvtk, &
           mtln_parsed, l%use_mtln_wires)
#else
            CALL launch_simulation (sgg,sggMtag,tag_numbers,sggMiNo, sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,&
            SINPML_fullsize,fullsize,finishedwithsuccess,Eps0,Mu0,tagtype, &
         !los del tipo l%             
            l%simu_devia,l%cfl, l%nEntradaRoot, l%finaltimestep, l%resume, l%saveall,l%makeholes, &
            l%connectendings, l%isolategroupgroups,l%stableradholland, l%flushsecondsFields,l%mtlnberenger, &
            l%flushsecondsData, l%layoutnumber, l%size, l%createmap, &
            l%inductance_model, l%inductance_order, l%wirethickness, l%maxCPUtime,time_desdelanzamiento, &
            l%nresumeable2, l%resume_fromold,l%groundwires,l%noSlantedcrecepelo,l%sgbc,l%sgbcDispersive,l%mibc,l%attfactorc,l%attfactorw,&
            l%alphamaxpar,l%alphaOrden,l%kappamaxpar,l%mur_second,l%MurAfterPML,l%MEDIOEXTRA,&
            l%singlefilewrite,maxSourceValue,l%NOcompomur,l%ade, &
            l%conformalskin,l%strictOLD,l%TAPARRABOS,l%wiresflavor,l%mindistwires,l%facesNF2FF,l%NF2FFDecim,l%vtkindex,&
            l%createh5bin,l%wirecrank, &
            l%opcionestotales,l%sgbcfreq,l%sgbcresol,l%sgbccrank,l%sgbcdepth,l%fatalerror,l%fieldtotl,l%permitscaling, &
            l%EpsMuTimeScale_input_parameters, &
            l%stochastic,l%mpidir,l%verbose,l%precision,l%hopf,l%ficherohopf,l%niapapostprocess,l%planewavecorr, &
            l%dontwritevtk,l%experimentalVideal,l%forceresampled,l%factorradius,l%factordelta,l%noconformalmapvtk,l%use_mtln_wires)
#endif
         deallocate (sggMiEx, sggMiEy, sggMiEz,sggMiHx, sggMiHy, sggMiHz,sggMiNo,sggMtag)
      else
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
      CALL get_secnds (l%time_out2)
      IF (l%layoutnumber == 0) THEN
          call print_credits(l)
          WRITE (dubuf,*) 'BEGUN '//trim (adjustl(l%nEntradaRoot)),' at ', time_comienzo%fecha(7:8), &
          & '/', time_comienzo%fecha(5:6), '/', time_comienzo%fecha(1:4),' , ',  &
          & time_comienzo%hora(1:2), ':', time_comienzo%hora(3:4)
          CALL print11 (l%layoutnumber, dubuf)
          WRITE (dubuf,*) 'ENDED '//trim (adjustl(l%nEntradaRoot)),' at ', l%time_out2%fecha(7:8), &
          & '/', l%time_out2%fecha(5:6), '/', l%time_out2%fecha(1:4),' , ',  &
          & l%time_out2%hora(1:2), ':', l%time_out2%hora(3:4)
          CALL print11 (l%layoutnumber, dubuf)
          WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
          CALL print11 (l%layoutnumber, dubuf)
          CALL print11 (l%layoutnumber, dubuf)
       ENDIF
         !!!!!!!        CALL CLOSEdxfFILE(l%layoutnumber,l%size)
         CALL CLOSEWARNINGFILE(l%layoutnumber,l%size,dummylog,l%stochastic,l%simu_devia) !aqui ya no se tiene en cuenta el l%fatalerror
#ifdef CompileWithMPI
         !wait until everything comes out
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
#ifdef CompileWithMPI
         CALL MPI_FINALIZE (l%ierr)
#endif
         stop
      endif
   END IF
   !
#ifdef CompileWithMPI
   !wait until everything comes out
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   !
   IF (l%layoutnumber == 0) THEN
      if (l%run) then
         OPEN (38, file='running')
         WRITE (38, '(a)') '!END'
         CLOSE (38,status='delete')
      endif
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
      WRITE (dubuf,*) 'DONE :  ', trim (adjustl(l%nEntradaRoot)), ' UNTIL n=', l%finaltimestep
      CALL print11 (l%layoutnumber, dubuf)
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
      call erasesignalingfiles(l%simu_devia)

   END IF

#ifdef CompileWithMPI
   !wait until everything comes out
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   !
   IF (l%deleteintermediates) THEN
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
      WRITE (dubuf,*) 'Attempting to delete all intermediate data files'
      CALL print11 (l%layoutnumber, dubuf)
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
      INQUIRE (file=trim(adjustl(l%nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt', EXIST=existe)
      IF (existe) THEN
         OPEN (19, file=trim(adjustl(l%nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt')
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
         IF (l%layoutnumber == 0) THEN
            OPEN (33, file=trim(adjustl(l%nEntradaRoot))//'_Outputlists.dat')
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
   if(l%input_conformal_flag)then
      call conf_sMesh%delete
      call conf_timeSteps%delete;
      call delete_conf_tools();
   end if
#endif
   !**************************************************************************************************
   !**************************************************************************************************
   !**************************************************************************************************

#ifdef CompileWithMPI
   call MPI_Barrier(SUBCOMM_MPI,l%ierr)
#endif
   CALL get_secnds (l%time_out2)
   IF (l%layoutnumber == 0) THEN
      call print_credits(l)
      WRITE (dubuf,*) 'BEGUN '//trim (adjustl(l%nEntradaRoot)),' at ', time_comienzo%fecha(7:8), &
      & '/', time_comienzo%fecha(5:6), '/', time_comienzo%fecha(1:4),' , ',  &
      & time_comienzo%hora(1:2), ':', time_comienzo%hora(3:4)
      CALL print11 (l%layoutnumber, dubuf)
      WRITE (dubuf,*) 'ENDED '//trim (adjustl(l%nEntradaRoot)),' at ', l%time_out2%fecha(7:8), &
      & '/', l%time_out2%fecha(5:6), '/', l%time_out2%fecha(1:4),' , ',  &
      & l%time_out2%hora(1:2), ':', l%time_out2%hora(3:4)
      CALL print11 (l%layoutnumber, dubuf)
      WRITE (dubuf,*) SEPARADOR // SEPARADOR // SEPARADOR
      CALL print11 (l%layoutnumber, dubuf)
      CALL print11 (l%layoutnumber, dubuf)
   ENDIF
   INQUIRE (file='relaunch', EXIST=l%relaunching)
#ifdef CompileWithMPI
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   ! Error reading check

#ifdef keeppause
   if (l%fatalerror) then
      fatalerror_aux=.true.
#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
      call MPI_AllReduce(fatalerror_aux, l%fatalerror, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, l%ierr)
#else
      l%fatalerror = fatalerror_aux
#endif
     if (l%fatalerror) l%relaunching=.true.
#ifdef CompileWithMPI
     CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
  endif
#endif

   IF (l%relaunching.and.(.not.finishedwithsuccess)) THEN
      IF (l%layoutnumber == 0) THEN
         CALL print11 (l%layoutnumber, SEPARADOR//SEPARADOR)
         CALL print11 (l%layoutnumber, 'Not finishing solicited either manually or by an error condition. Edit of create launch file and remove pause file ')
         CALL print11 (l%layoutnumber, SEPARADOR//SEPARADOR)
         OPEN (9, file='pause', FORM='formatted')
         write (9, '(a)') ' '
         CLOSE (9)
         OPEN (9, file='relaunch', FORM='formatted')
         write (9, '(a)') ' '
         CLOSE (9,status='delete')
      endif
      !!!!!
#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      IF (l%layoutnumber == 0) THEN
         CALL CloseReportingFiles
      endif
      GO TO 652
   END IF
!si ha acabado con exito sal borrando signal files
   IF (finishedwithsuccess) THEN
      IF (l%layoutnumber == 0) THEN
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
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif

   IF (l%layoutnumber == 0) THEN
      CALL CloseReportingFiles
   endif
   !**************************************************************************************************

#ifdef CompileWithMPI
   CALL MPI_FINALIZE (l%ierr)
#endif
   STOP
   !

contains
!END PROGRAM SEMBA_FDTD_launcher
!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!


#ifdef CompilePrivateVersion 
subroutine cargaNFDE(local_nfde,local_parser)
   CHARACTER (LEN=BUFSIZE) :: local_nfde
   TYPE (Parseador), POINTER :: local_parser
   INTEGER (KIND=8) :: numero,i8,troncho,longitud
   integer (kind=4) :: mpi_t_linea_t,longitud4
   IF (l%existeNFDE) THEN
       WRITE (dubuf,*) 'INIT Reading file '//trim (adjustl(whoami))//' ', trim (adjustl(local_nfde))
       CALL print11 (l%layoutnumber, dubuf)
!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
      if (l%layoutnumber==0) then
           NFDE_FILE => cargar_NFDE_FILE (local_nfde)
      !!!ya se allocatea dentro
      else
           ALLOCATE (NFDE_FILE)
      endif
      !
      write(dubuf,*) '[OK]';  call print11(l%layoutnumber,dubuf)
      !--->
      WRITE (dubuf,*) 'INIT Sharing file through MPI'; CALL print11 (l%layoutnumber, dubuf)
      !
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
      !
      numero=NFDE_FILE%numero
      call MPI_BCAST(numero, 1_4, MPI_INTEGER8, 0_4, SUBCOMM_MPI, l%ierr)      
      if (l%layoutnumber/=0) then
          NFDE_FILE%targ = 1
          NFDE_FILE%numero=numero
          ALLOCATE (NFDE_FILE%lineas(NFDE_FILE%numero))
      endif
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
      !CREAMOS EL DERIVED TYPE y lo enviamos !para evitar el error de Marconi asociado a PSM2_MQ_RECVREQS_MAX 100617

      CALL build_derived_t_linea(mpi_t_linea_t)

      !problema del limite de mandar mas de 2^29 bytes con MPI!!!  Los soluciono partiendo en maxmpibytes (2^27) (algo menos por prudencia)! 040716
      troncho=ceiling(maxmpibytes*1.0_8/(BUFSIZE*1.0_8+8.0_8),8)
     !!! print *,'numero,troncho ',numero,troncho
      do i8=1,numero,troncho
          longitud=min(troncho,numero-i8+1)
          CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
          if ((longitud>huge(1_4)).or.(longitud>maxmpibytes)) then
              print *,'Stop. Buggy error: MPI longitud greater that greatest integer*4'
              stop
          else
              longitud4=int(longitud,4)
          endif
          call MPI_BCAST(NFDE_FILE%lineas(i8),longitud4,mpi_t_linea_t,0_4,SUBCOMM_MPI,l%ierr)    
          CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
       !!!  if (l%layoutnumber==1) print *,'l%layoutnumber-->',l%layoutnumber, i8,i8+longitud-1 
       !!!  if (l%layoutnumber==1) print *,NFDE_FILE%lineas(i8)%len,' ',trim(adjustl(NFDE_FILE%lineas(i8)%dato)) 
       !!!  if (l%layoutnumber==1) print *,NFDE_FILE%lineas(i8+longitud-1)%len,' ',trim(adjustl(NFDE_FILE%lineas(i8+longitud-1)%dato))
          CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
         !      do i=1,numero
    !          call MPI_BCAST(NFDE_FILE%lineas(i)%len, 1_4, MPI_INTEGER4, 0_4, SUBCOMM_MPI, l%ierr)
    !          call MPI_BCAST(NFDE_FILE%lineas(i)%dato, BUFSIZE, MPI_CHARACTER, 0_4, SUBCOMM_MPI, l%ierr)
    !          CALL MPI_Barrier (SUBCOMM_MPI, l%ierr) !para evitar el error de Marconi asociado a PSM2_MQ_RECVREQS_MAX 100617
         !      end do
      end do
      !solo para debugeo
            !!!open(6729,file='comprob_'//trim(adjustl(dubuf))//'.nfde',form='formatted')
            !!!write(6729,'(2i12)') NFDE_FILE%numero,NFDE_FILE%targ
            !!!do i=1,numero
            !!!   write(6729,'(i6,a)') NFDE_FILE%lineas(i)%len,trim(adjustl(NFDE_FILE%lineas(i)%dato))
            !!!end do
            !!!close (6729)
      !!!!!!
#else
    NFDE_FILE => cargar_NFDE_FILE (local_nfde)
#endif
      write(dubuf,*) '[OK]';  call print11(l%layoutnumber,dubuf)
      !--->
   END IF    
   NFDE_FILE%mpidir=l%mpidir
!!!!!!!!!!!!!!!!!!!
   WRITE (dubuf,*) 'INIT interpreting geometrical data from ', trim (adjustl(local_nfde))
   CALL print11 (l%layoutnumber, dubuf)
!!!!!!!!!!
   if(newrotate) then
       verdadero_mpidir=NFDE_FILE%mpidir
       NFDE_FILE%mpidir=3     !no lo rota el parseador antiguo
   endif
   local_parser => newparser (NFDE_FILE)         
#ifdef CompileWithMPI            
   CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   if(newrotate) then      
       NFDE_FILE%mpidir=verdadero_mpidir   !restorealo
       call nfde_rotate (local_parser,NFDE_FILE%mpidir)   !lo rota el parseador nuevo  
#ifdef CompileWithMPI            
       CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   endif
   l%thereare_stoch=NFDE_FILE%thereare_stoch
   l%mpidir=NFDE_FILE%mpidir !bug 100419
!!!!!!!!!!!                             
  ! write(dubuf,*) '[OK]';  call print11(l%layoutnumber,dubuf)
   write(dubuf,*) '[OK] '//trim(adjustl(whoami))//' newparser (NFDE_FILE)';  call print11(l%layoutnumber,dubuf)       
#ifdef CompileWithMPI            
       CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
   return

end subroutine cargaNFDE
#endif

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

!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine NFDE2sgg     
!!!!!!!!!      
      real (kind=rkind) :: dt,finaldt
      logical fatalerror
      ! parser now holds all the .nfde info
      !first read the limits
#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      CALL read_limits_nogeom (l%layoutnumber,l%size, sgg, fullsize, SINPML_fullsize, parser,l%MurAfterPML,l%mur_exist)
    
      dtantesdecorregir=sgg%dt
      !!!!!corrige el delta de t si es necesario !sgg15 310715 bug distintos sgg%dt !!!!!!!!!!

      dxmin=minval(sgg%DX)
      dymin=minval(sgg%DY)
      dzmin=minval(sgg%DZ)
      !!!
      dtlay=(1.0_RKIND/(cluz*sqrt(((1.0_RKIND / dxmin)**2.0_RKIND )+((1.0_RKIND / dymin)**2.0_RKIND )+((1.0_RKIND / dzmin)**2.0_RKIND ))))
      dt=dtlay
#ifdef CompileWithMPI
      call MPIupdateMin(dtlay,dt)
#endif

      !!!write(dubuf,*) SEPARADOR//separador//separador
      !!!call print11(l%layoutnumber,dubuf)
      !!!write(dubuf,*) '--->dt,dxmin,dymin,dzmin,sgg%dt  ',dt,dxmin,dymin,dzmin,sgg%dt
      !!!call print11(l%layoutnumber,dubuf)
      !!!write(dubuf,*) SEPARADOR//separador//separador
      !!!call print11(l%layoutnumber,dubuf)

      if (l%forcecfl) then
         sgg%dt=dt*l%cfl
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(l%layoutnumber,dubuf)
         write(dubuf,*) 'Correcting sgg%dt with -l%cfl switch. New time step: ',sgg%dt
         call print11(l%layoutnumber,dubuf)
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(l%layoutnumber,dubuf)
      else
          if (sgg%dt > dt*heurCFL) then
             write(dubuf,*) SEPARADOR//separador//separador
             call print11(l%layoutnumber,dubuf)
             write(dubuf,*) 'Automatically correcting dt for stability reasons: '
             call print11(l%layoutnumber,dubuf)
             write(dubuf,*) 'Original dt: ',sgg%dt
             call print11(l%layoutnumber,dubuf)
             sgg%dt=dt*heurCFL
             write(dubuf,*) 'New dt: ',sgg%dt
             call print11(l%layoutnumber,dubuf)
             write(dubuf,*) SEPARADOR//separador//separador
             call print11(l%layoutnumber,dubuf)
          endif
      endif
      !!!!!!!!!!!!No es preciso re-sincronizar pero lo hago !!!!!!!!!!!!!!!!!!!!!!!!!!
      finaldt=sgg%dt
#ifdef CompileWithMPI
      call MPIupdateMin(real(sgg%dt,RKIND),finaldt)
#endif
      !!!!!!!!!!!!!!
      l%cfl=sgg%dt/dtlay
      write(dubuf,*) SEPARADOR//separador//separador
      call print11(l%layoutnumber,dubuf)
      write(dubuf,*) 'CFLN= ',l%cfl
      call print11(l%layoutnumber,dubuf)
      write(dubuf,*) SEPARADOR//separador//separador
      call print11(l%layoutnumber,dubuf)

      write(dubuf,*) SEPARADOR//separador//separador
      call print11(l%layoutnumber,dubuf)
      write(dubuf,*) 'Deltat= ',sgg%dt
      if (l%layoutnumber==0) call print11(l%layoutnumber,dubuf)
#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      write(dubuf,*) SEPARADOR//separador//separador
      call print11(l%layoutnumber,dubuf)
      if (l%mur_exist.and.l%mur_first) then
         l%mur_second=.false.
      else
         l%mur_second=.false. !arreglar cuando se arregle el bug de las mur second
         l%mur_first=.true. !arreglar cuando se arregle el bug de las mur second
      endif
#ifdef CompileWithMPI
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      !LATER OVERRRIDEN BY MPI
      !ALLOCATED ONE MORE TO KEEP PMC INFO FOR THE HX,HY,HZ FIELDS
      sgg%Alloc(1:6)%XI = fullsize(1:6)%XI - 1
      sgg%Alloc(1:6)%XE = fullsize(1:6)%XE + 1
      sgg%Alloc(1:6)%YI = fullsize(1:6)%YI - 1
      sgg%Alloc(1:6)%YE = fullsize(1:6)%YE + 1
      !REDUCE THE SWEEP AREA BY 1
      sgg%Sweep(1:6)%XI = fullsize(1:6)%XI
      sgg%Sweep(1:6)%XE = fullsize(1:6)%XE
      sgg%Sweep(1:6)%YI = fullsize(1:6)%YI
      sgg%Sweep(1:6)%YE = fullsize(1:6)%YE
      !
      IF (l%size == 1) THEN
         sgg%Alloc(1:6)%ZI = fullsize(1:6)%ZI - 1
         sgg%Alloc(1:6)%ZE = fullsize(1:6)%ZE + 1
         !REDUCE THE SWEEP AREA BY 1
         sgg%Sweep(1:6)%ZI = fullsize(1:6)%ZI
         sgg%Sweep(1:6)%ZE = fullsize(1:6)%ZE
         !!incluido aqui pq se precisa para clip 16/07/15
         DO field = iEx, iHz
            sgg%SINPMLSweep(field)%XI = Max (SINPML_fullsize(field)%XI, sgg%Sweep(field)%XI)
            sgg%SINPMLSweep(field)%XE = Min (SINPML_fullsize(field)%XE, sgg%Sweep(field)%XE)
            sgg%SINPMLSweep(field)%YI = Max (SINPML_fullsize(field)%YI, sgg%Sweep(field)%YI)
            sgg%SINPMLSweep(field)%YE = Min (SINPML_fullsize(field)%YE, sgg%Sweep(field)%YE)
            sgg%SINPMLSweep(field)%ZI = Max (SINPML_fullsize(field)%ZI, sgg%Sweep(field)%ZI)
            sgg%SINPMLSweep(field)%ZE = Min (SINPML_fullsize(field)%ZE, sgg%Sweep(field)%ZE)
         END DO
         !!fin 16/07/15
         WRITE (dubuf,*) 'INIT NFDE --------> GEOM'
         CALL print11 (l%layoutnumber, dubuf)
         CALL read_geomData (sgg,sggMtag,tag_numbers, sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, l%fichin, l%layoutnumber, l%size, SINPML_fullsize, fullsize, parser, &
         l%groundwires,l%attfactorc,l%mibc,l%sgbc,l%sgbcDispersive,l%MEDIOEXTRA,maxSourceValue,l%skindepthpre,l%createmapvtk,l%input_conformal_flag,l%CLIPREGION,l%boundwireradius,l%maxwireradius,l%updateshared,l%run_with_dmma, &
         eps0,mu0,.false.,l%hay_slanted_wires,l%verbose,l%ignoresamplingerrors,tagtype,l%wiresflavor)
#ifdef CompileWithMTLN
         if (trim(adjustl(l%extension))=='.json')  mtln_parsed = parser%mtln
         ! if (trim(adjustl(l%extension))=='.json')  mtln_solver = mtlnCtor(parser%mtln)   
#endif
         WRITE (dubuf,*) '[OK] ENDED NFDE --------> GEOM'
         CALL print11 (l%layoutnumber, dubuf)
         !writing
         slices = '!SLICES'
         WRITE (buff, '(i7)') sgg%Sweep(iHz)%ZE - sgg%Sweep(iHz)%ZI
         slices = trim (adjustl(slices)) // '_' // trim (adjustl(buff))
         IF (l%resume .AND. (slices /= l%slicesoriginales)) THEN
            buff='Different resumed/original MPI slices: '//trim(adjustl(slices))//' '//&
            & trim(adjustl(l%slicesoriginales))
            CALL stoponerror (l%layoutnumber, l%size, buff)
         END IF
         CALL print11 (l%layoutnumber, trim(adjustl(slices)))
         !end writing
         WRITE (buff, '(a,i7,a,i7)') '_________Spanning from z=', sgg%Sweep(iHz)%ZI, ' to z=', sgg%Sweep(iHz)%ZE
         CALL print11 (l%layoutnumber, trim(adjustl(buff)))
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#ifdef CompileWithStochastic
         if (l%stochastic) then
            buff='l%stochastic uncompatible with MPI l%size smaller than 2'
            CALL stoponerror (l%layoutnumber, l%size, buff)
         endif
#endif
#endif
      ELSE !del l%size==1       
#ifdef CompileWithMPI
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#ifdef CompileWithStochastic
         if (l%stochastic) then
            call HalvesStochasticMPI(l%layoutnumber,l%size,l%simu_devia)
         endif
#endif
                   
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)   
!!!ahora divide el espacio computacional
         CALL MPIdivide (sgg, fullsize, SINPML_fullsize, l%layoutnumber, l%size, l%forcing, l%forced, l%slicesoriginales, l%resume,l%fatalerror)
         !
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)   
         if (l%fatalerror) then
!intenta recuperarte
            return
         endif
     
         ! if the layout is pure PML then take at least a line of non PML to build the PML data insider read_geomDAta
         ! Uses extra memory but later matrix sggm is deallocated in favor of smaller sggMIEX, etc
         DO field = iEx, iHz
            tempalloc(field)%ZE = sgg%Alloc(field)%ZE
            tempalloc(field)%ZI = sgg%Alloc(field)%ZI
            sgg%Alloc(field)%ZE = Max (sgg%Alloc(field)%ZE, SINPML_fullsize(field)%ZI+1)
            sgg%Alloc(field)%ZI = Min (sgg%Alloc(field)%ZI, SINPML_fullsize(field)%ZE-1)
         END DO
         !   
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)  
         !!incluido aqui pq se precisa para clip 16/07/15
         DO field = iEx, iHz
            sgg%SINPMLSweep(field)%XI = Max (SINPML_fullsize(field)%XI, sgg%Sweep(field)%XI)
            sgg%SINPMLSweep(field)%XE = Min (SINPML_fullsize(field)%XE, sgg%Sweep(field)%XE)
            sgg%SINPMLSweep(field)%YI = Max (SINPML_fullsize(field)%YI, sgg%Sweep(field)%YI)
            sgg%SINPMLSweep(field)%YE = Min (SINPML_fullsize(field)%YE, sgg%Sweep(field)%YE)
            sgg%SINPMLSweep(field)%ZI = Max (SINPML_fullsize(field)%ZI, sgg%Sweep(field)%ZI)
            sgg%SINPMLSweep(field)%ZE = Min (SINPML_fullsize(field)%ZE, sgg%Sweep(field)%ZE)
         END DO
         !!fin 16/07/15
         WRITE (dubuf,*) 'INIT NFDE --------> GEOM'
         CALL print11 (l%layoutnumber, dubuf)           

         CALL read_geomData (sgg,sggMtag,tag_numbers, sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, l%fichin, l%layoutnumber, l%size, SINPML_fullsize, fullsize, parser, &
         l%groundwires,l%attfactorc,l%mibc,l%sgbc,l%sgbcDispersive,l%MEDIOEXTRA,maxSourceValue,l%skindepthpre,l%createmapvtk,l%input_conformal_flag,l%CLIPREGION,l%boundwireradius,l%maxwireradius,l%updateshared,l%run_with_dmma, &
         eps0,mu0,l%simu_devia,l%hay_slanted_wires,l%verbose,l%ignoresamplingerrors,tagtype,l%wiresflavor)


#ifdef CompileWithMPI
         !wait until everything comes out
         CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
         WRITE (dubuf,*) '[OK] ENDED NFDE --------> GEOM'
         CALL print11 (l%layoutnumber, dubuf)
         !restore back the indexes
         DO field = iEx, iHz
            sgg%Alloc(field)%ZE = tempalloc(field)%ZE
            sgg%Alloc(field)%ZI = tempalloc(field)%ZI
         END DO
#endif
         CONTINUE
      END IF !del l%size==1
      !
#ifdef CompileWithMPI
      !wait until everything comes out
      CALL MPI_Barrier (SUBCOMM_MPI, l%ierr)
#endif
      !!!!!!!!!!!!!lo dejo aqui debajo tambien aunque ya se ha calculado antes para lo del clipping
      DO field = iEx, iHz
         sgg%SINPMLSweep(field)%XI = Max (SINPML_fullsize(field)%XI, sgg%Sweep(field)%XI)
         sgg%SINPMLSweep(field)%XE = Min (SINPML_fullsize(field)%XE, sgg%Sweep(field)%XE)
         sgg%SINPMLSweep(field)%YI = Max (SINPML_fullsize(field)%YI, sgg%Sweep(field)%YI)
         sgg%SINPMLSweep(field)%YE = Min (SINPML_fullsize(field)%YE, sgg%Sweep(field)%YE)
         sgg%SINPMLSweep(field)%ZI = Max (SINPML_fullsize(field)%ZI, sgg%Sweep(field)%ZI)
         sgg%SINPMLSweep(field)%ZE = Min (SINPML_fullsize(field)%ZE, sgg%Sweep(field)%ZE)
      END DO
      return
   end subroutine
   !

   
  END PROGRAM SEMBA_FDTD_launcher
!
