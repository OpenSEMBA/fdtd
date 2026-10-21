module Report_m
   use FDETYPES_m

   implicit none
   private

   type :: time_t
      real(kind = 8) :: segundos
      character(LEN=BUFSIZE) :: hora
      character(LEN=BUFSIZE) :: fecha
   end type


   logical, save :: printea,calledStoponerrroonlyprint=.false.,warningfileIsOpen=.false.,verbose,file10isopen,file11isopen
   character(len=BUFSIZE), save :: warningFile = ' '
   character(len=BUFSIZE), save  :: whoami

   integer, save :: thefile !for mpi file management
   logical, save :: ignoreerrors
   !

   character(len=BUFSIZE), save :: myNInputRoot

   !!!logical, SAVE :: dxfFileIsOpen=.false.
   logical, save :: fatalerror=.false.
   integer, save :: CONTADORDEMENSAJES
   integer, save :: thefile2 !for mpi file management
   !!public StopOnError_OnlyPrint

   public StopOnError,CloseReportingFiles, &
   print11,Onprint,Offprint,file10isopen,file11isopen
   public WarnErrReport,INITWARNINGFILE,CLOSEWARNINGFILE,get_secnds,openfile_mpi,writefile_mpi, &
          closefile_mpi,reportmedia,erasesignalingfiles,openclosedelete,openclose
   public isFatalError, resetFatalError

   !part of the dxf
   !!!public dxfwrite,INITdxfFILE,CLOSEdxfFILE,writemmdxf,TRIMNULLCHAR
   public TRIMNULLCHAR,time_t




contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine OnPrint
      printea=.true.
      return
   end subroutine
   subroutine OffPrint
      printea=.false.
      return
   end subroutine
   !!!!!!!!!!!!!!!!!
   subroutine StopOnError(layoutnumber,num_procs,message,calledfrommain)
      character(len=BUFSIZE) :: lockFile
      logical , optional  :: calledfrommain
      character(len=*), intent(in) :: message
      integer(kind=4), optional  :: layoutnumber,num_procs
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif
      character(len=BUFSIZE) :: whoami

      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',num_procs,') '


      call print11(layoutnumber,trim(adjustl(whoami))//' ERROR: '//trim(adjustl(message)),.true.)

      !19/12/14 bug OLD1812. A stoponerror created by an embedded nodal source arrives here in MPI. The closewarn... does a barrier and prevents the process from dying.

      !the stoponerror calls must be reviewed and made more elegant. For now I abort crudely, commenting out without closing either warning or dxf (To do)

      !call CLOSEWARNINGFILE(layoutnumber,num_procs)
      !!!!call CLOSEdxfFILE(layoutnumber,num_procs)



      !!!!the return is only for the cluster
#ifdef keeppause
      if (present(calledfrommain)) then
         if (calledfrommain) then
            if (layoutnumber == 0) then
               open(38, FILE='pause')
               write(38, '(a)') '!END'
               close (38)
            end if
            call print11(layoutnumber,'Trying to relaunch. Correct error, create launch, and remove pause/warning '// &
                                      'file (or kill the process)',.true.)
!!!            call CloseReportingFiles !sgg 240817 the reporting files must not be closed
            return
         end if
      else
         if (layoutnumber == 0) then
            open(38, FILE='pause')
            write(38, '(a)') '!END'
            close (38)
         end if
         call print11(layoutnumber,'Stopping, but creating the signal file pause to prevent queuing losses!!! '// & '
                                   '(correct error and remove to continue)',.true.)
      end if
#else
      if (layoutnumber == 0) then
         lockFile='running'
          call openclosedelete(lockFile)
         !
         lockFile='pause'
          call openclosedelete(lockFile)
         !
         lockFile='relaunch'
          call openclosedelete(lockFile)
      end if
#endif

#ifdef CompileWithMPI
      call print11(layoutnumber,'Trying to kill all MPI processes (may fail!)...',.true.)
      call MPI_Abort(SUBCOMM_MPI, -1, ierr)
      call MPI_FINALIZE(ierr)
#endif
      call CloseReportingFiles
 
      stop 1

      return

   end subroutine




   subroutine CloseReportingFiles
      if (file10isopen) then
         close (10) !energy file
         file10isopen=.false.
      end if
      if (file11isopen) then
         close (11) !reporting file
         file11isopen=.false.
      end if
      return

   end subroutine



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---------------------------------------------------->
   subroutine get_secnds(time_out2)
      type(time_t), intent(out) :: time_out2
      integer(kind=8) :: clockCount, clockRate

      ! monotonic clock: callers only use differences of segundos, and wall-clock dates jump at month/DST changes
      call system_clock(clockCount, clockRate)
      time_out2%segundos = real(clockCount, 8) / real(clockRate, 8)
      call date_and_time(date=time_out2%fecha, time=time_out2%hora)
   end subroutine get_secnds





   subroutine INITWARNINGFILE(layoutnumber,num_procs,nInputRoot,verbosete,ignoreErrors1)
      character(len=*) :: nInputRoot
      integer(kind=4), intent(in) :: layoutnumber,num_procs
      !file management
      character(len=BUFSIZE) :: whoamishort
#ifdef CompileWithMPI
      integer(kind=MPI_OFFSET_KIND) disp
      integer(kind=4) :: ierr
#endif
      logical verbosete,ignoreerrors1       , itsopen2
      integer :: my_iostat
      character(len=BUFSIZE) :: lockFile
      verbose=verbosete

      ignoreerrors=ignoreerrors1

      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',num_procs,') '
      write(whoamishort,'(i5)') layoutnumber+1

      if (layoutnumber == 0) then          

      !!!inquire(unit=17, opened=itsopen2)
      !!!if (itsopen2) print *,'----------->17 open!!!'
        lockFile=trim(adjustl(nInputRoot))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
        call openclosedelete(lockFile)
      end if

      !!!#ifdef CompileWithMPI
      !!!if (SIZE/=0) then
      !!!    call MPI_Barrier(SUBCOMM_MPI,ierr)
      !!!    call MPI_FILE_open (SUBCOMM_MPI, trim(adjustl(nEntradaRoot))//'_tmpWarnings.txt', &
      !!!                           MPI_MODE_WRONLY + MPI_MODE_CREATE, &
      !!!                           MPI_INFO_NULL, thefile, ierr)
      !!!    disp = (layoutnumber+1) * BUFSIZE * maxmessages !I do not think more than 2000 messages per layout will occur
      !!!
      !!!    call MPI_FILE_SET_VIEW(thefile, disp, MPI_CHARACTER, &
      !!!                               MPI_CHARACTER, 'native', &
      !!!                               MPI_INFO_NULL, ierr)
      !!!ELSE
      !!!    open (17,file=trim(adjustl(nEntradaRoot))//'_tmpWarnings.txt',form='formatted')
      !!!end if
      !!!#else
      
      inquire(unit=17, opened=itsopen2)
      !!!if (itsopen2) print *,'----------->17 open!!!'
      lockFile=trim(adjustl(nInputRoot))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
      call opensolo(17,lockFile)
      !!!#endif

      warningfileIsOpen=.true.
      warningfile=nInputRoot
      fatalerror = .false.
      CONTADORDEMENSAJES=0

   end subroutine INITWARNINGFILE


   subroutine WarnErrReport(bufff,error)
      use iso_fortran_env, only : error_unit
      logical :: itsopen
      logical, optional :: error
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif

      character(len=*), intent(in) :: bufff
      character(len=BUFSIZE) :: buff2,buff3

      if (present(error)) then
         fatalerror = error .or. fatalerror
      end if

      !

      buff3=trim(adjustl(whoami))//' '//trim(adjustl(bufff))
      call trimnullchar(buff3)

      buff2=CHAR(13)//CHAR(10)//trim(adjustl(buff3))
      call trimnullchar(buff2)

      inquire(unit=17, opened=itsopen)
      if (itsopen) write (17,'(a)',err=154) trim(adjustl(buff3))
      write (error_unit,'(a)') trim(adjustl(buff3))

      goto 155
154   inquire(unit=17, opened=itsopen) 

      print *,itsopen,'- Cannot write into warning file the message: ',trim(adjustl(buff3))
155   return

   end subroutine WarnErrReport

   function isFatalError() result(res)
      logical :: res
      res = fatalerror
   end function

   subroutine resetFatalError()
      fatalerror = .false.
   end subroutine



   subroutine CLOSEWARNINGFILE(layoutnumber,num_procs,fatalerror_final,stoch_undivided,simu_devia)
      integer(kind=4), intent(in) :: layoutnumber,num_procs
      integer(kind=4) :: ierr,posic,i
      character(len=BUFSIZE) :: buf2
      character(len=BUFSIZE) :: dubuf
      logical :: fatalerror_final , lexis,stoch_undivided,simu_devia        , itsopen2
      character(LEN=BUFSIZE) :: whoamishort,whoami,chinstant
      integer :: my_iostat,file87
      character(len=BUFSIZE) :: lockFile

      if (.not.WarningFileIsOpen) return

      !!!#ifdef CompileWithMPI
      !!!if (SIZE/=0) then
      !!!    call MPI_FILE_close (thefile, ierr)
      !!!ELSE
      !!!    close (17)
      !!!end if
      !!!#else
      close (17)
      !!!#endif

#ifdef CompileWithMPI
      !wait until everything is closed
      call MPI_Barrier (SUBCOMM_MPI, ierr)
#endif

      !fixes the NULs
      if ((layoutnumber==0).or.((layoutnumber == num_procs/2).and.stoch_undivided)) then
         open (88,file=trim(adjustl(WarningFile))//'_Warnings.txt',form='formatted')
         posic=0
         do i=0,num_procs-1
            if (stoch_undivided) then
                write(whoamishort,'(i5)') i+1
            else
               if (simu_devia) then
                   write(whoamishort,'(i5)') num_procs+i+1
               else
                   write(whoamishort,'(i5)') i+1
               end if
            end if 
            inquire(file=trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt',exist=lexis)
            if (lexis) then         
      !!!inquire(unit=87, opened=itsopen2)
      !!!if (itsopen2) print *,'----------->87 open!!!'
               lockFile=trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
               call opensolo(87,lockFile)
               !
875            read(87,'(a)',end=876,err=876) buf2
               call trimnullchar(buf2)
               buf2=trim(adjustl(buf2))
               if ((buf2(1:1) /= ' ').and.(buf2(1:1) /=char(0))) then
                  write(88,'(a)') trim(adjustl(buf2))
                  posic=posic+1
               end if
               goto 875
               !
876            continue
               call closesolo(87)
#ifndef CorregirBugBorrado
!!               my_iostat=0
!!3467           if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',3467,'.',layoutnumber,&
!!!                         trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
!!               open (newunit=file87,file=trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt',&
!!                           err=3467,iostat=my_iostat,status='new',action='write')
               lockFile=trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
               call openclosedelete(lockFile)
#endif
            end if
         end do

         !
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(layoutnumber,dubuf)
         write(dubuf,*) 'Closing warning file. Number of messages: ',posic
         call print11(layoutnumber,dubuf)
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(layoutnumber,dubuf)

         close (88)
      end if

      warningfileIsOpen=.false.
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
      call MPI_AllReduce(fatalerror, fatalerror_final, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
#else
      fatalerror_final = fatalerror
#endif

      if (fatalerror_final.and.(ignoreErrors)) then
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(layoutnumber,dubuf)
         write(dubuf,*) 'There are ERRORS: The simulation will continue but Revise *Warnings file '
         call print11(layoutnumber,dubuf)
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(layoutnumber,dubuf)
      end if

      fatalerror_final=(fatalerror_final).and.(.not.ignoreErrors)

      return
   end subroutine

   subroutine trimnullchar(string)
      character(len=*), intent(inout) :: string
      integer :: i,longi
      integer :: ind
      longi=len(string)
      do i=1,bufsize
         ind = scan(string(i:longi), char(0))
         if (ind/=0) string(ind:ind) = ' '
      end do
      string=trim(adjustl(string))
   end subroutine trimnullchar



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine print11(layoutnumber,message,forceprint2)
      character(len=*), intent(in) :: message
      integer(kind=4), intent(in) :: layoutnumber
      logical :: soyImpresor,forceprint
      logical , optional  :: forceprint2

      forceprint=.false.
      if (present(forceprint2)) forceprint=forceprint2

      soyImpresor=((layoutnumber == 0).or.forceprint).and.printea
      if (message(1:1)=='&') then !respects spaces, the & is actually a space
#ifndef NoVerbose
        ! if (soyImpresor) write (*,'(a)') ' '//trim(message(2:))
         if (soyImpresor) write (*,'(a)') trim(message(2:))
#endif
        ! if (layoutnumber == 0)   write (11,'(a)',err=111) ' '//trim(message(2:))
         if (layoutnumber == 0)   write (11,'(a)',err=111) trim(message(2:))
      else !left-adjusts without respecting spaces
#ifndef NoVerbose
         if (soyImpresor) write (*,'(a)') trim(adjustl(message))
#endif
         if (layoutnumber == 0)    then
             write (11,'(a)',err=111) trim(adjustl(message))
           !!!  print *,'----11--> ', trim(adjustl(message))
         end if
         
      end if
      goto 112
111   continue 
      !sometimes 2 of those that fail in writing restarting fields try to write fort.11 at once. 
      !so it ignores it and continues
112   return

   end subroutine

   !part of the DXF



   !!!subroutine INITdxfFILE(layoutnumber,num_procs,nEntradaRoot)
   !!!character(len=*) :: nEntradaRoot
   !!!integer(kind=4) :: layoutnumber,num_procs
   !!!!file management
   !!!character(len=BUFSIZE) whoamishort
   !!!#ifdef CompileWithMPI
   !!!integer(kind=MPI_OFFSET_KIND) disp
   !!!integer(kind=4) :: ierr
   !!!#endif
   !!!write(whoamishort,'(i5)') layoutnumber+1
   !!!mynEntradaRoot=trim(adjustl(nentradaRoot))
   !!!dxfFileIsOpen=.true.
   !!!
   !!!if (layoutnumber == 0) then
   !!!    open (87,file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf',form='formatted')
   !!!    write(87,*) '!END'
   !!!    close (87,status='delete')
   !!!end if
   !!!
   !!!!!!#ifdef CompileWithMPI
   !!!!!!call MPI_Barrier(SUBCOMM_MPI,ierr)
   !!!!!!call MPI_FILE_open (SUBCOMM_MPI, trim(adjustl(mynEntradaRoot))//'.tmpdxf', &
   !!!!!!                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
   !!!!!!                       MPI_INFO_NULL, thefile2, ierr)
   !!!!!!disp = (layoutnumber+1) * dxflinesize * maxdxf !I do not think more than 2000 messages per layout will occur
   !!!!!!
   !!!!!!call MPI_FILE_SET_VIEW(thefile2, disp, MPI_CHARACTER, &
   !!!!!!                           MPI_CHARACTER, 'native', &
   !!!!!!                           MPI_INFO_NULL, ierr)
   !!!!!!#else
   !!!    open (97,file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf',form='formatted')
   !!!!!!#endif
   !!!
   !!!if (layoutnumber == 0) then
   !!!    write(dxfbuff,'(a)') '999'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)') 'SEMBA_FDTD'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '0'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'SECTION'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '2'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'HEADER'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '9'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '$ACADVER'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '1'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'AC1009'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '9'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '$SHADEDGE'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '70'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '1'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '0'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'ENDSEC'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '0'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'SECTION'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '2'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'ENTITIES'
   !!!    call DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '0'
   !!!    call DXFWRITE(DXFBUFF)
   !!!end if
   !!!
   !!!end subroutine INITdxfFILE
   !!!
   !!!
   !!!subroutine dxfwrite(bufff)
   !!!!
   !!!#ifdef CompileWithMPI
   !!!integer(kind=4) :: ierr
   !!!#endif
   !!!
   !!!character(len=dxflinesize) :: bufff,buff2
   !!!
   !!!buff2=CHAR(13)//CHAR(10)//bufff(1:dxflinesize-4)//CHAR(13)//CHAR(10)
   !!!!
   !!!!!#ifdef CompileWithMPI
   !!!!!        call MPI_FILE_WRITE(thefile2, buff2  , dxflinesize, MPI_CHARACTER, &
   !!!!!                                MPI_STATUS_IGNORE, ierr)
   !!!!!#else
   !!!        write (97,'(a)') trim(adjustl(bufff))
   !!!!!#endif
   !!!!
   !!!return
   !!!end subroutine dxfwrite
   !!!
   !!!
   !!!
   !!!subroutine CLOSEdxfFILE(layoutnumber,num_procs)
   !!!integer(kind=4), intent(in) :: layoutnumber,num_procs
   !!!integer(kind=4) :: ierr,i
   !!!integer(kind=8) :: posic
   !!!character(len=dxflinesize) :: buf2
   !!!character(len=BUFSIZE) :: dubuf
   !!!character(len=BUFSIZE) :: whoamishort
   !!!logical lexis
   !!!
   !!!if (.not.dxfFileIsOpen) return
   !!!
   !!!
   !!!
   !!!!!#ifdef CompileWithMPI
   !!!!!call MPI_FILE_close (thefile2, ierr)
   !!!!!#else
   !!!close (97)
   !!!!!#endif
   !!!
   !!!#ifdef CompileWithMPI
   !!!        call MPI_Barrier(SUBCOMM_MPI,ierr)
   !!!#endif
   !!!
   !!!!fixes the NULs
   !!!if (layoutnumber == 0) then
   !!!    open (988,file=trim(adjustl(mynEntradaRoot))//'.dxf',form='formatted')
   !!!    posic=0
   !!!    do i=0,num_procs-1
   !!!        write(whoamishort,'(i5)') i+1
   !!!        inquire(file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf',exist=lexis)
   !!!        if (lexis) then
   !!!            open (987,file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf',form='formatted')
   !!!!
   !!!875         read(987,'(a)',end=876,err=876) buf2
   !!!            call trimnullchar(buf2)
   !!!            buf2=trim(adjustl(buf2))
   !!!            if ((buf2(1:1) /= ' ').and.(buf2(1:1) /=char(0))) then
   !!!                write(988,'(a)') trim(adjustl(buf2))
   !!!                posic=posic+1
   !!!            end if
   !!!            goto 875
   !!!!
   !!!876         close (987)
   !!!!
   !!!            open (987,file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf')
   !!!            write(987,*) '!END'
   !!!            close (987,status='delete')
   !!!        end if
   !!!    end do
   !!!!end the file
   !!!    write(988,'(a)')  'ENDSEC'
   !!!    write(988,'(a)')  '0'
   !!!    write(988,'(a)')  'EOF'
   !!!    close (988)
   !!!
   !!!
   !!!
   !!!    write(dubuf,*)SEPARADOR//separador//separador
   !!!    call print11(layoutnumber,dubuf)
   !!!    write(dubuf,*) 'Closing dxf file. Number of lines: ',posic
   !!!    call print11(layoutnumber,dubuf)
   !!!end if
   !!!
   !!!continue
   !!!
   !!!dxfFileIsOpen=.false.
   !!!return
   !!!end subroutine

   !!!
   !!!subroutine writemmdxf(layoutnumber,sgg,sggMiHx,sggMiHy,sggMiHz)
   !!!integer i,j,k,layoutnumber
   !!!type(SGGFDTDINFO_t), intent(in) :: sgg
   !!!integer(kind=INTEGERSIZEOFMEDIAMATRICES), intent(in) :: &
   !!!           sggMiHx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE, &
   !!!                   sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE, &
   !!!                   sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE), &
   !!!           sggMiHy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE, &
   !!!                   sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE, &
   !!!                   sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE), &
   !!!           sggMiHz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE, &
   !!!                   sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE, &
   !!!                   sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE)
   !!!
   !!!!write ONLY PECS
   !!!      Do k=sgg%SINPMLSweep(iHx)%ZI,sgg%SINPMLSweep(iHx)%ZE
   !!!          Do j=sgg%SINPMLSweep(iHx)%YI,sgg%SINPMLSweep(iHx)%YE
   !!!              Do i=sgg%SINPMLSweep(iHx)%XI,sgg%SINPMLSweep(iHx)%XE
   !!!                if ((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec))  then
   !!!                    write(dxfbuff,'(a)') '3DFACE'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '8'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHx(i,j,k)+20
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '62'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHx(i,j,k)+20
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '10'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '20'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '30'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '11'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '21'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '31'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    call DXFWRITE(DXFBUFF)
   !!!
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '12'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '22'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '32'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '13'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '23'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '33'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '0'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                end if
   !!!            end do
   !!!        end do
   !!!    end do
   !!!!
   !!!      Do k=sgg%SINPMLSweep(iHy)%ZI,sgg%SINPMLSweep(iHy)%ZE
   !!!          Do j=sgg%SINPMLSweep(iHy)%YI,sgg%SINPMLSweep(iHy)%YE
   !!!              Do i=sgg%SINPMLSweep(iHy)%XI,sgg%SINPMLSweep(iHy)%XE
   !!!                if ((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec))  then
   !!!                    write(dxfbuff,'(a)') '3DFACE'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '8'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHy(i,j,k)+20
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '62'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHy(i,j,k)+20
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '10'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '20'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '30'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '11'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '21'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '31'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '12'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '22'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '32'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '13'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '23'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '33'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '0'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                end if
   !!!            end do
   !!!        end do
   !!!    end do
   !!!!
   !!!      Do k=sgg%SINPMLSweep(iHz)%ZI,sgg%SINPMLSweep(iHz)%ZE
   !!!          Do j=sgg%SINPMLSweep(iHz)%YI,sgg%SINPMLSweep(iHz)%YE
   !!!              Do i=sgg%SINPMLSweep(iHz)%XI,sgg%SINPMLSweep(iHz)%XE
   !!!                if ((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec))  then
   !!!                    write(dxfbuff,'(a)') '3DFACE'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '8'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHz(i,j,k)+20
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '62'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHz(i,j,k)+20
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '10'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '20'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '30'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '11'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '21'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '31'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    call DXFWRITE(DXFBUFF)
   !!!
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '12'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '22'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '32'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '13'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '23'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j+1
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '33'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '0'
   !!!                    call DXFWRITE(DXFBUFF)
   !!!                end if
   !!!            end do
   !!!        end do
   !!!    end do
   !!!    return
   !!!end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function openfile_mpi(layoutnumber, nombrefich) result(thefile8)
    integer(4), intent(in) :: layoutnumber
    character(len=BUFSIZE), intent(in) :: nombrefich
    integer(4) :: thefile8, iter, ios
    character(len=BUFSIZE) :: whoamishort
    logical :: file_exists
#ifdef CompileWithMPI
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer(4) :: ierr
#endif

    write(whoamishort, '(i5)') layoutnumber + 1
    do iter = 1, 10
        ios = 0
        open(newunit=thefile8, file=trim(adjustl(nombrefich)) // trim(adjustl(whoamishort)) // '_tmp', iostat=ios)
        if (ios /= 0) print *, 'Error opening temporary file: ', trim(adjustl(nombrefich)) // trim(adjustl(whoamishort)) // '_tmp'
        if ((layoutnumber == 0) .and. (ios == 0)) then
            call sleep(2)
            write(thefile8, '(a)') '!END'
            call sleep(2)
            ios = 0
            close(thefile8, status='delete',iostat=ios)
            if (ios /= 0) print *, 'Error deleting temporary file: ', trim(adjustl(nombrefich)) // trim(adjustl(whoamishort)) // '_tmp'
        end if
        call sleep(2)
        if (ios == 0) exit
    end do

#ifdef CompileWithMPI
    call MPI_Barrier(SUBCOMM_MPI, ierr)
#endif

end function openfile_mpi

   subroutine writefile_mpi(layoutnumber, thefile8,buff2)

      integer(kind=4) :: thefile8
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif
      integer(kind=4), intent(in) :: layoutnumber
      character(len=BUFSIZE) :: buff2,buff3

      buff3=CHAR(13)//CHAR(10)//trim(adjustl(buff2))
      !!!#ifdef CompileWithMPI
      !!!                call MPI_FILE_WRITE(thefile8, buff3  , 1024_4, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      !!!#else
      write(thefile8,'(a)') buff2
      !!!#endif

   end subroutine writefile_mpi

   subroutine closefile_mpi(layoutnumber,num_procs,nombrefich,thefile8)

      integer(kind=4) :: thefile8,thefile19
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif
      integer(kind=4), intent(in) :: layoutnumber,num_procs
      character(len=BUFSIZE) :: buff2
      character(len=BUFSIZE) :: nombrefich
      integer(kind=4) :: conta,i
      character(len=BUFSIZE) :: whoamishort
      logical lexis
      character(len=BUFSIZE) :: lockFile


      !!!!#ifdef CompileWithMPI
      !!!!call MPI_Barrier(SUBCOMM_MPI,ierr)
      !!!!call MPI_FILE_CLOSE(thefile8, ierr)
      !!!!#else
      close (thefile8)
      !!!#endif


#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      !fixes the NULs
      if (layoutnumber == 0) then
         open (newunit=thefile8,file=trim(adjustl(nombrefich)),form='formatted')
         do i=0,num_procs-1
            write(whoamishort,'(i5)') i+1
            inquire(file=trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp',exist=lexis)
            if (lexis) then
               open (newunit=thefile19,file=trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp',form='formatted')
               conta=0
875            read(thefile19,'(a)',end=876,err=876) buff2
               call trimnullchar(buff2)
               buff2=trim(adjustl(buff2))
               if ((buff2(1:1) /= ' ').and.(buff2(1:1) /=char(0))) then
                  write(thefile8,'(a)') trim(adjustl(buff2))
                  conta=conta+1
               end if
               goto 875
               !
876            continue
               close (thefile19)
               lockFile=trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp'
               call openclosedelete(lockFile)
            end if
         end do
         close (thefile8)

         !
      end if
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      !!!!end gnuplot
      return
   end subroutine closefile_mpi

   subroutine reportmedia(sgg)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      integer(kind=4) :: j
      character(len=BUFSIZE) :: buff

      do j = 0, sgg%NumMedia
         write(buff,*) '_____________________________'
         call WarnErrReport(Trim(buff))
         write(buff,*) 'MEDIO :  ', j
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Priority ', sgg%Med(j)%Priority
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Epr ', sgg%Med(j)%Epr
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Sigma ', sgg%Med(j)%Sigma
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Mur ', sgg%Med(j)%Mur
         call WarnErrReport(Trim(buff))
         write(buff,*) 'SigmaM ', sgg%Med(j)%SigmaM
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is PML ', sgg%Med(j)%Is%PML
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is PEC ', sgg%Med(j)%Is%PEC
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is ThinWIRE ', sgg%Med(j)%Is%ThinWire
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is MULTIWIRE ', sgg%Med(j)%Is%Multiwire
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is SlantedWIRE ', sgg%Med(j)%Is%SlantedWire
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is EDispersive ', sgg%Med(j)%Is%EDispersive
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is EDispersiveaANIS', sgg%Med(j)%Is%EDispersiveANIS
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is MDispersive ', sgg%Med(j)%Is%MDispersive
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is MDispersiveANIS ', sgg%Med(j)%Is%MDispersiveANIS
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is ThinSlot ', sgg%Med(j)%Is%ThinSlot
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is SGBC ', sgg%Med(j)%Is%SGBC
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is Lossy ', sgg%Med(j)%Is%Lossy
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is Multiport ', sgg%Med(j)%Is%multiport
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is AnisMultiport ', sgg%Med(j)%Is%anismultiport
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is MultiportPadding ', sgg%Med(j)%Is%multiportpadding
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is Dielectric ', sgg%Med(j)%Is%DIELECTRIC
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is ThinSlot ', sgg%Med(j)%Is%ThinSlot
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is Anisotropic ', sgg%Med(j)%Is%Anisotropic
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is Needed ', sgg%Med(j)%Is%Needed
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is already_YEEadvanced_byconformal ', sgg%Med(j)%Is%already_YEEadvanced_byconformal
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Iss split_and_useless ', sgg%Med(j)%Is%split_and_useless
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is Volume ', sgg%Med(j)%Is%Volume
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is Surface ', sgg%Med(j)%Is%Surface
         call WarnErrReport(Trim(buff))
         write(buff,*) 'Is Line ', sgg%Med(j)%Is%Line
         call WarnErrReport(Trim(buff))
         write(buff,*) '_____________________________'
         call WarnErrReport(Trim(buff))
      end do
      return
   end subroutine reportmedia

   subroutine erasesignalingfiles(simu_devia)
      logical :: simu_devia
      character(len=BUFSIZE) :: lockFile
      if (.not.simu_devia) then
          !!force erasing the signaling files
          lockFile='stop'
          call openclosedelete(lockFile)
          lockFile='stopflushing'
          call openclosedelete(lockFile)
          lockFile='flush'
          call openclosedelete(lockFile)
          lockFile='flushdata'
          call openclosedelete(lockFile)
          lockFile='unpack'
          call openclosedelete(lockFile)
          lockFile='postprocess'
          call openclosedelete(lockFile)
          lockFile='flushxdmf'
          call openclosedelete(lockFile)
          lockFile='flushvtk'
          call openclosedelete(lockFile)
          lockFile='snap'
          call openclosedelete(lockFile)         !
          !
          lockFile='stop_only'
          call openclosedelete(lockFile)
          lockFile='stopflushing_only'
          call openclosedelete(lockFile)
          lockFile='flush_only'
          call openclosedelete(lockFile)
          lockFile='flushdata_only'
          call openclosedelete(lockFile)
          !
          !
          lockFile='stop_dontwritevtk'
          call openclosedelete(lockFile)
          lockFile='stop_only_dontwritevtk'
          call openclosedelete(lockFile)
          lockFile='stopflushing_dontwritevtk'
          call openclosedelete(lockFile)
          lockFile='stopflushing_only_dontwritevtk'
          call openclosedelete(lockFile)
          lockFile='flush_dontwritevtk'
          call openclosedelete(lockFile)
          lockFile='flush_only_dontwritevtk'
          call openclosedelete(lockFile)
      !     
          lockFile='unpack'
          call openclosedelete(lockFile)
          lockFile='postprocess'
          call openclosedelete(lockFile)
          lockFile='flushxdmf'
          call openclosedelete(lockFile)
          lockFile='flushvtk'
          call openclosedelete(lockFile)
      end if
      return
        
   
  end subroutine erasesignalingfiles

   
   subroutine openclosedelete(ficherin)
      character(len=*) :: ficherin
      integer(kind=4) :: my_iostat, myunit
      my_iostat = 0
4216  if(my_iostat /= 0) write(*,FMT='(a)',advance='no'), '.'
      open  (newunit=myunit,file=trim(adjustl(ficherin)),form='formatted',err=4216,iostat=my_iostat)
      write (myunit,*) '!END'
      close (myunit,status='delete')
      return
   end subroutine openclosedelete
   
   subroutine openclose(ficherin)
      character(len=*) :: ficherin
      integer(kind=4) :: my_iostat, myunit
      my_iostat = 0
5216  if(my_iostat /= 0) write(*,FMT='(a)',advance='no'), '.'
      open  (newunit=myunit,file=trim(adjustl(ficherin)),form='formatted',err=5216,iostat=my_iostat)
      write (myunit,*) '!END'
      close (myunit)
      return
   end subroutine openclose
   
   
   subroutine opensolo(myunit,ficherin)
      character(len=*) :: ficherin
      integer(kind=4) :: my_iostat, myunit
      my_iostat = 0
6216  if(my_iostat /= 0) write(*,FMT='(a)',advance='no'), '.'
      open  (unit=myunit,file=trim(adjustl(ficherin)),form='formatted',err=6216,iostat=my_iostat)
      return
   end subroutine opensolo
   
   subroutine closesolo(myunit)
      integer(kind=4) :: my_iostat, myunit
      my_iostat = 0
7216  if(my_iostat /= 0) write(*,FMT='(a)',advance='no'), '.'
      close  (unit=myunit,err=7216,iostat=my_iostat)
      return
   end subroutine closesolo
   
end module
