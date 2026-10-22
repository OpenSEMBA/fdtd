module Report_m
   use FDETYPES_m

   implicit none
   private

   type :: time_t
      real(kind = 8) :: segundos
      character(LEN=BUFSIZE) :: hora
      character(LEN=BUFSIZE) :: fecha
   end type


   logical, save :: printea,warningfileIsOpen=.false.,verbose,file10isopen,file11isopen
   character(len=BUFSIZE), save :: warningFile = ' '
   character(len=BUFSIZE), save  :: whoami = ' '

   logical, save :: ignoreerrors
   logical, save :: fatalerror=.false.

   public StopOnError,CloseReportingFiles, &
   print11,Onprint,Offprint,file10isopen,file11isopen
   public WarnErrReport,INITWARNINGFILE,CLOSEWARNINGFILE,get_secnds, &
          erasesignalingfiles,openclosedelete,openclose
   public isFatalError, resetFatalError
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
