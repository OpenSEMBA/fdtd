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
 
Module Report
   use FDETYPES

   use snapxdmf

   implicit none
   private

   integer (kind=4), parameter  ::  reportingseconds=60
   type :: tiempo_t
      REAL ( kind = 8)  ::  segundos
      character ( LEN=BUFSIZE)  ::  hora
      character ( LEN=BUFSIZE)  ::  fecha
   end type


   !For timing
   REAL (KIND=8), SAVE ::  time_begin, time_end, time_begin2,time_begin3,time_begin_absoluto, time_end2,time_desdelanzamiento
   REAL (KIND=RKIND), SAVE ::  megaceldas,megaceldastotales,speedInst, speedGlobInst,speedAvg,speedGlobAvg
   REAL (KIND=RKIND), SAVE  ::  energy,energyTotal,oldenergyTotal,snapLevel
   type (tiempo_t), SAVE  ::  time_out2
   !
   character (LEN=BUFSIZE), SAVE :: charmeg
   integer (kind=4), SAVE   ::  reportedinstant,snapStep,snapHowMany,countersnap
   logical, SAVE :: printea,calledStoponerrroonlyprint=.false.,warningfileIsOpen=.false.,verbose,file10isopen,file11isopen
   CHARACTER (LEN=BUFSIZE), SAVE :: warningFile = ' '
   character (LEN=BUFSIZE), save  ::  whoami

   integer, save ::  thefile !for mpi file management
   logical, save :: ignoreerrors
   !
   type (coorsxyzP) , save  ::  Punto

   CHARACTER (LEN=BUFSIZE), SAVE :: mynEntradaRoot

   !!!logical, SAVE :: dxfFileIsOpen=.false.
   logical, SAVE :: fatalerror=.false.
   integer, SAVE :: CONTADORDEMENSAJES
   integer, save ::  thefile2 !for mpi file management
   !!public StopOnError_OnlyPrint

   public StopOnError,InitReporting,ReportExistence,InitTiming,Timing,CloseReportingFiles, &
   print11,Onprint,Offprint,file10isopen,file11isopen
   public WarnErrReport,INITWARNINGFILE,CLOSEWARNINGFILE,get_secnds,openfile_mpi,writefile_mpi, &
          closefile_mpi,reportmedia,erasesignalingfiles,openclosedelete,openclose

   !part of the dxf
   !!!public dxfwrite,INITdxfFILE,CLOSEdxfFILE,writemmdxf,TRIMNULLCHAR
   public TRIMNULLCHAR,tiempo_t




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
   subroutine StopOnError(layoutnumber,size,message,calledfrommain)
      character (LEN=BUFSIZE) :: ficherito
      logical , optional  ::  calledfrommain
      character (len=*), intent( IN) :: message
      integer (kind=4), optional  ::  layoutnumber,size
#ifdef CompileWithMPI
      integer (kind=4)  ::  ierr
#endif
      character (LEN=BUFSIZE)  ::  whoami

      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '


      call print11(layoutnumber,trim(adjustl(whoami))//' ERROR: '//trim(adjustl(message)),.true.)

      !19/12/14 bug OLD1812. Un stoponerror creado por un nodal source embebido llega aqui en MPI. El closewarn... hace un barrier e impide morir al proceso.

      !hay que revisar los stoponerror y hacerlos mas elegantes. De momento aborto a lo bestia comentanod sin cerrar ni warning ni dxf (To do)

      !CALL CLOSEWARNINGFILE(layoutnumber,size)
      !!!!CALL CLOSEdxfFILE(layoutnumber,size)



      !!!!lo de volver es solo para el cluster
#ifdef keeppause
      if (present(calledfrommain)) then
         if (calledfrommain) then
            IF (layoutnumber == 0) THEN
               OPEN (38, FILE='pause')
               WRITE (38, '(a)') '!END'
               CLOSE (38)
            endif
            call print11(layoutnumber,'Trying to relaunch. Correct error, create launch, and remove pause/warning '// &
                                      'file (or kill the process)',.true.)
!!!            call CloseReportingFiles !sgg 240817 no se deben cerrar los reporting files
            return
         endif
      else
         IF (layoutnumber == 0) THEN
            OPEN (38, FILE='pause')
            WRITE (38, '(a)') '!END'
            CLOSE (38)
         endif
         call print11(layoutnumber,'Stopping, but creating the signal file pause to prevent queuing losses!!! '// & '
                                   '(correct error and remove to continue)',.true.)
      endif
#else
      IF (layoutnumber == 0) THEN
         ficherito='running'
          call openclosedelete(ficherito)
         !
         ficherito='pause'
          call openclosedelete(ficherito)
         !
         ficherito='relaunch'
          call openclosedelete(ficherito)
      endif
#endif

#ifdef CompileWithMPI
      call print11(layoutnumber,'Trying to kill all MPI processes (may fail!)...',.true.)
      Call MPI_Abort(SUBCOMM_MPI, -1, ierr)
      call MPI_FINALIZE(ierr)
#endif
      call CloseReportingFiles



      STOP


      return

   end subroutine




   subroutine CloseReportingFiles
      if (file10isopen) then
         close (10) !energy file
         file10isopen=.false.
      endif
      if (file11isopen) then
         close (11) !reporting file
         file11isopen=.false.
      endif
      return

   end subroutine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine InitReporting(sgg,nEntradaRoot,resume,layoutnumber,size,nresumeable2,resume_fromold)
      !!!!!!!PML params!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer (kind=4), intent(in) :: layoutnumber,size
      type (SGGFDTDINFO), intent(INout)         ::  sgg
      !!!!!!!
      character (len=*), INTENT(IN)  ::  nEntradaRoot
      logical  ::  resume,resume_fromold
#ifdef CompileWithMPI
      integer (kind=4)  ::  ierr
#endif
      character (len=*), INTENT(IN)  ::  nresumeable2
      Logical  ::  errnofile
      character(len=BUFSIZE) :: buff

      character (LEN=BUFSIZE)  ::  whoami

      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      !
      Punto=Creapuntos(sgg) !crea coordenadas fisicas

      !

      IF (layoutnumber == 0) THEN  !only the master
         if (resume) then
            open (10,file=trim(adjustl(nEntradaRoot))//'_Energy.dat',form='formatted',position='append')
         else
            open (10,file=trim(adjustl(nEntradaRoot))//'_Energy.dat',form='formatted')
         endif
      endif
      file10isopen=.true.

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      call get_secnds(time_out2)
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      call print11(layoutnumber,SEPARADOR//separador//separador)


#ifndef CompileWithInt4
#define CompileWithInt4
#endif
      !

      if (resume) then
         errnofile=.false.
         if (resume_fromold) then
            INQUIRE (FILE=trim(adjustl(nresumeable2))//'.old', EXIST=errnofile)
         else
            INQUIRE (FILE=trim(adjustl(nresumeable2)), EXIST=errnofile)
         endif
         if (.not.errnofile) then
            if (resume_fromold) then
               buff='FILE '//trim(adjustl(nresumeable2))//'.old DOES NOT EXIST'
               call StopOnError(layoutnumber,size,buff)
            else
               buff='FILE '//trim(adjustl(nresumeable2))//' DOES NOT EXIST'
               call StopOnError(layoutnumber,size,buff)
            endif
         endif
         call print11(layoutnumber,SEPARADOR//SEPARADOR//SEPARADOR)
         call print11(layoutnumber,' ')
         if (resume_fromold) then
            call print11(layoutnumber,'Reading resuming data from '//trim(adjustl(nresumeable2))//'.old etc.')
         else
            call print11(layoutnumber,'Reading resuming data from '//trim(adjustl(nresumeable2))//' etc.')
         endif
         call print11(layoutnumber,SEPARADOR//sEPARADOR//SEPARADOR)
      endif


#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif


   end subroutine InitReporting


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ReportExistence(sgg,layoutnumber,size,thereare,mur_second,MurAfterPML)
      logical :: mur_second,MurAfterPML
      type (SGGFDTDINFO), intent(IN)  ::  sgg
      !
      type (logic_control), intent(in)  ::  thereare
      integer (kind=4), intent(in) :: layoutnumber,size
#ifdef CompileWithMPI
      integer (kind=4)  ::  ierr
#endif
      character (LEN=BUFSIZE)  ::  whoami
      character(len=BUFSIZE) :: buff

      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif


      call print11(layoutnumber,SEPARADOR//sEPARADOR//SEPARADOR)
      !!!
      !
      IF ((thereare%Multiports).or.(thereare%AnisMultiports))     then
#ifdef CompileWithNIBC
         continue
#else
         buff=trim(adjustl(whoami))//' MIBC unsupported. Recompile'
         call stoponerror(layoutnumber,size,buff)
#endif
      endif
      !
      !
      IF (thereare%ThinSlot)     then
#ifdef CompileWithDMMA
         continue
#else
         buff=trim(adjustl(whoami))//' Thin slots unsupported. Recompile'
         call stoponerror(layoutnumber,size,buff)
#endif
      endif
      !
      IF (thereare%EDispersives.or.thereare%MDispersives)     then
#ifdef CompileWithEDispersives
         continue
#else
         buff=trim(adjustl(whoami))//' Dispersive materials unsupported. Recompile'
         call stoponerror(layoutnumber,size,buff)
#endif
      endif
      !
      If (thereare%Wires)            then
#ifdef CompileWithWires
         continue
#else
#ifdef CompileWithBerengerWires
         continue
#else
#ifdef CompileWithSlantedWires
         continue
#else
         buff=trim(adjustl(whoami))//'  WIREs unsupported. Recompile'
         call stoponerror(layoutnumber,size,buff)
#endif
#endif
#endif
      endif
      !
      !!!!!!!!!!!!!
      if (thereAre%MagneticMedia) then
         buff=' has special H-media'
         call warnerrreport(buff)
      endif
      if (thereAre%PMLMagneticMedia) then
         buff=' has special PML H-media'
         call warnerrreport(buff)
      endif
      if ((thereare%NodalE).or.(thereare%NodalH))    then
         buff=' has Nodal sources'
         call warnerrreport(buff)
      endif
      if (thereare%Observation)    then
         buff=   ' has probes'
         call warnerrreport(buff)
      endif
      if (thereare%FarFields)    then
         buff=   ' has Far Field probes'
         call warnerrreport(buff)
      endif
      If (thereare%PlaneWaveBoxes) then
         buff=   ' has planewaves'
         call warnerrreport(buff)
      endif
      IF (thereare%Multiports)     then
         buff=   ' has MIBC Multiports'
         call warnerrreport(buff)
      endif
      IF (thereare%AnisMultiports)     then
         buff=   ' has MIBC Anisotropic Multiports'
         call warnerrreport(buff)
      endif
      IF (thereare%SGBCs)     then
         buff=   ' has Thin metal Materials'
         call warnerrreport(buff)
      endif
      IF ((thereare%Anisotropic).and.(.not.thereare%ThinSlot))     then
         buff=   ' has pure anisotropic media'
         call warnerrreport(buff)
      endif
#ifdef CompileWithDMMA
      IF (thereare%ThinSlot)     then
         buff=   ' has Thin Slots'
         call warnerrreport(buff)
      endif
#endif
      !
#ifdef CompileWithEDispersives
      IF (thereare%EDispersives)     then
         buff=   ' has electric dispersives'
         call warnerrreport(buff)
      endif
      IF (thereare%MDispersives)     then
         buff=   ' has magnetic dispersives'
         call warnerrreport(buff)
      endif
#endif
#ifdef CompileWithWires
      If (thereare%Wires)            then
         buff=   ' has Holland WIREs'
         call warnerrreport(buff)
      endif
#endif
#ifdef CompileWithBerengerWires
      If (thereare%Wires)            then
         buff=   ' has Multi-WIREs'
         call warnerrreport(buff)
      endif
#endif
#ifdef CompileWithSlantedWires
      If (thereare%Wires)            then
         buff=   ' has Slanted WIREs'
         call warnerrreport(buff)
      endif
#endif
      If (thereare%PMLBorders)      then
         if (sgg%Border%IsUpPML.or.sgg%Border%IsDownPML) then
            buff=   ' has PML regions inside Z'
            call warnerrreport(buff)
         endif
      endif
      If (thereare%MURBorders)      then
         if (sgg%Border%IsUpMUR.or.sgg%Border%IsDownMUR) then
            if (mur_second) then
               buff=   ' has MUR2 regions inside Z'
            else
               buff=   ' has MUR1 regions inside Z'
            endif
            call warnerrreport(buff)
         endif
      endif
      If (murAfterPML)      then
         if (mur_second) then
            buff=   ' CPML are backed by MUR1'
         else
            buff=   ' CPML are backed by MUR2'
         endif
         call warnerrreport(buff)
      endif
      If (thereare%PMCBorders)      then
         buff=   ' has PMC borders'
         call warnerrreport(buff)
      endif
      If (thereare%PECBorders)      then
         buff=   ' has PEC borders'
         call warnerrreport(buff)
      endif

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      call print11(layoutnumber,SEPARADOR//sEPARADOR//SEPARADOR)
      call warnerrreport(buff)

   end subroutine ReportExistence

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitTiming(layoutnumber,size,maxCPUtime,timedummy_desdelanzamiento, &
                         flushsecondsFields,flushsecondsData,initialtimestep,finaltimestep,c,maxSourceValue,sgg)
      type (SGGFDTDINFO), intent(IN)       :: sgg
      TYPE (tiempo_t) :: time_out2,time_comienzo
      type (XYZlimit_t), dimension(1:6)  ::  c
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4)  ::  maxCPUtime
#ifdef CompileWithMPI
      integer (kind=4)  ::  ierr
#endif
      integer (kind=4), intent(in)           ::  flushsecondsFields,flushsecondsData,initialtimestep,finaltimestep

      REAL (KIND=RKIND)   ::  MaxSourceValue
      REAL (KIND=8)   :: timedummy_desdelanzamiento

      character (LEN=BUFSIZE)  ::  whoami
      character (LEN=BUFSIZE)     ::  dubuf
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '
      
      time_desdelanzamiento=timedummy_desdelanzamiento
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      snapLevel=1.0e25_RKIND !*maxSourceValue
      snapStep=1
      snapHowMany=1
      countersnap=0
      !

      megaceldas=(1.0_RKIND*C(iEx)%ZE-1.0_RKIND*C(iEx)%ZI)*(1.0_RKIND*C(iEx)%YE-1.0_RKIND*C(iEx)%YI)* &
                 (1.0_RKIND*C(iEy)%XE-1.0_RKIND*C(iEy)%XI)/1.0e6_RKIND


#ifdef CompileWithMPI
      call MPI_BARRIER(SUBCOMM_MPI,ierr)
      call MPI_AllReduce( megaceldas, megaceldastotales, 1_4, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
#else
      megaceldastotales=megaceldas
#endif


      write(dubuf,*)  'Total Mcells: ',megaceldastotales
      call print11(layoutnumber,dubuf)


      IF (flushsecondsFIELDS/=0) then
         write(dubuf,*)  'Flushing restarting FIELDS every ',int(flushsecondsFIELDS/60.0_RKIND),' minutes'
         call print11(layoutnumber,dubuf)
      else
         if (maxCPUtime == topCPUtime) then
            call print11(layoutnumber,'NO flushing of restarting FIELDS scheduled')
         else
            write(dubuf,*)  'Flushing of restarting FIELDS at the end (mins) :',maxCPUtime
            call print11(layoutnumber,dubuf)
         endif
      endif
      IF (flushsecondsDATA/=0) then
         write(dubuf,*)  'Flushing observation DATA every  ',int(flushsecondsDATA/60.0_RKIND),' minutes and every ', &
                          BuffObse,' steps'
         call print11(layoutnumber,dubuf)
      else
         call print11(layoutnumber,'WARNING: NO flushing of observation DATA scheduled')
      endif
      write(dubuf,*)  'Reporting simulation info every  ',int(reportingseconds/60.0_RKIND),' minutes '
      call print11(layoutnumber,dubuf)

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      call get_secnds(time_out2)
      call print11(layoutnumber,SEPARADOR//separador//separador)
      write(dubuf,'(a,i7,a,e19.9e3,a,i9,a,e19.9e3)')  'Simulation from n=',initialtimestep,', t=',sgg%tiempo(initialtimestep),&
                                                      ' to n=',finaltimestep,', t=',sgg%tiempo(finaltimestep)
      call print11(layoutnumber,dubuf)
      write(dubuf,*)  'Date/time ', time_out2%fecha( 7: 8),'/',time_out2%fecha( 5: 6),'/',time_out2%fecha(1:4),'   ', &
                                    time_out2%hora( 1: 2), ':',time_out2%hora( 3: 4),':',time_out2%hora( 5: 6)
      call print11(layoutnumber,dubuf)
      time_begin_absoluto = time_out2%segundos
      time_begin = time_begin_absoluto
      time_begin2 = time_begin_absoluto
      time_begin3 = time_begin_absoluto

      reportedinstant=initialtimestep

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

   end subroutine InitTiming



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !---------------------------------------------------->
   subroutine get_secnds(time_out2)
      integer (kind=4), dimension(0:11)  ::  diasen, diasenbisiesto
      data diasen /31,59 ,90 ,120,151,181,212,243,273,304,334,365/
      data diasenbisiesto /31,60 ,91 ,121,152,182,213,244,274,305,335,366/
      type (tiempo_t), intent(out)  ::  time_out2

      REAL ( kind = 8)  ::  time_out

      !--->
      REAL ( kind = 8)  ::  s
      REAL ( kind = 8), save  ::  t_0 = 0.0_RKIND
      !--->
      integer (kind=4)  ::  h, m,month,day,year,cent
      !integer (kind=4)  ::  cent,year

      !--->

      !--->
      character ( LEN=BUFSIZE)  ::  caux
      character ( LEN=BUFSIZE)  ::  caux2
      character(5)  :: zone
      integer (kind=4) ,dimension(8) :: values
      !------------------------------------------------->
      !   if( esprimeravez) then
      !        call date_and_time( time = caux,date=caux2)
      !        read( caux2( 7: 8), '(i2)') day
      !!        read( caux2( 3: 4), '(i2)') year
      !        read( caux2( 5: 6), '(i2)') month
      !   if ((mod(year,4)==0).and.(year/=00)) then
      !        t_0 = (month-1) * diasenbisiesto(month-1) * 86400 + day * 86400 ! + (year-2000.) * 365 * 86400
      !   else
      !        t_0 = (month-1) * diasen(month-1) * 86400 + day * 86400 ! + (year-2000.) * 365 * 86400
      !   endif
      !        esprimeravez = .FALSE.
      !   endif
      call date_and_time(date=caux2,time = caux,zone=zone,values=values)
      read( caux( 1: 2), '(i2)') h
      read( caux( 3: 4), '(i2)') m
      read( caux( 5: 10), '(f6.3)') s
      read( caux2( 1: 2), '(i2)') cent
      read( caux2( 3: 4), '(i2)') year
      read( caux2( 5: 6), '(i2)') month
      read( caux2( 7: 8), '(i2)') day 

      if ((mod(year,4)==0).and.(year/=00)) then
         time_out = diasenbisiesto(month-1) * 86400 + (day-1) * 86400 + 3600.0 * h + 60.0 * m + s - t_0  + &
                    (year-2000.) * 365 * 86400.
      else
         time_out = diasen(month-1) * 86400 + (day-1) * 86400 + 3600.0 * h + 60.0 * m + s - t_0  + (year-2000.) * 365 * 86400.
      endif
      time_out2%segundos=time_out !seconds from year 2000
      time_out2%hora=caux !hora
      time_out2%fecha=caux2 !fecha


      !!if ((month >=4 ).and.(month <= 10)) then
      !!   h=h+1 !dst aproximado
      !!   write( time_out2%hora( 1: 2), '(i2)') h
      !!endif
      return
   endsubroutine get_secnds



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !**************************************************************************************************
   subroutine Timing(sgg, b, n, n_info, layoutnumber, size, maxCPUtime,flushsecondsFields, flushsecondsData, initialtimestep, &
   finaltimestep, performflushFields, performflushData,performUnpack, performpostprocess,performflushXdmf,performflushVTK, &
   parar, forcetiming,Ex,Ey,Ez,everflushed, nentradaroot,maxSourceValue,opcionestotales,simu_devia,dontwritevtk,permitscaling)
   
      logical :: simu_devia,dontwritevtk,stopdontwritevtk,stopflushingdontwritevtk,flushdontwritevtk,stoponlydontwritevtk
      !---------------------------> inputs <----------------------------------------------------------
      type (SGGFDTDINFO), intent(IN)              :: sgg              ! Simulation data.
      type( bounds_t), intent( IN)  ::  b
      CHARACTER (LEN=BUFSIZE), intent(in) :: opcionestotales
      integer( kind = 4), intent( IN)  ::  layoutnumber, size, n,maxCPUtime
      integer( kind = 4), intent( IN)  ::  flushsecondsFields, flushsecondsData, initialtimestep, finaltimestep
      !--->
      real (kind = RKIND), dimension( 0: b%Ex%NX-1, 0: b%Ex%NY-1, 0: b%Ex%NZ-1), intent( IN)  ::  Ex
      real (kind = RKIND), dimension( 0: b%Ey%NX-1, 0: b%Ey%NY-1, 0: b%Ey%NZ-1), intent( IN)  ::  Ey
      real (kind = RKIND), dimension( 0: b%Ez%NX-1, 0: b%Ez%NY-1, 0: b%Ez%NZ-1), intent( IN)  ::  Ez
      !--->
      !!!
      integer :: my_iostat
      !--->
      logical, intent( IN)  ::  forcetiming, everflushed,permitscaling
      !
      character (LEN=BUFSIZE) :: fichsnap,minmax,quien_es
      character( len = *), intent( IN)  ::  nEntradaRoot
      !---------------------------> input/output <----------------------------------------------------
      integer (kind=4), intent( INOUT)  ::  n_info
      logical, intent( INOUT)  ::  parar
      !---------------------------> outputS <---------------------------------------------------------
      logical, intent( OUT)  ::  performflushFIELDS, performflushDATA,performUnpack,performpostprocess,&
                                 performflushXdmf,performflushVTK
      !---------------------------> variables locales <-----------------------------------------------
      real (kind=rKIND)  ::  valor,maxSourceValue,LA,LV,LB
      logical  ::  hay_timing, l_aux, hay_flushFIELDS, hay_flushDATA, mustflushFIELDS, mustflushDATA,mustUnpack, &
                   mustPostprocess,mustflushXdmf , mustflushVTK ,   &
      pararflushing, pararNOflushing, stoponNaN , stoponNaN_aux,mustSnap,stop_only,stopflushing_only,flush_only,flushdata_only
      logical :: stopflushingonlydontwritevtk,flushonlydontwritevtk,flushdataonlydontwritevtk,flushdatadontwritevtk
      integer( kind = 4)  ::  in_aux, ini_i, fin_i, ini_j, fin_j, ini_k, fin_k, i, j, k
      character( LEN=BUFSIZE)  ::  whoamishort,whoami,chinstant
      character (LEN=BUFSIZE)     ::  dubuf
      character (LEN=BUFSIZE)     ::  dondex,dondey,dondez
      real (kind=rKIND), dimension (1:size) :: NEWlmaxval,NEWlmaxval_x,NEWlmaxval_y,NEWlmaxval_z
      integer( kind = 4), dimension (1:size) :: NEWlmaxval_i,NEWlmaxval_j,NEWlmaxval_k
      real (kind=rKIND), dimension (1:size) :: lmaxval,lmaxval_x,lmaxval_y,lmaxval_z
      integer( kind = 4), dimension (1:size) :: lmaxval_i,lmaxval_j,lmaxval_k
      real (kind=rKIND) :: qmaxval , qmaxval_x,qmaxval_y,qmaxval_z
      integer( kind = 4)  :: qmaxval_i,qmaxval_j,qmaxval_k,thefilenoflu
      real (kind=rKIND), dimension (1:size) :: NEWlminval,NEWlminval_x,NEWlminval_y,NEWlminval_z
      integer( kind = 4), dimension (1:size) :: NEWlminval_i,NEWlminval_j,NEWlminval_k
      real (kind=rKIND), dimension (1:size) :: lminval,lminval_x,lminval_y,lminval_z
      integer( kind = 4), dimension (1:size) :: lminval_i,lminval_j,lminval_k
      real (kind=rKIND) :: qminval , qminval_x,qminval_y,qminval_z
      integer( kind = 4)  :: qminval_i,qminval_j,qminval_k,dimxsnap,dimysnap,dimzsnap,veces,i1,j1,k1
      integer( kind = 4)  :: ini_ibox,fin_ibox,ini_jbox,fin_jbox,ini_kbox,fin_kbox
      integer( kind = 4)  :: ini_iboxsin,fin_iboxsin,ini_jboxsin,fin_jboxsin,ini_kboxsin,fin_kboxsin

      character (LEN=BUFSIZE) :: ficherito

      real (kind=4), dimension(:,:,:,:), allocatable  ::  snap


!!#ifdef PreventCrayBug
!!#ifdef CompileWithReal4
!!      logical IsNaNf
!!#else
!!      logical IsNaNd
!!#endif
!!#endif


#ifdef CompileWithMPI
      integer (kind=4)  ::  ierr
#endif
      !!!
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '
      write(whoamishort,'(i5)') layoutnumber+1

      !---------------------------> empieza Timing <--------------------------------------------------
#ifdef CompileWithMPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      call get_secnds( time_out2)
      time_end = time_out2%segundos
      !--->
      l_aux = (time_end - time_begin  >  reportingseconds) .or. forcetiming
#ifdef CompileWithMPI
      !print *,'layoutnumber+1,l_aux, hay_timing pre',layoutnumber+1,l_aux, hay_timing
      call MPI_AllReduce( l_aux, hay_timing, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr) !TODOS STOCH Y NO STOCH 050619
      !print *,'layoutnumber+1,l_aux, hay_timing post',layoutnumber+1,l_aux, hay_timing
#else
      hay_timing = l_aux
#endif
      n_info = n + 5 !after replaced depending on the speed
      !--->
      mustflushFIELDS = .FALSE.
      performflushFIELDS = .FALSE.
      mustflushDATA = .FALSE.
      performflushDATA = .FALSE.
      mustUnpack = .false.
      performUnpack = .FALSE.
      mustPostprocess = .false.
      performpostprocess = .false.
      mustflushXdmf=.false.
      mustflushVTK=.false.
      performflushXdmf=.false.
      performflushVTK=.false.
      energy=0.0_RKIND
      !--->
      if ( hay_timing) then !no calculation of time until at least 300 seconds lapse
         performflushFIELDS = .FALSE.
         performflushDATA = .FALSE.
         if (abs(time_end - time_begin_absoluto) < 1.0_RKIND) time_end = 60.0_RKIND+time_begin_absoluto
         if (abs(time_end - time_begin         ) < 1.0_RKIND) time_end = 60.0_RKIND+time_begin
         speedInst = ((N - reportedinstant + 1) * megaceldas / (time_end - time_begin))
         speedAvg = ((N - INITIALtimeSTEP + 1) * megaceldas / (time_end - time_begin_absoluto))
         if (speedAvg == 0) speedAvg=100.0_RKIND
         if (speedInst == 0) speedInst=speedAvg
#ifdef CompileWithMPI
         !print *,'layoutnumber+1,speedInst, speedGlobInst,speedAvg, speedGlobAvg pre', &
         !         layoutnumber+1,speedInst, speedGlobInst,speedAvg, speedGlobAvg
         call MPI_AllReduce( speedInst, speedGlobInst, 1_4, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( speedAvg, speedGlobAvg, 1_4, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         !print *,'layoutnumber+1,speedInst, speedGlobInst,speedAvg, speedGlobAvg post', &
         !         layoutnumber+1,speedInst, speedGlobInst,speedAvg, speedGlobAvg
#else
         speedGlobInst = speedInst
         speedGlobAvg = speedAvg
#endif
         !
         in_aux = n + max(int((reportingseconds / (megaceldastotales / speedGlobInst))) + 1,1)
#ifdef CompileWithMPI
         !print *,'layoutnumber+1,in_aux, n_info pre',layoutnumber+1,in_aux, n_info
         call MPI_AllReduce( in_aux, n_info, 1_4, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
         !print *,'layoutnumber+1,in_aux, n_info post',layoutnumber+1,in_aux, n_info
#else
         n_info = in_aux
#endif
         inquire( FILE = 'stop', EXIST = pararNOflushing)
         inquire( FILE = 'stop_only', EXIST = stop_only)
         inquire( FILE = 'stop_dontwritevtk', EXIST = stopdontwritevtk)
         inquire( FILE = 'stop_only_dontwritevtk', EXIST = stoponlydontwritevtk)
         if (pararnoflushing) then
             dontwritevtk=.false.
         endif
         if (stopdontwritevtk) then
             pararnoflushing=.true.
             dontwritevtk=.true.
         endif
         if (stoponlydontwritevtk) then
             stop_only=.true.
             dontwritevtk=.true.
         endif
         if (stop_only) then
             open(newunit=thefilenoflu,FILE = 'stop_only',action="read")
             read(thefilenoflu,*) quien_es
             if (trim(adjustl(quien_es))==trim(adjustl(nentradaroot))) then
                 pararNOflushing=.true.
             else
                 pararNOflushing=.false.
             endif
             close(thefilenoflu)
         endif
         !
         inquire( FILE = 'stopflushing', EXIST = pararflushing)
         inquire( FILE = 'stopflushing_only', EXIST = stopflushing_only)
         inquire( FILE = 'stopflushing_dontwritevtk', EXIST = stopflushingdontwritevtk)
         inquire( FILE = 'stopflushing_only_dontwritevtk', EXIST = stopflushingonlydontwritevtk)
         if (pararflushing) then
             dontwritevtk=.false.
         endif
         if (stopflushingdontwritevtk) then
             pararflushing=.true.
             dontwritevtk=.true.
         endif
         if (stopflushingonlydontwritevtk) then
             stopflushing_only=.true.
             dontwritevtk=.true.
         endif
         if (stopflushing_only) then
             open(newunit=thefilenoflu,FILE = 'stopflushing_only',action="read")
             read(thefilenoflu,*) quien_es
             if (trim(adjustl(quien_es))==trim(adjustl(nentradaroot))) then
                 pararflushing=.true.
             else
                 pararflushing=.false.
             endif
             close(thefilenoflu)
         endif
         !        
         inquire( FILE = 'flush', EXIST = mustflushFIELDS)
         inquire( FILE = 'flush_only', EXIST = flush_only)
         inquire( FILE = 'flush_dontwritevtk', EXIST = flushdontwritevtk)
         inquire( FILE = 'flush_only_dontwritevtk', EXIST = flushonlydontwritevtk)
         if (mustflushFIELDS) then
             dontwritevtk=.false.
         endif
         if (flushdontwritevtk) then
             mustflushFIELDS=.true.
             dontwritevtk=.true.
         endif
         if (flushdontwritevtk) then
             flush_only=.true.
             dontwritevtk=.true.
         endif
         if (flush_only) then
             open(newunit=thefilenoflu,FILE = 'flush_only',action="read")
             read(thefilenoflu,*) quien_es
             if (trim(adjustl(quien_es))==trim(adjustl(nentradaroot))) then
                 mustflushFIELDS=.true.
             else
                 mustflushFIELDS=.false.
             endif
             close(thefilenoflu)
         endif
         !
         inquire( FILE = 'flushdata', EXIST = mustflushdata)
         inquire( FILE = 'flushdata_only', EXIST = flushdata_only)
         inquire( FILE = 'flushdata_dontwritevtk', EXIST = flushdatadontwritevtk)
         inquire( FILE = 'flushdata_only_dontwritevtk', EXIST = flushdataonlydontwritevtk)
         if (mustflushdata) then
             dontwritevtk=.false.
         endif
         if (flushdatadontwritevtk) then
             mustflushdata=.true.
             dontwritevtk=.true.
         endif
         if (flushdataonlydontwritevtk) then
             flushdata_only=.true.
             dontwritevtk=.true.
         endif
         if (flushdata_only) then
             open(newunit=thefilenoflu,FILE = 'flushdata_only',action="read")
             read(thefilenoflu,*) quien_es
             if (trim(adjustl(quien_es))==trim(adjustl(nentradaroot))) then
                 mustflushdata=.true.
             else
                 mustflushdata=.false.
             endif
             close(thefilenoflu)
         endif
         !
         inquire( FILE = 'unpack', EXIST = mustUnpack)
         inquire( FILE = 'postprocess', EXIST = mustPostprocess)
         inquire( FILE = 'flushxdmf', EXIST = mustflushXdmf)
         inquire( FILE = 'flushvtk', EXIST = mustflushVTK)
         pararflushing = pararflushing .or.  &
         (ceiling((time_end-time_desdelanzamiento)/60.0_RKIND) >= maxCPUtime)
#ifdef CompileWithMPI
         l_aux=dontwritevtk
         call MPI_AllReduce( l_aux, dontwritevtk, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         !
         l_aux=pararnoflushing
         call MPI_AllReduce( l_aux, pararnoflushing, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         !
         l_aux=pararflushing
         call MPI_AllReduce( l_aux, pararflushing, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         !
         l_aux=mustUnpack
         call MPI_AllReduce( l_aux, mustUnpack, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         !
         l_aux=mustPostprocess
         call MPI_AllReduce( l_aux, mustPostprocess, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         !
         l_aux=mustflushXdmf
         call MPI_AllReduce( l_aux, mustflushXdmf, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         !
         l_aux=mustflushVTK
         call MPI_AllReduce( l_aux, mustflushVTK, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         !
         l_aux=mustflushFIELDS
         call MPI_AllReduce( l_aux, mustflushFIELDS, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         !
         l_aux=mustflushdata
         call MPI_AllReduce( l_aux, mustflushdata, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         call MPI_Barrier(MPI_COMM_WORLD,ierr) !050619 cambiado subcomm a mpi_comm_world
#endif
         mustflushfields = pararflushing .or. mustflushfields
         parar = pararNOflushing .or. pararflushing
         mustflushdata = mustflushdata .or. mustflushfields.or.parar !data is enforced to flush if fields are flushed
         !
         !to prevent duplicate writes on resuming
         !--->
         lmaxval  (1:size) =  0.0_RKIND
         lmaxval_i(1:size) =  0
         lmaxval_j(1:size) =  0
         lmaxval_k(1:size) =  0
         lmaxval_x(1:size) =  0.0_RKIND
         lmaxval_y(1:size) =  0.0_RKIND
         lmaxval_z(1:size) =  0.0_RKIND
         lminval  (1:size) =  0.0_RKIND
         lminval_i(1:size) =  0
         lminval_j(1:size) =  0
         lminval_k(1:size) =  0
         lminval_x(1:size) =  1e+20
         lminval_y(1:size) =  1e+20
         lminval_z(1:size) =  1e+20
         !

         valor = 0.0_RKIND
         !--->
         ini_i = b%sweepSINPMLEx%XI - b%Ex%XI
         fin_i = b%sweepSINPMLEx%XE - b%Ex%XI
         ini_j = b%sweepSINPMLEx%YI - b%Ex%YI
         fin_j = b%sweepSINPMLEx%YE - b%Ex%YI
         ini_k = b%sweepSINPMLEx%ZI - b%Ex%ZI
         fin_k = b%sweepSINPMLEx%ZE - b%Ex%ZI

         do k = ini_k, fin_k
            do j = ini_j, fin_j
               do i = ini_i, fin_i
                  valor = valor + Ex( i, j, k) * Ex( i, j, k)
               enddo
            enddo
         enddo
         !--->
         ini_i = b%sweepSINPMLEy%XI - b%Ey%XI
         fin_i = b%sweepSINPMLEy%XE - b%Ey%XI
         ini_j = b%sweepSINPMLEy%YI - b%Ey%YI
         fin_j = b%sweepSINPMLEy%YE - b%Ey%YI
         ini_k = b%sweepSINPMLEy%ZI - b%Ey%ZI
         fin_k = b%sweepSINPMLEy%ZE - b%Ey%ZI
         do k = ini_k, fin_k
            do j = ini_j, fin_j
               do i = ini_i, fin_i
                  valor = valor + Ey( i, j, k) * Ey( i, j, k)
               enddo
            enddo
         enddo
         !--->
         ini_i = b%sweepSINPMLEz%XI - b%Ez%XI
         fin_i = b%sweepSINPMLEz%XE - b%Ez%XI
         ini_j = b%sweepSINPMLEz%YI - b%Ez%YI
         fin_j = b%sweepSINPMLEz%YE - b%Ez%YI
         ini_k = b%sweepSINPMLEz%ZI - b%Ez%ZI
         fin_k = b%sweepSINPMLEz%ZE - b%Ez%ZI
         do k = ini_k, fin_k
            do j = ini_j, fin_j
               do i = ini_i, fin_i
                  valor = valor + Ez( i, j, k) * Ez( i, j, k)
               enddo
            enddo
         enddo
         !
         !--->
         energy = valor !!! quitado 241018 para evitar pasar el eps0----> 0.5_RKIND * Eps0 * valor
         !--->
         energytotal = energy
#ifdef CompileWithMPI
         !print *,'layoutnumber+1,energy,energytotal pre',layoutnumber+1,energy,energytotal
         call MPI_AllReduce( energy, energyTotal, 1_4, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         !print *,'layoutnumber+1,energy,energytotal post ',layoutnumber+1,energy,energytotal
#else
         energytotal = energy
#endif

         ini_iboxsin = max(max(b%sweepSINPMLEx%XI - b%Ex%XI, b%sweepSINPMLEy%XI - b%Ey%XI),b%sweepSINPMLEz%XI - b%Ez%XI)
         fin_iboxsin = min(min(b%sweepSINPMLEx%XE - b%Ex%XI, b%sweepSINPMLEy%XE - b%Ey%XI),b%sweepSINPMLEz%XE - b%Ez%XI)
         ini_jboxsin = max(max(b%sweepSINPMLEx%YI - b%Ex%YI, b%sweepSINPMLEy%YI - b%Ey%YI),b%sweepSINPMLEz%YI - b%Ez%YI)
         fin_jboxsin = min(min(b%sweepSINPMLEx%YE - b%Ex%YI, b%sweepSINPMLEy%YE - b%Ey%YI),b%sweepSINPMLEz%YE - b%Ez%YI)
         ini_kboxsin = max(max(b%sweepSINPMLEx%ZI - b%Ex%ZI, b%sweepSINPMLEy%ZI - b%Ey%ZI),b%sweepSINPMLEz%ZI - b%Ez%ZI)
         fin_kboxsin = min(min(b%sweepSINPMLEx%ZE - b%Ex%ZI, b%sweepSINPMLEy%ZE - b%Ey%ZI),b%sweepSINPMLEz%ZE - b%Ez%ZI)
         ini_ibox = max(max(b%sweepEx%XI - b%Ex%XI, b%sweepEy%XI - b%Ey%XI),b%sweepEz%XI - b%Ez%XI)
         fin_ibox = min(min(b%sweepEx%XE - b%Ex%XI, b%sweepEy%XE - b%Ey%XI),b%sweepEz%XE - b%Ez%XI)
         ini_jbox = max(max(b%sweepEx%YI - b%Ex%YI, b%sweepEy%YI - b%Ey%YI),b%sweepEz%YI - b%Ez%YI)
         fin_jbox = min(min(b%sweepEx%YE - b%Ex%YI, b%sweepEy%YE - b%Ey%YI),b%sweepEz%YE - b%Ez%YI)
         ini_kbox = max(max(b%sweepEx%ZI - b%Ex%ZI, b%sweepEy%ZI - b%Ey%ZI),b%sweepEz%ZI - b%Ez%ZI)
         fin_kbox = min(min(b%sweepEx%ZE - b%Ex%ZI, b%sweepEy%ZE - b%Ey%ZI),b%sweepEz%ZE - b%Ez%ZI)
         do k = ini_kbox, fin_kbox
            do j = ini_jbox, fin_jbox
               do i = ini_ibox, fin_ibox
                  valor = sqrt(Ex( i, j, k) * Ex( i, j, k) + Ey( i, j, k) * Ey( i, j, k)+Ez( i, j, k) * Ez( i, j, k))
                  if (lmaxval  (layoutnumber+1)< valor) then
                     lmaxval  (layoutnumber+1)= valor
                     lmaxval_i(layoutnumber+1)=i+b%Hx%XI
                     lmaxval_j(layoutnumber+1)=j+b%Hy%YI
                     lmaxval_k(layoutnumber+1)=k+b%Hz%ZI
                     lmaxval_x(layoutnumber+1)=Punto%PhysCoor(iHx)%x(lmaxval_i(layoutnumber+1))
                     lmaxval_y(layoutnumber+1)=Punto%PhysCoor(iHy)%y(lmaxval_j(layoutnumber+1))
                     lmaxval_z(layoutnumber+1)=Punto%PhysCoor(iHz)%z(lmaxval_k(layoutnumber+1))
                  endif
                  if (lminval  (layoutnumber+1)> valor) then
                     lminval  (layoutnumber+1)= valor
                     lminval_i(layoutnumber+1)=i+b%Hx%XI
                     lminval_j(layoutnumber+1)=j+b%Hy%YI
                     lminval_k(layoutnumber+1)=k+b%Hz%ZI
                     lminval_x(layoutnumber+1)=Punto%PhysCoor(iHx)%x(lminval_i(layoutnumber+1))
                     lminval_y(layoutnumber+1)=Punto%PhysCoor(iHy)%y(lminval_j(layoutnumber+1))
                     lminval_z(layoutnumber+1)=Punto%PhysCoor(iHz)%z(lminval_k(layoutnumber+1))
                  endif
               enddo
            enddo
         enddo

         !
         NEWlmaxval  (1:size) =0.0_RKIND
         NEWlmaxval_i(1:size) =0
         NEWlmaxval_j(1:size) =0
         NEWlmaxval_k(1:size) =0
#ifdef CompileWithMPI
         call MPI_AllReduce( LMAXVAL, NEWlmaxval  , size, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LMAXVAL_i, NEWlmaxval_I, size, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LMAXVAL_j, NEWlmaxval_J, size, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LMAXVAL_k, NEWlmaxval_K, size, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LMAXVAL_x, NEWlmaxval_x, size, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LMAXVAL_y, NEWlmaxval_y, size, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LMAXVAL_z, NEWlmaxval_z, size, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
#else
         NEWlmaxval   = LMAXVAL
         NEWlmaxval_i = LMAXVAL_I
         NEWlmaxval_j = LMAXVAL_J
         NEWlmaxval_k = LMAXVAL_K
         NEWlmaxval_x = LMAXVAL_x
         NEWlmaxval_y = LMAXVAL_y
         NEWlmaxval_z = LMAXVAL_z
#endif
         qmaxval=0.0_RKIND
         qmaxval_i = 0
         qmaxval_j = 0
         qmaxval_k = 0
         qmaxval_x = 0.0_RKIND
         qmaxval_y = 0.0_RKIND
         qmaxval_z = 0.0_RKIND
         do i=1,size
            if (abs(NEWlmaxval(i)) > qmaxval) then
               qmaxval   = abs(NEWlmaxval(i))
               qmaxval_i = newlmaxval_i(i)
               qmaxval_j = newlmaxval_j(i)
               qmaxval_k = newlmaxval_k(i)
               qmaxval_x = newlmaxval_x(i)
               qmaxval_y = newlmaxval_y(i)
               qmaxval_z = newlmaxval_z(i)
            endif
         end do

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         !
         NEWlminval  (1:size) =0.0_RKIND
         NEWlminval_i(1:size) =0
         NEWlminval_j(1:size) =0
         NEWlminval_k(1:size) =0
#ifdef CompileWithMPI
         call MPI_AllReduce( LminVAL, NEWlminval  , size, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LminVAL_i, NEWlminval_I, size, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LminVAL_j, NEWlminval_J, size, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LminVAL_k, NEWlminval_K, size, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LminVAL_x, NEWlminval_x, size, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LminVAL_y, NEWlminval_y, size, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( LminVAL_z, NEWlminval_z, size, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
#else
         NEWlminval   = LminVAL
         NEWlminval_i = LminVAL_I
         NEWlminval_j = LminVAL_J
         NEWlminval_k = LminVAL_K
         NEWlminval_x = LminVAL_x
         NEWlminval_y = LminVAL_y
         NEWlminval_z = LminVAL_z
#endif
         qminval=0.0_RKIND
         qminval_i = 0
         qminval_j = 0
         qminval_k = 0
         qminval_x = 0.0_RKIND
         qminval_y = 0.0_RKIND
         qminval_z = 0.0_RKIND
         do i=1,size
            if (abs(NEWlminval(i)) > qminval) then
               qminval   = abs(NEWlminval(i))
               qminval_i = newlminval_i(i)
               qminval_j = newlminval_j(i)
               qminval_k = newlminval_k(i)
               qminval_x = newlminval_x(i)
               qminval_y = newlminval_y(i)
               qminval_z = newlminval_z(i)
            endif
         end do

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

         !!!!!!!!!!!!!!!!!!!!!!!
         !escritura del fichero snap a voluntad o cuando se pase un umbral, cada minuto
         IF (layoutnumber == 0) THEN 
            inquire( FILE = 'snap', EXIST = mustSnap)
            if (mustsnap) then
               !  Clear the flushing signaling file
               open (35,FILE='snap')
               read(35,*,end=1153,err=1153) snapLevel
               read(35,*,end=1153,err=1153) snapStep
               read(35,*,end=1153,err=1153) snapHowMany
1153           close (35)
               call erasesignalingfiles(simu_devia)
            endif
         endif
#ifdef CompileWithMPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         call MPI_BCast( mustSnap, 1_4, MPI_LOGICAL, 0_4, MPI_COMM_WORLD, ierr)
         call MPI_BCast( snapLevel, 1_4, REALSIZE, 0_4, MPI_COMM_WORLD, ierr)
         call MPI_BCast( snapStep, 1_4, MPI_INTEGER, 0_4, MPI_COMM_WORLD, ierr)
         call MPI_BCast( snapHowMany, 1_4, MPI_INTEGER, 0_4, MPI_COMM_WORLD, ierr)
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

         if ((mustsnap.and.(lmaxval  (layoutnumber+1)> snapLevel)).or.(countersnap > 0)) then
            countersnap=countersnap + 1
            !
            dimxsnap=int((fin_ibox-ini_ibox)/snapstep) + 1
            dimysnap=int((fin_jbox-ini_jbox)/snapstep) + 1
            dimzsnap=int((fin_kbox-ini_kbox)/snapstep) + 1
            if (.not.allocated(snap)) allocate (snap(ini_ibox:ini_ibox+dimxsnap, &
                                                     ini_jbox:ini_jbox+dimysnap,ini_kbox:ini_kbox+dimzsnap,1))
            snap=0.0_RKIND
            !
            !!!!             k = ini_kbox - snapStep
            !!!!             do while (k < fin_kbox )
            !!!!                k = min (k  + snapStep , fin_kbox)
            !!!!                j = ini_jbox - snapStep
            !!!!                do while (j < fin_jbox )
            !!!!                    j = min(j  + snapStep , fin_jbox)
            !!!!                    i = ini_ibox - snapStep
            !!!!                    do while (i < fin_ibox )
            !!!!                        i = min (i  + snapStep , fin_ibox)
            !!!!                        veces=0
            !!!!                        valor=0.0_RKIND
            !!!!                        do k1=0,snapstep-1
            !!!!                        do j1=0,snapstep-1
            !!!!                        do i1=0,snapstep-1
            !!!!                        if ((i+i1 <= fin_ibox).and.(j+j1 <= fin_jbox).and.(k+k1 <= fin_kbox)) then
            !!!!                            valor = valor+sqrt(Ex(i+i1, j+j1, k+k1) * Ex( i+i1, j+j1, k+k1) + &
            !!!!                                               Ey( i+i1, j+j1, k+k1) * Ey(i+i1, j+j1, k+k1)+ &
            !!!!                                            Ez(i+i1, j+j1, k+k1) * Ez( i+i1, j+j1, k+k1))
            !!!!                            veces=veces+1
            !!!!                        endif
            !!!!                        end do
            !!!!                        end do
            !!!!                        end do
            !!!!                        snap(ini_ibox+int((i-ini_ibox)/snapstep),ini_jbox+int((j-ini_jbox)/snapstep), &
            !!!!                            ini_kbox+int((k-ini_kbox)/snapstep),1) = valor/veces
            !!!!                    enddo
            !!!!                enddo
            !!!!             enddo

            do k = ini_kbox, fin_kbox , snapStep
               do  j = ini_jbox, fin_jbox , snapStep
                  do i = ini_ibox, fin_ibox , snapStep
                     veces=0
                     valor=0.0_RKIND
                     do k1=0,snapstep-1
                        do j1=0,snapstep-1
                           do i1=0,snapstep-1
                              if ((i+i1 <= fin_ibox).and.(j+j1 <= fin_jbox).and.(k+k1 <= fin_kbox)) then
                                 valor = valor+sqrt(Ex(i+i1, j+j1, k+k1) * Ex( i+i1, j+j1, k+k1) + &
                                                    Ey( i+i1, j+j1, k+k1) * Ey(i+i1, j+j1, k+k1)+ &
                                 Ez(i+i1, j+j1, k+k1) * Ez( i+i1, j+j1, k+k1))
                                 veces=veces+1
                              endif
                           end do
                        end do
                     end do
                     snap(ini_ibox+int((i-ini_ibox)/snapstep),ini_jbox+int((j-ini_jbox)/snapstep), &
                          ini_kbox+int((k-ini_kbox)/snapstep),1) = valor/veces
                  enddo
               enddo
            enddo


            write(chinstant,'(i8)') n
            write(minmax,'(a,e15.4e3,a,e15.4e3,a)')  '_',lminval(layoutnumber+1),'_',lmaxval(layoutnumber+1),'_'
            fichsnap=trim(adjustl(nEntradaRoot))//'_snap_'//trim(adjustl(chinstant))//'_'// &
            trim(adjustl(whoamishort))

#ifdef CompileWithHDF
            ficherito=trim(adjustl(fichsnap))//'.h5' 
            call openclosedelete(ficherito)
            !
            call write_xdmfsnap(n,fichsnap,ini_ibox + b%Ex%XI , ini_ibox + dimxsnap + b%Ex%XI , &
            ini_jbox + b%Ex%YI , ini_jbox + dimysnap + b%Ex%YI , &
            ini_kbox + b%Ex%ZI , ini_kbox + dimzsnap + b%Ex%ZI , snap)
#endif
            !             open (35,file=trim(adjustl(fichsnap))//'.bin')
            !             write (35,*) '!END'
            !             close (35,status='delete')
            !             open (35,file=trim(adjustl(fichsnap))//'.bin',form='unformatted',status='new',action='write')
            !             write (35) n,lminval  (layoutnumber+1),lmaxval  (layoutnumber+1)
            !             write (35) ini_ibox,fin_ibox,ini_jbox,fin_jbox,ini_kbox,fin_kbox
            !             write (35) ini_iboxsin,fin_iboxsin,ini_jboxsin,fin_jboxsin,ini_kboxsin,fin_kboxsin
            !             write (35) (Punto%PhysCoor(iHx)%x(i+b%Hx%XI),i = ini_iboxsin, fin_iboxsin)
            !             write (35) (Punto%PhysCoor(iHy)%y(j+b%Hy%YI),j = ini_jboxsin, fin_jboxsin)
            !             write (35) (Punto%PhysCoor(iHz)%z(k+b%Hz%ZI),k = ini_kboxsin, fin_kboxsin)
            !#ifndef CompileWithHDF
            !             do k = ini_kbox, fin_kbox
            !                do j = ini_jbox, fin_jbox
            !                       write (35) (snap(i,j,k,1),i = ini_ibox, fin_ibox)
            !                enddo
            !             enddo
            !#endif
            !             close (35)

            write(dubuf,*)     trim(adjustl(whoami))//' Written Snap file at n= ',n,' max field over ',maxval(snap),'>',snapLevel
            deallocate (snap)
            call print11(layoutnumber,dubuf,.true.)
            if (countersnap >= snapHowMany) then
               mustsnap=.false.
               countersnap=0
            endif
         endif

#ifdef CompileWithMPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr) !TODOS STOCH O NO 060619
#endif
         !!!!!!!!!!!!!!!!!!!!!!!

         !
         !
         IF (layoutnumber == 0) THEN
            !
            write(dubuf,*) SEPARADOR,trim(adjustl(nentradaroot)),separador
            call print11(layoutnumber,dubuf)
            write(dubuf,*) 'Switches: '//trim(adjustl(opcionestotales))
            call print11(layoutnumber,dubuf)
            !if (size/=1) then
                write(dubuf,*) 'MPI Processes: ',size
                call print11(layoutnumber,dubuf)
            !endif
            !
            write(dubuf,*) 'Date/Time ', time_out2%fecha( 7: 8),'/', &
            time_out2%fecha( 5: 6),'/', &
            &               time_out2%fecha(1:4),'   ',time_out2%hora( 1: 2), &
            ':',time_out2%hora( 3: 4),':',time_out2%hora( 5: 6)
            call print11(layoutnumber,dubuf)
            !
            write(dubuf,*) 'Simulated:',n         ,'/',finaltimestep,' steps'
            call print11(layoutnumber,dubuf)
            !
            if (permitscaling) then
                write(dubuf,'(a,e19.9e3,a,e19.9e3,a,e19.9e3)') 'Time= ',sgg%tiempo(n),', dt0 (original)= ',dt0, &
                                                               ', dt(pscaled)= ',sgg%dt
            else
                write(dubuf,'(a,e19.9e3,a,e19.9e3,a,e19.9e3)') 'Time= ',sgg%tiempo(n),', dt0 = ',sgg%dt
            endif
            call print11(layoutnumber,dubuf)
            !
            if (energytotal > oldenergytotal) then
               write(dubuf,*)     'Total Energy (inc) :',energytotal
               if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
               oldenergytotal=energytotal
            else
               write(dubuf,*)     'Total Energy (dec) :',energytotal
               if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
               oldenergytotal=energytotal
            endif
            call print11(layoutnumber,dubuf)
            !
            if (qmaxval_x<-1e19) then
               write (dondex,'(a)') ' PML '
            else
               write (dondex,'(e15.4e3)') qmaxval_x
            endif
            if (qmaxval_y<-1e19) then
               write (dondey,'(a)') ' PML '
            else
               write (dondey,'(e15.4e3)') qmaxval_y
            endif
            if (qmaxval_z<-1e19) then
               write (dondez,'(a)') ' PML '
            else
               write (dondez,'(e15.4e3)') qmaxval_z
            endif

            write(dubuf,'(a,e13.4e3,a,3i5,a)')     'Max field: ',qmaxval,' at (',qmaxval_i,qmaxval_j,qmaxval_k,')=('// &
            ' '//trim(adjustl(dondex))//','//' '//trim(adjustl(dondey))//','//' '//trim(adjustl(dondez))//')'

            if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
            call print11(layoutnumber,dubuf)
            do i=1,size
               if (newlmaxval_x(i)<-1e19) then
                  write (dondex,'(a)') ' PML '
               else
                  write (dondex,'(e15.4e3)') newlmaxval_x(i)
               endif
               if (newlmaxval_y(i)<-1e19) then
                  write (dondey,'(a)') ' PML '
               else
                  write (dondey,'(e15.4e3)') newlmaxval_y(i)
               endif
               if (newlmaxval_z(i)<-1e19) then
                  write (dondez,'(a)') ' PML '
               else
                  write (dondez,'(e15.4e3)') newlmaxval_z(i)
               endif

               write(dubuf,'(a,i5,e15.4e3,a,e15.4e3,a,3i5,a)')  'Max field slice: ',i,NEWlmaxval(i),'/',maxSourceValue, &
               ' at (',newlmaxval_i(i),newlmaxval_j(i),newlmaxval_k(i),')=('// &
            ' '//trim(adjustl(dondex))//','//' '//trim(adjustl(dondey))//','//' '//trim(adjustl(dondez))//')'
               ! call print11(layoutnumber,dubuf) !comentado para que la salida sea menos verbose
            end do
            !

            !
            write(dubuf,*)     'Mins. since start  :',ceiling((time_end-time_desdelanzamiento)/60.0_RKIND)
            call print11(layoutnumber,dubuf)
            !
            write(dubuf,*)     'Mins. until end    :', &
            min(ceiling(((finaltimestep-n)*megaceldastotales)/speedGlobAvg/60.0_RKIND), &
            maxCPUtime-ceiling((time_end-time_desdelanzamiento)/60.0_RKIND))
            call print11(layoutnumber,dubuf)
            !
            if (everflushed) then
               write(dubuf,*) 'Mins. past flushing:',ceiling((time_end-time_begin2)/60.0_RKIND)
               call print11(layoutnumber,dubuf)
            else
               write(dubuf,*) 'Never flushed resuming fields.'
               call print11(layoutnumber,dubuf)
            endif
            if (flushsecondsFields /= 0) then
               write(dubuf,*) 'Mins. next flushing:',min(ceiling((flushsecondsFIELDS-(time_end-time_begin2))/60.0_RKIND), &
               ceiling(((finaltimestep-n)*megaceldastotales)/speedGlobAvg/60.0_RKIND))
               call print11(layoutnumber,dubuf)
            else
               if (maxCPUtime == topCPUtime) then
                  write(dubuf,*) 'Will Never flush resuming fields.'
                  call print11(layoutnumber,dubuf)
               else
                  write(dubuf,*) 'Flushing of restarting DATA at the end.'
                  call print11(layoutnumber,dubuf)
               endif
            endif
            !
            write(dubuf,*) 'Next info at step: ',n_info
            call print11(layoutnumber,dubuf)
            !
            write(dubuf,*) 'Total Mcells:',megaceldastotales
            call print11(layoutnumber,dubuf)
            !
            if (reportedinstant < n) then
               write(dubuf,*) 'Mcells/sec  :',speedGlobInst,'(',reportedinstant, ' to ',N,')'
               if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
               call print11(layoutnumber,dubuf)
            endif
            !
            write(dubuf,*) 'Mcells/sec  :',speedGlobAvg ,'(',INITIALtimeSTEP, ' to ',N,')'
            if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
            call print11(layoutnumber,dubuf)
            !
            write(dubuf,*) SEPARADOR//separador//separador
            call print11(layoutnumber,dubuf)
            !
            write (10,*)      sgg%tiempo(n),energytotal
            !write(67,'(i5)') nint(100.0_RKIND * n/finaltimestep) !percentage
            call flush(11)
            call flush(10)
         endif
         !
#ifdef CompileWithMPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         call get_secnds(time_out2)
         time_begin = time_out2%segundos !restart timing
         reportedinstant=n+1
      endif !every reporting seconds
      !

      !stop if this probe blows up
      stoponNaN_aux=.false.
!!!#ifdef ARCHITECTURE_SUN
!!#ifdef CompileWithReal4
!!      if (.false.) then
!!#else
!!      if (.false.) then
!!#endif
!!!#else
!!!#ifdef PreventCrayBug
!!#ifdef CompileWithReal4
!!      if (IsNaNf(energy)) then
!!#else
!!      if (IsNaNd(energy)) then
!!#endif
!!!#else
     !! if (IsNaN (energy)) then !quitado a mano para que PGI no se queje a 150623 !fm
!!!#endif
!!!#endif
         !
     !!    stoponNaN_aux=.true.
     !! endif
#ifdef CompileWithMPI
      call MPI_AllReduce( stoponNaN_aux, stoponNaN, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#else
      stoponNaN = stoponNaN_aux
#endif
      if (stoponNaN) then
         IF (layoutnumber == 0) THEN
            write(dubuf,*) 'ERROR, ABORTING: UNSTABILITIES. Possible reasons and fixes: '
            call print11(layoutnumber,dubuf)
            write(dubuf,*) '     In case of single wires: reduce WIRE radii and/or reduce sgg%dt'
            call print11(layoutnumber,dubuf)
            write(dubuf,*) '     In case of surface IBCs: reduce -att factor'
            call print11(layoutnumber,dubuf)
         endif
#ifdef CompileWithMPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         parar=.true.
         !            call StopOnError(layoutnumber,size,' Aborting')
      endif
      !
      l_aux = ( ((time_end-time_begin2) > flushsecondsFIELDS).AND. &
      (flushsecondsFIELDS/=0                       )).or. &
      mustflushFIELDS
#ifdef CompileWithMPI
      !print *,'layoutnumber+1,l_aux, hay_flushFIELDSl pre',layoutnumber+1,l_aux, hay_flushFIELDS
      call MPI_AllReduce( l_aux, hay_flushFIELDS, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      !print *,'layoutnumber+1,l_aux, hay_flushFIELDSl post',layoutnumber+1,l_aux, hay_flushFIELDS
#else
      hay_flushFIELDS = l_aux
#endif
      l_aux = ( ((time_end-time_begin3) > flushsecondsDATA).AND. &
      (flushsecondsDATA/=0                       )).or. &
      mustflushDATA
#ifdef CompileWithMPI
      !print *,'layoutnumber+1,l_aux, hay_flushDATA pre',layoutnumber+1,l_aux, hay_flushDATA
      call MPI_AllReduce( l_aux, hay_flushDATA, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
      !print *,'layoutnumber+1,l_aux, hay_flushDATA post',layoutnumber+1,l_aux, hay_flushDATA
#else
      hay_flushDATA = l_aux
#endif
      if (hay_flushFIELDS) then
         mustflushFIELDS=.false.
#ifdef CompileWithMPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         call get_secnds(time_out2)
         time_begin2=time_out2%segundos
         performflushFIELDS=.true.
         !  Clear the flushing signaling file
         IF (layoutnumber == 0) THEN !only the master proc mush erase this
             call erasesignalingfiles(simu_devia)
         endif
      endif
      if (hay_flushDATA) then
         !
         mustflushDATA=.false.
#ifdef CompileWithMPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         call get_secnds(time_out2)
         time_begin3=time_out2%segundos
         performflushDATA=.true.
         !  Clear the flushing signaling file
         IF (layoutnumber == 0) THEN !only the master proc mush erase this
             call erasesignalingfiles(simu_devia)
         endif
      endif
      if (mustunpack) then
         !
         mustunpack=.false.
         performunpack=.true.
         !  Clear the flushing signaling file
         IF (layoutnumber == 0) THEN !only the master proc mush erase this
             call erasesignalingfiles(simu_devia)
         endif
      endif
      if (mustpostprocess) then
         !
         mustpostprocess=.false.
         performpostprocess=.true.
         !  Clear the flushing signaling file
         IF (layoutnumber == 0) THEN !only the master proc mush erase this
             call erasesignalingfiles(simu_devia)
         endif
      endif
      if (mustflushXdmf) then
         !
         mustflushXdmf=.false.
         performflushXdmf=.true.
         !  Clear the flushing signaling file
         IF (layoutnumber == 0) THEN !only the master proc mush erase this
             call erasesignalingfiles(simu_devia)
         endif
      endif
      if (mustflushVTK) then
         !
         mustflushVTK=.false.
         performflushVTK=.true.
         IF (layoutnumber == 0) THEN !only the master proc mush erase this
             call erasesignalingfiles(simu_devia)
         endif
      endif
      !---------------------------> acaba Timing <----------------------------------------------------
      return
   endsubroutine Timing





   SUBROUTINE INITWARNINGFILE(layoutnumber,size,nEntradaRoot,verbosete,ignoreErrors1)
      character(len=*) :: nEntradaRoot
      integer (kind=4), intent(in) :: layoutnumber,size
      !file management
      character (LEN=BUFSIZE) ::  whoamishort
#ifdef CompileWithMPI
      integer(kind=MPI_OFFSET_KIND) disp
      integer (kind=4) :: ierr
#endif
      logical verbosete,ignoreerrors1       , itsopen2
      integer :: my_iostat
      character (LEN=BUFSIZE) :: ficherito
      verbose=verbosete

      ignoreerrors=ignoreerrors1

      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '
      write(whoamishort,'(i5)') layoutnumber+1

      IF (layoutnumber == 0) THEN          

      !!!inquire(unit=17, opened=itsopen2)
      !!!if (itsopen2) print *,'----------->17 open!!!'
        ficherito=trim(adjustl(nEntradaRoot))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
        call openclosedelete(ficherito)
      endif

      !!!#ifdef CompileWithMPI
      !!!IF (SIZE/=0) THEN
      !!!    call MPI_Barrier(SUBCOMM_MPI,ierr)
      !!!    call MPI_FILE_open (SUBCOMM_MPI, trim(adjustl(nEntradaRoot))//'_tmpWarnings.txt', &
      !!!                           MPI_MODE_WRONLY + MPI_MODE_CREATE, &
      !!!                           MPI_INFO_NULL, thefile, ierr)
      !!!    disp = (layoutnumber+1) * BUFSIZE * maxmessages !no creo que se den mas de 2000 mensajes por layout
      !!!
      !!!    call MPI_FILE_SET_VIEW(thefile, disp, MPI_CHARACTER, &
      !!!                               MPI_CHARACTER, 'native', &
      !!!                               MPI_INFO_NULL, ierr)
      !!!ELSE
      !!!    open (17,file=trim(adjustl(nEntradaRoot))//'_tmpWarnings.txt',form='formatted')
      !!!ENDIF
      !!!#else
      
      inquire(unit=17, opened=itsopen2)
      !!!if (itsopen2) print *,'----------->17 open!!!'
      ficherito=trim(adjustl(nEntradaRoot))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
      call opensolo(17,ficherito)
      !!!#endif

      warningfileIsOpen=.true.
      warningfile=nEntradaRoot
      fatalerror = .false.
      CONTADORDEMENSAJES=0

   END SUBROUTINE INITWARNINGFILE


   subroutine WarnErrReport(bufff,error)
      !
      logical :: itsopen
      logical, optional :: error
#ifdef CompileWithMPI
      integer (kind=4) :: ierr
#endif

      character(len=*), intent(in) :: bufff
      character(len=BUFSIZE) :: buff2,buff3


      if (present(error)) fatalerror=error .or. fatalerror

      !

      buff3=trim(adjustl(whoami))//' '//trim(adjustl(bufff))
      call trimnullchar(buff3)

      !buff2=CHAR(13)//CHAR(10)//trim(adjustl(bufff(1:bufsize-24)))//CHAR(13)//CHAR(10)
      buff2=CHAR(13)//CHAR(10)//trim(adjustl(buff3))
      call trimnullchar(buff2)
      !!!#ifdef CompileWithMPI
      !!!    CONTADORDEMENSAJES = CONTADORDEMENSAJES +1
      !!!    IF (CONTADORDEMENSAJES > maxmessages) call StopOnError(0,0,'ERROR: Relaunch with -maxmessages ', &
      !!!               CONTADORDEMENSAJES*10 )
      !!!    call MPI_FILE_WRITE(thefile, buff2  , BUFSIZE, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      !!!#else
      write (17,'(a)',err=154) trim(adjustl(buff3))
      !!!#endif
      !
      !!! elimino lo de sacarlo por pantalla 06/03/15
      !!!!!! if ((verbose).or.((present(error)).and.(error))) write (*,'(a)')  trim(adjustl(buff3))
      goto 155
154   inquire(unit=17, opened=itsopen) 
      print *,itsopen,'- Cannot write into warning file the message: ',trim(adjustl(buff3))
155   return
   end subroutine WarnErrReport



   SUBROUTINE CLOSEWARNINGFILE(layoutnumber,size,fatalerror_final,stoch_undivided,simu_devia)
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4) :: ierr,posic,i
      character (len=BUFSIZE) :: buf2
      character (LEN=BUFSIZE)     ::  dubuf
      logical :: fatalerror_final , lexis,stoch_undivided,simu_devia        , itsopen2
      character( LEN=BUFSIZE)  ::  whoamishort,whoami,chinstant
      integer :: my_iostat,file87
      character (LEN=BUFSIZE) :: ficherito

      if (.not.WarningFileIsOpen) return

      !!!#ifdef CompileWithMPI
      !!!IF (SIZE/=0) THEN
      !!!    call MPI_FILE_close (thefile, ierr)
      !!!ELSE
      !!!    close (17)
      !!!ENDIF
      !!!#else
      close (17)
      !!!#endif

#ifdef CompileWithMPI
      !wait until everything is closed
      CALL MPI_Barrier (SUBCOMM_MPI, ierr)
#endif

      !arregla los NUL
      IF ((layoutnumber==0).or.((layoutnumber == size/2).and.stoch_undivided)) THEN
         open (88,file=trim(adjustl(WarningFile))//'_Warnings.txt',form='formatted')
         posic=0
         do i=0,size-1
            if (stoch_undivided) then
                write(whoamishort,'(i5)') i+1
            else
               if (simu_devia) then
                   write(whoamishort,'(i5)') size+i+1
               else
                   write(whoamishort,'(i5)') i+1
               endif
            endif 
            inquire(file=trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt',exist=lexis)
            if (lexis) then         
      !!!inquire(unit=87, opened=itsopen2)
      !!!if (itsopen2) print *,'----------->87 open!!!'
               ficherito=trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
               call opensolo(87,ficherito)
               !
875            read(87,'(a)',end=876,err=876) buf2
               call trimnullchar(buf2)
               buf2=trim(adjustl(buf2))
               if ((buf2(1:1) /= ' ').and.(buf2(1:1) /=char(0))) then
                  write(88,'(a)') trim(adjustl(buf2))
                  posic=posic+1
               endif
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
               ficherito=trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
               call openclosedelete(ficherito)
#endif
            endif
         end do

         !
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(layoutnumber,dubuf)
         write(dubuf,*) 'Closing warning file. Number of messages: ',posic
         call print11(layoutnumber,dubuf)
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(layoutnumber,dubuf)

         close (88)
      endif

      warningfileIsOpen=.false.
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
      call MPI_AllReduce( fatalerror, fatalerror_final, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
#else
      fatalerror_final = fatalerror
#endif

      if (fatalerror_final.and.(ignoreErrors)) then
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(layoutnumber,dubuf)
         write(dubuf,*) 'There are ERRORS: The simulation will CONTINUE but Revise *Warnings file '
         call print11(layoutnumber,dubuf)
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(layoutnumber,dubuf)
      endif

      fatalerror_final=(fatalerror_final).and.(.not.ignoreErrors)

      RETURN
   END SUBROUTINE

   subroutine trimnullchar(string)
      character (len=*), intent(inout) :: string
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
      character (len=*), intent( IN) :: message
      integer (kind=4), intent( IN)  ::  layoutnumber
      logical :: soyImpresor,forceprint
      logical , optional  ::  forceprint2

      forceprint=.false.
      if (present(forceprint2)) forceprint=forceprint2

      soyImpresor=((layoutnumber == 0).or.forceprint).and.printea
      if (message(1:1)=='&') then !respeta los espacios, el & es un espacio en realidad
#ifndef NoVerbose
        ! if (soyImpresor) write (*,'(a)') ' '//trim(message(2:))
         if (soyImpresor) write (*,'(a)') trim(message(2:))
#endif
        ! IF (layoutnumber == 0)   write (11,'(a)',err=111) ' '//trim(message(2:))
         IF (layoutnumber == 0)   write (11,'(a)',err=111) trim(message(2:))
      else !ajusta a izquierda sin respetar espacios
#ifndef NoVerbose
         if (soyImpresor) write (*,'(a)') trim(adjustl(message))
#endif
         IF (layoutnumber == 0)    then
             write (11,'(a)',err=111) trim(adjustl(message))
           !!!  print *,'----11--> ', trim(adjustl(message))
         endif
         
      endif
      goto 112
111   continue 
      !fort.11 a veces lo intentan escribir 2 a la vez de los que dan fallos en writing restarting fields. 
      !asi que ignora y continua
112   return

   end subroutine

   !part of the DXF



   !!!SUBROUTINE INITdxfFILE(layoutnumber,size,nEntradaRoot)
   !!!character(len=*) :: nEntradaRoot
   !!!integer (kind=4) :: layoutnumber,size
   !!!!file management
   !!!character (LEN=BUFSIZE) whoamishort
   !!!#ifdef CompileWithMPI
   !!!integer(kind=MPI_OFFSET_KIND) disp
   !!!integer (kind=4) :: ierr
   !!!#endif
   !!!write(whoamishort,'(i5)') layoutnumber+1
   !!!mynEntradaRoot=trim(adjustl(nentradaRoot))
   !!!dxfFileIsOpen=.true.
   !!!
   !!!IF (layoutnumber == 0) THEN
   !!!    open (87,file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf',form='formatted')
   !!!    write(87,*) '!END'
   !!!    close (87,status='delete')
   !!!endif
   !!!
   !!!!!!#ifdef CompileWithMPI
   !!!!!!call MPI_Barrier(SUBCOMM_MPI,ierr)
   !!!!!!call MPI_FILE_open (SUBCOMM_MPI, trim(adjustl(mynEntradaRoot))//'.tmpdxf', &
   !!!!!!                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
   !!!!!!                       MPI_INFO_NULL, thefile2, ierr)
   !!!!!!disp = (layoutnumber+1) * dxflinesize * maxdxf !no creo que se den mas de 2000 mensajes por layout
   !!!!!!
   !!!!!!call MPI_FILE_SET_VIEW(thefile2, disp, MPI_CHARACTER, &
   !!!!!!                           MPI_CHARACTER, 'native', &
   !!!!!!                           MPI_INFO_NULL, ierr)
   !!!!!!#else
   !!!    open (97,file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf',form='formatted')
   !!!!!!#endif
   !!!
   !!!IF (layoutnumber == 0) THEN
   !!!    write(dxfbuff,'(a)') '999'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)') 'SEMBA_FDTD'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '0'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'SECTION'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '2'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'HEADER'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '9'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '$ACADVER'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '1'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'AC1009'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '9'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '$SHADEDGE'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '70'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '1'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '0'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'ENDSEC'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '0'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'SECTION'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '2'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  'ENTITIES'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!    write(dxfbuff,'(a)')  '0'
   !!!    CALL DXFWRITE(DXFBUFF)
   !!!endif
   !!!
   !!!END SUBROUTINE INITdxfFILE
   !!!
   !!!
   !!!subroutine dxfwrite(bufff)
   !!!!
   !!!#ifdef CompileWithMPI
   !!!integer (kind=4) :: ierr
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
   !!!SUBROUTINE CLOSEdxfFILE(layoutnumber,size)
   !!!integer (kind=4), intent(in) :: layoutnumber,size
   !!!integer (kind=4) :: ierr,i
   !!!integer (kind=8) :: posic
   !!!character (len=dxflinesize) :: buf2
   !!!character (LEN=BUFSIZE)     ::  dubuf
   !!!character (LEN=BUFSIZE) ::  whoamishort
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
   !!!!arregla los NUL
   !!!IF (layoutnumber == 0) THEN
   !!!    open (988,file=trim(adjustl(mynEntradaRoot))//'.dxf',form='formatted')
   !!!    posic=0
   !!!    do i=0,size-1
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
   !!!            endif
   !!!            goto 875
   !!!!
   !!!876         close (987)
   !!!!
   !!!            open (987,file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf')
   !!!            write(987,*) '!END'
   !!!            close (987,status='delete')
   !!!        endif
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
   !!!endif
   !!!
   !!!continue
   !!!
   !!!dxfFileIsOpen=.false.
   !!!RETURN
   !!!END SUBROUTINE

   !!!
   !!!subroutine writemmdxf(layoutnumber,sgg,sggMiHx,sggMiHy,sggMiHz)
   !!!integer i,j,k,layoutnumber
   !!!type (SGGFDTDINFO), intent(IN)     ::  sgg
   !!!INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES), intent(in)   :: &
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
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '8'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHx(i,j,k)+20
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '62'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHx(i,j,k)+20
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '10'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '20'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '30'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '11'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '21'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '31'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '12'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '22'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '32'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '13'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '23'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '33'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '0'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                endif
   !!!            end do
   !!!        end do
   !!!    end do
   !!!!
   !!!      Do k=sgg%SINPMLSweep(iHy)%ZI,sgg%SINPMLSweep(iHy)%ZE
   !!!          Do j=sgg%SINPMLSweep(iHy)%YI,sgg%SINPMLSweep(iHy)%YE
   !!!              Do i=sgg%SINPMLSweep(iHy)%XI,sgg%SINPMLSweep(iHy)%XE
   !!!                if ((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec))  then
   !!!                    write(dxfbuff,'(a)') '3DFACE'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '8'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHy(i,j,k)+20
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '62'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHy(i,j,k)+20
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '10'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '20'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '30'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '11'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '21'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '31'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '12'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '22'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '32'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '13'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '23'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '33'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '0'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                endif
   !!!            end do
   !!!        end do
   !!!    end do
   !!!!
   !!!      Do k=sgg%SINPMLSweep(iHz)%ZI,sgg%SINPMLSweep(iHz)%ZE
   !!!          Do j=sgg%SINPMLSweep(iHz)%YI,sgg%SINPMLSweep(iHz)%YE
   !!!              Do i=sgg%SINPMLSweep(iHz)%XI,sgg%SINPMLSweep(iHz)%XE
   !!!                if ((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec))  then
   !!!                    write(dxfbuff,'(a)') '3DFACE'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '8'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHz(i,j,k)+20
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '62'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') sggMiHz(i,j,k)+20
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '10'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '20'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '30'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '11'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '21'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '31'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '12'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '22'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '32'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '13'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') i
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '23'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') j+1
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(a)') '33'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    write(dxfbuff,'(i5)') k
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                    !
   !!!                    write(dxfbuff,'(a)') '0'
   !!!                    CALL DXFWRITE(DXFBUFF)
   !!!                endif
   !!!            end do
   !!!        end do
   !!!    end do
   !!!    return
   !!!end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   function openfile_mpi(layoutnumber,nombrefich) result (thefile8)

    logical :: borratedeunaputavez1,borratedeunaputavez2
    
      integer (kind=4)  ::   thefile8 !for mpi file management
#ifdef CompileWithMPI
      integer(kind=MPI_OFFSET_KIND) disp
      integer (kind=4)  ::  ierr !for mpi file management
#endif
      integer (kind=4), intent(in) :: layoutnumber
      character (LEN=BUFSIZE), intent(in) :: nombrefich
      character (LEN=BUFSIZE) whoamishort
      write(whoamishort,'(i5)') layoutnumber+1

!      IF (layoutnumber == 0) THEN
669      open (newunit=thefile8,file=trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp',err=667 )
         goto 668
667      print *,' '//trim(adjustl(whoamishort))//' ','--> no hay cojones con open para borrar unidad ',thefile8
         call sleep(2)
         goto 669
668     continue        
        call sleep(2) 
         write (thefile8,'(a)') '!END'
        call sleep(2)
         close (thefile8,status='delete')
!      endif
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      !!!#ifdef CompileWithMPI
      !!!call MPI_FILE_open (SUBCOMM_MPI, trim(adjustl(nombrefich))//'_tmp', &
      !!!                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
      !!!                       MPI_INFO_NULL, thefile8, ierr)
      !!!disp = (layoutnumber+1) * 1024 * 20000 !no creo que se den mas de 20000 escrituras por layout
      !!!
      !!!call MPI_FILE_SET_VIEW(thefile8, disp, MPI_CHARACTER, &
      !!!                           MPI_CHARACTER, 'native', &
      !!!                           MPI_INFO_NULL, ierr)
      !!!#else
!pgi    666     inquire (unit=thefile8,exist=borratedeunaputavez1)
        inquire (file=trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp' ,exist=borratedeunaputavez2)
!pgi    if (borratedeunaputavez1) then
!pgi         print *,' '//trim(adjustl(whoamishort))//' ','--> no hay cojones con inquire unidad ',thefile8
!pgi         call sleep(2)
!pgi         goto 666
!pgi    endif
        if (borratedeunaputavez2) then
        print *,' '//trim(adjustl(whoamishort))//' ','--> no hay cojones con inquire file fichero ',& 
                      trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp'
        call sleep(2)
        goto 669
      endif
      open (newunit=thefile8,file=trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp',err=767 )
      goto 768
767   print *,' '//trim(adjustl(whoamishort))//' ','--> no hay cojones con open definitivo unidad ',thefile8
      call sleep(2)
      goto 669
768   continue         
      !!!#endif
      return
   end function openfile_mpi



   subroutine writefile_mpi(layoutnumber, thefile8,buff2)

      integer (kind=4)  ::  thefile8
#ifdef CompileWithMPI
      integer (kind=4)  ::  ierr
#endif
      integer (kind=4), intent(in) :: layoutnumber
      character (LEN=BUFSIZE) :: buff2,buff3

      buff3=CHAR(13)//CHAR(10)//trim(adjustl(buff2))
      !!!#ifdef CompileWithMPI
      !!!                call MPI_FILE_WRITE(thefile8, buff3  , 1024_4, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      !!!#else
      write(thefile8,'(a)') buff2
      !!!#endif

   end subroutine writefile_mpi



   subroutine closefile_mpi(layoutnumber,size,nombrefich,thefile8)


      integer (kind=4)  ::  thefile8,thefile19
#ifdef CompileWithMPI
      integer (kind=4)  ::  ierr
#endif
      integer (kind=4), intent(in) :: layoutnumber,size
      character (LEN=BUFSIZE) :: buff2
      character (LEN=BUFSIZE) :: nombrefich
      integer (kind=4) :: conta,i
      character (LEN=BUFSIZE) ::  whoamishort
      logical lexis
      character (LEN=BUFSIZE) :: ficherito


      !!!!#ifdef CompileWithMPI
      !!!!call MPI_Barrier(SUBCOMM_MPI,ierr)
      !!!!call MPI_FILE_CLOSE(thefile8, ierr)
      !!!!#else
      close (thefile8)
      !!!#endif


#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      !arregla los NUL
      IF (layoutnumber == 0) THEN
         open (newunit=thefile8,file=trim(adjustl(nombrefich)),form='formatted' )
         do i=0,size-1
            write(whoamishort,'(i5)') i+1
            inquire(file=trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp',exist=lexis)
            if (lexis) then
               open (newunit=thefile19,file=trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp',form='formatted' )
               conta=0
875            read(thefile19,'(a)',end=876,err=876) buff2
               call trimnullchar(buff2)
               buff2=trim(adjustl(buff2))
               if ((buff2(1:1) /= ' ').and.(buff2(1:1) /=char(0))) then
                  write(thefile8,'(a)') trim(adjustl(buff2))
                  conta=conta+1
               endif
               goto 875
               !
876            continue
               close (thefile19)
               ficherito=trim(adjustl(nombrefich))//trim(adjustl(whoamishort))//'_tmp'
               call openclosedelete(ficherito)
            endif
         end do
         close (thefile8)

         !
      endif
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      !!!!end gnuplot
      return
   end subroutine closefile_mpi

   function creaPuntos(sgg)  result(punto) !crea coordenadas fisicas
      !
      type (SGGFDTDINFO), intent(INout)         ::  sgg
      type (coorsxyzP)  ::  Punto
      integer (Kind=4) :: i,j,k,field


      do field=iEx,iHz
         allocate (Punto%PhysCoor(field)%x(sgg%Sweep(field)%XI-1 : sgg%Sweep(field)%XE+1), &
         Punto%PhysCoor(field)%y(sgg%Sweep(field)%YI-1 : sgg%Sweep(field)%YE+1), &
         Punto%PhysCoor(field)%z(sgg%Sweep(field)%ZI-1 : sgg%Sweep(field)%ZE+1))
         Punto%PhysCoor(field)%x(sgg%Sweep(field)%XI-1 : sgg%Sweep(field)%XE+1) = -1e20
         Punto%PhysCoor(field)%y(sgg%Sweep(field)%YI-1 : sgg%Sweep(field)%YE+1) = -1e20
         Punto%PhysCoor(field)%z(sgg%Sweep(field)%ZI-1 : sgg%Sweep(field)%ZE+1) = -1e20
      end do

      !

      !
      field=iEx
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         Punto%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         Punto%PhysCoor(field)%y(j)=sgg%LineY(j)
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         Punto%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=iEy
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         Punto%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         Punto%PhysCoor(field)%y(j)=(sgg%LineY(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         Punto%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=iEz
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         Punto%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         Punto%PhysCoor(field)%y(j)=sgg%LineY(j)
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         Punto%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=iHx
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         Punto%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         Punto%PhysCoor(field)%y(j)=(sgg%LineY(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         Punto%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=iHy
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         Punto%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         Punto%PhysCoor(field)%y(j)=sgg%LineY(j)
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         Punto%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=iHz
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         Punto%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         Punto%PhysCoor(field)%y(j)=(sgg%LineY(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         Punto%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      !
      sgg%Punto = Punto
      return
   end function

   subroutine reportmedia(sgg)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      integer (kind=4) :: j
      character(len=BUFSIZE) :: buff

      DO j = 0, sgg%NumMedia
         WRITE (buff,*) '_____________________________'
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'MEDIO :  ', j
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Priority ', sgg%Med(j)%Priority
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Epr ', sgg%Med(j)%Epr
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Sigma ', sgg%Med(j)%Sigma
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Mur ', sgg%Med(j)%Mur
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'SigmaM ', sgg%Med(j)%SigmaM
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is PML ', sgg%Med(j)%Is%PML
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is PEC ', sgg%Med(j)%Is%PEC
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is ThinWIRE ', sgg%Med(j)%Is%ThinWire
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is SlantedWIRE ', sgg%Med(j)%Is%SlantedWire
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is EDispersive ', sgg%Med(j)%Is%EDispersive
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is EDispersiveaANIS', sgg%Med(j)%Is%EDispersiveANIS
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is MDispersive ', sgg%Med(j)%Is%MDispersive
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is MDispersiveANIS ', sgg%Med(j)%Is%MDispersiveANIS
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is ThinSlot ', sgg%Med(j)%Is%ThinSlot
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is SGBC ', sgg%Med(j)%Is%SGBC
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is Lossy ', sgg%Med(j)%Is%Lossy
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is Multiport ', sgg%Med(j)%Is%multiport
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is AnisMultiport ', sgg%Med(j)%Is%anismultiport
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is MultiportPadding ', sgg%Med(j)%Is%multiportpadding
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is Dielectric ', sgg%Med(j)%Is%dielectric
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is ThinSlot ', sgg%Med(j)%Is%ThinSlot
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is Anisotropic ', sgg%Med(j)%Is%Anisotropic
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is Needed ', sgg%Med(j)%Is%Needed
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is already_YEEadvanced_byconformal ', sgg%Med(j)%Is%already_YEEadvanced_byconformal
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Iss split_and_useless ', sgg%Med(j)%Is%split_and_useless
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is Volume ', sgg%Med(j)%Is%Volume
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is Surface ', sgg%Med(j)%Is%Surface
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) 'Is Line ', sgg%Med(j)%Is%Line
         call WarnErrReport(Trim(buff))
         WRITE (buff,*) '_____________________________'
         call WarnErrReport(Trim(buff))
      END DO
      return
   end subroutine reportmedia

   subroutine erasesignalingfiles(simu_devia)
      logical :: simu_devia
      CHARACTER (LEN=BUFSIZE) :: ficherito
      if (.not.simu_devia) then
          !!force erasing the signaling files
          ficherito='stop'
          call openclosedelete(ficherito)
          ficherito='stopflushing'
          call openclosedelete(ficherito)
          ficherito='flush'
          call openclosedelete(ficherito)
          ficherito='flushdata'
          call openclosedelete(ficherito)
          ficherito='unpack'
          call openclosedelete(ficherito)
          ficherito='postprocess'
          call openclosedelete(ficherito)
          ficherito='flushxdmf'
          call openclosedelete(ficherito)
          ficherito='flushvtk'
          call openclosedelete(ficherito)
          ficherito='snap'
          call openclosedelete(ficherito)         !
          !
          ficherito='stop_only'
          call openclosedelete(ficherito)
          ficherito='stopflushing_only'
          call openclosedelete(ficherito)
          ficherito='flush_only'
          call openclosedelete(ficherito)
          ficherito='flushdata_only'
          call openclosedelete(ficherito)
          !
          !
          ficherito='stop_dontwritevtk'
          call openclosedelete(ficherito)
          ficherito='stop_only_dontwritevtk'
          call openclosedelete(ficherito)
          ficherito='stopflushing_dontwritevtk'
          call openclosedelete(ficherito)
          ficherito='stopflushing_only_dontwritevtk'
          call openclosedelete(ficherito)
          ficherito='flush_dontwritevtk'
          call openclosedelete(ficherito)
          ficherito='flush_only_dontwritevtk'
          call openclosedelete(ficherito)
      !     
          ficherito='unpack'
          call openclosedelete(ficherito)
          ficherito='postprocess'
          call openclosedelete(ficherito)
          ficherito='flushxdmf'
          call openclosedelete(ficherito)
          ficherito='flushvtk'
          call openclosedelete(ficherito)
      endif
      return
        
   
  end subroutine erasesignalingfiles

   
   subroutine openclosedelete(ficherin)
      CHARACTER (LEN=*) :: ficherin
      integer (kind=4) :: my_iostat, myunit
      my_iostat = 0
4216  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.'
      open  (newunit=myunit,file=trim(adjustl(ficherin)),form='formatted',err=4216,iostat=my_iostat)
      write (myunit,*) '!END'
      close (myunit,status='delete')
      return
   end subroutine openclosedelete
   
   subroutine openclose(ficherin)
      CHARACTER (LEN=*) :: ficherin
      integer (kind=4) :: my_iostat, myunit
      my_iostat = 0
5216  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.'
      open  (newunit=myunit,file=trim(adjustl(ficherin)),form='formatted',err=5216,iostat=my_iostat)
      write (myunit,*) '!END'
      close (myunit)
      return
   end subroutine openclose
   
   
   subroutine opensolo(myunit,ficherin)
      CHARACTER (LEN=*) :: ficherin
      integer (kind=4) :: my_iostat, myunit
      my_iostat = 0
6216  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.'
      open  (unit=myunit,file=trim(adjustl(ficherin)),form='formatted',err=6216,iostat=my_iostat)
      return
   end subroutine opensolo
   
   subroutine closesolo(myunit)
      integer (kind=4) :: my_iostat, myunit
      my_iostat = 0
7216  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.'
      close  (unit=myunit,err=7216,iostat=my_iostat)
      return
   end subroutine closesolo
   
end module
