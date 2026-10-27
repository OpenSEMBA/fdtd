! Periodic progress reporting during time-stepping: prints speed, energy and
! max-field info, writes on-demand snapshots, and checks the signal files
! (stop, flush, flushdata, unpack, postprocess, ...) to decide what to perform.
module progress_m
   use FDETYPES_m
   use snapxdmf_m
   use Report_m, only: print11, get_secnds, time_t, erasesignalingfiles, &
                       openclosedelete, file10isopen, log_msg, LOG_ERROR, LOG_FATAL

   implicit none
   private

   integer(kind=4), parameter  :: REPORTINGSECONDS=60

   real(kind=8), save :: time_begin, time_end, time_begin2,time_begin3,time_begin_absolute, time_end2,time_launch
   real(kind=RKIND), save :: mega_cells,total_mega_cells,speedInst, speedGlobInst,speedAvg,speedGlobAvg
   real(kind=RKIND), save  :: energy,energyTotal,oldenergyTotal,snapLevel
   type(time_t), save  :: time_out2
   integer(kind=4), save   :: reportedinstant,snapStep,snapHowMany,countersnap
   type(coorsxyzP_t) , save  :: gridPoint

   public InitReporting, InitTiming, Timing

contains

   subroutine InitReporting(sgg,c)
      type(SGGFDTDINFO_t), intent(inout) :: sgg
      type(sim_control_t) :: c
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif
      logical  :: errnofile
      character(len=BUFSIZE) :: buff, whoami

      write(whoami,'(a,i5,a,i5,a)') '(',c%layoutnumber+1,'/',c%num_procs,') '
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      gridPoint=createPoints(sgg)

      if (c%layoutnumber == 0) then
         if (c%resume) then
            open (10,file=trim(adjustl(c%nInputRoot))//'_Energy.dat',form='formatted',position='append')
         else
            open (10,file=trim(adjustl(c%nInputRoot))//'_Energy.dat',form='formatted')
         end if
      end if
      file10isopen=.true.

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      call get_secnds(time_out2)
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      call print11(c%layoutnumber,SEPARADOR//SEPARADOR//SEPARADOR)

#ifndef CompileWithInt4
#define CompileWithInt4
#endif

      if (c%resume) then
         errnofile=.false.
         if (c%resume_fromold) then
            inquire(FILE=trim(adjustl(c%nresumeable2))//'.old', EXIST=errnofile)
         else
            inquire(FILE=trim(adjustl(c%nresumeable2)), EXIST=errnofile)
         end if
         if (.not.errnofile) then
            if (c%resume_fromold) then
               buff='FILE '//trim(adjustl(c%nresumeable2))//'.old DOES NOT EXIST'
               call log_msg(LOG_FATAL, 'PROGRESS', buff)
            else
               buff='FILE '//trim(adjustl(c%nresumeable2))//' DOES NOT EXIST'
               call log_msg(LOG_FATAL, 'PROGRESS', buff)
            end if
         end if
         call print11(c%layoutnumber,SEPARADOR//SEPARADOR//SEPARADOR)
         call print11(c%layoutnumber,' ')
         if (c%resume_fromold) then
            call print11(c%layoutnumber,'Reading resuming data from '//trim(adjustl(c%nresumeable2))//'.old etc.')
         else
            call print11(c%layoutnumber,'Reading resuming data from '//trim(adjustl(c%nresumeable2))//' etc.')
         end if
         call print11(c%layoutnumber,SEPARADOR//SEPARADOR//SEPARADOR)
      end if

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

   end subroutine InitReporting

   ! observationFlushSteps: step period of the probe buffer flush, 0 when there are no probes.
   subroutine InitTiming(sgg, control, initialtimestep, observationFlushSteps)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      type(sim_control_t) , intent(in):: control
      integer(kind=4), intent(in) :: initialtimestep
      integer(kind=4), intent(in) :: observationFlushSteps
      type(time_t) :: time_out2
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif
      character(len=BUFSIZE) :: whoami
      character(len=BUFSIZE) :: dubuf
      write(whoami,'(a,i5,a,i5,a)') '(',control%layoutnumber+1,'/',control%num_procs,') '

      time_launch=control%time_desdelanzamiento
      snapLevel=1.0e25_RKIND
      snapStep=1
      snapHowMany=1
      countersnap=0

      mega_cells = (1.0_RKIND*sgg%sweep(IEX)%ZE-1.0_RKIND*sgg%sweep(IEX)%ZI)* &
                   (1.0_RKIND*sgg%sweep(IEX)%YE-1.0_RKIND*sgg%sweep(IEX)%YI)* &
                   (1.0_RKIND*sgg%sweep(IEY)%XE-1.0_RKIND*sgg%sweep(IEY)%XI)/1.0e6_RKIND

#ifdef CompileWithMPI
      call MPI_BARRIER(SUBCOMM_MPI,ierr)
      call MPI_AllReduce(mega_cells, total_mega_cells, 1_4, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
#else
      total_mega_cells=mega_cells
#endif

      write(dubuf,*)  'Total Mcells: ',total_mega_cells
      call print11(control%layoutnumber,dubuf)

      if (control%flushsecondsFIELDS/=0) then
         write(dubuf,'(a,i0,a)')  'Restart FIELDS flush: every ',control%flushsecondsFIELDS/60,' minutes'
      else if (control%maxCPUtime == TOPCPUTIME) then
         dubuf = 'Restart FIELDS flush: none scheduled'
      else
         write(dubuf,'(a,i0,a)')  'Restart FIELDS flush: at the end, after ',control%maxCPUtime,' minutes'
      end if
      call print11(control%layoutnumber,dubuf)
      if (observationFlushSteps == 0) then
         dubuf = 'Observation DATA flush: no probes'
      else if (control%flushsecondsDATA == 0 .or. control%flushsecondsDATA >= TOPCPUTIME*60) then
         write(dubuf,'(a,i0,a)')  'Observation DATA flush: every ',observationFlushSteps,' steps; time-based: none'
      else
         write(dubuf,'(a,i0,a,i0,a)')  'Observation DATA flush: every ',observationFlushSteps, &
                                       ' steps; time-based: every ',control%flushsecondsDATA/60,' minutes'
      end if
      call print11(control%layoutnumber,dubuf)
      write(dubuf,'(a,i0,a)')  'Reporting simulation info every ',REPORTINGSECONDS/60,' minutes'
      call print11(control%layoutnumber,dubuf)

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      call get_secnds(time_out2)
      call print11(control%layoutnumber,SEPARADOR//SEPARADOR//SEPARADOR)
      write(dubuf,'(a,i7,a,e19.9e3,a,i9,a,e19.9e3)')  'Simulation from n=',initialtimestep,', t=',sgg%time(initialtimestep),&
                                                      ' to n=',control%finaltimestep,', t=',sgg%time(control%finaltimestep)
      call print11(control%layoutnumber,dubuf)
      write(dubuf,*)  'Date/time ', time_out2%fecha(7: 8),'/',time_out2%fecha(5: 6),'/',time_out2%fecha(1:4),'   ', &
                                    time_out2%hora(1: 2), ':',time_out2%hora(3: 4),':',time_out2%hora(5: 6)
      call print11(control%layoutnumber,dubuf)
      time_begin_absolute = time_out2%segundos
      time_begin = time_begin_absolute
      time_begin2 = time_begin_absolute
      time_begin3 = time_begin_absolute

      reportedinstant=initialtimestep

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

   end subroutine InitTiming

   subroutine Timing(sgg, b, n, n_info, layoutnumber, num_procs, maxCPUtime,flushsecondsFields, flushsecondsData, initialtimestep, &
   finaltimestep, perform, &
   stop_flag, forcetiming,Ex,Ey,Ez,everflushed, nInputRoot,maxSourceValue,all_options,simu_devia,dontwritevtk,permitscaling)

      logical :: simu_devia,dontwritevtk,stopdontwritevtk,stopflushingdontwritevtk,flushdontwritevtk,stoponlydontwritevtk
      type(SGGFDTDINFO_t), intent(in)              :: sgg
      type(bounds_t), intent(in) :: b
      character(len=BUFSIZE), intent(in) :: all_options
      integer(kind = 4), intent(in) :: layoutnumber, num_procs, n,maxCPUtime
      integer(kind = 4), intent(in) :: flushsecondsFields, flushsecondsData, initialtimestep, finaltimestep
      real(kind = RKIND), dimension(0: b%Ex%NX-1, 0: b%Ex%NY-1, 0: b%Ex%NZ-1), intent(in) :: Ex
      real(kind = RKIND), dimension(0: b%Ey%NX-1, 0: b%Ey%NY-1, 0: b%Ey%NZ-1), intent(in) :: Ey
      real(kind = RKIND), dimension(0: b%Ez%NX-1, 0: b%Ez%NY-1, 0: b%Ez%NZ-1), intent(in) :: Ez
      integer :: my_iostat
      logical, intent(in) :: forcetiming, everflushed,permitscaling
      character(len=BUFSIZE) :: snap_file,minmax,target_name
      character(len = *), intent(in) :: nInputRoot
      integer(kind=4), intent(inout) :: n_info
      logical, intent(inout) :: stop_flag
      type(perform_t), intent(out) :: perform
      real(kind=rKIND) :: valueField,maxSourceValue,LA,LV,LB
      logical  :: do_timing, l_aux, do_flushFIELDS, do_flushDATA, mustflushFIELDS, mustflushDATA,mustUnpack, &
                   mustPostprocess,mustflushXdmf , mustflushVTK ,   &
      stopFlushing, stopNoFlushing, stoponNaN , stoponNaN_aux,mustSnap,stop_only,stopflushing_only,flush_only,flushdata_only
      logical :: stopflushingonlydontwritevtk,flushonlydontwritevtk,flushdataonlydontwritevtk,flushdatadontwritevtk
      integer(kind = 4) :: in_aux, ini_i, fin_i, ini_j, fin_j, ini_k, fin_k, i, j, k
      character(LEN=BUFSIZE) :: whoamishort,whoami,chinstant
      character(len=BUFSIZE) :: dubuf
      character(len=BUFSIZE) :: where_x,where_y,where_z
      real(kind=rKIND), dimension(1:num_procs) :: NEWlmaxval,NEWlmaxval_x,NEWlmaxval_y,NEWlmaxval_z
      integer(kind = 4), dimension(1:num_procs) :: NEWlmaxval_i,NEWlmaxval_j,NEWlmaxval_k
      real(kind=rKIND), dimension(1:num_procs) :: lmaxval,lmaxval_x,lmaxval_y,lmaxval_z
      integer(kind = 4), dimension(1:num_procs) :: lmaxval_i,lmaxval_j,lmaxval_k
      real(kind=rKIND) :: qmaxval , qmaxval_x,qmaxval_y,qmaxval_z
      integer(kind = 4) :: qmaxval_i,qmaxval_j,qmaxval_k,thefilenoflu
      real(kind=rKIND), dimension(1:num_procs) :: NEWlminval,NEWlminval_x,NEWlminval_y,NEWlminval_z
      integer(kind = 4), dimension(1:num_procs) :: NEWlminval_i,NEWlminval_j,NEWlminval_k
      real(kind=rKIND), dimension(1:num_procs) :: lminval,lminval_x,lminval_y,lminval_z
      integer(kind = 4), dimension(1:num_procs) :: lminval_i,lminval_j,lminval_k
      real(kind=rKIND) :: qminval , qminval_x,qminval_y,qminval_z
      integer(kind = 4) :: qminval_i,qminval_j,qminval_k,dimxsnap,dimysnap,dimzsnap,n_samples,i1,j1,k1
      integer(kind = 4) :: ini_ibox,fin_ibox,ini_jbox,fin_jbox,ini_kbox,fin_kbox
      integer(kind = 4) :: ini_iboxsin,fin_iboxsin,ini_jboxsin,fin_jboxsin,ini_kboxsin,fin_kboxsin

      character(len=BUFSIZE) :: lockFile

      real(kind=4), dimension(:,:,:,:), allocatable  :: snap

#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',num_procs,') '
      write(whoamishort,'(i5)') layoutnumber+1

#ifdef CompileWithMPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
      call get_secnds(time_out2)
      time_end = time_out2%segundos
      l_aux = (time_end - time_begin  >  REPORTINGSECONDS) .or. forcetiming
#ifdef CompileWithMPI
      call MPI_AllReduce(l_aux, do_timing, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#else
      do_timing = l_aux
#endif
      n_info = n + 5
      call perform%reset()
      mustflushFIELDS = .FALSE.
      mustflushDATA = .FALSE.
      mustUnpack = .false.
      mustPostprocess = .false.
      mustflushXdmf=.false.
      mustflushVTK=.false.
      energy=0.0_RKIND
      if (do_timing) then
         perform%flushFIELDS = .FALSE.
         perform%flushDATA = .FALSE.
         if (abs(time_end - time_begin_absolute) < 1.0_RKIND) time_end = 60.0_RKIND+time_begin_absolute
         if (abs(time_end - time_begin) < 1.0_RKIND) time_end = 60.0_RKIND+time_begin
         speedInst = ((N - reportedinstant + 1) * mega_cells / (time_end - time_begin))
         speedAvg = ((N - INITIALtimeSTEP + 1) * mega_cells / (time_end - time_begin_absolute))
         if (speedAvg == 0) speedAvg=100.0_RKIND
         if (speedInst == 0) speedInst=speedAvg
#ifdef CompileWithMPI
         call MPI_AllReduce(speedInst, speedGlobInst, 1_4, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce(speedAvg, speedGlobAvg, 1_4, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#else
         speedGlobInst = speedInst
         speedGlobAvg = speedAvg
#endif
         in_aux = n + max(int((REPORTINGSECONDS / (total_mega_cells / speedGlobInst))) + 1,1)
#ifdef CompileWithMPI
         call MPI_AllReduce(in_aux, n_info, 1_4, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
#else
         n_info = in_aux
#endif
         inquire(FILE = 'stop', EXIST = stopNoFlushing)
         inquire(FILE = 'stop_only', EXIST = stop_only)
         inquire(FILE = 'stop_dontwritevtk', EXIST = stopdontwritevtk)
         inquire(FILE = 'stop_only_dontwritevtk', EXIST = stoponlydontwritevtk)
         if (stopNoFlushing) then
             dontwritevtk=.false.
         end if
         if (stopdontwritevtk) then
             stopNoFlushing=.true.
             dontwritevtk=.true.
         end if
         if (stoponlydontwritevtk) then
             stop_only=.true.
             dontwritevtk=.true.
         end if
         if (stop_only) then
             open(newunit=thefilenoflu,FILE = 'stop_only',action="read")
             read(thefilenoflu,*) target_name
             if (trim(adjustl(target_name))==trim(adjustl(nInputRoot))) then
                 stopNoFlushing=.true.
             else
                 stopNoFlushing=.false.
             end if
             close(thefilenoflu)
         end if
         inquire(FILE = 'stopflushing', EXIST = stopFlushing)
         inquire(FILE = 'stopflushing_only', EXIST = stopflushing_only)
         inquire(FILE = 'stopflushing_dontwritevtk', EXIST = stopflushingdontwritevtk)
         inquire(FILE = 'stopflushing_only_dontwritevtk', EXIST = stopflushingonlydontwritevtk)
         if (stopFlushing) then
             dontwritevtk=.false.
         end if
         if (stopflushingdontwritevtk) then
             stopFlushing=.true.
             dontwritevtk=.true.
         end if
         if (stopflushingonlydontwritevtk) then
             stopflushing_only=.true.
             dontwritevtk=.true.
         end if
         if (stopflushing_only) then
             open(newunit=thefilenoflu,FILE = 'stopflushing_only',action="read")
             read(thefilenoflu,*) target_name
             if (trim(adjustl(target_name))==trim(adjustl(nInputRoot))) then
                 stopFlushing=.true.
             else
                 stopFlushing=.false.
             end if
             close(thefilenoflu)
         end if
         inquire(FILE = 'flush', EXIST = mustflushFIELDS)
         inquire(FILE = 'flush_only', EXIST = flush_only)
         inquire(FILE = 'flush_dontwritevtk', EXIST = flushdontwritevtk)
         inquire(FILE = 'flush_only_dontwritevtk', EXIST = flushonlydontwritevtk)
         if (mustflushFIELDS) then
             dontwritevtk=.false.
         end if
         if (flushdontwritevtk) then
             mustflushFIELDS=.true.
             dontwritevtk=.true.
         end if
         if (flushonlydontwritevtk) then
             flush_only=.true.
             dontwritevtk=.true.
         end if
         if (flush_only) then
             open(newunit=thefilenoflu,FILE = 'flush_only',action="read")
             read(thefilenoflu,*) target_name
             if (trim(adjustl(target_name))==trim(adjustl(nInputRoot))) then
                 mustflushFIELDS=.true.
             else
                 mustflushFIELDS=.false.
             end if
             close(thefilenoflu)
         end if
         inquire(FILE = 'flushdata', EXIST = mustflushdata)
         inquire(FILE = 'flushdata_only', EXIST = flushdata_only)
         inquire(FILE = 'flushdata_dontwritevtk', EXIST = flushdatadontwritevtk)
         inquire(FILE = 'flushdata_only_dontwritevtk', EXIST = flushdataonlydontwritevtk)
         if (mustflushdata) then
             dontwritevtk=.false.
         end if
         if (flushdatadontwritevtk) then
             mustflushdata=.true.
             dontwritevtk=.true.
         end if
         if (flushdataonlydontwritevtk) then
             flushdata_only=.true.
             dontwritevtk=.true.
         end if
         if (flushdata_only) then
             open(newunit=thefilenoflu,FILE = 'flushdata_only',action="read")
             read(thefilenoflu,*) target_name
             if (trim(adjustl(target_name))==trim(adjustl(nInputRoot))) then
                 mustflushdata=.true.
             else
                 mustflushdata=.false.
             end if
             close(thefilenoflu)
         end if
         inquire(FILE = 'unpack', EXIST = mustUnpack)
         inquire(FILE = 'postprocess', EXIST = mustPostprocess)
         inquire(FILE = 'flushxdmf', EXIST = mustflushXdmf)
         inquire(FILE = 'flushvtk', EXIST = mustflushVTK)
         stopFlushing = stopFlushing .or.  &
         (ceiling((time_end-time_launch)/60.0_RKIND) >= maxCPUtime)
#ifdef CompileWithMPI
         l_aux=dontwritevtk
         call MPI_AllReduce(l_aux, dontwritevtk, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         l_aux=stopNoFlushing
         call MPI_AllReduce(l_aux, stopNoFlushing, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         l_aux=stopFlushing
         call MPI_AllReduce(l_aux, stopFlushing, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         l_aux=mustUnpack
         call MPI_AllReduce(l_aux, mustUnpack, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         l_aux=mustPostprocess
         call MPI_AllReduce(l_aux, mustPostprocess, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         l_aux=mustflushXdmf
         call MPI_AllReduce(l_aux, mustflushXdmf, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         l_aux=mustflushVTK
         call MPI_AllReduce(l_aux, mustflushVTK, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         l_aux=mustflushFIELDS
         call MPI_AllReduce(l_aux, mustflushFIELDS, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         l_aux=mustflushdata
         call MPI_AllReduce(l_aux, mustflushdata, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         mustflushfields = stopFlushing .or. mustflushfields
         stop_flag = stopNoFlushing .or. stopFlushing
         mustflushdata = mustflushdata .or. mustflushfields.or.stop_flag
         lmaxval  (1:num_procs) =  0.0_RKIND
         lmaxval_i(1:num_procs) =  0
         lmaxval_j(1:num_procs) =  0
         lmaxval_k(1:num_procs) =  0
         lmaxval_x(1:num_procs) =  0.0_RKIND
         lmaxval_y(1:num_procs) =  0.0_RKIND
         lmaxval_z(1:num_procs) =  0.0_RKIND
         lminval  (1:num_procs) =  0.0_RKIND
         lminval_i(1:num_procs) =  0
         lminval_j(1:num_procs) =  0
         lminval_k(1:num_procs) =  0
         lminval_x(1:num_procs) =  1e+20
         lminval_y(1:num_procs) =  1e+20
         lminval_z(1:num_procs) =  1e+20

         valueField = 0.0_RKIND
         ini_i = b%sweepSINPMLEx%XI - b%Ex%XI
         fin_i = b%sweepSINPMLEx%XE - b%Ex%XI
         ini_j = b%sweepSINPMLEx%YI - b%Ex%YI
         fin_j = b%sweepSINPMLEx%YE - b%Ex%YI
         ini_k = b%sweepSINPMLEx%ZI - b%Ex%ZI
         fin_k = b%sweepSINPMLEx%ZE - b%Ex%ZI

         do k = ini_k, fin_k
            do j = ini_j, fin_j
               do i = ini_i, fin_i
                  valueField = valueField + Ex(i, j, k) * Ex(i, j, k)
               end do
            end do
         end do
         ini_i = b%sweepSINPMLEy%XI - b%Ey%XI
         fin_i = b%sweepSINPMLEy%XE - b%Ey%XI
         ini_j = b%sweepSINPMLEy%YI - b%Ey%YI
         fin_j = b%sweepSINPMLEy%YE - b%Ey%YI
         ini_k = b%sweepSINPMLEy%ZI - b%Ey%ZI
         fin_k = b%sweepSINPMLEy%ZE - b%Ey%ZI
         do k = ini_k, fin_k
            do j = ini_j, fin_j
               do i = ini_i, fin_i
                  valueField = valueField + Ey(i, j, k) * Ey(i, j, k)
               end do
            end do
         end do
         ini_i = b%sweepSINPMLEz%XI - b%Ez%XI
         fin_i = b%sweepSINPMLEz%XE - b%Ez%XI
         ini_j = b%sweepSINPMLEz%YI - b%Ez%YI
         fin_j = b%sweepSINPMLEz%YE - b%Ez%YI
         ini_k = b%sweepSINPMLEz%ZI - b%Ez%ZI
         fin_k = b%sweepSINPMLEz%ZE - b%Ez%ZI
         do k = ini_k, fin_k
            do j = ini_j, fin_j
               do i = ini_i, fin_i
                  valueField = valueField + Ez(i, j, k) * Ez(i, j, k)
               end do
            end do
         end do
         energy = valueField
         energytotal = energy
#ifdef CompileWithMPI
         call MPI_AllReduce(energy, energyTotal, 1_4, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
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
                  valueField = sqrt(Ex(i, j, k) * Ex(i, j, k) + Ey(i, j, k) * Ey(i, j, k)+Ez(i, j, k) * Ez(i, j, k))
                  if (lmaxval  (layoutnumber+1)< valueField) then
                     lmaxval  (layoutnumber+1)= valueField
                     lmaxval_i(layoutnumber+1)=i+b%Hx%XI
                     lmaxval_j(layoutnumber+1)=j+b%Hy%YI
                     lmaxval_k(layoutnumber+1)=k+b%Hz%ZI
                     lmaxval_x(layoutnumber+1)=gridPoint%PhysCoor(IHX)%x(lmaxval_i(layoutnumber+1))
                     lmaxval_y(layoutnumber+1)=gridPoint%PhysCoor(IHY)%y(lmaxval_j(layoutnumber+1))
                     lmaxval_z(layoutnumber+1)=gridPoint%PhysCoor(IHZ)%z(lmaxval_k(layoutnumber+1))
                  end if
                  if (lminval  (layoutnumber+1)> valueField) then
                     lminval  (layoutnumber+1)= valueField
                     lminval_i(layoutnumber+1)=i+b%Hx%XI
                     lminval_j(layoutnumber+1)=j+b%Hy%YI
                     lminval_k(layoutnumber+1)=k+b%Hz%ZI
                     lminval_x(layoutnumber+1)=gridPoint%PhysCoor(IHX)%x(lminval_i(layoutnumber+1))
                     lminval_y(layoutnumber+1)=gridPoint%PhysCoor(IHY)%y(lminval_j(layoutnumber+1))
                     lminval_z(layoutnumber+1)=gridPoint%PhysCoor(IHZ)%z(lminval_k(layoutnumber+1))
                  end if
               end do
            end do
         end do

         NEWlmaxval  (1:num_procs) =0.0_RKIND
         NEWlmaxval_i(1:num_procs) =0
         NEWlmaxval_j(1:num_procs) =0
         NEWlmaxval_k(1:num_procs) =0
#ifdef CompileWithMPI
         call MPI_AllReduce(LMAXVAL, NEWlmaxval  , num_procs, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LMAXVAL_i, NEWlmaxval_I, num_procs, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LMAXVAL_j, NEWlmaxval_J, num_procs, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LMAXVAL_k, NEWlmaxval_K, num_procs, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LMAXVAL_x, NEWlmaxval_x, num_procs, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LMAXVAL_y, NEWlmaxval_y, num_procs, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LMAXVAL_z, NEWlmaxval_z, num_procs, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
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
         do i=1,num_procs
            if (abs(NEWlmaxval(i)) > qmaxval) then
               qmaxval   = abs(NEWlmaxval(i))
               qmaxval_i = newlmaxval_i(i)
               qmaxval_j = newlmaxval_j(i)
               qmaxval_k = newlmaxval_k(i)
               qmaxval_x = newlmaxval_x(i)
               qmaxval_y = newlmaxval_y(i)
               qmaxval_z = newlmaxval_z(i)
            end if
         end do

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         NEWlminval  (1:num_procs) =0.0_RKIND
         NEWlminval_i(1:num_procs) =0
         NEWlminval_j(1:num_procs) =0
         NEWlminval_k(1:num_procs) =0
#ifdef CompileWithMPI
         call MPI_AllReduce(LminVAL, NEWlminval  , num_procs, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LminVAL_i, NEWlminval_I, num_procs, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LminVAL_j, NEWlminval_J, num_procs, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LminVAL_k, NEWlminval_K, num_procs, MPI_INTEGER, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LminVAL_x, NEWlminval_x, num_procs, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LminVAL_y, NEWlminval_y, num_procs, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
         call MPI_AllReduce(LminVAL_z, NEWlminval_z, num_procs, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr)
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
         do i=1,num_procs
            if (abs(NEWlminval(i)) > qminval) then
               qminval   = abs(NEWlminval(i))
               qminval_i = newlminval_i(i)
               qminval_j = newlminval_j(i)
               qminval_k = newlminval_k(i)
               qminval_x = newlminval_x(i)
               qminval_y = newlminval_y(i)
               qminval_z = newlminval_z(i)
            end if
         end do

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

         if (layoutnumber == 0) then
            inquire(FILE = 'snap', EXIST = mustSnap)
            if (mustsnap) then
               open (35,FILE='snap')
               read(35,*,end=1153,err=1153) snapLevel
               read(35,*,end=1153,err=1153) snapStep
               read(35,*,end=1153,err=1153) snapHowMany
1153           close (35)
               call erasesignalingfiles(simu_devia)
            end if
         end if
#ifdef CompileWithMPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         call MPI_BCast(mustSnap, 1_4, MPI_LOGICAL, 0_4, MPI_COMM_WORLD, ierr)
         call MPI_BCast(snapLevel, 1_4, REALSIZE, 0_4, MPI_COMM_WORLD, ierr)
         call MPI_BCast(snapStep, 1_4, MPI_INTEGER, 0_4, MPI_COMM_WORLD, ierr)
         call MPI_BCast(snapHowMany, 1_4, MPI_INTEGER, 0_4, MPI_COMM_WORLD, ierr)
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

         if ((mustsnap.and.(lmaxval  (layoutnumber+1)> snapLevel)).or.(countersnap > 0)) then
            countersnap=countersnap + 1
            dimxsnap=int((fin_ibox-ini_ibox)/snapstep) + 1
            dimysnap=int((fin_jbox-ini_jbox)/snapstep) + 1
            dimzsnap=int((fin_kbox-ini_kbox)/snapstep) + 1
            if (.not.allocated(snap)) allocate (snap(ini_ibox:ini_ibox+dimxsnap, &
                                                     ini_jbox:ini_jbox+dimysnap,ini_kbox:ini_kbox+dimzsnap,1))
            snap=0.0_RKIND

            do k = ini_kbox, fin_kbox , snapStep
               do  j = ini_jbox, fin_jbox , snapStep
                  do i = ini_ibox, fin_ibox , snapStep
                     n_samples=0
                     valueField=0.0_RKIND
                     do k1=0,snapstep-1
                        do j1=0,snapstep-1
                           do i1=0,snapstep-1
                              if ((i+i1 <= fin_ibox).and.(j+j1 <= fin_jbox).and.(k+k1 <= fin_kbox)) then
                                 valueField = valueField+sqrt(Ex(i+i1, j+j1, k+k1) * Ex(i+i1, j+j1, k+k1) + &
                                                    Ey(i+i1, j+j1, k+k1) * Ey(i+i1, j+j1, k+k1)+ &
                                 Ez(i+i1, j+j1, k+k1) * Ez(i+i1, j+j1, k+k1))
                                 n_samples=n_samples+1
                              end if
                           end do
                        end do
                     end do
                     snap(ini_ibox+int((i-ini_ibox)/snapstep),ini_jbox+int((j-ini_jbox)/snapstep), &
                          ini_kbox+int((k-ini_kbox)/snapstep),1) = valueField/n_samples
                  end do
               end do
            end do

            write(chinstant,'(i8)') n
            write(minmax,'(a,e15.4e3,a,e15.4e3,a)')  '_',lminval(layoutnumber+1),'_',lmaxval(layoutnumber+1),'_'
            snap_file=trim(adjustl(nInputRoot))//'_snap_'//trim(adjustl(chinstant))//'_'// &
            trim(adjustl(whoamishort))

            lockFile=trim(adjustl(snap_file))//'.h5'
            call openclosedelete(lockFile)

            call write_xdmfsnap(n,snap_file,ini_ibox + b%Ex%XI , ini_ibox + dimxsnap + b%Ex%XI , &
            ini_jbox + b%Ex%YI , ini_jbox + dimysnap + b%Ex%YI , &
            ini_kbox + b%Ex%ZI , ini_kbox + dimzsnap + b%Ex%ZI , snap)

            write(dubuf,*)     trim(adjustl(whoami))//' Written Snap file at n= ',n,' max field over ',maxval(snap),'>',snapLevel
            deallocate(snap)
            call print11(layoutnumber,dubuf,.true.)
            if (countersnap >= snapHowMany) then
               mustsnap=.false.
               countersnap=0
            end if
         end if

#ifdef CompileWithMPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

         if (layoutnumber == 0) then
            write(dubuf,*) SEPARADOR,trim(adjustl(nInputRoot)),SEPARADOR
            call print11(layoutnumber,dubuf)
            write(dubuf,*) 'Switches: '//trim(adjustl(all_options))
            call print11(layoutnumber,dubuf)
                write(dubuf,*) 'MPI Processes: ',num_procs
                call print11(layoutnumber,dubuf)
            write(dubuf,*) 'Date/Time ', time_out2%fecha(7: 8),'/', &
            time_out2%fecha(5: 6),'/', &
            &               time_out2%fecha(1:4),'   ',time_out2%hora(1: 2), &
            ':',time_out2%hora(3: 4),':',time_out2%hora(5: 6)
            call print11(layoutnumber,dubuf)
            write(dubuf,*) 'Simulated:',n         ,'/',finaltimestep,' steps'
            call print11(layoutnumber,dubuf)
            if (permitscaling) then
                write(dubuf,'(a,e19.9e3,a,e19.9e3,a,e19.9e3)') 'Time= ',sgg%time(n),', dt0 (original)= ',dt0, &
                                                               ', dt(pscaled)= ',sgg%dt
            else
                write(dubuf,'(a,e19.9e3,a,e19.9e3,a,e19.9e3)') 'Time= ',sgg%time(n),', dt0 = ',sgg%dt
            end if
            call print11(layoutnumber,dubuf)
            if (energytotal > oldenergytotal) then
               write(dubuf,*)     'Total Energy (inc) :',energytotal
               if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
               oldenergytotal=energytotal
            else
               write(dubuf,*)     'Total Energy (dec) :',energytotal
               if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
               oldenergytotal=energytotal
            end if
            call print11(layoutnumber,dubuf)
            if (qmaxval_x<-1e19) then
               write (where_x,'(a)') ' PML '
            else
               write (where_x,'(e15.4e3)') qmaxval_x
            end if
            if (qmaxval_y<-1e19) then
               write (where_y,'(a)') ' PML '
            else
               write (where_y,'(e15.4e3)') qmaxval_y
            end if
            if (qmaxval_z<-1e19) then
               write (where_z,'(a)') ' PML '
            else
               write (where_z,'(e15.4e3)') qmaxval_z
            end if

            write(dubuf,'(a,e13.4e3,a,3i5,a)')     'Max field: ',qmaxval,' at (',qmaxval_i,qmaxval_j,qmaxval_k,')=('// &
            ' '//trim(adjustl(where_x))//','//' '//trim(adjustl(where_y))//','//' '//trim(adjustl(where_z))//')'

            if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
            call print11(layoutnumber,dubuf)
            do i=1,num_procs
               if (newlmaxval_x(i)<-1e19) then
                  write (where_x,'(a)') ' PML '
               else
                  write (where_x,'(e15.4e3)') newlmaxval_x(i)
               end if
               if (newlmaxval_y(i)<-1e19) then
                  write (where_y,'(a)') ' PML '
               else
                  write (where_y,'(e15.4e3)') newlmaxval_y(i)
               end if
               if (newlmaxval_z(i)<-1e19) then
                  write (where_z,'(a)') ' PML '
               else
                  write (where_z,'(e15.4e3)') newlmaxval_z(i)
               end if

               write(dubuf,'(a,i5,e15.4e3,a,e15.4e3,a,3i5,a)')  'Max field slice: ',i,NEWlmaxval(i),'/',maxSourceValue, &
               ' at (',newlmaxval_i(i),newlmaxval_j(i),newlmaxval_k(i),')=('// &
            ' '//trim(adjustl(where_x))//','//' '//trim(adjustl(where_y))//','//' '//trim(adjustl(where_z))//')'
            end do

            write(dubuf,*)     'Mins. since start  :',ceiling((time_end-time_launch)/60.0_RKIND)
            call print11(layoutnumber,dubuf)
            write(dubuf,*)     'Mins. until end    :', &
            min(ceiling(((finaltimestep-n)*total_mega_cells)/speedGlobAvg/60.0_RKIND), &
            maxCPUtime-ceiling((time_end-time_launch)/60.0_RKIND))
            call print11(layoutnumber,dubuf)
            if (everflushed) then
               write(dubuf,*) 'Mins. past flushing:',ceiling((time_end-time_begin2)/60.0_RKIND)
               call print11(layoutnumber,dubuf)
            else
               write(dubuf,*) 'Never flushed resuming fields.'
               call print11(layoutnumber,dubuf)
            end if
            if (flushsecondsFields /= 0) then
               write(dubuf,*) 'Mins. next flushing:',min(ceiling((flushsecondsFIELDS-(time_end-time_begin2))/60.0_RKIND), &
               ceiling(((finaltimestep-n)*total_mega_cells)/speedGlobAvg/60.0_RKIND))
               call print11(layoutnumber,dubuf)
            else
               if (maxCPUtime == TOPCPUTIME) then
                  write(dubuf,*) 'Will Never flush resuming fields.'
                  call print11(layoutnumber,dubuf)
               else
                  write(dubuf,*) 'Flushing of restarting DATA at the end.'
                  call print11(layoutnumber,dubuf)
               end if
            end if
            write(dubuf,*) 'Next info at step: ',n_info
            call print11(layoutnumber,dubuf)
            write(dubuf,*) 'Total Mcells:',total_mega_cells
            call print11(layoutnumber,dubuf)
            if (reportedinstant < n) then
               write(dubuf,*) 'Mcells/sec  :',speedGlobInst,'(',reportedinstant, ' to ',N,')'
               if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
               call print11(layoutnumber,dubuf)
            end if
            write(dubuf,*) 'Mcells/sec  :',speedGlobAvg ,'(',INITIALtimeSTEP, ' to ',N,')'
            if (simu_devia) dubuf=trim(adjustl(dubuf))//' (Stoch)'
            call print11(layoutnumber,dubuf)
            write(dubuf,*) SEPARADOR//SEPARADOR//SEPARADOR
            call print11(layoutnumber,dubuf)
            write (10,*)      sgg%time(n),energytotal
            call flush(11)
            call flush(10)
         end if
#ifdef CompileWithMPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         call get_secnds(time_out2)
         time_begin = time_out2%segundos
         reportedinstant=n+1
      end if

      stoponNaN_aux=.false.
#ifdef CompileWithMPI
      call MPI_AllReduce(stoponNaN_aux, stoponNaN, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#else
      stoponNaN = stoponNaN_aux
#endif
      if (stoponNaN) then
         if (layoutnumber == 0) then
            call log_msg(LOG_ERROR, 'PROGRESS', 'Instabilities detected, aborting. Single wires: reduce WIRE radii '// &
                         'and/or sgg%dt. Surface IBCs: reduce -att factor')
         end if
#ifdef CompileWithMPI
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         stop_flag=.true.
      end if
      l_aux = ( ((time_end-time_begin2) > flushsecondsFIELDS).AND. &
      (flushsecondsFIELDS/=0                       )).or. &
      mustflushFIELDS
#ifdef CompileWithMPI
      call MPI_AllReduce(l_aux, do_flushFIELDS, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#else
      do_flushFIELDS = l_aux
#endif
      l_aux = ( ((time_end-time_begin3) > flushsecondsDATA).AND. &
      (flushsecondsDATA/=0                       )).or. &
      mustflushDATA
#ifdef CompileWithMPI
      call MPI_AllReduce(l_aux, do_flushDATA, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#else
      do_flushDATA = l_aux
#endif
      if (do_flushFIELDS) then
         mustflushFIELDS=.false.
#ifdef CompileWithMPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         call get_secnds(time_out2)
         time_begin2=time_out2%segundos
         perform%flushFIELDS=.true.
         if (layoutnumber == 0) then
             call erasesignalingfiles(simu_devia)
         end if
      end if
      if (do_flushDATA) then
         mustflushDATA=.false.
#ifdef CompileWithMPI
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif
         call get_secnds(time_out2)
         time_begin3=time_out2%segundos
         perform%flushDATA=.true.
         if (layoutnumber == 0) then
             call erasesignalingfiles(simu_devia)
         end if
      end if
      if (mustunpack) then
         mustunpack=.false.
         perform%unpackFlag=.true.
         if (layoutnumber == 0) then
             call erasesignalingfiles(simu_devia)
         end if
      end if
      if (mustpostprocess) then
         mustpostprocess=.false.
         perform%postprocess=.true.
         if (layoutnumber == 0) then
             call erasesignalingfiles(simu_devia)
         end if
      end if
      if (mustflushXdmf) then
         mustflushXdmf=.false.
         perform%flushXdmf=.true.
         if (layoutnumber == 0) then
             call erasesignalingfiles(simu_devia)
         end if
      end if
      if (mustflushVTK) then
         mustflushVTK=.false.
         perform%flushVTK=.true.
         if (layoutnumber == 0) then
             call erasesignalingfiles(simu_devia)
         end if
      end if
      return
   end subroutine Timing

   function createPoints(sgg)  result(gridPoint)
      type(SGGFDTDINFO_t), intent(inout) :: sgg
      type(coorsxyzP_t) :: gridPoint
      integer(kind=4) :: i,j,k,field

      do field=IEX,IHZ
         allocate (gridPoint%PhysCoor(field)%x(sgg%Sweep(field)%XI-1 : sgg%Sweep(field)%XE+1), &
         gridPoint%PhysCoor(field)%y(sgg%Sweep(field)%YI-1 : sgg%Sweep(field)%YE+1), &
         gridPoint%PhysCoor(field)%z(sgg%Sweep(field)%ZI-1 : sgg%Sweep(field)%ZE+1))
         gridPoint%PhysCoor(field)%x(sgg%Sweep(field)%XI-1 : sgg%Sweep(field)%XE+1) = -1e20
         gridPoint%PhysCoor(field)%y(sgg%Sweep(field)%YI-1 : sgg%Sweep(field)%YE+1) = -1e20
         gridPoint%PhysCoor(field)%z(sgg%Sweep(field)%ZI-1 : sgg%Sweep(field)%ZE+1) = -1e20
      end do

      field=IEX
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         gridPoint%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         gridPoint%PhysCoor(field)%y(j)=sgg%LineY(j)
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         gridPoint%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=IEY
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         gridPoint%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         gridPoint%PhysCoor(field)%y(j)=(sgg%LineY(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         gridPoint%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=IEZ
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         gridPoint%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         gridPoint%PhysCoor(field)%y(j)=sgg%LineY(j)
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         gridPoint%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=IHX
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         gridPoint%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         gridPoint%PhysCoor(field)%y(j)=(sgg%LineY(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         gridPoint%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=IHY
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         gridPoint%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         gridPoint%PhysCoor(field)%y(j)=sgg%LineY(j)
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         gridPoint%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=IHZ
      do i=sgg%SINPMLSweep(field)%XI-1,sgg%SINPMLSweep(field)%XE+1
         gridPoint%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%SINPMLSweep(field)%YI-1,sgg%SINPMLSweep(field)%YE+1
         gridPoint%PhysCoor(field)%y(j)=(sgg%LineY(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%SINPMLSweep(field)%ZI-1,sgg%SINPMLSweep(field)%ZE+1
         gridPoint%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      sgg%gridPoint = gridPoint
      return
   end function

end module progress_m
