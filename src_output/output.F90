module output
   use FDETYPES
   use Report
   use mod_domain
   use mod_outputUtils
   use mod_pointProbeOutput
   use mod_wireProbeOutput
   use mod_bulkProbeOutput
   use mod_movieProbeOutput
   use mod_frequencySliceProbeOutput
   use mod_farFieldOutput
   use mtln_solver_mod, only: mtln_solver_t => mtln_t
   use Wire_bundles_mtln_mod, only: GetSolverPtr
   use mod_mapVTKOutput
   

   implicit none
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: solver_output_t
   public :: GetOutputs
   public :: init_outputs
   public :: update_outputs
   public :: flush_outputs
   public :: close_outputs

   public :: POINT_PROBE_ID, WIRE_CURRENT_PROBE_ID, WIRE_CHARGE_PROBE_ID, BULK_PROBE_ID, VOLUMIC_CURRENT_PROBE_ID, &
             MOVIE_PROBE_ID, FREQUENCY_SLICE_PROBE_ID, FAR_FIELD_PROBE_ID
   !===========================

   !===========================
   !  Private interface summary
   !===========================
   private :: get_required_output_count
   !===========================

   integer(kind=SINGLE), parameter :: UNDEFINED_PROBE = -1, &
                                      POINT_PROBE_ID = 0, &
                                      WIRE_CURRENT_PROBE_ID = 1, &
                                      WIRE_CHARGE_PROBE_ID = 2, &
                                      BULK_PROBE_ID = 3, &
                                      VOLUMIC_CURRENT_PROBE_ID = 4, &
                                      MOVIE_PROBE_ID = 5, &
                                      FREQUENCY_SLICE_PROBE_ID = 6, &
                                      FAR_FIELD_PROBE_ID = 7, &
                                      MAPVTK_ID = 8

   REAL(KIND=RKIND), save           ::  eps0, mu0
   REAL(KIND=RKIND), pointer, dimension(:), save  ::  InvEps, InvMu
   type(solver_output_t), pointer, dimension(:), save  ::  outputs
   type(problem_info_t), save, target :: problemInfo

   interface init_solver_output
      module procedure &
         init_point_probe_output, &
         init_wire_current_probe_output, &
         init_wire_charge_probe_output, &
         init_bulk_probe_output, &
         init_movie_probe_output, &
         init_frequency_slice_probe_output, &
         init_farField_probe_output, &
         init_mapvtk_output
   end interface

   interface update_solver_output
      module procedure &
         update_point_probe_output, &
         update_wire_current_probe_output, &
         update_wire_charge_probe_output, &
         update_bulk_probe_output, &
         update_movie_probe_output, &
         update_frequency_slice_probe_output, &
         update_farField_probe_output
   end interface

   interface flush_solver_output
      module procedure &
         flush_point_probe_output, &
         flush_wire_current_probe_output, &
         flush_wire_charge_probe_output, &
         flush_bulk_probe_output, &
         flush_movie_probe_output, &
         flush_frequency_slice_probe_output, &
         flush_farField_probe_output

   end interface

   interface close_solver_output
      module procedure close_frequency_slice_probe_output
   end interface
contains

   function GetOutputs() result(r)
      type(solver_output_t), pointer, dimension(:)  ::  r
      r => outputs
      return
   end function

   function GetProblemInfo() result(r)
      type(problem_info_t), pointer ::  r
      r => problemInfo
      return
   end function

   subroutine init_outputs(sgg, media, sinpml_fullsize, materialTags, bounds, control, observationsExists, wiresExists)
      type(SGGFDTDINFO), intent(in) ::  sgg
      type(media_matrices_t), target, intent(in) :: media
      type(limit_t), dimension(:), target, intent(in)  ::  SINPML_fullsize
      type(bounds_t),intent(in), target :: bounds
      type(taglist_t),intent(in), target :: materialTags
      type(sim_control_t), intent(in) :: control
      logical, intent(inout) :: wiresExists
      logical, intent(out) :: observationsExists

      type(domain_t) :: domain
      type(spheric_domain_t) :: sphericRange
      type(cell_coordinate_t) :: lowerBound, upperBound
      integer(kind=SINGLE) :: i, ii, outputRequestType
      integer(kind=SINGLE) :: NODE
      integer(kind=SINGLE) :: outputCount
      integer(kind=SINGLE) :: requestedOutputs
      character(len=BUFSIZE) :: outputTypeExtension

#ifdef CompileWithMTLN
      logical :: thereAreMtlnObservations = .false.
#endif
      observationsExists = .false.
      requestedOutputs = get_required_output_count(sgg)

      problemInfo%geometryToMaterialData => media
      problemInfo%materialList => sgg%Med
      problemInfo%simulationBounds => bounds
      problemInfo%problemDimension => SINPML_fullsize
      problemInfo%materialTag => materialTags

      outputs => NULL()
      allocate (outputs(requestedOutputs))

      allocate (InvEps(0:sgg%NumMedia - 1), InvMu(0:sgg%NumMedia - 1))
      outputCount = 0

      InvEps(0:sgg%NumMedia - 1) = 1.0_RKIND/(Eps0*sgg%Med(0:sgg%NumMedia - 1)%Epr)
      InvMu(0:sgg%NumMedia - 1) = 1.0_RKIND/(Mu0*sgg%Med(0:sgg%NumMedia - 1)%Mur)

      !do ii = 1, sgg%NumberRequest
      !do i = 1, sgg%Observation(ii)%nP
      !   call eliminate_unnecesary_observation_points(sgg%Observation(ii)%P(i), output(ii)%item(i), &
      !     sgg%Sweep, sgg%SINPMLSweep, sgg%Observation(ii)%P(1)%ZI, sgg%Observation(ii)%P(1)%ZE, control%layoutnumber, control%size)
      !end do
      !end do

#ifdef CompileWithMTLN
      block
         type(mtln_solver_t), pointer :: mtln_solver
         integer :: i, j
         mtln_solver => GetSolverPtr()
         if (mtln_solver%number_of_bundles > 0) then
            do i = 1, ubound(mtln_solver%bundles, 1)
               if (ubound(mtln_solver%bundles(i)%probes, 1) /= 0) then
                  do j = 1, ubound(mtln_solver%bundles(i)%probes, 1)
                     if (mtln_solver%bundles(i)%probes(j)%in_layer) thereAreMtlnObservations = .true.
                  end do
               end if
            end do
         end if
      end block
#endif

      do ii = 1, sgg%NumberRequest
         domain = preprocess_domain(sgg%Observation(ii), sgg%tiempo, sgg%dt, control%finaltimestep)
         if (domain%domainType == UNDEFINED_DOMAIN) cycle
         do i = 1, sgg%Observation(ii)%nP
            lowerBound%x = sgg%observation(ii)%P(i)%XI
            lowerBound%y = sgg%observation(ii)%P(i)%YI
            lowerBound%z = sgg%observation(ii)%P(i)%ZI

            upperBound%x = sgg%observation(ii)%P(i)%XE
            upperBound%y = sgg%observation(ii)%P(i)%YE
            upperBound%z = sgg%observation(ii)%P(i)%ZE
            NODE = sgg%observation(ii)%P(i)%NODE

            outputTypeExtension = trim(adjustl(control%nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))

            outputRequestType = sgg%observation(ii)%P(i)%what
            select case (outputRequestType)
            case (mapvtk)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = MAPVTK_ID

               allocate (outputs(outputCount)%mapvtkOutput)
               call init_solver_output(outputs(outputCount)%mapvtkOutput, lowerBound, upperBound, outputRequestType, outputTypeExtension, control%mpidir, problemInfo)
               call create_geometry_simulation_vtu(outputs(outputCount)%mapvtkOutput, control, sgg%LineX, sgg%LineY, sgg%LineZ)
               !! call adjust_computation_range --- Required due to issues in mpi region edges

            case (iEx, iEy, iEz, iHx, iHy, iHz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = POINT_PROBE_ID

               allocate (outputs(outputCount)%pointProbe)
               call init_solver_output(outputs(outputCount)%pointProbe, lowerBound, outputRequestType, domain, outputTypeExtension, control%mpidir, sgg%dt)
            case (iJx, iJy, iJz)
               if (wiresExists) then
                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = WIRE_CURRENT_PROBE_ID

                  allocate (outputs(outputCount)%wireCurrentProbe)
                  call init_solver_output(outputs(outputCount)%wireCurrentProbe, lowerBound, NODE, outputRequestType, domain, problemInfo%materialList, outputTypeExtension, control%mpidir, control%wiresflavor)
               end if

            case (iQx, iQy, iQz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = WIRE_CHARGE_PROBE_ID

               allocate (outputs(outputCount)%wireChargeProbe)
               call init_solver_output(outputs(outputCount)%wireChargeProbe, lowerBound, NODE, outputRequestType, domain, outputTypeExtension, control%mpidir, control%wiresflavor)

            case (iBloqueJx, iBloqueJy, iBloqueJz, iBloqueMx, iBloqueMy, iBloqueMz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = BULK_PROBE_ID

               allocate (outputs(outputCount)%bulkCurrentProbe)
               call init_solver_output(outputs(outputCount)%bulkCurrentProbe, lowerBound, upperBound, outputRequestType, domain, outputTypeExtension, control%mpidir)
               !! call adjust_computation_range --- Required due to issues in mpi region edges

            case (iCur, iMEC, iMHC, iCurX, iCurY, iCurZ, iExC, iEyC, iEzC, iHxC, iHyC, iHzC)
               call adjust_bound_range()

               if (domain%domainType == TIME_DOMAIN) then

                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = MOVIE_PROBE_ID
                  allocate (outputs(outputCount)%movieProbe)
                  call init_solver_output(outputs(outputCount)%movieProbe, lowerBound, upperBound, outputRequestType, domain, control, problemInfo, outputTypeExtension)
               else if (domain%domainType == FREQUENCY_DOMAIN) then

                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = FREQUENCY_SLICE_PROBE_ID
                  allocate (outputs(outputCount)%frequencySliceProbe)
                  call init_solver_output(outputs(outputCount)%frequencySliceProbe, lowerBound, upperBound, sgg%dt, outputRequestType, domain, outputTypeExtension, control, problemInfo)
               end if
            case (farfield)
               sphericRange = preprocess_polar_range(sgg%Observation(ii))

               outputCount = outputCount + 1
               outputs(outputCount)%outputID = FAR_FIELD_PROBE_ID
               allocate (outputs(outputCount)%farFieldOutput)
               call init_solver_output(outputs(outputCount)%farFieldOutput, sgg, lowerBound, upperBound, outputRequestType, domain, sphericRange, outputTypeExtension, sgg%Observation(ii)%FileNormalize, control, problemInfo, eps0, mu0)
            case default
               call stoponerror(0, 0, 'OutputRequestType type not implemented yet on new observations')
            end select
         end do
      end do
      if (outputCount /= requestedOutputs) then
         call remove_unused_outputs(outputs)
         outputCount = size(outputs)
      end if
      if (outputCount /= 0) observationsExists = .true.
#ifdef CompileWithMTLN
      observationsExists = observationsExists .or. thereAreMtlnObservations
#endif
      if (observationsExists) call registerOutputFiles(control, outputCount)
      return
   contains
      subroutine adjust_bound_range()
         select case (outputRequestType)
         case (iExC, iEyC, iHzC, iMhC)
            lowerBound%z = max(sgg%Sweep(fieldo(outputRequestType, 'Z'))%ZI, sgg%observation(ii)%P(i)%ZI)
            upperBound%z = min(sgg%Sweep(fieldo(outputRequestType, 'Z'))%ZE - 1, sgg%observation(ii)%P(i)%ZE)
         case (iEzC, iHxC, iHyC, iMeC)
            lowerBound%z = max(sgg%Sweep(fieldo(outputRequestType, 'Z'))%ZI, sgg%observation(ii)%P(i)%ZI)
            upperbound%z = min(sgg%Sweep(fieldo(outputRequestType, 'Z'))%ZE, sgg%observation(ii)%P(i)%ZE)
         case (iCur, iCurX, iCurY, iCurZ)
            lowerBound%z = max(sgg%Sweep(fieldo(outputRequestType, 'X'))%ZI, sgg%observation(ii)%P(i)%ZI) !ojo estaba sweep(iEz) para ser conservador...puede dar problemas!! 03/07/15
            upperbound%z = min(sgg%Sweep(fieldo(outputRequestType, 'X'))%ZE, sgg%observation(ii)%P(i)%ZE) !ojo estaba sweep(iEz) para ser conservador...puede dar problemas!! 03/07/15
         end select
      end subroutine
      function preprocess_domain(observation, timeArray, simulationTimeStep, finalStepIndex) result(newDomain)
         type(Obses_t), intent(in) :: observation
         real(kind=RKIND_tiempo), pointer, dimension(:), intent(in) :: timeArray
         real(kind=RKIND_tiempo), intent(in) :: simulationTimeStep
         integer(kind=4), intent(in) :: finalStepIndex
         type(domain_t) :: newDomain

         integer(kind=SINGLE) :: nFreq

         if (observation%TimeDomain) then
            newdomain = domain_t(real(observation%InitialTime, kind=RKIND_tiempo), &
                                 real(observation%FinalTime, kind=RKIND_tiempo), &
                                 real(observation%TimeStep, kind=RKIND_tiempo))

            newdomain%tstep = max(newdomain%tstep, simulationTimeStep)

            if (10.0_RKIND*(newdomain%tstop - newdomain%tstart)/min(simulationTimeStep, newdomain%tstep) >= huge(1_4)) then
               newdomain%tstop = newdomain%tstart + min(simulationTimeStep, newdomain%tstep)*huge(1_4)/10.0_RKIND
            end if

            if (newDomain%tstart < newDomain%tstep) then
               newDomain%tstart = 0.0_RKIND_tiempo
            end if

            if (newDomain%tstep > (newdomain%tstop - newdomain%tstart)) then
               newDomain%tstop = newDomain%tstart + newDomain%tstep
            end if

         elseif (observation%FreqDomain) then
            !Just linear progression for now. Need to bring logartihmic info to here
            nFreq = int((observation%FinalFreq - observation%InitialFreq)/observation%FreqStep, kind=SINGLE) + 1_SINGLE
            newdomain = domain_t(observation%InitialFreq, observation%FinalFreq, nFreq, logarithmicspacing=.false.)

            newDomain%fstep = min(newDomain%fstep, 2.0_RKIND/simulationTimeStep)
            if ((newDomain%fstep > newDomain%fstop - newDomain%fstart) .or. (newDomain%fstep == 0)) then
               newDomain%fstep = newDomain%fstop - newDomain%fstart
               newDomain%fstop = newDomain%fstart + newDomain%fstep
            end if

            newDomain%fnum = int((newDomain%fstop - newDomain%fstart)/newDomain%fstep, kind=SINGLE)

         else
            newDomain = domain_t()
         end if
         return
      end function preprocess_domain

      function preprocess_polar_range(observation) result(sphericDomain)
         type(spheric_domain_t) :: sphericDomain
         type(Obses_t), intent(in) :: observation

         sphericDomain%phiStart = observation%phiStart
         sphericDomain%phiStop = observation%phiStop
         sphericDomain%phiStep = observation%phiStep
         sphericDomain%thetaStart = observation%thetaStart
         sphericDomain%thetaStop = observation%thetaStop
         sphericDomain%thetaStep = observation%thetaStep
      end function preprocess_polar_range

   end subroutine init_outputs

   subroutine update_outputs(control, discreteTimeArray, timeIndx, fieldsReference)
      integer(kind=SINGLE), intent(in) :: timeIndx
      real(kind=RKIND_tiempo), dimension(:), intent(in) :: discreteTimeArray
      integer(kind=SINGLE) :: i, id
      type(sim_control_t), intent(in) :: control
      real(kind=RKIND), pointer, dimension(:, :, :) :: fieldComponent
      type(field_data_t) :: fieldReference
      type(fields_reference_t), intent(in) :: fieldsReference
      real(kind=RKIND_tiempo) :: discreteTime

      discreteTime = discreteTimeArray(timeIndx)

      do i = 1, size(outputs)
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID)
            fieldComponent => get_field_component(outputs(i)%pointProbe%component, fieldsReference) !Cada componente requiere de valores deiferentes pero estos valores no se como conseguirlos
            call update_solver_output(outputs(i)%pointProbe, discreteTime, fieldComponent)
         case (WIRE_CURRENT_PROBE_ID)
            call update_solver_output(outputs(i)%wireCurrentProbe, discreteTime, control, InvEps, InvMu)
         case (WIRE_CHARGE_PROBE_ID)
            call update_solver_output(outputs(i)%wireChargeProbe, discreteTime)
         case (BULK_PROBE_ID)
            fieldReference = get_field_reference(outputs(i)%bulkCurrentProbe%component, fieldsReference)
            call update_solver_output(outputs(i)%bulkCurrentProbe, discreteTime, fieldReference)
         case (MOVIE_PROBE_ID)
            call update_solver_output(outputs(i)%movieProbe, discreteTime, fieldsReference, control, problemInfo)
         case (FREQUENCY_SLICE_PROBE_ID)
            call update_solver_output(outputs(i)%frequencySliceProbe, discreteTime, fieldsReference, control, problemInfo)
         case (FAR_FIELD_PROBE_ID)
            call update_solver_output(outputs(i)%farFieldOutput, timeIndx, problemInfo%simulationBounds, fieldsReference)
         case(MAPVTK_ID)
         case default
            call stoponerror(0, 0, 'Output update not implemented')
         end select
      end do

   end subroutine update_outputs

   subroutine flush_outputs(simulationTimeArray, simulationTimeIndex, control, fields, bounds, farFieldFlushRequested)
      implicit none
      type(fields_reference_t), target :: fields
      type(fields_reference_t), pointer :: fieldsPtr
      type(sim_control_t), intent(in) :: control
      type(bounds_t), intent(in) :: bounds
      logical, intent(in) :: farFieldFlushRequested
      real(KIND=RKIND_tiempo), pointer, dimension(:), intent(in) :: simulationTimeArray
      integer, intent(in) :: simulationTimeIndex
      integer :: outIdx

      fieldsPtr => fields

      do outIdx = 1, size(outputs)
         select case (outputs(outIdx)%outputID)
         case (POINT_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%pointProbe)
         case (WIRE_CURRENT_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%wireCurrentProbe)
         case (WIRE_CHARGE_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%wireChargeProbe)
         case (BULK_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%bulkCurrentProbe)
         case (MOVIE_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%movieProbe)
         case (FREQUENCY_SLICE_PROBE_ID)
            call flush_solver_output(outputs(outIdx)%frequencySliceProbe)
         case (FAR_FIELD_PROBE_ID)
            if (farFieldFlushRequested) call flush_solver_output(outputs(outIdx)%farFieldOutput, simulationTimeArray, simulationTimeIndex, control, fieldsPtr, bounds)
         case default
         end select
      end do
   end subroutine flush_outputs

   subroutine remove_unused_outputs(output_list)
      implicit none
      type(solver_output_t), pointer, intent(inout) :: output_list(:)

      type(solver_output_t), allocatable :: tmp(:)
      integer :: i, n, k

      n = count(output_list%outputID /= UNDEFINED_PROBE)

      allocate (tmp(n))

      ! Copy valid elements
      k = 0
      do i = 1, size(output_list)
         if (output_list(i)%outputID /= UNDEFINED_PROBE) then
            k = k + 1
            tmp(k) = output_list(i)   ! deep copy of all allocatable components
         end if
      end do

      ! Replace the saved pointer target safely
      if (associated(output_list)) deallocate (output_list)
      allocate (output_list(n))
      output_list = tmp

   end subroutine remove_unused_outputs

   subroutine close_outputs()
      integer :: i
      do i = 1, size(outputs)
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID)
         case (WIRE_CURRENT_PROBE_ID)
         case (WIRE_CHARGE_PROBE_ID)
         case (BULK_PROBE_ID)
         case (VOLUMIC_CURRENT_PROBE_ID)
         case (MOVIE_PROBE_ID)
         case (FREQUENCY_SLICE_PROBE_ID)
            call close_solver_output(outputs(i)%frequencySliceProbe)
         end select
      end do
   end subroutine

   subroutine create_pvd(pvdPath)
      implicit none
      character(len=*), intent(out) :: pvdPath
      integer :: ios
      integer :: unit

      open (newunit=unit, file=trim(pvdPath), status="replace", action="write", iostat=ios)
      if (ios /= 0) stop "Error al crear archivo PVD"

      ! Escribimos encabezados XML
      write (unit, *) '<?xml version="1.0"?>'
      write (unit, *) '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
      write (unit, *) '  <Collection>'
      close (unit)
   end subroutine create_pvd

   subroutine close_pvd(pvdPath)
      implicit none
      character(len=*), intent(in) :: pvdPath
      integer :: unit
      integer :: ios
      open (newunit=unit, file=trim(pvdPath), status="old", action="write", iostat=ios)
      if (ios /= 0) stop "Error al abrir archivo PVD"
      write (unit, *) '  </Collection>'
      write (unit, *) '</VTKFile>'
      close (unit)
   end subroutine close_pvd

   function get_required_output_count(sgg) result(count)
      type(SGGFDTDINFO), intent(in) :: sgg
      integer(kind=SINGLE) ::i, count
      count = 0
      do i = 1, sgg%NumberRequest
         count = count + sgg%Observation(i)%nP
      end do
      return
   end function

   subroutine registerOutputFiles(control, outputCount)
      type(sim_control_t), intent(in) :: control
      integer, intent(in) :: outputCount

      character(LEN=BUFSIZE)  ::  whoami, whoamishort, outputRequestFile
      integer :: iostat, i, unit

      write (whoamishort, '(i5)') control%layoutnumber + 1
      write (whoami, '(a,i5,a,i5,a)') '(', control%layoutnumber + 1, '/', control%size, ') '
      write (outputRequestFile, *) trim(adjustl(control%nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'

      call create_file_with_path(outputRequestFile, iostat)
      if (iostat /= 0) call StopOnError(control%layoutnumber, control%size, 'Error while creating new outputrequestRegister file...')

      open (newunit=unit, file=trim(adjustl(outputRequestFile)), status='replace', action='write', position='append', iostat=iostat)
      do i = 1, outputCount
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID)
            if (any(outputs(i)%pointProbe%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
               write (unit, *) trim(adjustl(outputs(i)%pointProbe%filePathTime))
            end if
            if (any(outputs(i)%pointProbe%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
               write (unit, *) trim(adjustl(outputs(i)%pointProbe%filePathFreq))
            end if
         case (WIRE_CURRENT_PROBE_ID)
            write (unit, *) trim(adjustl(outputs(i)%wireCurrentProbe%filePathTime))
         case (WIRE_CHARGE_PROBE_ID)
            write (unit, *) trim(adjustl(outputs(i)%wireChargeProbe%filePathTime))
         case (BULK_PROBE_ID)
            write (unit, *) trim(adjustl(outputs(i)%bulkCurrentProbe%filePathTime))
         case (MOVIE_PROBE_ID)
            write (unit, *) trim(adjustl(outputs(i)%movieProbe%filePathTime))
         case (FREQUENCY_SLICE_PROBE_ID)
            write (unit, *) trim(adjustl(outputs(i)%frequencySliceProbe%filePathFreq))
         case (FAR_FIELD_PROBE_ID)
            write (unit, *) trim(adjustl(outputs(i)%farFieldOutput%filePathFreq))
         case (MAPVTK_ID)
            write (unit, *) trim(adjustl(outputs(i)%mapvtkOutput%path))
         case default
            call stoponerror(0, 0, 'Output update not implemented')
         end select
      end do

      write (unit, *) 'END!'
      close (unit)
   end subroutine

end module output
