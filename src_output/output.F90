module output
   use FDETYPES
   use Report
   use mod_domain
   use mod_outputUtils
   use mod_pointProbeOutput
   use mod_wireProbeOutput
   use mod_bulkProbeOutput
   use mod_volumicProbeOutput
   use mod_movieProbeOutput
   use mod_frequencySliceProbeOutput
   use mod_farFieldOutput

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

   integer(kind=SINGLE), parameter :: POINT_PROBE_ID = 0, &
                                      WIRE_CURRENT_PROBE_ID = 1, &
                                      WIRE_CHARGE_PROBE_ID = 2, &
                                      BULK_PROBE_ID = 3, &
                                      VOLUMIC_CURRENT_PROBE_ID = 4, &
                                      MOVIE_PROBE_ID = 5, &
                                      FREQUENCY_SLICE_PROBE_ID = 6, &
                                      FAR_FIELD_PROBE_ID = 7

   REAL(KIND=RKIND), save           ::  eps0, mu0
   REAL(KIND=RKIND), pointer, dimension(:), save  ::  InvEps, InvMu
   type(solver_output_t), pointer, dimension(:), save  ::  outputs
   type(problem_info_t), save :: problemInfo

   interface init_solver_output
      module procedure &
         init_point_probe_output, &
         init_wire_current_probe_output, &
         init_wire_charge_probe_output, &
         init_bulk_probe_output, &
         init_volumic_probe_output, &
         init_movie_probe_output, &
         init_frequency_slice_probe_output, &
         init_farField_probe_output
   end interface

   interface create_empty_files
      module procedure &
         create_point_probe_output_files, &
         create_wire_current_probe_output, &
         create_wire_charge_probe_output, &
         create_bulk_probe_output
   end interface

   interface update_solver_output
      module procedure &
         update_point_probe_output, &
         update_wire_current_probe_output, &
         update_wire_charge_probe_output, &
         update_bulk_probe_output, &
         update_volumic_probe_output, &
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

   subroutine init_outputs(sgg, media, sinpml_fullsize, bounds, control, observationsExists, wiresExists)
      type(SGGFDTDINFO), intent(in) ::  sgg
      type(media_matrices_t), intent(in) :: media
      type(limit_t), dimension(:), intent(in)  ::  SINPML_fullsize
      type(bounds_t) :: bounds
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

      observationsExists = .false.
      requestedOutputs = get_required_output_count(sgg)

      problemInfo%geometryToMaterialData => media
      problemInfo%materialList => sgg%Med
      problemInfo%simulationBounds => bounds
      problemInfo%problemDimension => SINPML_fullsize

      outputs => NULL()
      allocate (outputs(requestedOutputs))

      allocate (InvEps(0:sgg%NumMedia - 1), InvMu(0:sgg%NumMedia - 1))
      outputCount = 0

      InvEps(0:sgg%NumMedia - 1) = 1.0_RKIND/(Eps0*sgg%Med(0:sgg%NumMedia - 1)%Epr)
      InvMu(0:sgg%NumMedia - 1) = 1.0_RKIND/(Mu0*sgg%Med(0:sgg%NumMedia - 1)%Mur)

      do ii = 1, sgg%NumberRequest
      do i = 1, sgg%Observation(ii)%nP
         call eliminate_unnecesary_observation_points(sgg%Observation(ii)%P(i), output(ii)%item(i), &
           sgg%Sweep, sgg%SINPMLSweep, sgg%Observation(ii)%P(1)%ZI, sgg%Observation(ii)%P(1)%ZE, control%layoutnumber, control%size)
      end do
      end do

      do ii = 1, sgg%NumberRequest
         do i = 1, sgg%Observation(ii)%nP
            lowerBound%x = sgg%observation(ii)%P(i)%XI
            lowerBound%y = sgg%observation(ii)%P(i)%YI
            lowerBound%z = sgg%observation(ii)%P(i)%ZI

            upperBound%x = sgg%observation(ii)%P(i)%XE
            upperBound%y = sgg%observation(ii)%P(i)%YE
            upperBound%z = sgg%observation(ii)%P(i)%ZE
            NODE = sgg%observation(ii)%P(i)%NODE

            domain = preprocess_domain(sgg%Observation(ii), sgg%tiempo, sgg%dt, control%finaltimestep)
            outputTypeExtension = trim(adjustl(control%nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))

            outputRequestType = sgg%observation(ii)%P(i)%what
            select case (outputRequestType)
            case (iEx, iEy, iEz, iHx, iHy, iHz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = POINT_PROBE_ID

               allocate (outputs(outputCount)%pointProbe)
               call init_solver_output(outputs(outputCount)%pointProbe, lowerBound, outputRequestType, domain, outputTypeExtension, control, sgg%dt)
               call create_empty_files(outputs(outputCount)%pointProbe)
            case (iJx, iJy, iJz)
               if (wiresExists) then
                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = WIRE_CURRENT_PROBE_ID

                  allocate (outputs(outputCount)%wireCurrentProbe)
                  call init_solver_output(outputs(outputCount)%wireCurrentProbe, lowerBound, NODE, outputRequestType, domain, outputTypeExtension, control, problemInfo)
                  call create_empty_files(outputs(outputCount)%wireCurrentProbe)
               end if

            case (iQx, iQy, iQz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = WIRE_CHARGE_PROBE_ID

               allocate (outputs(outputCount)%wireChargeProbe)
               call init_solver_output(outputs(outputCount)%wireChargeProbe, lowerBound, NODE, outputRequestType, domain, outputTypeExtension, control)
               call create_empty_files(outputs(outputCount)%wireChargeProbe)

            case (iBloqueJx, iBloqueJy, iBloqueJz, iBloqueMx, iBloqueMy, iBloqueMz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = BULK_PROBE_ID

               allocate (outputs(outputCount)%bulkCurrentProbe)
               call init_solver_output(outputs(outputCount)%bulkCurrentProbe, lowerBound, upperBound, outputRequestType, domain, outputTypeExtension, control)
               call create_empty_files(outputs(outputCount)%bulkCurrentProbe)
               !! call adjust_computation_range --- Required due to issues in mpi region edges

            case (iCur, iMEC, iMHC, iCurX, iCurY, iCurZ, iExC, iEyC, iEyC, iHxC, iHyC, iHyC)
               call adjust_bound_range()

               if (domain%domainType == TIME_DOMAIN) then

                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = MOVIE_PROBE_ID
                  allocate (outputs(outputCount)%movieProbe)
                  call init_solver_output(outputs(outputCount)%movieProbe, lowerBound, upperBound, outputRequestType, domain, outputTypeExtension, control, problemInfo)
                  call create_pvd(outputs(outputCount)%movieProbe%path, outputs(outputCount)%movieProbe%PDVUnit)

               else if (domain%domainType == FREQUENCY_DOMAIN) then

                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = FREQUENCY_SLICE_PROBE_ID
                  allocate (outputs(outputCount)%frequencySliceProbe)
                  call init_solver_output(outputs(outputCount)%frequencySliceProbe, lowerBound, upperBound, sgg%dt, outputRequestType, domain, outputTypeExtension, control, problemInfo)
                  call create_pvd(outputs(outputCount)%frequencySliceProbe%path, outputs(outputCount)%frequencySliceProbe%PDVUnit)

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

      if (outputCount /= 0) observationsExists = .true.
      return
   contains
      subroutine adjust_bound_range()
         select case (outputRequestType)
         case (iExC, iEyC, iHzC, iMhC)
            lowerBound%z = max(sgg%Sweep(fieldo(field, 'Z'))%ZI, sgg%observation(ii)%P(i)%ZI)
            upperBound%z = min(sgg%Sweep(fieldo(field, 'Z'))%ZE - 1, sgg%observation(ii)%P(i)%ZE)
         case (iEzC, iHxC, iHyC, iMeC)
            lowerBound%z = max(sgg%Sweep(fieldo(field, 'Z'))%ZI, sgg%observation(ii)%P(i)%ZI)
            upperbound%z = min(sgg%Sweep(fieldo(field, 'Z'))%ZE, sgg%observation(ii)%P(i)%ZE)
         case (iCur, iCurX, iCurY, iCurZ)
            lowerBound%z = max(sgg%Sweep(fieldo(field, 'X'))%ZI, sgg%observation(ii)%P(i)%ZI) !ojo estaba sweep(iEz) para ser conservador...puede dar problemas!! 03/07/15
            upperbound%z = min(sgg%Sweep(fieldo(field, 'X'))%ZE, sgg%observation(ii)%P(i)%ZE) !ojo estaba sweep(iEz) para ser conservador...puede dar problemas!! 03/07/15
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
            call stoponerror(0, 0, 'No domain present')
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

   subroutine create_output_files()
      integer(kind=SINGLE) :: i
      do i = 1, size(outputs)
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID); call create_empty_files(outputs(i)%pointProbe)
         end select
      end do
   end subroutine create_output_files

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
            fieldComponent => get_field_component(outputs(i)%pointProbe%fieldComponent, fieldsReference) !Cada componente requiere de valores deiferentes pero estos valores no se como conseguirlos
            call update_solver_output(outputs(i)%pointProbe, discreteTime, fieldComponent)
         case (WIRE_CURRENT_PROBE_ID)
         call update_solver_output(outputs(i)%wireCurrentProbe, discreteTime, contorl, InvEps, InvMu)
         case (WIRE_CHARGE_PROBE_ID)
            call update_solver_output(outputs(i)%wireChargeProbe, discreteTime)
         case (BULK_PROBE_ID)
            fieldReference = get_field_reference(outputs(i)%bulkCurrentProbe%fieldComponent, fieldsReference)
            call update_solver_output(outputs(i)%bulkCurrentProbe, discreteTime, fieldReference)
         case (MOVIE_PROBE_ID)
       call update_solver_output(outputs(i)%movieProbe, discreteTime, problemInfo, fieldsReference)
         case (FREQUENCY_SLICE_PROBE_ID)
            call update_solver_output(outputs(i)%frequencySliceProbe, discreteTime, problemInfo, fieldsReference)
         case (FAR_FIELD_PROBE_ID)
            call update_solver_output(outputs(i)%farFieldOutput, timeIndx, problemInfo, fieldsReference)
         case default
            call stoponerror(0, 0, 'Output update not implemented')
         end select
      end do

   end subroutine update_outputs

   subroutine flush_outputs(simulationTimeArray, simulationTimeIndex, control, fields, bounds, farFieldFlushRequested)
      type(fields_reference_t), target :: fields
      type(fields_reference_t), pointer :: fieldsPtr
      type(sim_control_t), intent(in) :: control
      type(bounds_t), intent(in) :: bounds
      logical, intent(in) :: farFieldFlushRequested
      real(KIND=RKIND_tiempo), pointer, dimension(:), intent(in) :: simulationTimeArray
      integer, intent(in) :: simulationTimeIndex
      integer :: i

      fieldsPtr => fields

      do i = 1, size(outputs)
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID)
            call flush_solver_output(outputs(i)%pointProbe)
         case (WIRE_CURRENT_PROBE_ID)
            call flush_solver_output(outputs(i)%wireCurrentProbe)
         case (WIRE_CHARGE_PROBE_ID)
            call flush_solver_output(outputs(i)%wireChargeProbe)
         case (BULK_PROBE_ID)
            call flush_solver_output(outputs(i)%bulkCurrentProbe)
         case (MOVIE_PROBE_ID)
            call flush_solver_output(outputs(i)%movieProbe)
         case (FREQUENCY_SLICE_PROBE_ID)
            call flush_solver_output(outputs(i)%frequencySliceProbe)
         case (FAR_FIELD_PROBE_ID)
            if (farFieldFlushRequested) call flush_solver_output(outputs(i)%farFieldOutput, simulationTimeArray, simulationTimeIndex, control, fieldsPtr, bounds)
         case default
         end select
      end do
   end subroutine flush_outputs

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
            call close_pvd(outputs(i)%movieProbe%PDVUnit)
         case (FREQUENCY_SLICE_PROBE_ID)
         end select
      end do
   end subroutine

   subroutine create_pvd(pdvPath, unitPVD)
      implicit none
      character(len=*), intent(in) :: pdvPath
      integer, intent(out) :: unitPVD
      integer :: ios

      open (newunit=unitPVD, file=trim(pdvPath)//".pvd", status="replace", action="write", iostat=ios)
      if (ios /= 0) stop "Error al crear archivo PVD"

      ! Escribimos encabezados XML
      write (unitPVD, *) '<?xml version="1.0"?>'
      write (unitPVD, *) '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">'
      write (unitPVD, *) '  <Collection>'
   end subroutine create_pvd

   subroutine close_pvd(unitPVD)
      implicit none
      integer, intent(in) :: unitPVD

      write (unitPVD, *) '  </Collection>'
      write (unitPVD, *) '</VTKFile>'
      close (unitPVD)
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

 subroutine eliminate_unnecessary_observation_points(observation_probe, output_item, sweep, SINPMLSweep, ZI, ZE, layoutnumber, size)
      type(item_t), intent(inout) :: output_item
      type(observable_t), intent(inout) :: observation_probe
      type(XYZlimit_t), dimension(1:6), intent(in) :: sweep, SINPMLSweep
      integer(kind=4), intent(in) :: ZI, ZE, layoutnumber, size
      integer(kind=4) :: field

      ! Initialize output_item trancos
      output_item%Xtrancos = observation_probe%Xtrancos
      output_item%Ytrancos = observation_probe%Ytrancos
      output_item%Ztrancos = observation_probe%Ztrancos

      output_item%XItrancos = ceiling(real(observation_probe%XI)/real(output_item%Xtrancos))
      output_item%YItrancos = ceiling(real(observation_probe%YI)/real(output_item%Ytrancos))
      output_item%ZItrancos = ceiling(real(observation_probe%ZI)/real(output_item%Ztrancos))

      output_item%XEtrancos = int(observation_probe%XE/output_item%Xtrancos)
      output_item%YEtrancos = int(observation_probe%YE/output_item%Ytrancos)
      output_item%ZEtrancos = int(observation_probe%ZE/output_item%Ztrancos)

#ifdef CompileWithMPI
      output_item%MPISubComm = -1
#endif

      field = observation_probe%What

      select case (field)
      case (iBloqueJx, iBloqueJy, iBloqueMx, iBloqueMy, iExC, iEyC, iHzC, iMhC, iEzC, iHxC, iHyC, iMeC)
         call eliminate_observation_block(observation_probe, output_item, sweep, field, layoutnumber, size)
      case (iEx, iVx, iEy, iVy, iHz, iBloqueMz, iJx, iJy, iQx, iQy)
         call eliminate_observation_range(observation_probe, sweep, field, layoutnumber, size, lower_inclusive=.false.)
      case (iEz, iVz, iJz, iQz, iBloqueJz, iHx, iHy)
         call eliminate_observation_range(observation_probe, sweep, field, layoutnumber, size, lower_inclusive=.true.)
      case (iCur, iCurX, iCurY, iCurZ, mapvtk)
         call eliminate_observation_current(observation_probe, output_item, sweep, field, layoutnumber, size)
      case (FarField)
         call eliminate_observation_farfield(observation_probe, output_item, SINPMLSweep, ZI, ZE, layoutnumber, size)
      end select
   end subroutine

! Generic subroutine for block observations
   subroutine eliminate_observation_block(obs, out, sweep, field, layoutnumber, size)
      type(observable_t), intent(inout) :: obs
      type(item_t), intent(inout) :: out
      type(XYZlimit_t), dimension(1:6), intent(in) :: sweep
      integer, intent(in) :: field, layoutnumber, size

      call eliminate_observation_range_generic(obs, out, sweep(fieldo(field, 'Z'))%ZI, &
                                               sweep(fieldo(field, 'Z'))%ZE, layoutnumber, size)
   end subroutine

! Generic Z-range check with optional inclusive lower bound
   subroutine eliminate_observation_range(obs, sweep, field, layoutnumber, size, lower_inclusive)
      type(observable_t), intent(inout) :: obs
      type(XYZlimit_t), dimension(1:6), intent(in) :: sweep
      integer, intent(in) :: field, layoutnumber, size
      logical, intent(in) :: lower_inclusive

      if (lower_inclusive) then
         if ((obs%ZI > sweep(fieldo(field, 'Z'))%ZE) .or. (obs%ZI < sweep(fieldo(field, 'Z'))%ZI)) obs%What = nothing
      else
        if ((obs%ZI >= sweep(fieldo(field,'Z'))%ZE) .and. (layoutnumber /= size-1) .or. (obs%ZI < sweep(fieldo(field,'Z'))%ZI)) obs%What = nothing
      end if
   end subroutine

! Generic subroutine for currents
   subroutine eliminate_observation_current(obs, out, sweep, field, layoutnumber, size)
      type(observable_t), intent(inout) :: obs
      type(item_t), intent(inout) :: out
      type(XYZlimit_t), dimension(1:6), intent(in) :: sweep
      integer, intent(in) :: field, layoutnumber, size

  call eliminate_observation_range_generic(obs, out, sweep(fieldo(field, 'Z'))%ZI, sweep(fieldo(field, 'Z'))%ZE, layoutnumber, size)
      if ((field == iCur .or. field == iCurX .or. field == iCurY .or. field == mapvtk)) then
         obs%ZE = min(obs%ZE, sweep(iHx)%ZE)
      end if
   end subroutine

! Far field specialized
   subroutine eliminate_observation_farfield(obs, out, sweep, ZI, ZE, layoutnumber, size)
      type(observable_t), intent(inout) :: obs
      type(item_t), intent(inout) :: out
      type(XYZlimit_t), dimension(1:6), intent(in) :: sweep
      integer(kind=4), intent(in) :: ZI, ZE, layoutnumber, size

      call eliminate_observation_range_generic(obs, out, sweep(iHz)%ZI, sweep(iHz)%ZE, layoutnumber, size, ZI, ZE)
   end subroutine

! The ultimate generic routine for MPI and Z-limits
   subroutine eliminate_observation_range_generic(obs, out, Z_lower, Z_upper, layoutnumber, size, Zstart, Zend)
      type(observable_t), intent(inout) :: obs
      type(item_t), intent(inout) :: out
      integer, intent(in) :: Z_lower, Z_upper, layoutnumber, size
      integer, optional, intent(in) :: Zstart, Zend

      integer :: zi_local, ze_local
      zi_local = merge(Zstart, obs%ZI, present(Zstart))
      ze_local = merge(Zend, obs%ZE, present(Zend))

      if ((zi_local > Z_upper) .or. (ze_local < Z_lower)) then
         obs%What = nothing
#ifdef CompileWithMPI
         out%MPISubComm = -1
      else
         out%MPISubComm = 1
      end if
      out%MPIRoot = 0
      if ((obs%ZI >= Z_lower) .and. (obs%ZI <= Z_upper)) out%MPIRoot = layoutnumber
      call MPIinitSubcomm(layoutnumber, size, out%MPISubComm, out%MPIRoot, out%MPIGroupIndex)
#endif
      end if
      end subroutine

   end module output
