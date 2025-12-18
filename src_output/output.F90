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

   type solver_output_t
      integer(kind=SINGLE) :: outputID
      type(point_probe_output_t), allocatable :: pointProbe !iEx, iEy, iEz, iHx, iHy, iHz
      type(wire_current_probe_output_t), allocatable :: wireCurrentProbe !Jx, Jy, Jz
      type(wire_charge_probe_output_t), allocatable :: wireChargeProbe !Qx, Qy, Qz
      type(bulk_current_probe_output_t), allocatable :: bulkCurrentProbe !BloqueXJ, BloqueYJ, BloqueZJ, BloqueXM, BloqueYM, BloqueZM
      type(volumic_current_probe_t), allocatable :: volumicCurrentProbe !icurX, icurY, icurZ
      type(volumic_field_probe_output_t), allocatable :: volumicFieldProbe
      type(line_integral_probe_output_t), allocatable :: lineIntegralProbe
      type(movie_probe_output_t), allocatable :: movieProbe !iCur if timeDomain
      type(frequency_slice_probe_output_t), allocatable :: frequencySliceProbe !iCur if freqDomain
      type(far_field_probe_output_t), allocatable :: farFieldOutput !farfield
   end type solver_output_t

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

   subroutine init_outputs(sgg, media, sinpml_fullsize, control, outputs, ThereAreWires)
      type(SGGFDTDINFO), intent(in) ::  sgg
      type(media_matrices_t), pointer, intent(in) :: media
      type(limit_t), pointer, dimension(:), intent(in)  ::  SINPML_fullsize
      type(sim_control_t), intent(inout) :: control
      type(solver_output_t), dimension(:), allocatable, intent(out) :: outputs
      logical :: ThereAreWires

      type(domain_t) :: domain
      type(spheric_domain_t) :: sphericRange
      type(cell_coordinate_t) :: lowerBound, upperBound
      integer(kind=SINGLE) :: i, ii, outputRequestType
      integer(kind=SINGLE) :: NODE
      integer(kind=SINGLE) :: outputCount
      integer(kind=SINGLE) :: requestedOutputs
      character(len=BUFSIZE) :: outputTypeExtension

      OutputRequested = .false.
      requestedOutputs = get_required_output_count(sgg)

      outputs => NULL()
      allocate (outputs(requestedOutputs))

      allocate (InvEps(0:sgg%NumMedia - 1), InvMu(0:sgg%NumMedia - 1))
      outputCount = 0

      InvEps(0:sgg%NumMedia - 1) = 1.0_RKIND/(Eps0*sgg%Med(0:sgg%NumMedia - 1)%Epr)
      InvMu(0:sgg%NumMedia - 1) = 1.0_RKIND/(Mu0*sgg%Med(0:sgg%NumMedia - 1)%Mur)

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
               call init_solver_output(outputs(outputCount)%pointProbe, lowerBound, outputRequestType, domain, outputTypeExtension, control%mpidir, sgg%dt)
               call create_empty_files(outputs(outputCount)%pointProbe)
            case (iJx, iJy, iJz)
               if (ThereAreWires) then
                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = WIRE_CURRENT_PROBE_ID

                  allocate (outputs(outputCount)%wireCurrentProbe)
                  call init_solver_output(outputs(outputCount)%wireCurrentProbe, lowerBound, NODE, outputRequestType, domain, sgg%Med, outputTypeExtension, control%mpidir, control%wiresflavor)
                  call create_empty_files(outputs(outputCount)%wireCurrentProbe)
               end if

            case (iQx, iQy, iQz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = WIRE_CHARGE_PROBE_ID

               allocate (outputs(outputCount)%wireChargeProbe)
               call init_solver_output(outputs(outputCount)%wireChargeProbe, lowerBound, NODE, outputRequestType, domain, outputTypeExtension, control%mpidir, control%wiresflavor)
               call create_empty_files(outputs(outputCount)%wireChargeProbe)

            case (iBloqueJx, iBloqueJy, iBloqueJz, iBloqueMx, iBloqueMy, iBloqueMz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = BULK_PROBE_ID

               allocate (outputs(outputCount)%bulkCurrentProbe)
               call init_solver_output(outputs(outputCount)%bulkCurrentProbe, lowerBound, upperBound, outputRequestType, domain, outputTypeExtension, control%mpidir)
               call create_empty_files(outputs(outputCount)%bulkCurrentProbe)
               !! call adjust_computation_range --- Required due to issues in mpi region edges

            case (iCurX, iCurY, iCurZ)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = VOLUMIC_CURRENT_PROBE_ID

               allocate (outputs(outputCount)%volumicCurrentProbe)
               call init_solver_output(outputs(outputCount)%volumicCurrentProbe, lowerBound, upperBound, outputRequestType, domain, media, sgg%Med, sinpml_fullsize, outputTypeExtension, control%mpidir, sgg%dt)
               
            case (iCur)
               if (domain%domainType == TIME_DOMAIN) then

                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = MOVIE_PROBE_ID
                  allocate (outputs(outputCount)%movieProbe)
                  call init_solver_output(outputs(outputCount)%movieProbe, lowerBound, upperBound, outputRequestType, domain, media, sgg%Med, SINPML_fullsize, outputTypeExtension, control%mpidir)
                  call create_pvd(outputs(outputCount)%movieProbe%path, outputs(outputCount)%movieProbe%PDVUnit)

               else if ( domain%domainType == FREQUENCY_DOMAIN ) then

                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = FREQUENCY_SLICE_PROBE_ID
                  allocate (outputs(outputCount)%frequencySliceProbe)
                  call init_solver_output(outputs(outputCount)%frequencySliceProbe, lowerBound, upperBound, outputRequestType, domain, media, sgg%Med, SINPML_fullsize, outputTypeExtension, control%mpidir, sgg%dt)
                  call create_pvd(outputs(outputCount)%frequencySliceProbe%path, outputs(outputCount)%frequencySliceProbe%PDVUnit)

               end if
            case (farfield)
               sphericRange = preprocess_polar_range(sgg%Observation(ii))

               outputCount = outputCount + 1
               outputs(outputCount)%outputID = FAR_FIELD_PROBE_ID
               allocate (outputs(outputCount)%farFieldOutput)
               call init_solver_output(outputs(outputCount)%farFieldOutput, sgg, lowerBound, upperBound,outputRequestType, domain, sphericRange, control, outputTypeExtension, sgg%Observation(ii)%FileNormalize, eps0, mu0, media, SINPML_fullsize, bounds)
            case default
               call stoponerror(0, 0, 'OutputRequestType type not implemented yet on new observations')
            end select
         end do
      end do

      if (outputCount /= 0) OutputRequested = .true.
      return
   contains
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

   subroutine create_output_files(outputs)
      type(solver_output_t), dimension(:), intent(inout) :: outputs
      integer(kind=SINGLE) :: i
      do i = 1, size(outputs)
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID); call create_empty_files(outputs(i)%pointProbe)
         end select
      end do
   end subroutine create_output_files

   subroutine update_outputs(outputs, geometryMedia, materialList, SINPML_fullsize , control, step, fields)
      type(solver_output_t), dimension(:), intent(inout) :: outputs
      real(kind=RKIND_tiempo) :: step
      integer(kind=SINGLE) :: i, id
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t),dimension(:), pointer :: materialList
      type(limit_t), pointer, dimension(:), intent(in)  ::  SINPML_fullsize
      type(sim_control_t), intent(in) :: control
      real(kind=RKIND), pointer, dimension(:, :, :) :: fieldComponent
      type(field_data_t), pointer :: fieldReference
      type(fields_reference_t), target :: fields
      type(fields_reference_t), pointer :: fieldsPtr

      fieldsPtr => fields

      do i = 1, size(outputs)
         select case (outputs(i)%outputID)
         case (POINT_PROBE_ID)
            fieldComponent => get_field_component(outputs(i)%pointProbe%fieldComponent, fields) !Cada componente requiere de valores deiferentes pero estos valores no se como conseguirlos
            call update_solver_output(outputs(i)%pointProbe, step, fieldComponent)
         case (WIRE_CURRENT_PROBE_ID)
            call update_solver_output(outputs(i)%wireCurrentProbe, step, control%wiresflavor, control%wirecrank, InvEps, InvMu)
         case (WIRE_CHARGE_PROBE_ID)
            call update_solver_output(outputs(i)%wireChargeProbe, step)
         case (BULK_PROBE_ID)
            fieldReference => get_field_reference(outputs(i)%bulkCurrentProbe%fieldComponent, fields)
            call update_solver_output(outputs(i)%bulkCurrentProbe, step, fieldReference)
         case (MOVIE_PROBE_ID)
            call update_solver_output(outputs(i)%movieProbe, step, geometryMedia, materialList, SINPML_fullsize, fieldsPtr)
         case(FREQUENCY_SLICE_PROBE_ID)
            call update_solver_output(outputs(i)%frequencySliceProbe, step, geometryMedia, materialList, SINPML_fullsize, fieldsPtr)
         case default
            call stoponerror(0, 0, 'Output update not implemented')
         end select
      end do

   contains
      function get_field_component(fieldId, fieldsReference) result(field)
         integer(kind=SINGLE), intent(in) :: fieldId
         type(fields_reference_t), intent(in) :: fieldsReference
         real(kind=RKIND), pointer, dimension(:, :, :) :: field
         select case (fieldId)
         case (iEx); field => fieldsReference%E%x
         case (iEy); field => fieldsReference%E%y
         case (iEz); field => fieldsReference%E%z
         case (iHx); field => fieldsReference%H%x
         case (iHy); field => fieldsReference%H%y
         case (iHz); field => fieldsReference%H%z
         end select
      end function get_field_component

      function get_field_reference(fieldId, fieldsReference) result(field)
         integer(kind=SINGLE), intent(in) :: fieldId
         type(fields_reference_t), intent(in) :: fieldsReference
         type(field_data_t), pointer :: field
         select case (fieldId)
         case (iBloqueJx, iBloqueJy, iBloqueJz)
            field%x => fieldsReference%E%x
            field%y => fieldsReference%E%y
            field%z => fieldsReference%E%z

            field%deltaX => fieldsReference%E%deltax
            field%deltaY => fieldsReference%E%deltay
            field%deltaZ => fieldsReference%E%deltaz
         case (iBloqueMx, iBloqueMy, iBloqueMz)
            field%x => fieldsReference%H%x
            field%y => fieldsReference%H%y
            field%z => fieldsReference%H%z

            field%deltaX => fieldsReference%H%deltax
            field%deltaY => fieldsReference%H%deltay
            field%deltaZ => fieldsReference%H%deltaz
         end select
      end function get_field_reference

   end subroutine update_outputs

   subroutine flush_outputs(outputs)
      type(solver_output_t), dimension(:), intent(inout) :: outputs
      integer :: i
      do i = 1, size(outputs)
         select case(outputs(i)%outputID)
         case(POINT_PROBE_ID)
            call flush_point_probe_output(outputs(i)%pointProbe)
         case(WIRE_CURRENT_PROBE_ID)
         case(WIRE_CHARGE_PROBE_ID)
         case(BULK_PROBE_ID)
         case(VOLUMIC_CURRENT_PROBE_ID)
         case(MOVIE_PROBE_ID)
            call flush_solver_output(outputs(i)%movieProbe)
         case(FREQUENCY_SLICE_PROBE_ID)
            call flush_solver_output(outputs(i)%frequencySliceProbe)
         end select
      end do
   end subroutine flush_outputs

   subroutine close_outputs(outputs)
      type(solver_output_t), dimension(:), intent(inout) :: outputs
      integer :: i
      do i = 1, size(outputs)
         select case(outputs(i)%outputID)
         case(POINT_PROBE_ID)
         case(WIRE_CURRENT_PROBE_ID)
         case(WIRE_CHARGE_PROBE_ID)
         case(BULK_PROBE_ID)
         case(VOLUMIC_CURRENT_PROBE_ID)
         case(MOVIE_PROBE_ID)
            call close_pvd(outputs(i)%movieProbe%PDVUnit)
         case(FREQUENCY_SLICE_PROBE_ID)
         end select
      end do
   end subroutine


   subroutine create_pvd(pdvPath, unitPVD)
      implicit none
      character(len=*), intent(in) :: pdvPath
      integer, intent(out) :: unitPVD
      integer :: ios

      open(newunit=unitPVD, file=trim(pdvPath)//".pvd", status="replace", action="write", iostat=ios)
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

end module output
