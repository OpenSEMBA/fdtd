module output
   use FDETYPES
   use mod_domain
   use mod_outputUtils
   use mod_pointProbeOutput
   use mod_wireCurrentProbeOutput
   use mod_wireChargeProbeOutput
   use mod_bulkProbe
   use mod_volumicProbe

   implicit none

   integer(kind=SINGLE), parameter :: POINT_PROBE_ID = 0, &
                                      WIRE_CURRENT_PROBE_ID = 1, &
                                      WIRE_CHARGE_PROBE_ID = 2, &
                                      BULK_PROBE_ID = 3, &
                                      VOLUMIC_CURRENT_PROBE_ID = 4

   REAL(KIND=RKIND), save           ::  eps0, mu0
   REAL(KIND=RKIND), pointer, dimension(:), save  ::  InvEps, InvMu

   type solver_output_t
      integer(kind=SINGLE) :: outputID
      type(point_probe_output_t), allocatable :: pointProbe
      type(wire_current_probe_output_t), allocatable :: wireCurrentProbe
      type(wire_charge_probe_output_t), allocatable :: wireChargeProbe
      type(bulk_probe_output_t), allocatable :: bulkProbe
      type(volumic_current_probe_t), allocatable :: volumicCurrentProbe
      !type(volumic_field_probe_t), allocatable :: volumicFieldProbe
      !type(bulk_current_probe_output_t), allocatable :: bulkCurrentProbe
      !type(far_field_t), allocatable :: farField
      !type(time_movie_output_t), allocatable :: timeMovie
      !type(frequency_slice_output_t), allocatable :: frequencySlice
   end type solver_output_t

   interface init_solver_output
      module procedure &
         init_point_probe_output, &
         init_wire_current_probe_output, &
         init_wire_charge_probe_output, &
         init_bulk_probe_output, &
         init_volumic_probe_output
      !init_far_field, &
      !initime_movie_output, &
      !init_frequency_slice_output
   end interface

   interface update_solver_output
      module procedure &
         update_point_probe_output, &
         update_wire_current_probe_output, &
         update_wire_charge_probe_output, &
         update_bulk_probe_output
      !update_bulk_current_probe_output, &
      !update_far_field, &
      !updateime_movie_output, &
      !update_frequency_slice_output
   end interface

   interface flush_solver_output
      module procedure &
         flush_point_probe_output
      !flush_wire_probe_output, &
      !flush_bulk_current_probe_output, &
      !flush_far_field, &
      !flushime_movie_output, &
      !flush_frequency_slice_output
   end interface

   interface delete_solver_output
      module procedure &
         delete_point_probe_output
      !delete_wire_probe_output, &
      !delete_bulk_current_probe_output, &
      !delete_far_field, &
      !deleteime_movie_output, &
      !delete_frequency_slice_output
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
      integer(kind=SINGLE) :: i, ii, outputRequestType
      integer(kind=SINGLE) :: I1, J1, K1, I2, J2, K2, NODE
      integer(kind=SINGLE) :: outputCount
      character(len=BUFSIZE) :: outputTypeExtension

      allocate (outputs(sgg%NumberRequest))

      allocate (InvEps(0:sgg%NumMedia), InvMu(0:sgg%NumMedia))
      outputCount = 0

      InvEps(0:sgg%NumMedia) = 1.0_RKIND/(Eps0*sgg%Med(0:sgg%NumMedia)%Epr)
      InvMu(0:sgg%NumMedia) = 1.0_RKIND/(Mu0*sgg%Med(0:sgg%NumMedia)%Mur)

      do ii = 1, sgg%NumberRequest
         do i = 1, sgg%Observation(ii)%nP
            I1 = sgg%observation(ii)%P(i)%XI
            J1 = sgg%observation(ii)%P(i)%YI
            K1 = sgg%observation(ii)%P(i)%ZI
            I2 = sgg%observation(ii)%P(i)%XE
            J2 = sgg%observation(ii)%P(i)%YE
            K2 = sgg%observation(ii)%P(i)%ZE
            NODE = sgg%observation(ii)%P(i)%NODE

            domain = preprocess_domain(sgg%Observation(ii), sgg%tiempo, sgg%dt, control%finaltimestep)
            outputTypeExtension = trim(adjustl(control%nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))

            outputRequestType = sgg%observation(ii)%P(i)%what
            select case (outputRequestType)
            case (iEx, iEy, iEz, iHx, iHy, iHz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = POINT_PROBE_ID

               allocate (outputs(outputCount)%pointProbe)
               call init_solver_output(outputs(outputCount)%pointProbe, I1, J1, K1, outputRequestType, domain, outputTypeExtension, control%mpidir)

            case (iJx, iJy, iJz)
               if (ThereAreWires) then
                  outputCount = outputCount + 1
                  outputs(outputCount)%outputID = WIRE_CURRENT_PROBE_ID

                  allocate (outputs(outputCount)%wireCurrentProbe)
                  call init_solver_output(outputs(outputCount)%wireCurrentProbe, I1, J1, K1, NODE, outputRequestType, domain, sgg%Med, outputTypeExtension, control%mpidir, control%wiresflavor)
               end if

            case (iQx, iQy, iQz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = WIRE_CHARGE_PROBE_ID

               allocate (outputs(outputCount)%wireChargeProbe)
               call init_solver_output(outputs(outputCount)%wireChargeProbe, I1, J1, K1, NODE, outputRequestType, domain, outputTypeExtension, control%mpidir, control%wiresflavor)
            case (iBloqueJx, iBloqueJy, iBloqueJz, iBloqueMx, iBloqueMy, iBloqueMz)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = BULK_PROBE_ID

               allocate (outputs(outputCount)%bulkProbe)
               call init_solver_output(outputs(outputCount)%bulkProbe, I1, J1, K1, I2, J2, K2, outputRequestType, domain, outputTypeExtension, control%mpidir)
               !! call adjust_computation_range --- Required due to issues in mpi region edges

            case (iCur, iCurX, iCurY, iCurZ)
               outputCount = outputCount + 1
               outputs(outputCount)%outputID = VOLUMIC_CURRENT_PROBE_ID

               allocate (outputs(outputCount)%volumicCurrentProbe)
               call init_solver_output(outputs(outputCount)%volumicCurrentProbe, I1, J1, K1, I2, J2, K2, outputRequestType, domain, media, sgg%Med, sinpml_fullsize, outputTypeExtension, control%mpidir)

            case default
               call stoponerror(0, 0, 'OutputRequestType type not implemented yet on new observations')
            end select
         end do
      end do
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
            nFreq = int((observation%FinalFreq - observation%InitialFreq)/observation%FreqStep, kind=SINGLE)
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

   end subroutine init_outputs

   subroutine update_outputs(outputs, control, step, fields)
      type(solver_output_t), dimension(:), intent(inout) :: outputs
      real(kind=RKIND_tiempo) :: step
      integer(kind=SINGLE) :: i, id
      type(sim_control_t), intent(in) :: control
      real(kind=RKIND), pointer, dimension(:, :, :) :: fieldComponent
      type(field_data_t), pointer :: fieldReference
      type(fields_reference_t) :: fields

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
            fieldReference => get_field_reference(outputs(i)%bulkProbe%fieldComponent, fields)
            call update_solver_output(outputs(i)%bulkProbe, step, fieldReference)
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

end module output
