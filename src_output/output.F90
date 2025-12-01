module output
   use FDETYPES
   use mod_domain
   use mod_outputUtils
   use mod_pointProbeOutput
   use mod_wireCurrentProbeOutput

   implicit none
   
   

   integer(kind=SINGLE), parameter :: POINT_PROBE_ID = 0, &
                                      WIRE_CURRENT_PROBE_ID = 0

   type solver_output_t
      integer(kind=SINGLE) :: outputID
      type(point_probe_output_t), allocatable :: pointProbe
      type(wire_current_probe_output_t), allocatable :: wireCurrentProbe
      !type(bulk_current_probe_output_t), allocatable :: bulkCurrentProbe
      !type(far_field_t), allocatable :: farField
      !type(time_movie_output_t), allocatable :: timeMovie
      !type(frequency_slice_output_t), allocatable :: frequencySlice
   end type solver_output_t


   interface init_solver_output
      module procedure &
         init_point_probe_output, &
         init_wire_current_probe_output
      !init_bulk_current_probe_output, &
      !init_far_field, &
      !initime_movie_output, &
      !init_frequency_slice_output
   end interface

   interface update_solver_output
      module procedure &
         update_point_probe_output
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

   subroutine init_outputs(sgg, control, outputs, ThereAreWires)
      type(SGGFDTDINFO), intent(in) ::  sgg
      type(sim_control_t), intent(inout) :: control
      type(solver_output_t), dimension(:), allocatable, intent(out) :: outputs
      logical :: ThereAreWires

      type(domain_t) :: domain
      integer(kind=SINGLE) :: i, ii, outputRequestType
      integer(kind=SINGLE) :: I1, J1, K1, NODE
      integer(kind=SINGLE) :: outputCount = 0
      character(len=BUFSIZE) :: outputTypeExtension
      allocate (outputs(sgg%NumberRequest))

      call retrive_wires_data()

      do ii = 1, sgg%NumberRequest
         do i = 1, sgg%Observation(ii)%nP
            I1 = sgg%observation(ii)%P(i)%XI
            J1 = sgg%observation(ii)%P(i)%YI
            K1 = sgg%observation(ii)%P(i)%ZI
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
            case default
               call stoponerror('OutputRequestType type not implemented yet on new observations')
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
            call stoponerror('No domain present')
         end if
         return
      end function preprocess_domain

   end subroutine init_outputs

   subroutine update_outputs(outputs, step, Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh, alloc)
      type(solver_output_t), dimension(:), intent(inout) :: outputs
      real(kind=RKIND_tiempo) :: step
      integer(kind=SINGLE) :: i, id
      type(XYZlimit_t), dimension(1:6), intent(in) :: alloc
      real(kind=RKIND), pointer, dimension(:, :, :) :: fieldPointer

      real(KIND=RKIND), intent(in), target     :: &
         Ex(alloc(iEx)%XI:alloc(iEx)%XE, alloc(iEx)%YI:alloc(iEx)%YE, alloc(iEx)%ZI:alloc(iEx)%ZE), &
         Ey(alloc(iEy)%XI:alloc(iEy)%XE, alloc(iEy)%YI:alloc(iEy)%YE, alloc(iEy)%ZI:alloc(iEy)%ZE), &
         Ez(alloc(iEz)%XI:alloc(iEz)%XE, alloc(iEz)%YI:alloc(iEz)%YE, alloc(iEz)%ZI:alloc(iEz)%ZE), &
         Hx(alloc(iHx)%XI:alloc(iHx)%XE, alloc(iHx)%YI:alloc(iHx)%YE, alloc(iHx)%ZI:alloc(iHx)%ZE), &
         Hy(alloc(iHy)%XI:alloc(iHy)%XE, alloc(iHy)%YI:alloc(iHy)%YE, alloc(iHy)%ZI:alloc(iHy)%ZE), &
         Hz(alloc(iHz)%XI:alloc(iHz)%XE, alloc(iHz)%YI:alloc(iHz)%YE, alloc(iHz)%ZI:alloc(iHz)%ZE)
      !--->
      real(KIND=RKIND), dimension(:), intent(in)   :: dxh(alloc(iEx)%XI:alloc(iEx)%XE), &
                                                      dyh(alloc(iEy)%YI:alloc(iEy)%YE), &
                                                      dzh(alloc(iEz)%ZI:alloc(iEz)%ZE), &
                                                      dxe(alloc(iHx)%XI:alloc(iHx)%XE), &
                                                      dye(alloc(iHy)%YI:alloc(iHy)%YE), &
                                                      dze(alloc(iHz)%ZI:alloc(iHz)%ZE)

      do i = 1, size(outputs)
         id = outputs(i)%outputID
         select case (id)
         case (POINT_PROBE_ID)
            fieldPointer => get_field_component(outputs(i)%pointProbe%fieldComponent) !Cada componente requiere de valores deiferentes pero estos valores no se como conseguirlos
            call update_solver_output(outputs(i)%pointProbe, step, fieldPointer)
         case default
            call stoponerror('Output update not implemented')
         end select
      end do

   contains
      function get_field_component(fieldId) result(field)
         integer(kind=SINGLE), intent(in) :: fieldId
         real(kind=RKIND), pointer, dimension(:, :, :) :: field
         select case (fieldId)
         case (iEx); field => Ex
         case (iEy); field => Ey
         case (iEz); field => Ez
         case (iHx); field => Hx
         case (iHy); field => Hy
         case (iHz); field => Hz
         end select
      end function get_field_component

   end subroutine update_outputs


end module output
