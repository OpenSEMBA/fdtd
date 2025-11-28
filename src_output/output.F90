module output
   use FDETYPES
   use mod_domain
   use mod_outputUtils
   use mod_pointProbeOutput
   use wiresHolland_constants
   use HollandWires
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

   
   type current_values_t
      real(kind=RKIND) :: current = 0.0_RKIND, deltaVoltage = 0.0_RKIND
      real(kind=RKIND) :: plusVoltage = 0.0_RKIND, minusVoltage = 0.0_RKIND, voltageDiference = 0.0_RKIND
   end type
   type wire_current_probe_output_t
      integer(kind=SINGLE) :: columnas = 6_SINGLE !reference, corriente, -e*dl, vplus, vminus, vplus-vminus
      type(domain_t) :: domain
      integer(kind=SINGLE) :: xCoord, yCoord, zCoord
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: currentComponent
      integer(kind=SINGLE) :: sign = +1
      type(CurrentSegments), pointer :: segment

      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(BuffObse) :: timeStep = 0.0_RKIND
      type(current_values_t), dimension(BuffObse) :: currentValues
   end type wire_current_probe_output_t

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
      !update_wire_probe_output, &
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

   

   subroutine init_wire_current_probe_output(this, iCoord, jCoord, kCoord, node, field, domain, media, outputTypeExtension, mpidir, wiresflavor)
      type(wire_current_probe_output_t), intent(out) :: this
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord, node
      integer(kind=SINGLE), intent(in) :: field, mpidir
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      character(len=*), intent(in) :: wiresflavor
      type(domain_t), intent(in) :: domain
      type(MediaData_t), pointer, dimension(:), intent(in) :: media

      

      type(Thinwires_t), pointer  ::  Hwireslocal
#ifdef CompileWithBerengerWires
      type(TWires), pointer  ::  Hwireslocal_Berenger
#endif
#ifdef CompileWithSlantedWires
      type(WiresData), pointer  ::  Hwireslocal_Slanted
#endif

      select case (trim(adjustl(wiresflavor)))
      case ('holland', 'transition'); Hwireslocal => GetHwires()
#ifdef CompileWithBerengerWires
      case ('berenger'); Hwireslocal_Berenger => GetHwires_Berenger()
#endif
#ifdef CompileWithSlantedWires
      case ('slanted', 'semistructured'); Hwireslocal_Slanted => GetHwires_Slanted()
#endif
      end select

      call find_segment()

      this%xCoord = iCoord
      this%yCoord = jCoord
      this%zCoord = kCoord

      this%currentComponent = field

      this%domain = domain
      this%path = get_output_path()

   contains
      subroutine find_segment()
         integer(kind=SINGLE) :: n, iwi, iwj, node2
         type(CurrentSegments), pointer :: currentSegment
         logical :: found = .false.
         character(len=BUFSIZE) :: buff

         select case (trim(adjustl(wiresflavor)))
         case ('holland', 'transition')
            this%segment => HWireslocal%NullSegment
            do n = 1, HWireslocal%NumCurrentSegments
               currentSegment => HWireslocal%CurrentSegment(n)
               if ((currentSegment%origindex == node) .and. &
                   (currentSegment%i == iCoord) .and. (currentSegment%j == jCoord) .and. (currentSegment%k == kCoord) .and. &
                   (currentSegment%tipofield*10 == field)) then
                  found = .true.
                  this%segment => currentSegment
                  if (currentSegment%orientadoalreves) this%sign = -1
               end if
            end do
#ifdef CompileWithBerengerWires
         case ('berenger')
            do n = 1, Hwireslocal_Berenger%NumSegments
               currentSegment => Hwireslocal_Berenger%Segments(n)
               if (currentSegment%IndexSegment == node) then
                  found = .true.
                  this%segmentBerenger => currentSegment
                  if (currentSegment%orientadoalreves) this%sign = -1
               end if
            end do
#endif
#ifdef CompileWithSlantedWires
         case ('slanted', 'semistructured')
            do n = 1, Hwireslocal_Slanted%NumSegments
               currentSegment => Hwireslocal_Slanted%Segments(n)
               if (currentSegment%ptr%Index == node) then
                  found = .true.
                  this%segmentSlanted => currentSegment%ptr
               end if
            end do
#endif
         end select

         if (.not. found) then
            select case (trim(adjustl(wiresflavor)))
            case ('holland', 'transition')
               buscarabono: do iwi = 1, Hwireslocal%NumDifferentWires
                  do iwj = 1, media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%numsegmentos
                     if ((node == media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex) .and. &
                         media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                        node2 = media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multiraboDE
                        do n = 1, HWireslocal%NumCurrentSegments
                           currentSegment => HWireslocal%CurrentSegment(n)
                           if (currentSegment%origindex == node2) then
                              found = .true.
                              this%segment => currentSegment
                              if (currentSegment%orientadoalreves) this%sign = -1
                           end if
                        end do
                        exit buscarabono
                     end if
                  end do
               end do buscarabono
#ifdef CompileWithSlantedWires
            case ('slanted', 'semistructured')
               do n = 1, Hwireslocal_Slanted%NumSegments
                  currentSegment => Hwireslocal_Slanted%Segments(n)
                  if (currentSegment%ptr%elotroindice == node) then
                     found = .true.
                     this%segmentSlanted => currentSegment%ptr
                  end if
               end do
#endif
            end select
         end if

         if (.not. found) then
            write (buff, '(a,4i7,a)') 'ERROR: WIRE probe ', node, iCoord, jCoord, kCoord, ' DOES NOT EXIST'
            CALL WarnErrReport(buff, .true.)
         end if
      end subroutine find_segment

      function get_output_path() result(outputPath)
         character(len=BUFSIZE) :: outputPath
         character(len=BUFSIZE)  ::  charNO
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension, prefixNodeExtension

         write (charNO, '(i7)') node
         prefixNodeExtension = 's'//trim(adjustl(charNO))
         probeBoundsExtension = get_probe_bounds_extension()
         prefixFieldExtension = get_prefix_extension(field, mpidir)

         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_' &
            //trim(adjustl(probeBoundsExtension))//'_'//trim(adjustl(prefixNodeExtension))
         return
      end function get_output_path

      function get_probe_bounds_extension() result(ext)
         character(len=BUFSIZE) :: ext
         character(len=BUFSIZE)  ::  chari, charj, chark

         write (chari, '(i7)') iCoord
         write (charj, '(i7)') jCoord
         write (chark, '(i7)') kCoord

#if CompileWithMPI
         if (mpidir == 3) then
            ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
         elseif (mpidir == 2) then
            ext = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))
         elseif (mpidir == 1) then
            ext = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))
         else
            call stoponerror('Buggy error in mpidir. ')
         end if
#else
         ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
#endif

         return
      end function get_probe_bounds_extension

   end subroutine init_wire_current_probe_output

   

end module output
