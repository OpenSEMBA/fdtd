module output
   use FDETYPES
   use mod_domain
   use mod_outputUtils
   implicit none
   character(len=4) :: datFileExtension = '.dat', timeExtension = 'tm', frequencyExtension = 'fq'
   integer(kind=SINGLE), parameter :: MAX_SERIALIZED_COUNT = 500, FILE_UNIT = 400

   integer(kind=SINGLE), parameter :: POINT_PROBE_ID = 0

   type solver_output_t
      integer(kind=SINGLE) :: outputID
      type(point_probe_output_t), allocatable :: pointProbe
      !type(wire_probe_output_t), allocatable :: wireProbe
      !type(bulk_current_probe_output_t), allocatable :: bulkCurrentProbe
      !type(far_field_t), allocatable :: farField
      !type(time_movie_output_t), allocatable :: timeMovie
      !type(frequency_slice_output_t), allocatable :: frequencySlice
   end type solver_output_t

   type point_probe_output_t
      integer(kind=SINGLE) :: columnas = 2_SINGLE
      type(domain_t) :: domain
      integer(kind=SINGLE) :: xCoord, yCoord, zCoord
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent
      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE, nFreq = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(MAX_SERIALIZED_COUNT) :: timeStep
      real(kind=RKIND), dimension(MAX_SERIALIZED_COUNT) :: valueForTime
      real(kind=RKIND), dimension(:), allocatable :: frequencySlice
      real(kind=CKIND), dimension(:), allocatable :: valueForFreq
   end type point_probe_output_t

   interface init_solver_output
      module procedure &
         init_point_probe_output
      !init_wire_probe_output, &
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

   subroutine init_outputs(sgg, control, outputs)
      type(SGGFDTDINFO), intent(in) ::  sgg
      type(sim_control_t), intent(inout) :: control
      type(solver_output_t), dimension(:), intent(out) :: outputs

      integer(kind=SINGLE) :: outputCount = 0
      allocate (outputs(sgg%NumberRequest))

      do ii = 1, sgg%NumberRequest
         do i = 1, sgg%Observation(ii)%nP
            I1 = sgg%observation(ii)%P(i)%XI
            J1 = sgg%observation(ii)%P(i)%YI
            K1 = sgg%observation(ii)%P(i)%ZI

            field = sgg%observation(ii)%P(i)%what
            select case (field)
            case (iEx, iEy, iEz, iVx, iVy, iVz, iJx, iJy, iJz, iQx, iQy, iQz, iHx, iHy, iHz, lineIntegral)
               outputCount = outputCount + 1

               outputs(outputCount)%outputID = POINT_PROBE_ID

               domain = preprocess_domain(sgg%Observation(ii), sgg%tiempo, sgg%dt, finaltimestep)

               outputTypeExtension = trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))

               allocate (outputs(outputCount)%pointProbe)
               init_solver_output(outputs(outputCount)%pointProbe, I1, J1, K1, field, domain, outputTypeExtension, mpidir)
            case default
               call stoponerror('Field type not implemented yet on new observations')
            end select
         end do
      end do
      return
   contains
      function preprocess_domain(observation, timeArray, timeStep, finalStepIndex) result(newDomain)
         type(Obses_t), intent(in) :: observation
         real(kind=RKIND_tiempo), pointer, dimension(:), intent(in) :: timeArray
         real(kind=RKIND_tiempo), intent(in) :: timeStep
         integer(kind=4), intent(in) :: finalStepIndex
         type(domain_t) :: newDomain

         integer(kind=SINGLE) :: nFreq

         if (observation%TimeDomain) then
            newdomain = domain_t(observation%InitialTime, observation%FinalTime, observation%TimeStep)

            newdomain%tstep = max(newdomain%tstep, timeStep)

            if (10.0_RKIND*(newdomain%tstop - newdomain%tstart)/min(timeStep, newdomain%tstep) >= huge(1_4)) then
               newdomain%tstop = newdomain%tstart + min(timeStep, newdomain%tstep)*huge(1_4)/10.0_RKIND
            end if

            if (newDomain%tstart < newDomain%tstep) then
               newDomain%tstart = 0.0_RKIND_tiempo
            end if

            if (newDomain%tstep > (newdomain%tstop - newdomain%tstart)) then
               newDomain%tstop = newDomain%tstart + newDomain%tstep
            end if

         elseif (observation%FreqDomain) then
            !Just linear progression for now. Need to bring logartihmic info to here
            nFreq = int((observation%FinalFreq - observation%InitialFreq) / observation%FreqStep, kind=SINGLE)
            newdomain = domain_t(observation%InitialFreq, observation%FinalFreq, nFreq, logarithmicspacing=.false.)

            newDomain%fstep = min(newDomain%fstep, 2.0_RKIND/dt)
            if ((newDomain%fstep > newDomain%fstop - newDomain%fstart) .or. (newDomain%fstep == 0)) then
               newDomain%fstep = newDomain%fstop - newDomain%fstart
               newDomain%fstop = newDomain%fstart + observation%fstep
            end if

            newDomain%fnum =  int((newDomain%fstop - newDomain%fstart) / newDomain%fstep, kind=SINGLE)

         else
            call stoponerror('No domain present')
         end if
         return
      end function preprocess_domain

   end subroutine init_observations

   subroutine update_outputs(outputs, step, Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh)
      type(solver_output_t), dimension(:), intent(inout) :: outputs
      real(kind=RKIND_tiempo) :: step
      integer(kind=SINGLE) :: i, id


      REAL(KIND=RKIND), intent(in), target     :: &
        Ex(sgg%alloc(iEx)%XI:sgg%alloc(iEx)%XE, sgg%alloc(iEx)%YI:sgg%alloc(iEx)%YE, sgg%alloc(iEx)%ZI:sgg%alloc(iEx)%ZE), &
        Ey(sgg%alloc(iEy)%XI:sgg%alloc(iEy)%XE, sgg%alloc(iEy)%YI:sgg%alloc(iEy)%YE, sgg%alloc(iEy)%ZI:sgg%alloc(iEy)%ZE), &
        Ez(sgg%alloc(iEz)%XI:sgg%alloc(iEz)%XE, sgg%alloc(iEz)%YI:sgg%alloc(iEz)%YE, sgg%alloc(iEz)%ZI:sgg%alloc(iEz)%ZE), &
        Hx(sgg%alloc(iHx)%XI:sgg%alloc(iHx)%XE, sgg%alloc(iHx)%YI:sgg%alloc(iHx)%YE, sgg%alloc(iHx)%ZI:sgg%alloc(iHx)%ZE), &
        Hy(sgg%alloc(iHy)%XI:sgg%alloc(iHy)%XE, sgg%alloc(iHy)%YI:sgg%alloc(iHy)%YE, sgg%alloc(iHy)%ZI:sgg%alloc(iHy)%ZE), &
        Hz(sgg%alloc(iHz)%XI:sgg%alloc(iHz)%XE, sgg%alloc(iHz)%YI:sgg%alloc(iHz)%YE, sgg%alloc(iHz)%ZI:sgg%alloc(iHz)%ZE)
      !--->
      REAL(KIND=RKIND), dimension(:), intent(in)   :: dxh(sgg%ALLOC(iEx)%XI:sgg%ALLOC(iEx)%XE), &
                                                      dyh(sgg%ALLOC(iEy)%YI:sgg%ALLOC(iEy)%YE), &
                                                      dzh(sgg%ALLOC(iEz)%ZI:sgg%ALLOC(iEz)%ZE), &
                                                      dxe(sgg%alloc(iHx)%XI:sgg%alloc(iHx)%XE), &
                                                      dye(sgg%alloc(iHy)%YI:sgg%alloc(iHy)%YE), &
                                                      dze(sgg%alloc(iHz)%ZI:sgg%alloc(iHz)%ZE)


      do i = 1, size(outputs)
         id = outputs(i)%outputID
         select case(id)
         case(POINT_PROBE_ID)
            field => get_field_component(outputs(i)%pointProbe%fieldComponent) !Cada componente requiere de valores deiferentes pero estos valores no se como conseguirlos
            update_solver_output(outputs(i)%pointProbe, step, field)
         case default
            call stoponerror('Output update not implemented')
         end select
      end do

      contains
      function get_field_component(fieldId) result(field)
      integer(kind=SINGLE), intent(in) :: fieldId
      select case(fieldId)
      case(iEx); field => Ex 
      case(iEy); field => Ey 
      case(iEz); field => Ez 
      case(iHx); field => Hx 
      case(iHy); field => Hy 
      case(iHz); field => Hz 
      end select
      end function get_field_component

   end subroutine update_outputs

   subroutine init_point_probe_output(this, iCoord, jCoord, kCoord, field, domain, outputTypeExtension, mpidir)
      type(point_probe_output_t), intent(out) :: this
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      type(domain_t), intent(in) :: domain

      character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
      integer(kind=SINGLE) :: i

      this%xCoord = iCoord
      this%yCoord = jCoord
      this%zCoord = kCoord

      this%domain = domain
      this%path = get_output_path()

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         this%nFreq = this%domain%fnum
         allocate (this%frequencySlice(this%domain%fnum))
         allocate (this%valueForFreq(this%domain%fnum))
         do i = 1, this%nFreq
            call init_frequency_slice(this%frequencySlice, this%domain)
         end do
         this%valueForFreq = (0.0_RKIND, 0.0_RKIND)
      end if

   contains
      function get_output_path() result(outputPath)
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_probe_bounds_extension()
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//trim(adjustl(probeBoundsExtension))
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
   end subroutine init_point_probe_output

   subroutine
      

   subroutine update_point_probe_output(this, step, field)
      type(point_probe_output_t), intent(inout) :: this
      real(kind=RKIND), pointer, dimension(:, :, :) :: field
      real(kind=RKIND_tiempo) :: step
      integer(kind=SINGLE) :: iter

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         this%serializedTimeSize = this%serializedTimeSize + 1
         this%timeStep(this%serializedTimeSize) = step
         this%valueForTime(this%serializedTimeSize) = field(this%xCoord, this%yCoord, this%zCoord)
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         do iter = 1, this%nFreq
            this%valueForFreq(iter) = &
      this%valueForFreq(iter) + field(this%xCoord, this%yCoord, this%zCoord)*get_auxExp(this%frequencySlice(iter), this%fieldComponent)
         end do
      end if
   end subroutine update_point_probe_output

   subroutine flush_point_probe_output(this)
      type(point_probe_output_t), intent(inout) :: this

      integer(kind=SINGLE) :: timeUnitFile, frequencyUnitFile, status
      character(len=BUFSIZE) :: timeFileName, frequencyFileName
      integer(kind=SINGLE) :: i

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         timeFileName = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
         timeUnitFile = FILE_UNIT + 1

         status = open_file(timeUnitFile, timeFileName)
         if (status /= 0) call stoponerror('Failed to open timeDomainFile. ')

         do i = 1, this%serializedTimeSize
            write (timeUnitFile, '(F12.4, 2X, F12.4)') this%timeStep(i), this%valueForTime(i)
         end do

         status = close_file(timeUnitFile)
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         frequencyFileName = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
         frequencyUnitFile = FILE_UNIT + 2

         OPEN (UNIT=frequencyUnitFile, FILE=frequencyFileName, STATUS='REPLACE', ACTION='WRITE', iostat=status)
         if (status /= 0) call stoponerror('Failed to open frequencyDomainFile. ')

         do i = 1, this%nFreq
            write (frequencyUnitFile, '(F12.4, 2X, F12.4)') this%frequencySlice(i), this%valueForFreq(i)
         end do

         status = close_file(frequencyUnitFile)
      end if
   end subroutine flush_point_probe_output

   subroutine delete_point_probe_output()
      !TODO
   end subroutine delete_point_probe_output

end module output
