module pointProbeOutput_m
   use FDETYPES_m
   use utils_m
   use allocationUtils_m, only: alloc_and_init
   use outputTypes_m
   use domain_m
   use outputUtils_m
   use ilumina_m, only: Incid
   use Report_m, only: log_msg, LOG_INFO, LOG_WARN, LOG_FATAL

   implicit none

   private

   public :: init_point_probe_output
   public :: update_point_probe_output
   public :: flush_point_probe_output
   public :: close_point_probe_output

contains
   subroutine init_point_probe_output(this, coordinates, field, domain, outputTypeExtension, mpidir, timeInterval, &
                                      hasIncident, append_to_file)
      type(point_probe_output_t), intent(out) :: this
      type(cell_coordinate_t) :: coordinates
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=*), intent(in) :: outputTypeExtension
      type(domain_t), intent(in) :: domain

      real(kind=RKIND_TIME), intent(in) :: timeInterval
      logical, intent(in), optional :: hasIncident
      logical, intent(in), optional :: append_to_file

      integer(kind=SINGLE) :: i
      integer :: artifact_kinds(2)
      character(len=BUFSIZE) :: artifact_paths(2)

      this%mainCoords = coordinates

      this%component = field

      this%domain = domain
      this%path = get_output_path()
      if (present(hasIncident)) this%hasIncident = hasIncident

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         call alloc_and_init(this%timeStep, time_buffer_capacity(this%domain), 0.0_RKIND_TIME)
         call alloc_and_init(this%valueForTime, time_buffer_capacity(this%domain), 0.0_RKIND)
         if (this%hasIncident) call alloc_and_init(this%incidentForTime, time_buffer_capacity(this%domain), 0.0_RKIND)
      end if
      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         this%nFreq = this%domain%fnum
         this%quadratureDt = timeInterval
         allocate (this%frequencySlice(this%domain%fnum))
         call alloc_and_init(this%valueForFreq, this%domain%fnum, (0.0_CKIND, 0.0_CKIND))
         if (this%hasIncident) call alloc_and_init(this%incidentForFreq, this%domain%fnum, (0.0_CKIND, 0.0_CKIND))
         call init_frequency_slice(this%frequencySlice, this%domain)
         this%valueForFreq = (0.0_RKIND, 0.0_RKIND)

         allocate (this%auxExp_E(this%nFreq))
         allocate (this%auxExp_H(this%nFreq))
         do i = 1, this%nFreq
            this%auxExp_E(i) = MCPI2*this%frequencySlice(i)
            this%auxExp_H(i) = this%auxExp_E(i)
         end do
         if (this%domain%transferFlag) call initialise_normalization_spectrum(this)
      end if

      if (this%domain%domainType == BOTH_DOMAIN) then
         allocate (this%artifacts(2))
         artifact_paths(1) = trim(this%path)//'_'//TIMEEXTENSION//DATFILEEXTENSION
         artifact_paths(2) = trim(this%path)//'_'//FREQUENCYEXTENSION//DATFILEEXTENSION
         artifact_kinds(:2) = OUTPUT_ARTIFACT_TEXT
         call declare_probe_artifacts(this%artifacts, artifact_paths(:2), artifact_kinds(:2))
         this%filePathTime = this%artifacts(1)%relative_path
         this%filePathFreq = this%artifacts(2)%relative_path
         call open_data_file(this%filePathTime, this%path, TIMEEXTENSION, DATFILEEXTENSION, &
                             time_header(), this%fileUnitTime, append=append_to_file)
         call open_data_file(this%filePathFreq, this%path, FREQUENCYEXTENSION, DATFILEEXTENSION, &
                             frequency_header(this), this%fileUnitFreq)
         call check_output_file(this%fileUnitTime, this%filePathTime)
         call check_output_file(this%fileUnitFreq, this%filePathFreq)
      else if (this%domain%domainType == TIME_DOMAIN) then
         allocate (this%artifacts(1))
         artifact_paths(1) = trim(this%path)//'_'//TIMEEXTENSION//DATFILEEXTENSION
         artifact_kinds(1) = OUTPUT_ARTIFACT_TEXT
         call declare_probe_artifacts(this%artifacts, artifact_paths(:1), artifact_kinds(:1))
         this%filePathTime = this%artifacts(1)%relative_path
         call open_data_file(this%filePathTime, this%path, TIMEEXTENSION, DATFILEEXTENSION, &
                             time_header(), this%fileUnitTime, append=append_to_file)
         call check_output_file(this%fileUnitTime, this%filePathTime)
      else if (this%domain%domainType == FREQUENCY_DOMAIN) then
         allocate (this%artifacts(1))
         artifact_paths(1) = trim(this%path)//'_'//FREQUENCYEXTENSION//DATFILEEXTENSION
         artifact_kinds(1) = OUTPUT_ARTIFACT_TEXT
         call declare_probe_artifacts(this%artifacts, artifact_paths(:1), artifact_kinds(:1))
         this%filePathFreq = this%artifacts(1)%relative_path
         call open_data_file(this%filePathFreq, this%path, FREQUENCYEXTENSION, DATFILEEXTENSION, &
                             frequency_header(this), this%fileUnitFreq)
         call check_output_file(this%fileUnitFreq, this%filePathFreq)
      end if

   contains
      subroutine check_output_file(unit, filePath)
         integer, intent(in) :: unit
         character(len=*), intent(in) :: filePath

         if (unit == 0) call log_msg(LOG_FATAL, 'OUTPUT', 'Unable to open probe output file: '//trim(filePath))
      end subroutine check_output_file

      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_coordinates_extension(this%mainCoords, mpidir)
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

      function time_header() result(header)
         character(len=16) :: header

         if (this%hasIncident) then
            header = 't field incident'
         else
            header = 't field'
         end if
      end function time_header

   end subroutine init_point_probe_output

   function frequency_header(this) result(header)
      type(point_probe_output_t), intent(in) :: this
      character(len=80) :: header
      if (this%hasIncident) then
         header = 'frequency magnitude phase incident_magnitude incident_phase'
      else
         header = 'frequency magnitude phase'
      end if
   end function frequency_header

   subroutine initialise_normalization_spectrum(this)
      type(point_probe_output_t), intent(inout) :: this
      integer :: unit, ioStatus, i
      real(kind=RKIND_TIME) :: t0, t1, t, dt
      real(kind=RKIND) :: scalarValue

      open(newunit=unit, file=trim(this%domain%normalizationFile), status='old', action='read', iostat=ioStatus)
      if (ioStatus /= 0) call log_msg(LOG_FATAL, 'OUTPUT', 'Unable to read point-probe normalization file')
      read(unit, *, iostat=ioStatus) t0, scalarValue
      read(unit, *, iostat=ioStatus) t1, scalarValue
      if (ioStatus /= 0) call log_msg(LOG_FATAL, 'OUTPUT', 'Normalization file needs at least two samples')
      dt = abs(t1 - t0)
      if (dt <= tiny(1.0_RKIND_TIME)) call log_msg(LOG_FATAL, 'OUTPUT', 'Normalization file has a zero sampling interval')
      rewind(unit)
      call alloc_and_init(this%normalizationForFreq, this%nFreq, (0.0_CKIND, 0.0_CKIND))
      do
         read(unit, *, iostat=ioStatus) t, scalarValue
         if (ioStatus < 0) exit
         if (ioStatus /= 0) call log_msg(LOG_FATAL, 'OUTPUT', 'Invalid normalization sample')
         do i = 1, this%nFreq
            this%normalizationForFreq(i) = this%normalizationForFreq(i) + dt*scalarValue*exp(MCPI2*this%frequencySlice(i)*t)
         end do
      end do
      close(unit)
      if (any(abs(this%normalizationForFreq) <= tiny(1.0_RKIND))) call log_msg(LOG_FATAL, 'OUTPUT', 'Zero normalization spectrum')
   end subroutine initialise_normalization_spectrum

   subroutine update_point_probe_output(this, step, field, sgg, saveTimeSample)
      type(point_probe_output_t), intent(inout) :: this
      real(kind=RKIND), pointer, dimension(:, :, :), intent(in) :: field
      real(kind=RKIND_TIME), intent(in) :: step
      type(SGGFDTDINFO_t), intent(in), optional :: sgg
      logical, intent(in), optional :: saveTimeSample

      integer(kind=SINGLE) :: iter
      logical :: recordTimeSample, still_planewave_time
      real(kind=RKIND) :: incidentValue

      recordTimeSample = .true.
      if (present(saveTimeSample)) recordTimeSample = saveTimeSample

      if (recordTimeSample .and. any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         this%nTime = this%nTime + 1
         this%timeStep(this%nTime) = step
         this%valueForTime(this%nTime) = field(this%mainCoords%x, this%mainCoords%y, this%mainCoords%z)
         if (this%hasIncident .and. present(sgg)) then
            still_planewave_time = .false.
            this%incidentForTime(this%nTime) = Incid(sgg, 1, this%component, real(step + sgg%dt, RKIND), &
                                                     this%mainCoords%x, this%mainCoords%y, this%mainCoords%z, &
                                                     still_planewave_time, .true.)
         end if
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         incidentValue = 0.0_RKIND
         if (this%hasIncident .and. present(sgg)) then
            still_planewave_time = .false.
            incidentValue = Incid(sgg, 1, this%component, real(step + sgg%dt, RKIND), &
                                  this%mainCoords%x, this%mainCoords%y, this%mainCoords%z, still_planewave_time, .true.)
         end if
         select case (this%component)
         case (IEX, IEY, IEZ)
            do iter = 1, this%nFreq
               this%valueForFreq(iter) = &
                  this%valueForFreq(iter) + field(this%mainCoords%x, this%mainCoords%y, this%mainCoords%z)* &
                  this%quadratureDt*exp(this%auxExp_E(iter)*step)
               if (this%hasIncident .and. present(sgg)) this%incidentForFreq(iter) = this%incidentForFreq(iter) + &
                  incidentValue*this%quadratureDt*exp(this%auxExp_E(iter)*step)
            end do
         case (IHX, IHY, IHZ)
            do iter = 1, this%nFreq
               this%valueForFreq(iter) = &
                  this%valueForFreq(iter) + field(this%mainCoords%x, this%mainCoords%y, this%mainCoords%z)* &
                  this%quadratureDt*exp(this%auxExp_H(iter)*(step + 0.5_RKIND_TIME*this%quadratureDt))
               if (this%hasIncident .and. present(sgg)) this%incidentForFreq(iter) = this%incidentForFreq(iter) + &
                  incidentValue*this%quadratureDt*exp(this%auxExp_H(iter)*(step + 0.5_RKIND_TIME*this%quadratureDt))
            end do
         end select

      end if
   end subroutine update_point_probe_output

   subroutine flush_point_probe_output(this)
      type(point_probe_output_t), intent(inout) :: this
      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         call flush_time_domain(this)
         call clear_time_data()
      end if
      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         call flush_frequency_domain(this)
      end if
   contains

      subroutine flush_time_domain(this)
         type(point_probe_output_t), intent(in) :: this
         integer :: i

         if (this%nTime <= 0) then
#ifdef CompileWithDebug
            call log_msg(LOG_INFO, 'OUTPUT', 'No data to write.')
#endif
             return
         end if
         if (this%fileUnitTime == 0) call log_msg(LOG_FATAL, 'OUTPUT', 'Point probe output is not open: '//trim(this%filePathTime))

         do i = 1, this%nTime
            if (this%hasIncident) then
               write (this%fileUnitTime, FMT) this%timeStep(i), this%valueForTime(i), this%incidentForTime(i)
            else
               write (this%fileUnitTime, FMT) this%timeStep(i), this%valueForTime(i)
            end if
         end do
         flush (this%fileUnitTime)
      end subroutine flush_time_domain

      subroutine flush_frequency_domain(this)
         type(point_probe_output_t), intent(in) :: this
         integer :: i
         complex(kind=CKIND) :: spectrum

         if (.not. allocated(this%frequencySlice) .or. .not. allocated(this%valueForFreq)) then
            call log_msg(LOG_WARN, 'OUTPUT', 'Frequency arrays not allocated; skipping write.')
            return
         end if

         if (this%nFreq <= 0) then
#ifdef CompileWithDebug
            call log_msg(LOG_INFO, 'OUTPUT', 'No data to write.')
#endif
             return
         end if
         if (this%fileUnitFreq == 0) call log_msg(LOG_FATAL, 'OUTPUT', 'Point probe output is not open: '//trim(this%filePathFreq))

         rewind (this%fileUnitFreq)
         write (this%fileUnitFreq, '(A)') frequency_header(this)

         do i = 1, this%nFreq
            spectrum = this%valueForFreq(i)
            if (this%domain%transferFlag) spectrum = spectrum/this%normalizationForFreq(i)
            if (this%hasIncident) then
               write (this%fileUnitFreq, FMT) this%frequencySlice(i), abs(spectrum), atan2(aimag(spectrum), real(spectrum)), &
                                              abs(this%incidentForFreq(i)), &
                                              atan2(aimag(this%incidentForFreq(i)), real(this%incidentForFreq(i)))
            else
               write (this%fileUnitFreq, FMT) this%frequencySlice(i), abs(spectrum), atan2(aimag(spectrum), real(spectrum))
            end if
         end do
         endfile (this%fileUnitFreq)
         flush (this%fileUnitFreq)
      end subroutine flush_frequency_domain

      subroutine clear_time_data()
         this%timeStep = 0.0_RKIND_TIME
         this%valueForTime = 0.0_RKIND
         if (this%hasIncident) this%incidentForTime = 0.0_RKIND

         this%nTime = 0
      end subroutine clear_time_data

   end subroutine flush_point_probe_output

   subroutine close_point_probe_output(this)
      type(point_probe_output_t), intent(inout) :: this

      if (this%fileUnitTime /= 0) then
         close (this%fileUnitTime)
         this%fileUnitTime = 0
      end if
      if (this%fileUnitFreq /= 0) then
         close (this%fileUnitFreq)
         this%fileUnitFreq = 0
      end if
   end subroutine close_point_probe_output
end module
