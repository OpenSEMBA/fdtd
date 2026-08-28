module pointProbeOutput_m
   use FDETYPES_m
   use utils_m
   use allocationUtils_m, only: alloc_and_init
   use outputTypes_m
   use domain_m
   use outputUtils_m
   use ilumina_m, only: Incid

   implicit none

   private

   public :: init_point_probe_output
   public :: update_point_probe_output
   public :: flush_point_probe_output

contains
   subroutine init_point_probe_output(this, coordinates, field, domain, outputTypeExtension, mpidir, timeInterval, hasIncident)
      type(point_probe_output_t), intent(out) :: this
      type(cell_coordinate_t) :: coordinates
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=*), intent(in) :: outputTypeExtension
      type(domain_t), intent(in) :: domain

      real(kind=RKIND_tiempo), intent(in) :: timeInterval
      logical, intent(in), optional :: hasIncident

      integer(kind=SINGLE) :: i
      integer :: artifact_kinds(2)
      character(len=BUFSIZE) :: artifact_paths(2)

      this%mainCoords = coordinates

      this%component = field

      this%domain = domain
      this%path = get_output_path()
      if (present(hasIncident)) this%hasIncident = hasIncident

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         call alloc_and_init(this%timeStep, OUTPUT_TIME_BUFFER_SIZE, 0.0_RKIND_tiempo)
         call alloc_and_init(this%valueForTime, OUTPUT_TIME_BUFFER_SIZE, 0.0_RKIND)
         if (this%hasIncident) call alloc_and_init(this%incidentForTime, OUTPUT_TIME_BUFFER_SIZE, 0.0_RKIND)
      end if
      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         this%nFreq = this%domain%fnum
         this%quadratureDt = timeInterval
         allocate (this%frequencySlice(this%domain%fnum))
         call alloc_and_init(this%valueForFreq, this%domain%fnum, (0.0_CKIND, 0.0_CKIND))
         call init_frequency_slice(this%frequencySlice, this%domain)
         this%valueForFreq = (0.0_RKIND, 0.0_RKIND)

         allocate (this%auxExp_E(this%nFreq))
         allocate (this%auxExp_H(this%nFreq))
         do i = 1, this%nFreq
            this%auxExp_E(i) = mcpi2*this%frequencySlice(i)
            this%auxExp_H(i) = this%auxExp_E(i)
         end do
      end if

      if (this%domain%domainType == BOTH_DOMAIN) then
         allocate (this%artifacts(2))
         artifact_paths(1) = trim(this%path)//'_'//timeExtension//datFileExtension
         artifact_paths(2) = trim(this%path)//'_'//frequencyExtension//datFileExtension
         artifact_kinds(:2) = OUTPUT_ARTIFACT_TEXT
         call declare_probe_artifacts(this%artifacts, artifact_paths(:2), artifact_kinds(:2))
         this%filePathTime = this%artifacts(1)%relative_path
         this%filePathFreq = this%artifacts(2)%relative_path
         call create_data_file(this%filePathTime, this%path, timeExtension, datFileExtension, time_header())
         call create_data_file(this%filePathFreq, this%path, frequencyExtension, datFileExtension, &
                               'frequency real imaginary')
      else if (this%domain%domainType == TIME_DOMAIN) then
         allocate (this%artifacts(1))
         artifact_paths(1) = trim(this%path)//'_'//timeExtension//datFileExtension
         artifact_kinds(1) = OUTPUT_ARTIFACT_TEXT
         call declare_probe_artifacts(this%artifacts, artifact_paths(:1), artifact_kinds(:1))
         this%filePathTime = this%artifacts(1)%relative_path
         call create_data_file(this%filePathTime, this%path, timeExtension, datFileExtension, time_header())
      else if (this%domain%domainType == FREQUENCY_DOMAIN) then
         allocate (this%artifacts(1))
         artifact_paths(1) = trim(this%path)//'_'//frequencyExtension//datFileExtension
         artifact_kinds(1) = OUTPUT_ARTIFACT_TEXT
         call declare_probe_artifacts(this%artifacts, artifact_paths(:1), artifact_kinds(:1))
         this%filePathFreq = this%artifacts(1)%relative_path
         call create_data_file(this%filePathFreq, this%path, frequencyExtension, datFileExtension, &
                               'frequency real imaginary')
      end if

   contains
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

   subroutine update_point_probe_output(this, step, field, sgg)
      type(point_probe_output_t), intent(inout) :: this
      real(kind=RKIND), pointer, dimension(:, :, :), intent(in) :: field
      real(kind=RKIND_tiempo), intent(in) :: step
      type(SGGFDTDINFO_t), intent(in), optional :: sgg

      integer(kind=SINGLE) :: iter
      logical :: still_planewave_time

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
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
         select case (this%component)
         case (iEx, iEy, iEz)
            do iter = 1, this%nFreq
               this%valueForFreq(iter) = &
                  this%valueForFreq(iter) + field(this%mainCoords%x, this%mainCoords%y, this%mainCoords%z)* &
                  this%quadratureDt*exp(this%auxExp_E(iter)*step)
            end do
         case (iHx, iHy, iHz)
            do iter = 1, this%nFreq
               this%valueForFreq(iter) = &
                  this%valueForFreq(iter) + field(this%mainCoords%x, this%mainCoords%y, this%mainCoords%z)* &
                  this%quadratureDt*exp(this%auxExp_H(iter)*(step + 0.5_RKIND_tiempo*this%quadratureDt))
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
         integer :: unit

         if (this%nTime <= 0) then
            print *, "No data to write."
            return
         end if
         open (newunit=unit, file=this%filePathTime, status="old", action="write", position="append")

         do i = 1, this%nTime
            if (this%hasIncident) then
               write (unit, fmt) this%timeStep(i), this%valueForTime(i), this%incidentForTime(i)
            else
               write (unit, fmt) this%timeStep(i), this%valueForTime(i)
            end if
         end do

         close (unit)
      end subroutine flush_time_domain

      subroutine flush_frequency_domain(this)
         type(point_probe_output_t), intent(in) :: this
         integer :: i
         integer :: unit

         if (.not. allocated(this%frequencySlice) .or. .not. allocated(this%valueForFreq)) then
            print *, "Error: arrays not allocated."
            return
         end if

         if (this%nFreq <= 0) then
            print *, "No data to write."
            return
         end if
         open (newunit=unit, file=this%filePathFreq, status="replace", action="write")
         write (unit, '(A)') 'frequency real imaginary'

         do i = 1, this%nFreq
            write (unit, fmt) this%frequencySlice(i), real(this%valueForFreq(i)), aimag(this%valueForFreq(i))
         end do

         close (unit)
      end subroutine flush_frequency_domain

      subroutine clear_time_data()
         this%timeStep = 0.0_RKIND_tiempo
         this%valueForTime = 0.0_RKIND
         if (this%hasIncident) this%incidentForTime = 0.0_RKIND

         this%nTime = 0
      end subroutine clear_time_data

   end subroutine flush_point_probe_output
end module
