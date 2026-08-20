module pointProbeOutput_m
   use FDETYPES_m
   use utils_m
   use allocationUtils_m, only: alloc_and_init
   use outputTypes_m
   use domain_m
     use outputUtils_m
     use allocationUtils_m, only: alloc_and_init
   use outputBinary_m, only: append_binary_real64, write_binary_complex_record64, BINARY_WRITER_SUCCESS
   use directoryUtils_m, only: create_file_with_path, get_last_component, join_path
   use, intrinsic :: iso_fortran_env, only: real64
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
        integer :: artifact_kinds(4)
        character(len=BUFSIZE) :: artifact_paths(4)

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
           artifact_paths(1) = trim(this%path)//'_'//timeExtension//datFileExtension
           artifact_paths(2) = trim(this%path)//'_'//timeExtension//binaryExtension
           artifact_paths(3) = trim(this%path)//'_'//frequencyExtension//datFileExtension
           artifact_paths(4) = trim(this%path)//'_'//frequencyExtension//binaryExtension
           artifact_kinds = [OUTPUT_ARTIFACT_TEXT, OUTPUT_ARTIFACT_BINARY, &
                             OUTPUT_ARTIFACT_TEXT, OUTPUT_ARTIFACT_BINARY]
           call declare_probe_artifacts(this%artifacts, artifact_paths, artifact_kinds)
           this%filePathTime = this%artifacts(1)%relative_path
           this%filePathFreq = this%artifacts(3)%relative_path
            call create_data_file(this%filePathTime, this%path, timeExtension, datFileExtension, time_header())
            call create_data_file(this%filePathFreq, this%path, frequencyExtension, datFileExtension, &
                                  'frequency real imaginary')
           call configure_time_binary_artifact(this%artifacts(2))
           call configure_binary_artifact(this%artifacts(4), 24_8, BINARY_COMPONENTS_SCALAR_FREQUENCY, &
                                          BINARY_COMPLEX_REAL_IMAG)
        else if (this%domain%domainType == TIME_DOMAIN) then
           artifact_paths(1) = trim(this%path)//'_'//timeExtension//datFileExtension
           artifact_paths(2) = trim(this%path)//'_'//timeExtension//binaryExtension
           artifact_kinds(:2) = [OUTPUT_ARTIFACT_TEXT, OUTPUT_ARTIFACT_BINARY]
           call declare_probe_artifacts(this%artifacts, artifact_paths(:2), artifact_kinds(:2))
           this%filePathTime = this%artifacts(1)%relative_path
            call create_data_file(this%filePathTime, this%path, timeExtension, datFileExtension, time_header())
           call configure_time_binary_artifact(this%artifacts(2))
        else if (this%domain%domainType == FREQUENCY_DOMAIN) then
           artifact_paths(1) = trim(this%path)//'_'//frequencyExtension//datFileExtension
           artifact_paths(2) = trim(this%path)//'_'//frequencyExtension//binaryExtension
           artifact_kinds(:2) = [OUTPUT_ARTIFACT_TEXT, OUTPUT_ARTIFACT_BINARY]
           call declare_probe_artifacts(this%artifacts, artifact_paths(:2), artifact_kinds(:2))
           this%filePathFreq = this%artifacts(1)%relative_path
            call create_data_file(this%filePathFreq, this%path, frequencyExtension, datFileExtension, &
                                  'frequency real imaginary')
          call configure_binary_artifact(this%artifacts(2), 24_8, BINARY_COMPONENTS_SCALAR_FREQUENCY, &
                                         BINARY_COMPLEX_REAL_IMAG)
        end if
        if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
           call create_file_with_path(this%artifacts(2)%relative_path, i)
        end if
        if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
           if (this%domain%domainType == BOTH_DOMAIN) then
              call create_file_with_path(this%artifacts(4)%relative_path, i)
           else
              call create_file_with_path(this%artifacts(2)%relative_path, i)
           end if
        end if

    contains
       function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_coordinates_extension(this%mainCoords, mpidir)
         prefixFieldExtension = get_prefix_extension(field, mpidir)
          outputPath = &
             trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
          outputPath = join_path(outputPath, get_last_component(outputPath))
         return
       end function get_output_path

       subroutine configure_binary_artifact(artifact, record_bytes, component_order, complex_representation)
          type(output_artifact_t), intent(inout) :: artifact
          integer(kind=8), intent(in) :: record_bytes
          character(len=*), intent(in) :: component_order
          integer, intent(in) :: complex_representation

          artifact%byte_order = BINARY_ENDIAN_LITTLE
          artifact%numeric_representation = BINARY_NUMERIC_REAL64
          artifact%complex_representation = complex_representation
          artifact%record_bytes = record_bytes
          artifact%component_order = component_order
       end subroutine configure_binary_artifact

       subroutine configure_time_binary_artifact(artifact)
          type(output_artifact_t), intent(inout) :: artifact

          if (this%hasIncident) then
             call configure_binary_artifact(artifact, 24_8, BINARY_COMPONENTS_SCALAR_TIME_INCIDENT, &
                                            BINARY_COMPLEX_UNSPECIFIED)
          else
             call configure_binary_artifact(artifact, 16_8, BINARY_COMPONENTS_SCALAR_TIME, &
                                            BINARY_COMPLEX_UNSPECIFIED)
          end if
       end subroutine configure_time_binary_artifact

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
         select case(this%component)
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
          call write_time_binary(this)
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
          call write_frequency_binary(this)
       end subroutine flush_frequency_domain

       subroutine write_time_binary(this)
          type(point_probe_output_t), intent(in) :: this
          real(real64), allocatable :: records(:)
          integer :: i, status

          if (this%hasIncident) then
             allocate(records(3 * this%nTime))
             do i = 1, this%nTime
                records(3 * i - 2:3 * i) = [real(this%timeStep(i), real64), real(this%valueForTime(i), real64), &
                                             real(this%incidentForTime(i), real64)]
             end do
          else
             allocate(records(2 * this%nTime))
             do i = 1, this%nTime
                records(2 * i - 1:2 * i) = [real(this%timeStep(i), real64), real(this%valueForTime(i), real64)]
             end do
          end if
          call append_binary_real64(this%artifacts(2)%relative_path, this%artifacts(2), records, status)
       end subroutine write_time_binary

       subroutine write_frequency_binary(this)
          type(point_probe_output_t), intent(in) :: this
          real(real64), allocatable :: records(:)
          integer :: i, artifact_index, status

          if (this%domain%domainType == BOTH_DOMAIN) then
             artifact_index = 4
          else
             artifact_index = 2
          end if
          allocate(records(3 * this%nFreq))
          do i = 1, this%nFreq
             records(3 * i - 2) = real(this%frequencySlice(i), real64)
             records(3 * i - 1) = real(this%valueForFreq(i), real64)
             records(3 * i) = real(aimag(this%valueForFreq(i)), real64)
          end do
          call write_binary_complex_record64(this%artifacts(artifact_index)%relative_path, &
                                             this%artifacts(artifact_index), records, status)
       end subroutine write_frequency_binary

      subroutine clear_time_data()
          this%timeStep = 0.0_RKIND_tiempo
          this%valueForTime = 0.0_RKIND
          if (this%hasIncident) this%incidentForTime = 0.0_RKIND

         this%nTime = 0
      end subroutine clear_time_data

   end subroutine flush_point_probe_output
end module
