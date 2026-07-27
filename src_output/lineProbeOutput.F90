module lineProbeOutput_m
   use FDETYPES_m, only: RKIND, RKIND_tiempo, BUFSIZE, direction_t, iEx, iEy, iEz
   use outputTypes_m, only: field_data_t, line_probe_output_t, domain_t, TIME_DOMAIN, &
                             OUTPUT_ARTIFACT_TEXT, OUTPUT_ARTIFACT_BINARY, BINARY_ENDIAN_LITTLE, &
                             BINARY_NUMERIC_REAL32, BINARY_COMPLEX_UNSPECIFIED, binaryExtension, &
                             datFileExtension, timeExtension, declare_probe_artifacts
   use allocationUtils_m, only: alloc_and_init
   use directoryUtils_m, only: create_file_with_path
   use outputBinary_m, only: append_binary_real32, BINARY_WRITER_SUCCESS
   use outputTransport_m, only: output_transport_t, reduce_scalar_batch, OUTPUT_TRANSPORT_SUCCESS
   implicit none
   private

   public :: calculate_line_integral, init_line_probe_output, update_line_probe_output, &
             flush_line_probe_output, complete_line_probe_sample, reduce_line_probe_sample

contains

   function calculate_line_integral(segments, electric_field) result(value)
      type(direction_t), intent(in) :: segments(:)
      type(field_data_t), intent(in) :: electric_field
      real(kind=RKIND) :: value
      integer :: segment_index, orientation

      value = 0.0_RKIND
      do segment_index = 1, size(segments)
         orientation = segments(segment_index)%orientation
         select case (abs(orientation))
         case (iEx)
            value = value + electric_field%x(segments(segment_index)%x, segments(segment_index)%y, &
                                              segments(segment_index)%z) * sign(1, orientation) * &
                            electric_field%deltaX(segments(segment_index)%x)
         case (iEy)
            value = value + electric_field%y(segments(segment_index)%x, segments(segment_index)%y, &
                                              segments(segment_index)%z) * sign(1, orientation) * &
                            electric_field%deltaY(segments(segment_index)%y)
         case (iEz)
            value = value + electric_field%z(segments(segment_index)%x, segments(segment_index)%y, &
                                              segments(segment_index)%z) * sign(1, orientation) * &
                            electric_field%deltaZ(segments(segment_index)%z)
         end select
      end do
   end function calculate_line_integral

   subroutine init_line_probe_output(this, segments, domain, output_path)
      type(line_probe_output_t), intent(out) :: this
      type(direction_t), intent(in) :: segments(:)
      type(domain_t), intent(in) :: domain
      character(len=*), intent(in) :: output_path
      character(len=BUFSIZE) :: artifact_paths(2)
      integer :: artifact_kinds(2), ios

      this%domain = domain
      this%path = output_path
      allocate(this%segments(size(segments)))
      this%segments = segments
      call alloc_and_init(this%timeStep, BUFSIZE, 0.0_RKIND_tiempo)
      call alloc_and_init(this%valueForTime, BUFSIZE, 0.0_RKIND)

      artifact_paths(1) = trim(this%path)//'_'//timeExtension//datFileExtension
      artifact_paths(2) = trim(this%path)//'_'//timeExtension//binaryExtension
      artifact_kinds = [OUTPUT_ARTIFACT_TEXT, OUTPUT_ARTIFACT_BINARY]
      call declare_probe_artifacts(this%artifacts, artifact_paths, artifact_kinds)
      this%artifacts(2)%byte_order = BINARY_ENDIAN_LITTLE
      this%artifacts(2)%numeric_representation = BINARY_NUMERIC_REAL32
      this%artifacts(2)%complex_representation = BINARY_COMPLEX_UNSPECIFIED
      this%artifacts(2)%record_bytes = 8
      this%artifacts(2)%component_order = 'time,line_integral'
      call create_file_with_path(this%artifacts(1)%relative_path, ios)
      call create_file_with_path(this%artifacts(2)%relative_path, ios)
   end subroutine init_line_probe_output

   subroutine update_line_probe_output(this, step, electric_field)
      type(line_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(field_data_t), intent(in) :: electric_field

      if (this%domain%domainType /= TIME_DOMAIN .or. size(this%segments) == 0) return
      this%nTime = this%nTime + 1
      this%timeStep(this%nTime) = step
      this%valueForTime(this%nTime) = calculate_line_integral(this%segments, electric_field)
   end subroutine update_line_probe_output

   subroutine complete_line_probe_sample(this)
      type(line_probe_output_t), intent(inout) :: this

      if (this%nTime == 0) return
      this%nTimesFlushed = this%nTimesFlushed + this%nTime
      this%nTime = 0
      this%timeStep = 0.0_RKIND_tiempo
      this%valueForTime = 0.0_RKIND
   end subroutine complete_line_probe_sample

   subroutine reduce_line_probe_sample(transport, local_value, canonical_value, status)
      type(output_transport_t), intent(in) :: transport
      real(kind=RKIND), intent(in) :: local_value
      real(kind=RKIND), intent(out) :: canonical_value
      integer, intent(out) :: status
      real(kind=RKIND), allocatable :: reduced(:)

      call reduce_scalar_batch(transport, [local_value], reduced, status)
      if (status == OUTPUT_TRANSPORT_SUCCESS) canonical_value = reduced(1)
   end subroutine reduce_line_probe_sample

   subroutine flush_line_probe_output(this)
      type(line_probe_output_t), intent(inout) :: this
      integer :: index, ios, unit
      real(kind=RKIND), allocatable :: records(:)
      real(kind=RKIND_tiempo) :: time_value

      if (this%nTime == 0) return
      open(newunit=unit, file=this%artifacts(1)%relative_path, status='old', action='write', position='append', iostat=ios)
      if (ios /= 0) return
      do index = 1, this%nTime
         write(unit, '(ES24.16E3,1X,ES24.16E3)', iostat=ios) this%timeStep(index), this%valueForTime(index)
         if (ios /= 0) exit
      end do
      close(unit)
      if (ios /= 0) return

      allocate(records(2 * this%nTime))
      do index = 1, this%nTime
         records(2 * index - 1) = real(this%timeStep(index), RKIND)
         records(2 * index) = this%valueForTime(index)
      end do
      call append_binary_real32(this%artifacts(2)%relative_path, this%artifacts(2), real(records, kind=4), ios)
      if (ios == BINARY_WRITER_SUCCESS) call complete_line_probe_sample(this)
   end subroutine flush_line_probe_output

end module lineProbeOutput_m
