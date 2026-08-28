module lineProbeOutput_m
   use FDETYPES_m, only: RKIND, RKIND_tiempo, SINGLE, BUFSIZE, direction_t, xyzlimit_t, iEx, iEy, iEz
   use outputTypes_m, only: field_data_t, line_probe_output_t, domain_t, TIME_DOMAIN, OUTPUT_TIME_BUFFER_SIZE, &
                            OUTPUT_ARTIFACT_TEXT, datFileExtension, timeExtension, declare_probe_artifacts
   use allocationUtils_m, only: alloc_and_init
   use directoryUtils_m, only: create_file_with_path
#ifdef CompileWithMPI
   use FDETYPES_m, only: REALSIZE
   use mpi
#endif
   implicit none
   private

   public :: calculate_line_integral, init_line_probe_output, update_line_probe_output, &
               flush_line_probe_output, complete_line_probe_sample, line_segment_is_local

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

   subroutine init_line_probe_output(this, segments, domain, output_path, sweeps, rank, rank_count)
      type(line_probe_output_t), intent(out) :: this
      type(direction_t), intent(in) :: segments(:)
      type(domain_t), intent(in) :: domain
      character(len=*), intent(in) :: output_path
      type(xyzlimit_t), intent(in), optional :: sweeps(:)
      integer(kind=SINGLE), intent(in), optional :: rank, rank_count
      character(len=BUFSIZE) :: artifact_paths(1)
       integer :: artifact_kinds(1), ios, unit, segment_index, local_count

      this%domain = domain
      this%path = output_path
#ifdef CompileWithMPI
      if (present(rank_count)) this%isDistributed = rank_count > 1
      if (this%isDistributed .and. present(rank)) this%isWriter = rank == 0
#endif
      local_count = size(segments)
      if (present(sweeps)) then
         local_count = 0
         do segment_index = 1, size(segments)
            if (line_segment_is_local(segments(segment_index), sweeps, rank, rank_count)) local_count = local_count + 1
         end do
      end if
      allocate(this%segments(local_count))
      if (present(sweeps)) then
         local_count = 0
         do segment_index = 1, size(segments)
            if (.not. line_segment_is_local(segments(segment_index), sweeps, rank, rank_count)) cycle
            local_count = local_count + 1
            this%segments(local_count) = segments(segment_index)
         end do
      else
         this%segments = segments
      end if
      call alloc_and_init(this%timeStep, OUTPUT_TIME_BUFFER_SIZE, 0.0_RKIND_tiempo)
      call alloc_and_init(this%valueForTime, OUTPUT_TIME_BUFFER_SIZE, 0.0_RKIND)

      artifact_paths(1) = trim(this%path)//'_'//timeExtension//datFileExtension
      artifact_kinds = OUTPUT_ARTIFACT_TEXT
      call declare_probe_artifacts(this%artifacts, artifact_paths, artifact_kinds)
       if (this%isWriter) then
          call create_file_with_path(this%artifacts(1)%relative_path, ios)
          if (ios /= 0) return
          open(newunit=unit, file=this%artifacts(1)%relative_path, status='old', action='write', position='append', iostat=ios)
          if (ios /= 0) return
          write(unit, '(A)', iostat=ios) 't lineIntegral'
          close(unit)
          if (ios /= 0) return
       end if
   end subroutine init_line_probe_output

   pure logical function line_segment_is_local(segment, sweeps, rank, rank_count)
      type(direction_t), intent(in) :: segment
      type(xyzlimit_t), intent(in) :: sweeps(:)
      integer, intent(in), optional :: rank, rank_count
      integer(kind=SINGLE) :: component, local_rank, local_rank_count, owned_upper_z

      local_rank = 0
      local_rank_count = 1
      if (present(rank)) local_rank = rank
      if (present(rank_count)) local_rank_count = rank_count
      component = abs(segment%orientation)
      if (component < lbound(sweeps, 1) .or. component > ubound(sweeps, 1)) then
         line_segment_is_local = .false.
         return
      end if
      owned_upper_z = sweeps(component)%ZE
      if (local_rank < local_rank_count - 1 .and. any(component == [iEx, iEy])) owned_upper_z = owned_upper_z - 1
      line_segment_is_local = segment%x >= sweeps(component)%XI .and. segment%x <= sweeps(component)%XE .and. &
                              segment%y >= sweeps(component)%YI .and. segment%y <= sweeps(component)%YE .and. &
                              segment%z >= sweeps(component)%ZI .and. segment%z <= owned_upper_z
   end function line_segment_is_local

   subroutine update_line_probe_output(this, step, electric_field)
      type(line_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(field_data_t), intent(in) :: electric_field
#ifdef CompileWithMPI
      integer :: ierr
      real(kind=RKIND) :: local_value
#endif

      if (this%domain%domainType /= TIME_DOMAIN) return
      if (size(this%segments) == 0 .and. .not. this%isDistributed) return
      this%nTime = this%nTime + 1
      this%timeStep(this%nTime) = step
      this%valueForTime(this%nTime) = calculate_line_integral(this%segments, electric_field)
#ifdef CompileWithMPI
      if (this%isDistributed) then
         local_value = this%valueForTime(this%nTime)
         call MPI_Allreduce(local_value, this%valueForTime(this%nTime), 1, REALSIZE, MPI_SUM, MPI_COMM_WORLD, ierr)
      end if
#endif
   end subroutine update_line_probe_output

   subroutine complete_line_probe_sample(this)
      type(line_probe_output_t), intent(inout) :: this

      if (this%nTime == 0) return
      this%nTimesFlushed = this%nTimesFlushed + this%nTime
      this%nTime = 0
      this%timeStep = 0.0_RKIND_tiempo
      this%valueForTime = 0.0_RKIND
   end subroutine complete_line_probe_sample

   subroutine flush_line_probe_output(this)
      type(line_probe_output_t), intent(inout) :: this
       integer :: index, ios, unit

      if (this%nTime == 0) return
      if (this%isWriter) then
         open(newunit=unit, file=this%artifacts(1)%relative_path, status='old', action='write', position='append', iostat=ios)
         if (ios /= 0) return
         do index = 1, this%nTime
            write(unit, '(ES24.16E3,1X,ES24.16E3)', iostat=ios) this%timeStep(index), this%valueForTime(index)
            if (ios /= 0) exit
         end do
         close(unit)
         if (ios /= 0) return
      end if
      call complete_line_probe_sample(this)
   end subroutine flush_line_probe_output

end module lineProbeOutput_m
