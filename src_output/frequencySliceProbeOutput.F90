module frequencySliceProbeOutput_m
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use FDETYPES_m
   use utils_m
   use allocationUtils_m, only: alloc_and_init
   use report_m
   use outputTypes_m
   use outputUtils_m
   use volumicProbeUtils_m
   use directoryUtils_m
   use outputBinary_m, only: validate_binary_layout, write_binary_complex_record64, BINARY_WRITER_SUCCESS
   use outputVisualisation_m, only: initialise_frequency_slice_visualisation, begin_visualisation_step, &
                                   write_visualisation_attribute, write_visualisation_attribute_hyperslab, end_visualisation_step, &
                                    close_visualisation, verify_volumetric_visualisation, VISUALISATION_SUCCESS, &
                                    VISUALISATION_ATTRIBUTE_X, VISUALISATION_ATTRIBUTE_Y, VISUALISATION_ATTRIBUTE_Z, &
                                    VISUALISATION_ATTRIBUTE_X_PHASE, VISUALISATION_ATTRIBUTE_Y_PHASE, &
                                    VISUALISATION_ATTRIBUTE_Z_PHASE, VISUALISATION_ATTRIBUTE_TAG, &
                                    VISUALISATION_ATTRIBUTE_MEDIA, VISUALISATION_ATTRIBUTE_TAG_X, &
                                    VISUALISATION_ATTRIBUTE_TAG_Y, VISUALISATION_ATTRIBUTE_TAG_Z, &
                                    VISUALISATION_ATTRIBUTE_MEDIA_X, VISUALISATION_ATTRIBUTE_MEDIA_Y, &
                                    VISUALISATION_ATTRIBUTE_MEDIA_Z
   use outputCollective_m, only: OUTPUT_PUBLICATION_COLLECTIVE, OUTPUT_PUBLICATION_ROOT_AGGREGATION
   use outputTransport_m, only: output_transport_t, init_output_transport, transfer_flush_batch, &
                                OUTPUT_TRANSPORT_SUCCESS
#ifdef CompileWithMPI
   use mpi
#endif
   implicit none
#if defined(CompileWithMPI) && defined(IFXCompiler)
   external :: MPI_Allgather, MPI_Allgatherv
#endif
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: init_frequency_slice_probe_output
   public :: update_frequency_slice_probe_output
   public :: flush_frequency_slice_probe_output
   public :: close_frequency_slice_probe_output
   !===========================

   !===========================
   !  Private interface summary
   !===========================
   private :: save_field
   private :: save_field_module
   private :: save_field_component
   private :: save_current
   private :: save_current_module
   private :: save_current_component
   private :: initialise_frequency_metadata
   !===========================

   !===========================

contains

   subroutine init_frequency_slice_probe_output(this, publication, timeInterval, field, domain, &
                                                outputTypeExtension, control, problemInfo)
      type(frequency_slice_probe_output_t), intent(out) :: this
      type(volumetric_publication_t), intent(in) :: publication
      real(kind=RKIND_tiempo), intent(in) :: timeInterval
      integer(kind=SINGLE), intent(in) :: field
      type(domain_t), intent(in) :: domain
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      type(sim_control_t), intent(in) :: control
      type(problem_info_t), intent(in) :: problemInfo

      integer :: i
      integer :: error
      character(len=BUFSIZE) :: filename
#ifdef CompileWithMPI
      integer :: mpi_error
#endif

      this%publication = publication
      this%mainCoords = publication%global_lower
      this%auxCoords = publication%global_upper
      if (publication%local_participates) call set_local_bounds(this)
      this%component = field !This can refer to electric, magnetic or currentDensity
      this%domain = domain
      this%nFreq = domain%fnum
      this%quadratureDt = timeInterval

      call alloc_and_init(this%frequencySlice, this%nFreq, 0.0_RKIND)
      call init_frequency_slice(this%frequencySlice, this%domain)

      this%nPoints = 0_SINGLE
      if (publication%local_participates) then
         call find_and_store_important_coords(this%mainCoords, this%auxCoords, this%component, &
                                              problemInfo, this%nPoints, this%coords)
      else
         allocate (this%coords(3, 0))
      end if
      allocate (this%tagNumber(3, this%nPoints), this%mediaType(3, this%nPoints))
      this%tagNumber = 0_IKINDMTAG
      this%mediaType = -1.0_RKIND
      if (publication%local_participates) call store_classification(this, problemInfo)

      call alloc_and_init(this%xValueForFreq, this%nFreq, this%nPoints, (0.0_CKIND, 0.0_CKIND))
      call alloc_and_init(this%yValueForFreq, this%nFreq, this%nPoints, (0.0_CKIND, 0.0_CKIND))
      call alloc_and_init(this%zValueForFreq, this%nFreq, this%nPoints, (0.0_CKIND, 0.0_CKIND))

      call alloc_and_init(this%auxExp_E, this%nFreq, (0.0_CKIND, 0.0_CKIND))
      call alloc_and_init(this%auxExp_H, this%nFreq, (0.0_CKIND, 0.0_CKIND))

      do i = 1, this%nFreq
         this%auxExp_E(i) = mcpi2*this%frequencySlice(i)
         this%auxExp_H(i) = this%auxExp_E(i)
      end do

      this%path = get_output_path_freq(this, outputTypeExtension, field, control)

      filename = get_last_component(this%path)
      this%filesPath = join_path(this%path, filename)

      call initialise_frequency_metadata(this, error, control%mpidir)
      if (error /= 0) call StopOnError(control%layoutnumber, control%num_procs, &
                                       'Unable to initialise frequency slice output metadata')

      if (.not. publication%local_participates) then
         this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_ACTIVE
         return
      end if

      call gather_global_coordinates(this, control)

      if (publication%local_is_owner) then
         call create_folder(this%path, error)
         if (error /= 0) call StopOnError(control%layoutnumber, control%num_procs, &
                                          'Unable to create frequency slice output directory')
         call create_bin_file(this%filesPath, error)
         if (error /= 0) call StopOnError(control%layoutnumber, control%num_procs, &
                                          'Unable to create frequency slice binary output')
      end if
#ifdef CompileWithMPI
      call MPI_Barrier(this%publication%communicator, mpi_error)
      if (mpi_error /= MPI_SUCCESS) call StopOnError(control%layoutnumber, control%num_procs, &
                                                     'Unable to synchronise frequency slice output participants')
#endif

      if (publication%mode == OUTPUT_PUBLICATION_COLLECTIVE .or. publication%local_is_owner) then
         call create_frequency_writer(this, problemInfo, error)
         if (error /= 0) call StopOnError(control%layoutnumber, control%num_procs, &
                                          'Unable to initialise frequency slice HDF5 output')
      end if

      this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_ACTIVE
   end subroutine init_frequency_slice_probe_output

   subroutine set_local_bounds(this)
      type(frequency_slice_probe_output_t), intent(inout) :: this

      this%mainCoords%x = this%publication%global_lower%x + int(this%publication%file_offset(1), SINGLE)
      this%mainCoords%y = this%publication%global_lower%y + int(this%publication%file_offset(2), SINGLE)
      this%mainCoords%z = this%publication%global_lower%z + int(this%publication%file_offset(3), SINGLE)
      this%auxCoords%x = this%mainCoords%x + int(this%publication%local_shape(1), SINGLE) - 1_SINGLE
      this%auxCoords%y = this%mainCoords%y + int(this%publication%local_shape(2), SINGLE) - 1_SINGLE
      this%auxCoords%z = this%mainCoords%z + int(this%publication%local_shape(3), SINGLE) - 1_SINGLE
   end subroutine set_local_bounds

   subroutine gather_global_coordinates(this, control)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(sim_control_t), intent(in) :: control
      integer, allocatable :: point_counts(:), point_displacements(:)
      integer, allocatable :: coordinate_counts(:), coordinate_displacements(:)
      integer :: i, rank_count
#ifdef CompileWithMPI
      integer :: mpi_error
#endif

      rank_count = this%publication%communicator_size
      allocate (point_counts(rank_count), point_displacements(rank_count))
      point_counts = 0
      point_displacements = 0
#ifdef CompileWithMPI
      call MPI_Allgather(this%nPoints, 1, MPI_INTEGER, point_counts, 1, MPI_INTEGER, &
                         this%publication%communicator, mpi_error)
      if (mpi_error /= MPI_SUCCESS) call StopOnError(control%layoutnumber, control%num_procs, &
                                                     'Unable to gather frequency slice point counts')
#else
      point_counts(1) = this%nPoints
#endif

      do i = 2, rank_count
         point_displacements(i) = point_displacements(i - 1) + point_counts(i - 1)
      end do
      this%publication%point_offset = int(point_displacements(this%publication%communicator_rank + 1), int64)
      this%publication%global_point_count = sum(int(point_counts, int64))
      if (this%publication%global_point_count <= 0_int64 .or. &
          this%publication%global_point_count > int(huge(1), int64)) then
         call StopOnError(control%layoutnumber, control%num_procs, &
                          'Invalid global frequency slice point count')
      end if

      allocate (this%globalCoords(3, int(this%publication%global_point_count)))
      allocate (coordinate_counts(rank_count), coordinate_displacements(rank_count))
      coordinate_counts = 3*point_counts
      coordinate_displacements = 3*point_displacements
#ifdef CompileWithMPI
      call MPI_Allgatherv(this%coords, 3*this%nPoints, MPI_INTEGER, this%globalCoords, &
                          coordinate_counts, coordinate_displacements, MPI_INTEGER, &
                          this%publication%communicator, mpi_error)
      if (mpi_error /= MPI_SUCCESS) call StopOnError(control%layoutnumber, control%num_procs, &
                                                     'Unable to gather frequency slice coordinates')
#else
      this%globalCoords = this%coords
#endif
   end subroutine gather_global_coordinates

   subroutine create_frequency_writer(this, problemInfo, error)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(out) :: error

      real(real64), allocatable :: points(:, :)
      character(len=BUFSIZE) :: attribute_names(14), diagnostic
      logical :: attribute_enabled(14)
      integer :: i, status

      error = 0
      attribute_names = ''
      attribute_enabled = .false.
      attribute_names(1:6) = [character(len=BUFSIZE) :: &
                              'xMagnitude', 'yMagnitude', 'zMagnitude', 'xPhase', 'yPhase', 'zPhase']
      attribute_enabled(1:6) = .true.
      if (any([iCur, iMEC, iMHC] == this%component)) then
         attribute_names(VISUALISATION_ATTRIBUTE_TAG_X) = 'tagnumber_x'
         attribute_names(VISUALISATION_ATTRIBUTE_TAG_Y) = 'tagnumber_y'
         attribute_names(VISUALISATION_ATTRIBUTE_TAG_Z) = 'tagnumber_z'
         attribute_names(VISUALISATION_ATTRIBUTE_MEDIA_X) = 'mediatype_x'
         attribute_names(VISUALISATION_ATTRIBUTE_MEDIA_Y) = 'mediatype_y'
         attribute_names(VISUALISATION_ATTRIBUTE_MEDIA_Z) = 'mediatype_z'
         attribute_enabled(VISUALISATION_ATTRIBUTE_TAG_X:VISUALISATION_ATTRIBUTE_MEDIA_Z) = .true.
      else
         attribute_names(VISUALISATION_ATTRIBUTE_TAG) = 'tagnumber'
         attribute_names(VISUALISATION_ATTRIBUTE_MEDIA) = 'mediatype'
         attribute_enabled(VISUALISATION_ATTRIBUTE_TAG) = .true.
         attribute_enabled(VISUALISATION_ATTRIBUTE_MEDIA) = .true.
      end if
      allocate (points(3, int(this%publication%global_point_count)))
      do i = 1, int(this%publication%global_point_count)
         if (associated(problemInfo%xSteps) .and. &
             associated(problemInfo%ySteps) .and. &
             associated(problemInfo%zSteps)) then
            if (this%globalCoords(1, i) >= lbound(problemInfo%xSteps, 1) .and. &
                this%globalCoords(1, i) <= ubound(problemInfo%xSteps, 1) .and. &
                this%globalCoords(2, i) >= lbound(problemInfo%ySteps, 1) .and. &
                this%globalCoords(2, i) <= ubound(problemInfo%ySteps, 1) .and. &
                this%globalCoords(3, i) >= lbound(problemInfo%zSteps, 1) .and. &
                this%globalCoords(3, i) <= ubound(problemInfo%zSteps, 1)) then
               points(:, i) = [real(problemInfo%xSteps(this%globalCoords(1, i)), real64), &
                               real(problemInfo%ySteps(this%globalCoords(2, i)), real64), &
                               real(problemInfo%zSteps(this%globalCoords(3, i)), real64)]
            else
               points(:, i) = real(this%globalCoords(:, i), real64)
            end if
         else
            points(:, i) = real(this%globalCoords(:, i), real64)
         end if
      end do

      call initialise_frequency_slice_visualisation(this%visualisation, trim(this%filesPath), points, &
                                       attribute_names, attribute_enabled, this%publication%mode == OUTPUT_PUBLICATION_COLLECTIVE, &
                                                    this%publication%communicator, status, diagnostic)
      if (status /= VISUALISATION_SUCCESS) then
         error = 1
         print *, trim(diagnostic)
      end if
      deallocate (points)
   end subroutine create_frequency_writer

   subroutine create_bin_file(filePath, error)
      character(len=*), intent(in) :: filePath
      integer, intent(out) :: error
      call create_file_with_path(add_extension(filePath, binaryExtension), error)
   end subroutine

   subroutine initialise_frequency_metadata(this, error, mpidir)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      integer, intent(out) :: error
      integer(kind=SINGLE), intent(in) :: mpidir
      character(len=BUFSIZE) :: base_name

      base_name = get_last_component(this%filesPath)
      this%metadata%probe_id = trim(base_name)
      this%metadata%quantity = get_prefix_extension(this%component, mpidir)
      this%metadata%lower_bound = this%publication%global_lower
      this%metadata%upper_bound = this%publication%global_upper
      this%metadata%domain_type = FREQUENCY_DOMAIN
      this%metadata%ownership%participant_ranks = this%publication%participant_ranks
      this%metadata%ownership%scalar_writer_rank = this%publication%owner_rank
      if (allocated(this%metadata%artifacts)) deallocate (this%metadata%artifacts)
      allocate (this%metadata%artifacts(3))
      this%metadata%artifacts(1)%kind = OUTPUT_ARTIFACT_BINARY
      this%metadata%artifacts(1)%relative_path = trim(base_name)//binaryExtension
      this%metadata%artifacts(1)%byte_order = BINARY_ENDIAN_LITTLE
      this%metadata%artifacts(1)%numeric_representation = BINARY_NUMERIC_REAL64
      this%metadata%artifacts(1)%complex_representation = BINARY_COMPLEX_REAL_IMAG
      this%metadata%artifacts(1)%record_bytes = 80
      this%metadata%artifacts(1)%component_order = &
         'frequency,x,y,z,Ex.real,Ex.imag,Ey.real,Ey.imag,Ez.real,Ez.imag'
      this%metadata%artifacts(2)%kind = OUTPUT_ARTIFACT_VISUALISATION_METADATA
      this%metadata%artifacts(2)%relative_path = trim(base_name)//'.xdmf'
      this%metadata%artifacts(3)%kind = OUTPUT_ARTIFACT_VISUALISATION_DATA
      this%metadata%artifacts(3)%relative_path = trim(base_name)//'.h5'

      call validate_binary_layout(this%metadata%artifacts(1), error)
      if (error /= BINARY_WRITER_SUCCESS) return
      error = 0
   end subroutine initialise_frequency_metadata

   subroutine write_visualisation(this)
      type(frequency_slice_probe_output_t), intent(inout) :: this

      type(output_transport_t) :: transport
      real(real64), allocatable :: xMagnitude(:), yMagnitude(:), zMagnitude(:)
      real(real64), allocatable :: xPhase(:), yPhase(:), zPhase(:)
      character(len=BUFSIZE) :: diagnostic
      integer :: f, i, status, transport_status

      allocate (xMagnitude(this%nPoints), yMagnitude(this%nPoints), &
                zMagnitude(this%nPoints), xPhase(this%nPoints), yPhase(this%nPoints), &
                zPhase(this%nPoints))
      if (this%publication%mode == OUTPUT_PUBLICATION_ROOT_AGGREGATION) then
         call init_output_transport(transport, 0, transport_status, this%publication%communicator)
         if (transport_status /= OUTPUT_TRANSPORT_SUCCESS) then
            call StopOnError(0, 0, 'Unable to initialise frequency slice output transport')
         end if
      end if
      call write_classification_attributes(this, transport)

      do f = 1, this%nFreq
         do i = 1, this%nPoints
            xMagnitude(i) = real(abs(this%xValueForFreq(f, i)), real64)
            yMagnitude(i) = real(abs(this%yValueForFreq(f, i)), real64)
            zMagnitude(i) = real(abs(this%zValueForFreq(f, i)), real64)
            xPhase(i) = atan2(real(aimag(this%xValueForFreq(f, i)), real64), &
                              real(this%xValueForFreq(f, i), real64))
            yPhase(i) = atan2(real(aimag(this%yValueForFreq(f, i)), real64), &
                              real(this%yValueForFreq(f, i), real64))
            zPhase(i) = atan2(real(aimag(this%zValueForFreq(f, i)), real64), &
                              real(this%zValueForFreq(f, i), real64))
         end do

         select case (this%publication%mode)
         case (OUTPUT_PUBLICATION_COLLECTIVE)
            call begin_visualisation_step(this%visualisation, real(this%frequencySlice(f), real64), &
                                          status, diagnostic)
            call require_visualisation_success(status, diagnostic)
            call write_collective_attribute(this, VISUALISATION_ATTRIBUTE_X, xMagnitude)
            call write_collective_attribute(this, VISUALISATION_ATTRIBUTE_Y, yMagnitude)
            call write_collective_attribute(this, VISUALISATION_ATTRIBUTE_Z, zMagnitude)
            call write_collective_attribute(this, VISUALISATION_ATTRIBUTE_X_PHASE, xPhase)
            call write_collective_attribute(this, VISUALISATION_ATTRIBUTE_Y_PHASE, yPhase)
            call write_collective_attribute(this, VISUALISATION_ATTRIBUTE_Z_PHASE, zPhase)
            call end_visualisation_step(this%visualisation, status, diagnostic)
            call require_visualisation_success(status, diagnostic)
         case (OUTPUT_PUBLICATION_ROOT_AGGREGATION)
            if (this%publication%local_is_owner) then
               call begin_visualisation_step(this%visualisation, real(this%frequencySlice(f), real64), &
                                             status, diagnostic)
               call require_visualisation_success(status, diagnostic)
            end if
            call gather_and_write_attribute(this, transport, VISUALISATION_ATTRIBUTE_X, xMagnitude)
            call gather_and_write_attribute(this, transport, VISUALISATION_ATTRIBUTE_Y, yMagnitude)
            call gather_and_write_attribute(this, transport, VISUALISATION_ATTRIBUTE_Z, zMagnitude)
            call gather_and_write_attribute(this, transport, VISUALISATION_ATTRIBUTE_X_PHASE, xPhase)
            call gather_and_write_attribute(this, transport, VISUALISATION_ATTRIBUTE_Y_PHASE, yPhase)
            call gather_and_write_attribute(this, transport, VISUALISATION_ATTRIBUTE_Z_PHASE, zPhase)
            if (this%publication%local_is_owner) then
               call end_visualisation_step(this%visualisation, status, diagnostic)
               call require_visualisation_success(status, diagnostic)
            end if
         case default
            call StopOnError(0, 0, 'Unsupported frequency slice publication mode')
         end select
      end do
      deallocate (xMagnitude, yMagnitude, zMagnitude, xPhase, yPhase, zPhase)
   end subroutine write_visualisation

   subroutine write_classification_attributes(this, transport)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(output_transport_t), intent(in) :: transport
      integer, parameter :: tag_attributes(3) = [VISUALISATION_ATTRIBUTE_TAG_X, &
                                                 VISUALISATION_ATTRIBUTE_TAG_Y, &
                                                 VISUALISATION_ATTRIBUTE_TAG_Z]
      integer, parameter :: media_attributes(3) = [VISUALISATION_ATTRIBUTE_MEDIA_X, &
                                                   VISUALISATION_ATTRIBUTE_MEDIA_Y, &
                                                   VISUALISATION_ATTRIBUTE_MEDIA_Z]
      integer :: axis

      if (any([iCur, iMEC, iMHC] == this%component)) then
         do axis = 1, 3
            call write_classification_attribute(this, transport, tag_attributes(axis), axis, .true.)
            call write_classification_attribute(this, transport, media_attributes(axis), axis, .false.)
         end do
      else
         do axis = 1, 3
            if (.not. classification_axis_enabled(this%component, axis)) cycle
            call write_classification_attribute(this, transport, VISUALISATION_ATTRIBUTE_TAG, axis, .true.)
            call write_classification_attribute(this, transport, VISUALISATION_ATTRIBUTE_MEDIA, axis, .false.)
         end do
      end if
   end subroutine write_classification_attributes

   subroutine write_classification_attribute(this, transport, attribute, axis, write_tag)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(output_transport_t), intent(in) :: transport
      integer, intent(in) :: attribute, axis
      logical, intent(in) :: write_tag
      real(real64), allocatable :: values(:)
      integer :: i

      allocate (values(this%nPoints))
      do i = 1, this%nPoints
         if (write_tag) then
            values(i) = real(this%tagNumber(axis, i), real64)
         else
            values(i) = real(this%mediaType(axis, i), real64)
         end if
      end do
      select case (this%publication%mode)
      case (OUTPUT_PUBLICATION_COLLECTIVE)
         call write_collective_attribute(this, attribute, values)
      case (OUTPUT_PUBLICATION_ROOT_AGGREGATION)
         call gather_and_write_attribute(this, transport, attribute, values)
      end select
   end subroutine write_classification_attribute

   subroutine write_collective_attribute(this, attribute, values)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      integer, intent(in) :: attribute
      real(real64), intent(in) :: values(:)
      character(len=BUFSIZE) :: diagnostic
      integer :: status

      call write_visualisation_attribute_hyperslab(this%visualisation, attribute, values, &
                                                   [this%publication%point_offset], [int(this%nPoints, int64)], status, diagnostic)
      call require_visualisation_success(status, diagnostic)
   end subroutine write_collective_attribute

   subroutine gather_and_write_attribute(this, transport, attribute, local_values)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(output_transport_t), intent(in) :: transport
      integer, intent(in) :: attribute
      real(real64), intent(in) :: local_values(:)
      real(real64), allocatable :: gathered_values(:)
      integer, allocatable :: counts(:), displacements(:)
      character(len=BUFSIZE) :: diagnostic
      integer :: status, transport_status

      call transfer_flush_batch(transport, local_values, gathered_values, counts, displacements, transport_status)
      if (transport_status /= OUTPUT_TRANSPORT_SUCCESS) then
         call StopOnError(0, 0, 'Unable to gather frequency slice visualisation values')
      end if
      if (.not. this%publication%local_is_owner) return
      if (size(gathered_values, kind=int64) /= this%publication%global_point_count) then
         call StopOnError(0, 0, 'Invalid gathered frequency slice visualisation size')
      end if
      call write_visualisation_attribute(this%visualisation, attribute, gathered_values, status, diagnostic)
      call require_visualisation_success(status, diagnostic)
   end subroutine gather_and_write_attribute

   subroutine require_visualisation_success(status, diagnostic)
      integer, intent(in) :: status
      character(len=*), intent(in) :: diagnostic

      if (status /= VISUALISATION_SUCCESS) then
         call StopOnError(0, 0, 'Frequency slice HDF5 write failed: '//trim(diagnostic))
      end if
   end subroutine require_visualisation_success

   subroutine store_classification(this, problemInfo)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(problem_info_t), intent(in) :: problemInfo
      integer :: axis, field, i

      do axis = 1, 3
         if (.not. classification_axis_enabled(this%component, axis)) cycle
         field = get_volumetric_classification_field(this%component, axis)
         do i = 1, this%nPoints
            if (.not. classification_point_is_valid(this%component, field, this%coords(:, i), problemInfo)) cycle
            this%tagNumber(axis, i) = get_output_tag_number(field, this%coords(:, i), problemInfo)
            this%mediaType(axis, i) = get_output_media_type(field, this%coords(:, i), problemInfo)
         end do
      end do
   end subroutine store_classification

   logical function classification_axis_enabled(component, axis)
      integer(kind=SINGLE), intent(in) :: component
      integer, intent(in) :: axis

      classification_axis_enabled = any([iCur, iMEC, iMHC] == component) .or. &
                                    (axis == 1 .and. any([iCurX, iExC, iHxC] == component)) .or. &
                                    (axis == 2 .and. any([iCurY, iEyC, iHyC] == component)) .or. &
                                    (axis == 3 .and. any([iCurZ, iEzC, iHzC] == component))
   end function classification_axis_enabled

   logical function classification_point_is_valid(component, field, position, problemInfo)
      integer(kind=SINGLE), intent(in) :: component
      integer, intent(in) :: field, position(3)
      type(problem_info_t), intent(in) :: problemInfo

      if (any([iCur, iCurX, iCurY, iCurZ] == component)) then
         classification_point_is_valid = isValidPointForCurrent(field, position(1), position(2), position(3), problemInfo)
      else
         classification_point_is_valid = isValidPointForField(field, position(1), position(2), position(3), problemInfo)
      end if
   end function classification_point_is_valid

   subroutine write_bin_file(this)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(output_transport_t) :: transport
      real(real64), allocatable :: local_records(:), gathered_records(:), canonical_records(:)
      integer, allocatable :: counts(:), displacements(:)
      integer :: i, f, record_index, status, transport_status, frequency_offset

      call init_output_transport(transport, 0, transport_status, this%publication%communicator)
      if (transport_status /= OUTPUT_TRANSPORT_SUCCESS) then
         call StopOnError(0, 0, 'Unable to initialise frequency slice binary transport')
      end if
      allocate (local_records(10*this%nPoints))
      if (this%publication%local_is_owner) then
         allocate (canonical_records(10*int(this%publication%global_point_count)*this%nFreq))
      end if

      do f = 1, this%nFreq
         record_index = 0
         do i = 1, this%nPoints
            local_records(record_index + 1:record_index + 10) = [real(this%frequencySlice(f), real64), &
                                real(this%coords(1, i), real64), real(this%coords(2, i), real64), real(this%coords(3, i), real64), &
                                            real(this%xValueForFreq(f, i), real64), real(aimag(this%xValueForFreq(f, i)), real64), &
                                            real(this%yValueForFreq(f, i), real64), real(aimag(this%yValueForFreq(f, i)), real64), &
                                              real(this%zValueForFreq(f, i), real64), real(aimag(this%zValueForFreq(f, i)), real64)]
            record_index = record_index + 10
         end do
         call transfer_flush_batch(transport, local_records, gathered_records, counts, displacements, transport_status)
         if (transport_status /= OUTPUT_TRANSPORT_SUCCESS) then
            call StopOnError(0, 0, 'Unable to gather frequency slice binary records')
         end if
         if (this%publication%local_is_owner) then
            if (size(gathered_records, kind=int64) /= 10_int64*this%publication%global_point_count) then
               call StopOnError(0, 0, 'Invalid gathered frequency slice binary size')
            end if
            frequency_offset = 10*int(this%publication%global_point_count)*(f - 1)
            canonical_records(frequency_offset + 1:frequency_offset + size(gathered_records)) = gathered_records
         end if
      end do
      if (this%publication%local_is_owner) then
         call write_binary_complex_record64(add_extension(this%filesPath, binaryExtension), &
                                            this%metadata%artifacts(1), canonical_records, status)
         if (status /= BINARY_WRITER_SUCCESS) then
            call StopOnError(0, 0, 'Unable to replace frequency slice binary output')
         end if
      end if
   end subroutine

   function get_output_path_freq(this, outputTypeExtension, field, control) result(outputPath)
      type(frequency_slice_probe_output_t), intent(in) :: this
      character(len=*), intent(in) :: outputTypeExtension
      integer(kind=SINGLE), intent(in) :: field
      type(sim_control_t), intent(in) :: control
      character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
      character(len=BUFSIZE) :: outputPath
      probeBoundsExtension = get_coordinates_extension(this%publication%global_lower, &
                                                       this%publication%global_upper, control%mpidir)
      prefixFieldExtension = get_prefix_extension(field, control%mpidir)
      outputPath = &
         trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
   end function get_output_path_freq

   subroutine update_frequency_slice_probe_output(this, step, fieldsReference, control, problemInfo)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(sim_control_t), intent(in) :: control
      type(problem_info_t), intent(in) :: problemInfo
      type(fields_reference_t), intent(in) :: fieldsReference

      integer(kind=4) :: request
      if (.not. this%publication%local_participates .or. &
          this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_COMPLETE .or. &
          this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED) return
      request = this%component

      if (any(VOLUMIC_M_MEASURE == request)) then
         select case (request)
         case (iCur); call save_current_module(this, fieldsReference, step, problemInfo)
         case (iMEC); call save_field_module(this, fieldsReference%E, step, request, problemInfo)
         case (iMHC); call save_field_module(this, fieldsReference%H, step, request, problemInfo)
         case default; call StopOnError(control%layoutnumber, control%num_procs, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_X_MEASURE == request)) then
         select case (request)
         case (iCurX); call save_current_component(this, this%xValueForFreq, fieldsReference, problemInfo, iEx, this%auxExp_E, this%nFreq, step)
         case (iExC); call save_field_component(this, this%xValueForFreq, fieldsReference%E%x, step, problemInfo, iEx)
         case (iHxC); call save_field_component(this, this%xValueForFreq, fieldsReference%H%x, step, problemInfo, iHx)
         case default; call StopOnError(control%layoutnumber, control%num_procs, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_Y_MEASURE == request)) then
         select case (request)
         case (iCurY); call save_current_component(this, this%yValueForFreq, fieldsReference, problemInfo, iEy, this%auxExp_E, this%nFreq, step)
         case (iEyC); call save_field_component(this, this%yValueForFreq, fieldsReference%E%y, step, problemInfo, iEy)
         case (iHyC); call save_field_component(this, this%yValueForFreq, fieldsReference%H%y, step, problemInfo, iHy)
         case default; call StopOnError(control%layoutnumber, control%num_procs, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_Z_MEASURE == request)) then
         select case (request)
         case (iCurZ); call save_current_component(this, this%zValueForFreq, fieldsReference, problemInfo, iEz, this%auxExp_E, this%nFreq, step)
         case (iEzC); call save_field_component(this, this%zValueForFreq, fieldsReference%E%z, step, problemInfo, iEz)
         case (iHzC); call save_field_component(this, this%zValueForFreq, fieldsReference%H%z, step, problemInfo, iHz)
         case default; call StopOnError(control%layoutnumber, control%num_procs, "Volumic measure not supported")
         end select
      end if
   end subroutine update_frequency_slice_probe_output

   subroutine save_current_module(this, fieldsReference, step, problemInfo)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(fields_reference_t), intent(in) :: fieldsReference
      type(problem_info_t), intent(in) :: problemInfo
      real(kind=RKIND_tiempo), intent(in) :: step

      integer :: i, j, k, coordIdx

      coordIdx = 0
      do k = this%mainCoords%z, this%auxCoords%z
      do j = this%mainCoords%y, this%auxCoords%y
      do i = this%mainCoords%x, this%auxCoords%x
         if (isValidPointForCurrent(iCur, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_current(this%xValueForFreq, iEx, coordIdx, i, j, k, fieldsReference, this%auxExp_E, &
                              this%quadratureDt, this%nFreq, step)
            call save_current(this%yValueForFreq, iEy, coordIdx, i, j, k, fieldsReference, this%auxExp_E, &
                              this%quadratureDt, this%nFreq, step)
            call save_current(this%zValueForFreq, iEz, coordIdx, i, j, k, fieldsReference, this%auxExp_E, &
                              this%quadratureDt, this%nFreq, step)
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_current_component(this, currentData, fieldsReference, problemInfo, fieldDir, auxExp, nFreq, step)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      complex(kind=CKIND), intent(inout) :: currentData(:, :)
      type(fields_reference_t), intent(in) :: fieldsReference
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: fieldDir, nFreq
      complex(kind=ckind), intent(in), dimension(:) :: auxExp
      real(kind=RKIND_tiempo), intent(in) :: step

      integer :: i, j, k, coordIdx

      coordIdx = 0
      do k = this%mainCoords%z, this%auxCoords%z
      do j = this%mainCoords%y, this%auxCoords%y
      do i = this%mainCoords%x, this%auxCoords%x
         if (isValidPointForCurrent(fieldDir, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_current(currentData, fieldDir, coordIdx, i, j, k, fieldsReference, auxExp, &
                              this%quadratureDt, nFreq, step)
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_current(valorComplex, direction, coordIdx, i, j, k, fieldsReference, auxExponential, &
                           quadratureDt, nFreq, step)
      integer, intent(in) :: direction
      complex(kind=CKIND), intent(inout) :: valorComplex(:, :)
      complex(kind=CKIND), intent(in) :: auxExponential(:)
      integer, intent(in) :: i, j, k, coordIdx, nFreq
      type(fields_reference_t), intent(in) :: fieldsReference
      real(kind=RKIND_tiempo), intent(in) :: quadratureDt
      real(kind=RKIND_tiempo), intent(in) :: step

      integer :: iter
      complex(kind=CKIND) :: z_cplx = (0.0_RKIND, 0.0_RKIND)
      real(kind=rkind) :: jdir

      jdir = computej(direction, i, j, k, fieldsReference)

      do iter = 1, nFreq
         valorComplex(iter, coordIdx) = valorComplex(iter, coordIdx) + quadratureDt* &
                                        exp(auxExponential(iter)*step)*jdir
      end do
   end subroutine

   subroutine save_field_module(this, fieldInfo, simTime, request, problemInfo)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(field_data_t), intent(in) :: fieldInfo
      real(kind=RKIND_tiempo), intent(in) :: simTime
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: request

      complex(kind=CKIND), dimension(this%nFreq) :: auxExponential
      integer :: i, j, k, coordIdx

      if (iMHC == request) auxExponential = this%quadratureDt*exp(this%auxExp_H*(simTime + 0.5_RKIND_tiempo*this%quadratureDt))
      if (iMEC == request) auxExponential = this%quadratureDt*exp(this%auxExp_E*simTime)

      coordIdx = 0
      do k = this%mainCoords%z, this%auxCoords%z
      do j = this%mainCoords%y, this%auxCoords%y
      do i = this%mainCoords%x, this%auxCoords%x
         if (isValidPointForField(request, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_field(this%xValueForFreq, auxExponential, fieldInfo%x(i, j, k), this%nFreq, coordIdx)
            call save_field(this%yValueForFreq, auxExponential, fieldInfo%y(i, j, k), this%nFreq, coordIdx)
            call save_field(this%zValueForFreq, auxExponential, fieldInfo%z(i, j, k), this%nFreq, coordIdx)
         end if
      end do
      end do
      end do

   end subroutine

   subroutine save_field_component(this, fieldData, fieldComponent, simTime, problemInfo, fieldDir)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      complex(kind=CKIND), intent(inout) :: fieldData(:, :)
      real(kind=RKIND), pointer, intent(in) :: fieldComponent(:, :, :)
      real(kind=RKIND_tiempo), intent(in) :: simTime
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: fieldDir

      complex(kind=CKIND), dimension(this%nFreq) :: auxExponential
      integer :: i, j, k, coordIdx

      if (any(MAGNETIC_FIELD_DIRECTION == fieldDir)) then
         auxExponential = this%quadratureDt*exp(this%auxExp_H*(simTime + 0.5_RKIND_tiempo*this%quadratureDt))
      end if
      if (any(ELECTRIC_FIELD_DIRECTION == fieldDir)) auxExponential = this%quadratureDt*exp(this%auxExp_E*simTime)

      coordIdx = 0
      do k = this%mainCoords%z, this%auxCoords%z
      do j = this%mainCoords%y, this%auxCoords%y
      do i = this%mainCoords%x, this%auxCoords%x
         if (isValidPointForField(fieldDir, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_field(fieldData, auxExponential, fieldComponent(i, j, k), this%nFreq, coordIdx)
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_field(valorComplex, auxExp, fieldValue, nFreq, coordIdx)
      complex(kind=CKIND), intent(inout) :: valorComplex(:, :)
      complex(kind=CKIND), intent(in) :: auxExp(:)
      real(KIND=RKIND), intent(in) :: fieldValue
      integer(KIND=SINGLE), intent(in) :: nFreq, coordIdx

      integer :: freq

      do freq = 1, nFreq
         valorComplex(freq, coordIdx) = valorComplex(freq, coordIdx) + auxExp(freq)*fieldValue
      end do
   end subroutine

   subroutine flush_frequency_slice_probe_output(this)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      if (.not. this%publication%local_participates .or. &
          this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_COMPLETE .or. &
          this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED) return
      call write_bin_file(this)
   end subroutine flush_frequency_slice_probe_output

   subroutine close_frequency_slice_probe_output(this)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      integer :: error
      character(len=BUFSIZE) :: diagnostic
#ifdef CompileWithMPI
      integer :: mpi_error
#endif

      if (.not. this%publication%local_participates) then
         this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
         this%metadata%lifecycle%diagnostic = ''
         return
      end if
      if ((this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_COMPLETE .or. &
           this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED) .and. &
          .not. this%publication%owns_communicator) return
      if (this%metadata%lifecycle%state /= OUTPUT_LIFECYCLE_COMPLETE .and. &
          this%metadata%lifecycle%state /= OUTPUT_LIFECYCLE_FAILED) then
         call flush_frequency_slice_probe_output(this)
         call write_visualisation(this)
      end if
      call close_visualisation(this%visualisation, error, diagnostic)
      if (error /= VISUALISATION_SUCCESS) then
         call StopOnError(0, 0, 'Unable to close frequency slice HDF5 output: '//trim(diagnostic))
      end if
#ifdef CompileWithMPI
      call MPI_Barrier(this%publication%communicator, mpi_error)
      if (mpi_error /= MPI_SUCCESS) then
         call StopOnError(0, 0, 'Unable to synchronise closing frequency slice participants')
      end if
#endif
      if (this%publication%local_is_owner) then
         call verify_volumetric_visualisation(this%filesPath, error)
         if (file_exists(add_extension(this%filesPath, binaryExtension)) .and. &
             error == VISUALISATION_SUCCESS) then
            this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
            this%metadata%lifecycle%diagnostic = ''
         else
            this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_FAILED
            this%metadata%lifecycle%diagnostic = 'Required frequency slice output artifacts are incomplete'
         end if
         if (this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED) then
            call StopOnError(0, 0, trim(this%metadata%lifecycle%diagnostic))
         end if
      else
         this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
         this%metadata%lifecycle%diagnostic = ''
      end if
#ifdef CompileWithMPI
      call MPI_Barrier(this%publication%communicator, mpi_error)
      if (mpi_error /= MPI_SUCCESS) then
         call StopOnError(0, 0, 'Unable to complete frequency slice publication')
      end if
      if (this%publication%owns_communicator) then
         call MPI_Comm_free(this%publication%communicator, mpi_error)
         if (mpi_error /= MPI_SUCCESS) then
            call StopOnError(0, 0, 'Unable to free frequency slice output communicator')
         end if
         this%publication%owns_communicator = .false.
         this%publication%communicator = MPI_COMM_NULL
      end if
#endif
   end subroutine

end module frequencySliceProbeOutput_m
