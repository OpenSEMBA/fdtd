module movieProbeOutput_m
   use FDETYPES_m
   use utils_m
   use allocationUtils_m, only: alloc_and_init
   use report_m
   use outputTypes_m
   use outputUtils_m
   use volumicProbeUtils_m
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use outputCollective_m, only: OUTPUT_PUBLICATION_COLLECTIVE, OUTPUT_PUBLICATION_ROOT_AGGREGATION
   use outputTransport_m, only: output_transport_t, init_output_transport, transfer_flush_batch, &
                                OUTPUT_TRANSPORT_SUCCESS, OUTPUT_TRANSPORT_INVALID_CONTEXT
   use outputBinary_m, only: validate_binary_layout, append_binary_real64, BINARY_WRITER_SUCCESS
   use outputVisualisation_m, only: initialise_movie_visualisation, begin_visualisation_step, &
                                   write_visualisation_attribute, write_visualisation_attribute_hyperslab, end_visualisation_step, &
                                 flush_visualisation, close_visualisation, visualisation_is_open, verify_volumetric_visualisation, &
                                    VISUALISATION_SUCCESS, VISUALISATION_ATTRIBUTE_X, VISUALISATION_ATTRIBUTE_Y, &
                                    VISUALISATION_ATTRIBUTE_Z, VISUALISATION_ATTRIBUTE_TAG, VISUALISATION_ATTRIBUTE_MEDIA, &
                                    VISUALISATION_ATTRIBUTE_TAG_X, VISUALISATION_ATTRIBUTE_TAG_Y, &
                                    VISUALISATION_ATTRIBUTE_TAG_Z, VISUALISATION_ATTRIBUTE_MEDIA_X, &
                                    VISUALISATION_ATTRIBUTE_MEDIA_Y, VISUALISATION_ATTRIBUTE_MEDIA_Z
   use mapVTKOutput_m, only: write_geometry_companion
   use directoryUtils_m, only: add_extension, create_file_with_path, create_folder, file_exists, &
                               get_last_component, join_path
#ifdef CompileWithMPI
   use mpi
#endif
   implicit none(type, external)
   private

   !===========================
   ! Public interface
   !===========================
   public :: init_movie_probe_output
   public :: update_movie_probe_output
   public :: flush_movie_probe_output
   public :: close_movie_probe_output
   !===========================

   !===========================
   ! Private helpers
   !===========================
   ! Output & File Management
   private :: clear_memory_data

contains

   !===========================
   ! Public routines
   !===========================

   subroutine init_movie_probe_output(this, publication, field, domain, control, problemInfo, outputTypeExtension)
      type(movie_probe_output_t), intent(out) :: this
      type(volumetric_publication_t), intent(in) :: publication
      integer(kind=SINGLE), intent(in)        :: field
      type(domain_t), intent(in)              :: domain
      type(sim_control_t), intent(in)         :: control
      type(problem_info_t), intent(in)        :: problemInfo
      character(len=BUFSIZE), intent(in)      :: outputTypeExtension

      integer :: error
      type(output_transport_t) :: transport
      character(len=BUFSIZE) :: filename
      integer :: geometry_status
      character(len=BUFSIZE) :: geometry_diagnostic

      this%publication = publication
      this%mainCoords = publication%global_lower
      this%auxCoords = publication%global_upper
      if (publication%local_participates) then
         this%mainCoords%x = publication%global_lower%x + int(publication%file_offset(1), kind=SINGLE)
         this%mainCoords%y = publication%global_lower%y + int(publication%file_offset(2), kind=SINGLE)
         this%mainCoords%z = publication%global_lower%z + int(publication%file_offset(3), kind=SINGLE)
         this%auxCoords%x = this%mainCoords%x + int(publication%local_shape(1), kind=SINGLE) - 1_SINGLE
         this%auxCoords%y = this%mainCoords%y + int(publication%local_shape(2), kind=SINGLE) - 1_SINGLE
         this%auxCoords%z = this%mainCoords%z + int(publication%local_shape(3), kind=SINGLE) - 1_SINGLE
      end if
      this%component = field
      this%domain = domain

      this%path = get_output_path(publication%global_lower, publication%global_upper, &
                                  outputTypeExtension, field, control%mpidir)
      filename = get_last_component(this%path)
      this%filesPath = join_path(this%path, filename)
      call initialise_movie_metadata(this, error, control%mpidir)
      if (error /= 0) call StopOnError(control%layoutnumber, control%num_procs, &
                                       'Unable to initialise movie output metadata')
      this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_ACTIVE

      if (.not. publication%local_participates) return

      call find_and_store_important_coords(this%mainCoords, this%auxCoords, this%component, &
                                           problemInfo, this%nPoints, this%coords)
      allocate (this%tagNumber(3, this%nPoints), this%mediaType(3, this%nPoints))
      this%tagNumber = 0_IKINDMTAG
      this%mediaType = -1.0_RKIND
      call store_classification(this, problemInfo)
      call alloc_and_init(this%timeStep, OUTPUT_TIME_BUFFER_SIZE, 0.0_RKIND_tiempo)
      call alloc_and_init(this%xValueForTime, OUTPUT_TIME_BUFFER_SIZE, this%nPoints, 0.0_RKIND)
      call alloc_and_init(this%yValueForTime, OUTPUT_TIME_BUFFER_SIZE, this%nPoints, 0.0_RKIND)
      call alloc_and_init(this%zValueForTime, OUTPUT_TIME_BUFFER_SIZE, this%nPoints, 0.0_RKIND)

       if (publication%local_is_owner) then
          call create_folder(this%path, error)
          if (error /= 0) call StopOnError(control%layoutnumber, control%num_procs, &
                                           'Unable to create movie output directory')
          call create_bin_file(this%filesPath, error)
          if (error /= 0) call StopOnError(control%layoutnumber, control%num_procs, &
                                           'Unable to create movie binary output')
       end if
       call synchronise_movie_participants(this, control)
       call write_geometry_companion(this%filesPath, this%mainCoords, this%auxCoords, problemInfo, &
                                     this%publication%communicator, geometry_status, geometry_diagnostic)
       if (this%publication%local_is_owner .and. geometry_status /= VISUALISATION_SUCCESS) then
          call StopOnError(control%layoutnumber, control%num_procs, &
                           'Unable to create movie geometry: '//trim(geometry_diagnostic))
       end if

       if (publication%mode == OUTPUT_PUBLICATION_COLLECTIVE .or. publication%local_is_owner) then
         call create_movie_files(this, problemInfo, error)
         if (error /= 0) call StopOnError(control%layoutnumber, control%num_procs, &
                                          'Unable to initialise movie HDF5 output')
      end if
      call init_output_transport(transport, root_rank=0, status=error, &
                                 communicator=this%publication%communicator)
      if (error /= OUTPUT_TRANSPORT_SUCCESS) then
         call StopOnError(control%layoutnumber, control%num_procs, &
                          'Unable to initialise movie classification transport')
      end if
      call write_classification_attributes(this, transport)

   end subroutine init_movie_probe_output

   subroutine create_bin_file(filePath, error)
      character(len=*), intent(in) :: filePath
      integer, intent(out) :: error
      call create_file_with_path(add_extension(filePath, binaryExtension), error)
   end subroutine

   subroutine synchronise_movie_participants(this, control)
      type(movie_probe_output_t), intent(in) :: this
      type(sim_control_t), intent(in) :: control
#ifdef CompileWithMPI
      integer :: error

      if (this%publication%communicator_size > 1) then
         call MPI_Barrier(this%publication%communicator, error)
         if (error /= MPI_SUCCESS) call StopOnError(control%layoutnumber, control%num_procs, &
                                                    'Unable to synchronise movie output participants')
      end if
#endif
   end subroutine synchronise_movie_participants

   subroutine initialise_movie_metadata(this, error, mpidir)
      type(movie_probe_output_t), intent(inout) :: this
      integer, intent(out) :: error
      integer(kind=SINGLE), intent(in) :: mpidir
      character(len=BUFSIZE) :: base_name

      base_name = get_last_component(this%filesPath)
      this%metadata%probe_id = trim(base_name)
      this%metadata%quantity = get_prefix_extension(this%component, mpidir)
      this%metadata%lower_bound = this%publication%global_lower
      this%metadata%upper_bound = this%publication%global_upper
      this%metadata%domain_type = TIME_DOMAIN
      this%metadata%ownership%participant_ranks = this%publication%participant_ranks
      this%metadata%ownership%scalar_writer_rank = this%publication%owner_rank
      if (allocated(this%metadata%artifacts)) deallocate (this%metadata%artifacts)
      allocate (this%metadata%artifacts(5))
      this%metadata%artifacts(1)%kind = OUTPUT_ARTIFACT_BINARY
      this%metadata%artifacts(1)%relative_path = trim(base_name)//binaryExtension
      this%metadata%artifacts(1)%byte_order = BINARY_ENDIAN_LITTLE
      this%metadata%artifacts(1)%numeric_representation = BINARY_NUMERIC_REAL64
      this%metadata%artifacts(1)%record_bytes = 56
      this%metadata%artifacts(1)%component_order = 'time,x,y,z,Ex,Ey,Ez'
      this%metadata%artifacts(2)%kind = OUTPUT_ARTIFACT_VISUALISATION_METADATA
      this%metadata%artifacts(2)%relative_path = trim(base_name)//'.xdmf'
      this%metadata%artifacts(3)%kind = OUTPUT_ARTIFACT_VISUALISATION_DATA
      this%metadata%artifacts(3)%relative_path = trim(base_name)//'.h5'
      this%metadata%artifacts(4)%kind = OUTPUT_ARTIFACT_GEOMETRY
      this%metadata%artifacts(4)%relative_path = trim(base_name)//'_geometry.xdmf'
      this%metadata%artifacts(5)%kind = OUTPUT_ARTIFACT_VISUALISATION_DATA
      this%metadata%artifacts(5)%relative_path = trim(base_name)//'_geometry.h5'

      call validate_binary_layout(this%metadata%artifacts(1), error)
      if (error /= BINARY_WRITER_SUCCESS) return
      error = 0
   end subroutine initialise_movie_metadata

   subroutine create_movie_files(this, problemInfo, error)
      type(movie_probe_output_t), intent(inout) :: this
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(out) :: error

      character(len=BUFSIZE) :: attribute_names(14), diagnostic
      character(len=BUFSIZE) :: attributeBaseName
      logical :: attribute_enabled(14)
      integer :: status

      error = 0
      attribute_names = ''
      attribute_enabled = .false.
      select case (this%component)
      case (iCur, iMEC, iMHC)
         attributeBaseName = 'CurrenDensity'
         if (this%component == iMEC) attributeBaseName = 'ElectricField'
         if (this%component == iMHC) attributeBaseName = 'MagneticField'
         attribute_names(VISUALISATION_ATTRIBUTE_X) = trim(attributeBaseName)//'X'
         attribute_names(VISUALISATION_ATTRIBUTE_Y) = trim(attributeBaseName)//'Y'
         attribute_names(VISUALISATION_ATTRIBUTE_Z) = trim(attributeBaseName)//'Z'
         attribute_enabled(1:3) = .true.
         attribute_names(VISUALISATION_ATTRIBUTE_TAG_X) = 'tagnumber_x'
         attribute_names(VISUALISATION_ATTRIBUTE_TAG_Y) = 'tagnumber_y'
         attribute_names(VISUALISATION_ATTRIBUTE_TAG_Z) = 'tagnumber_z'
         attribute_names(VISUALISATION_ATTRIBUTE_MEDIA_X) = 'mediatype_x'
         attribute_names(VISUALISATION_ATTRIBUTE_MEDIA_Y) = 'mediatype_y'
         attribute_names(VISUALISATION_ATTRIBUTE_MEDIA_Z) = 'mediatype_z'
         attribute_enabled(VISUALISATION_ATTRIBUTE_TAG_X:VISUALISATION_ATTRIBUTE_MEDIA_Z) = .true.
      case (iCurX, iEXC, iHXC)
         attributeBaseName = 'CurrenDensity'
         if (this%component == iEXC) attributeBaseName = 'ElectricField'
         if (this%component == iHXC) attributeBaseName = 'MagneticField'
         attribute_names(VISUALISATION_ATTRIBUTE_X) = trim(attributeBaseName)//'X'
         attribute_enabled(VISUALISATION_ATTRIBUTE_X) = .true.
      case (iCurY, iEyC, iHyC)
         attributeBaseName = 'CurrenDensity'
         if (this%component == iEyC) attributeBaseName = 'ElectricField'
         if (this%component == iHyC) attributeBaseName = 'MagneticField'
         attribute_names(VISUALISATION_ATTRIBUTE_Y) = trim(attributeBaseName)//'Y'
         attribute_enabled(VISUALISATION_ATTRIBUTE_Y) = .true.
      case (iCurZ, iEZC, iHzC)
         attributeBaseName = 'CurrenDensity'
         if (this%component == iEZC) attributeBaseName = 'ElectricField'
         if (this%component == iHzC) attributeBaseName = 'MagneticField'
         attribute_names(VISUALISATION_ATTRIBUTE_Z) = trim(attributeBaseName)//'Z'
         attribute_enabled(VISUALISATION_ATTRIBUTE_Z) = .true.
      end select
      if (.not. any([iCur, iMEC, iMHC] == this%component)) then
         attribute_names(VISUALISATION_ATTRIBUTE_TAG) = 'tagnumber'
         attribute_names(VISUALISATION_ATTRIBUTE_MEDIA) = 'mediatype'
         attribute_enabled(VISUALISATION_ATTRIBUTE_TAG) = .true.
         attribute_enabled(VISUALISATION_ATTRIBUTE_MEDIA) = .true.
      end if
      call initialise_movie_visualisation(this%visualisation, trim(this%filesPath), &
                                real(problemInfo%xSteps(this%publication%global_lower%x:this%publication%global_upper%x), real64), &
                                real(problemInfo%ySteps(this%publication%global_lower%y:this%publication%global_upper%y), real64), &
                                real(problemInfo%zSteps(this%publication%global_lower%z:this%publication%global_upper%z), real64), &
                                       attribute_names, attribute_enabled, this%publication%mode == OUTPUT_PUBLICATION_COLLECTIVE, &
                                          this%publication%communicator, status, diagnostic)
      if (status /= VISUALISATION_SUCCESS) then
         error = 1
         print *, trim(diagnostic)
      end if
   end subroutine create_movie_files

   subroutine update_movie_probe_output(this, step, fieldsReference, control, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in)       :: step
      type(fields_reference_t), intent(in)      :: fieldsReference
      type(sim_control_t), intent(in)           :: control
      type(problem_info_t), intent(in)          :: problemInfo

      integer(kind=4) :: request
      if (.not. this%publication%local_participates .or. &
          this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_COMPLETE .or. &
          this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED) return

      request = this%component
      this%nTime = this%nTime + 1

      ! Determine which save routine to call
      if (any(VOLUMIC_M_MEASURE == request)) then
         select case (request)
         case (iCur)
            call save_current_module(this, fieldsReference, step, problemInfo)
         case (iMEC)
            call save_field_module(this, fieldsReference%E, request, step, problemInfo)
         case (iMHC)
            call save_field_module(this, fieldsReference%H, request, step, problemInfo)
         case default
            call StopOnError(control%layoutnumber, control%num_procs, "Volumic measure not supported")
         end select
      else if (any(VOLUMIC_X_MEASURE == request)) then
         select case (request)
         case (iCurX)
            call save_current_component(this, this%xValueForTime, fieldsReference, step, problemInfo, iEx)
         case (iExC)
            call save_field_component(this, this%xValueForTime, fieldsReference%E%x, step, problemInfo, iEx)
         case (iHxC)
            call save_field_component(this, this%xValueForTime, fieldsReference%H%x, step, problemInfo, iHx)
         case default
            call StopOnError(control%layoutnumber, control%num_procs, "Volumic measure not supported")
         end select
      else if (any(VOLUMIC_Y_MEASURE == request)) then
         select case (request)
         case (iCurY)
            call save_current_component(this, this%yValueForTime, fieldsReference, step, problemInfo, iEy)
         case (iEyC)
            call save_field_component(this, this%yValueForTime, fieldsReference%E%y, step, problemInfo, iEy)
         case (iHyC)
            call save_field_component(this, this%yValueForTime, fieldsReference%H%y, step, problemInfo, iHy)
         case default
            call StopOnError(control%layoutnumber, control%num_procs, "Volumic measure not supported")
         end select
      else if (any(VOLUMIC_Z_MEASURE == request)) then
         select case (request)
         case (iCurZ)
            call save_current_component(this, this%zValueForTime, fieldsReference, step, problemInfo, iEz)
         case (iEzC)
            call save_field_component(this, this%zValueForTime, fieldsReference%E%z, step, problemInfo, iEz)
         case (iHzC)
            call save_field_component(this, this%zValueForTime, fieldsReference%H%z, step, problemInfo, iHz)
         case default
            call StopOnError(control%layoutnumber, control%num_procs, "Volumic measure not supported")
         end select
      end if
   end subroutine update_movie_probe_output

   subroutine flush_movie_probe_output(this)
      type(movie_probe_output_t), intent(inout) :: this
      type(output_transport_t) :: transport
      integer :: error
      character(len=BUFSIZE) :: diagnostic

      if (.not. this%publication%local_participates .or. this%nTime == 0 .or. &
          this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_COMPLETE .or. &
          this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED) return
      call init_output_transport(transport, root_rank=0, status=error, &
                                 communicator=this%publication%communicator)
      if (error /= OUTPUT_TRANSPORT_SUCCESS) then
         call StopOnError(0, 0, 'Unable to initialise movie output transport')
      end if

      call write_bin_file(this, transport)
      call write_visualisation(this, transport)
      call flush_visualisation(this%visualisation, error, diagnostic)
      if (error /= VISUALISATION_SUCCESS) call StopOnError(0, 0, &
                                                           'Unable to flush movie HDF5 output: '//trim(diagnostic))
      call clear_memory_data(this)
   end subroutine flush_movie_probe_output

   subroutine close_movie_probe_output(this)
      type(movie_probe_output_t), intent(inout) :: this
      integer :: error
      logical :: lifecycle_is_terminal
      character(len=BUFSIZE) :: diagnostic

      lifecycle_is_terminal = this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_COMPLETE .or. &
                              this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED
      if (this%publication%local_participates .and. .not. lifecycle_is_terminal) then
         call flush_movie_probe_output(this)
      end if
      call close_visualisation(this%visualisation, error, diagnostic)
      if (error /= VISUALISATION_SUCCESS) then
         call StopOnError(0, 0, 'Unable to close movie HDF5 output: '//trim(diagnostic))
      end if
#ifdef CompileWithMPI
      if (this%publication%owns_communicator .and. this%publication%local_participates) then
         call MPI_Comm_free(this%publication%communicator, error)
         if (error /= MPI_SUCCESS) call StopOnError(0, 0, 'Unable to free movie output communicator')
         this%publication%owns_communicator = .false.
         this%publication%communicator = MPI_COMM_NULL
      end if
#endif

      if (lifecycle_is_terminal) return
      if (.not. this%publication%local_is_owner) then
         this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
         this%metadata%lifecycle%diagnostic = ''
         return
      end if

      call verify_volumetric_visualisation(this%filesPath, error)
      if (file_exists(add_extension(this%filesPath, binaryExtension)) .and. &
          file_exists(trim(this%filesPath)//'_geometry.xdmf') .and. &
          file_exists(trim(this%filesPath)//'_geometry.h5') .and. &
          error == VISUALISATION_SUCCESS) then
         this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
         this%metadata%lifecycle%diagnostic = ''
      else
         this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_FAILED
         this%metadata%lifecycle%diagnostic = 'Required movie output artifacts are incomplete'
      end if
      if (this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED) then
         call StopOnError(0, 0, trim(this%metadata%lifecycle%diagnostic))
      end if
   end subroutine close_movie_probe_output

   !===========================
   ! Private routines
   !===========================

   subroutine write_bin_file(this, transport)
      type(movie_probe_output_t), intent(inout) :: this
      type(output_transport_t), intent(in) :: transport
      integer :: counts_status, i, record_index, status, time_index
      integer, allocatable :: counts(:), displacements(:)
      real(real64), allocatable :: gathered_records(:), local_records(:)

      allocate (local_records(7*this%nPoints))
      do time_index = 1, this%nTime
         record_index = 0
         do i = 1, this%nPoints
            local_records(record_index + 1:record_index + 7) = [real(this%timeStep(time_index), real64), &
                                                                real(this%coords(1, i), real64), real(this%coords(2, i), real64), &
                                                 real(this%coords(3, i), real64), real(this%xValueForTime(time_index, i), real64), &
                                                                real(this%yValueForTime(time_index, i), real64), &
                                                                real(this%zValueForTime(time_index, i), real64)]
            record_index = record_index + 7
         end do
         call transfer_flush_batch(transport, local_records, gathered_records, counts, displacements, counts_status)
         if (counts_status /= OUTPUT_TRANSPORT_SUCCESS) then
            call StopOnError(0, 0, 'Unable to gather movie binary records')
         end if
         if (this%publication%local_is_owner) then
            call append_binary_real64(add_extension(this%filesPath, binaryExtension), &
                                      this%metadata%artifacts(1), gathered_records, status)
            if (status /= BINARY_WRITER_SUCCESS) then
               call StopOnError(0, 0, 'Unable to append movie binary records')
            end if
         end if
      end do
   end subroutine write_bin_file

   subroutine write_visualisation(this, transport)
      type(movie_probe_output_t), intent(inout) :: this
      type(output_transport_t), intent(in) :: transport

      character(len=BUFSIZE) :: diagnostic
      integer :: status, time_index

      call write_classification_attributes(this, transport)
      do time_index = 1, this%nTime
         if (visualisation_is_open(this%visualisation)) then
            call begin_visualisation_step(this%visualisation, real(this%timeStep(time_index), real64), &
                                          status, diagnostic)
            call require_visualisation_success(status, diagnostic, 'Movie HDF5 write failed: ')
         end if
         if (any([iCur, iMEC, iMHC, iCurX, iExC, iHxC] == this%component)) then
            call write_external_attribute(this, VISUALISATION_ATTRIBUTE_X, this%xValueForTime, &
                                          time_index, transport)
         end if
         if (any([iCur, iMEC, iMHC, iCurY, iEyC, iHyC] == this%component)) then
            call write_external_attribute(this, VISUALISATION_ATTRIBUTE_Y, this%yValueForTime, &
                                          time_index, transport)
         end if
         if (any([iCur, iMEC, iMHC, iCurZ, iEzC, iHzC] == this%component)) then
            call write_external_attribute(this, VISUALISATION_ATTRIBUTE_Z, this%zValueForTime, &
                                          time_index, transport)
         end if
         if (visualisation_is_open(this%visualisation)) then
            call end_visualisation_step(this%visualisation, status, diagnostic)
            call require_visualisation_success(status, diagnostic, 'Movie HDF5 write failed: ')
         end if
      end do
   end subroutine write_visualisation

   subroutine write_external_attribute(this, attribute, values, time_index, transport)
      type(movie_probe_output_t), intent(inout) :: this
      integer, intent(in) :: attribute
      real(RKIND), intent(in) :: values(:, :)
      integer, intent(in) :: time_index
      type(output_transport_t), intent(in) :: transport

      integer :: i, local_index, nx, ny
      real(real64), allocatable :: field(:)

      nx = int(this%publication%local_shape(1))
      ny = int(this%publication%local_shape(2))
      allocate (field(int(product(this%publication%local_shape))))
      field = 0.0_real64
      do i = 1, this%nPoints
         local_index = this%coords(1, i) - this%mainCoords%x + 1 + &
                       (this%coords(2, i) - this%mainCoords%y)*nx + &
                       (this%coords(3, i) - this%mainCoords%z)*nx*ny
         field(local_index) = real(values(time_index, i), real64)
      end do
      call publish_dense_attribute(this, attribute, field, transport)
   end subroutine write_external_attribute

   subroutine write_classification_attributes(this, transport)
      type(movie_probe_output_t), intent(inout) :: this
      type(output_transport_t), intent(in) :: transport
      integer, parameter :: tag_attributes(3) = [VISUALISATION_ATTRIBUTE_TAG_X, &
                                                 VISUALISATION_ATTRIBUTE_TAG_Y, &
                                                 VISUALISATION_ATTRIBUTE_TAG_Z]
      integer, parameter :: media_attributes(3) = [VISUALISATION_ATTRIBUTE_MEDIA_X, &
                                                   VISUALISATION_ATTRIBUTE_MEDIA_Y, &
                                                   VISUALISATION_ATTRIBUTE_MEDIA_Z]
      integer :: axis

      if (this%classificationWritten) return
      if (any([iCur, iMEC, iMHC] == this%component)) then
         do axis = 1, 3
            call write_external_classification_attribute(this, tag_attributes(axis), axis, .true., transport)
            call write_external_classification_attribute(this, media_attributes(axis), axis, .false., transport)
         end do
      else
         do axis = 1, 3
            if (.not. classification_axis_enabled(this%component, axis)) cycle
            call write_external_classification_attribute(this, VISUALISATION_ATTRIBUTE_TAG, axis, .true., transport)
            call write_external_classification_attribute(this, VISUALISATION_ATTRIBUTE_MEDIA, axis, .false., transport)
         end do
      end if
      this%classificationWritten = .true.
   end subroutine write_classification_attributes

   subroutine write_external_classification_attribute(this, attribute, axis, write_tag, transport)
      type(movie_probe_output_t), intent(inout) :: this
      integer, intent(in) :: attribute, axis
      logical, intent(in) :: write_tag
      type(output_transport_t), intent(in) :: transport

      integer :: i, local_index, nx, ny
      real(real64), allocatable :: values(:)

      nx = int(this%publication%local_shape(1))
      ny = int(this%publication%local_shape(2))
      allocate (values(int(product(this%publication%local_shape))))
      values = 0.0_real64
      do i = 1, this%nPoints
         local_index = this%coords(1, i) - this%mainCoords%x + 1 + &
                       (this%coords(2, i) - this%mainCoords%y)*nx + &
                       (this%coords(3, i) - this%mainCoords%z)*nx*ny
         if (write_tag) then
            values(local_index) = real(this%tagNumber(axis, i), real64)
         else
            values(local_index) = real(this%mediaType(axis, i), real64)
         end if
      end do
      call publish_dense_attribute(this, attribute, values, transport)
   end subroutine write_external_classification_attribute

   subroutine publish_dense_attribute(this, attribute, local_values, transport)
      type(movie_probe_output_t), intent(inout) :: this
      integer, intent(in) :: attribute
      real(real64), intent(in) :: local_values(:)
      type(output_transport_t), intent(in) :: transport
      character(len=BUFSIZE) :: diagnostic
      real(real64), allocatable :: global_values(:)
      integer :: status, transport_status

      select case (this%publication%mode)
      case (OUTPUT_PUBLICATION_COLLECTIVE)
         call write_visualisation_attribute_hyperslab(this%visualisation, attribute, local_values, &
                                                     this%publication%file_offset, this%publication%local_shape, status, diagnostic)
         call require_visualisation_success(status, diagnostic, 'Movie HDF5 hyperslab write failed: ')
      case (OUTPUT_PUBLICATION_ROOT_AGGREGATION)
         call gather_global_dense(this, local_values, transport, global_values, transport_status)
         if (transport_status /= OUTPUT_TRANSPORT_SUCCESS) then
            call StopOnError(0, 0, 'Unable to aggregate movie HDF5 values')
         end if
         if (this%publication%local_is_owner) then
            call write_visualisation_attribute(this%visualisation, attribute, global_values, status, diagnostic)
            call require_visualisation_success(status, diagnostic, 'Movie HDF5 write failed: ')
         end if
      case default
         call StopOnError(0, 0, 'Unsupported movie output publication mode')
      end select
   end subroutine publish_dense_attribute

   subroutine require_visualisation_success(status, diagnostic, context)
      integer, intent(in) :: status
      character(len=*), intent(in) :: diagnostic, context

      if (status /= VISUALISATION_SUCCESS) then
         call StopOnError(0, 0, context//trim(diagnostic))
      end if
   end subroutine require_visualisation_success

   subroutine gather_global_dense(this, local_values, transport, global_values, status)
      type(movie_probe_output_t), intent(in) :: this
      real(real64), intent(in) :: local_values(:)
      type(output_transport_t), intent(in) :: transport
      real(real64), allocatable, intent(out) :: global_values(:)
      integer, intent(out) :: status

      integer, allocatable :: counts(:), displacements(:)
      real(real64), allocatable :: gathered_values(:), local_batch(:)
      integer(int64) :: global_shape(3), offset(3), shape(3)
      integer :: global_index, i, j, k, local_index, rank_index, value_start

      allocate (local_batch(6 + size(local_values)))
      local_batch(1:3) = real(this%publication%file_offset, real64)
      local_batch(4:6) = real(this%publication%local_shape, real64)
      local_batch(7:) = local_values
      call transfer_flush_batch(transport, local_batch, gathered_values, counts, displacements, status)
      if (status /= OUTPUT_TRANSPORT_SUCCESS) return

      if (transport%rank /= transport%root_rank) then
         allocate (global_values(0))
         return
      end if

      global_shape = [ &
                     int(this%publication%global_upper%x, int64) - int(this%publication%global_lower%x, int64) + 1_int64, &
                     int(this%publication%global_upper%y, int64) - int(this%publication%global_lower%y, int64) + 1_int64, &
                     int(this%publication%global_upper%z, int64) - int(this%publication%global_lower%z, int64) + 1_int64]
      allocate (global_values(int(product(global_shape))))
      global_values = 0.0_real64
      do rank_index = 1, transport%rank_count
         if (counts(rank_index) < 6) then
            status = OUTPUT_TRANSPORT_INVALID_CONTEXT
            return
         end if
         value_start = displacements(rank_index)
         offset = nint(gathered_values(value_start + 1:value_start + 3), kind=int64)
         shape = nint(gathered_values(value_start + 4:value_start + 6), kind=int64)
         if (any(shape <= 0_int64) .or. any(offset < 0_int64) .or. &
             any(offset + shape > global_shape) .or. &
             int(counts(rank_index) - 6, int64) /= product(shape)) then
            status = OUTPUT_TRANSPORT_INVALID_CONTEXT
            return
         end if
         do k = 1, int(shape(3))
         do j = 1, int(shape(2))
         do i = 1, int(shape(1))
            local_index = i + (j - 1)*int(shape(1)) + &
                          (k - 1)*int(shape(1)*shape(2))
            global_index = int(offset(1)) + i + &
                           (int(offset(2)) + j - 1)*int(global_shape(1)) + &
                           (int(offset(3)) + k - 1)*int(global_shape(1)*global_shape(2))
            global_values(global_index) = gathered_values(value_start + 6 + local_index)
         end do
         end do
         end do
      end do
   end subroutine gather_global_dense

   subroutine store_classification(this, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      type(problem_info_t), intent(in) :: problemInfo

      integer :: axis, i, field

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

   function get_output_path(global_lower, global_upper, outputTypeExtension, field, mpidir) result(path)
      type(cell_coordinate_t), intent(in) :: global_lower, global_upper
      character(len=*), intent(in) :: outputTypeExtension
      integer(kind=SINGLE), intent(in) :: field, mpidir
      character(len=BUFSIZE) :: path, probeBoundsExtension, prefixFieldExtension

      probeBoundsExtension = get_coordinates_extension(global_lower, global_upper, mpidir)
      prefixFieldExtension = get_prefix_extension(field, mpidir)
      path = trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
   end function get_output_path

   subroutine save_current_module(this, fieldsReference, simTime, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      type(fields_reference_t), intent(in)      :: fieldsReference
      real(kind=RKIND_tiempo), intent(in)       :: simTime
      type(problem_info_t), intent(in)          :: problemInfo

      integer :: i, j, k, coordIdx
      this%timeStep(this%nTime) = simTime
      coordIdx = 0
      do k = this%mainCoords%z, this%auxCoords%z
      do j = this%mainCoords%y, this%auxCoords%y
      do i = this%mainCoords%x, this%auxCoords%x
         if (isValidPointForCurrent(iCur, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_current(this%xValueForTime, this%nTime, coordIdx, iEx, i, j, k, fieldsReference)
            call save_current(this%yValueForTime, this%nTime, coordIdx, iEy, i, j, k, fieldsReference)
            call save_current(this%zValueForTime, this%nTime, coordIdx, iEz, i, j, k, fieldsReference)
         end if
      end do
      end do
      end do
   end subroutine save_current_module

   subroutine save_current_component(this, currentData, fieldsReference, simTime, problemInfo, fieldDir)
      type(movie_probe_output_t), intent(inout) :: this
      real(kind=RKIND), intent(inout)           :: currentData(:, :)
      type(fields_reference_t), intent(in)      :: fieldsReference
      real(kind=RKIND_tiempo), intent(in)       :: simTime
      type(problem_info_t), intent(in)          :: problemInfo
      integer, intent(in)                       :: fieldDir

      integer :: i, j, k, coordIdx
      this%timeStep(this%nTime) = simTime
      coordIdx = 0
      do k = this%mainCoords%z, this%auxCoords%z
      do j = this%mainCoords%y, this%auxCoords%y
      do i = this%mainCoords%x, this%auxCoords%x
         if (isValidPointForCurrent(fieldDir, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_current(currentData, this%nTime, coordIdx, fieldDir, i, j, k, fieldsReference)
         end if
      end do
      end do
      end do
   end subroutine save_current_component

   subroutine save_current(currentData, timeIdx, coordIdx, field, i, j, k, fieldsReference)
      real(kind=RKIND), intent(inout)        :: currentData(:, :)
      integer(kind=SINGLE), intent(in)       :: timeIdx, coordIdx, field, i, j, k
      type(fields_reference_t), intent(in)   :: fieldsReference

      currentData(timeIdx, coordIdx) = computeJ(field, i, j, k, fieldsReference)
   end subroutine save_current

   subroutine save_field_module(this, field, request, simTime, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      type(field_data_t), intent(in)            :: field
      real(kind=RKIND_tiempo), intent(in)       :: simTime
      type(problem_info_t), intent(in)          :: problemInfo
      integer, intent(in)                       :: request

      integer :: i, j, k, coordIdx
      this%timeStep(this%nTime) = simTime
      coordIdx = 0
      do k = this%mainCoords%z, this%auxCoords%z
      do j = this%mainCoords%y, this%auxCoords%y
      do i = this%mainCoords%x, this%auxCoords%x
         if (isValidPointForField(request, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_field(this%xValueForTime, this%nTime, coordIdx, field%x(i, j, k))
            call save_field(this%yValueForTime, this%nTime, coordIdx, field%y(i, j, k))
            call save_field(this%zValueForTime, this%nTime, coordIdx, field%z(i, j, k))
         end if
      end do
      end do
      end do
   end subroutine save_field_module

   subroutine save_field_component(this, fieldData, fieldComponent, simTime, problemInfo, fieldDir)
      type(movie_probe_output_t), intent(inout) :: this
      real(kind=RKIND), intent(inout)           :: fieldData(:, :)
      real(kind=RKIND), intent(in)              :: fieldComponent(:, :, :)
      real(kind=RKIND_tiempo), intent(in)       :: simTime
      type(problem_info_t), intent(in)          :: problemInfo
      integer, intent(in)                       :: fieldDir

      integer :: i, j, k, coordIdx
      this%timeStep(this%nTime) = simTime
      coordIdx = 0
      do k = this%mainCoords%z, this%auxCoords%z
      do j = this%mainCoords%y, this%auxCoords%y
      do i = this%mainCoords%x, this%auxCoords%x
         if (isValidPointForField(fieldDir, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_field(fieldData, this%nTime, coordIdx, fieldComponent(i, j, k))
         end if
      end do
      end do
      end do
   end subroutine save_field_component

   subroutine save_field(fieldData, timeIdx, coordIdx, fieldValue)
      real(kind=RKIND), intent(inout)  :: fieldData(:, :)
      integer(kind=SINGLE), intent(in) :: timeIdx, coordIdx
      real(kind=RKIND), intent(in)     :: fieldValue

      fieldData(timeIdx, coordIdx) = fieldValue
   end subroutine save_field

   subroutine clear_memory_data(this)
      type(movie_probe_output_t), intent(inout) :: this
      this%nTimesFlushed = this%nTimesFlushed + this%nTime
      this%nTime = 0
      this%timeStep = 0.0_RKIND
      if (any(VOLUMIC_M_MEASURE == this%component)) then
         this%xValueForTime = 0.0_RKIND
         this%yValueForTime = 0.0_RKIND
         this%zValueForTime = 0.0_RKIND
      else if (any(VOLUMIC_X_MEASURE == this%component)) then
         this%xValueForTime = 0.0_RKIND
      else if (any(VOLUMIC_Y_MEASURE == this%component)) then
         this%yValueForTime = 0.0_RKIND
      else if (any(VOLUMIC_Z_MEASURE == this%component)) then
         this%zValueForTime = 0.0_RKIND
      end if
   end subroutine clear_memory_data

end module movieProbeOutput_m
