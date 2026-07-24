module movieProbeOutput_m
   use FDETYPES_m
   use utils_m
   use allocationUtils_m, only: alloc_and_init
   use report_m
   use outputTypes_m
   use outputUtils_m
   use allocationUtils_m, only: alloc_and_init
   use volumicProbeUtils_m
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use xdmf_hdf5_m, only: xdmf_options_t, xdmf_status_t, &
      xdmf_attribute_id_t, XDMF_SERIES_TIME, XDMF_CENTER_NODE, &
      XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64
   use outputBinary_m, only: validate_binary_layout, open_binary_append, BINARY_WRITER_SUCCESS
   use outputMetadata_m, only: publish_initial_probe_metadata, publish_final_probe_metadata, OUTPUT_METADATA_SUCCESS
   use outputVisualisation_m, only: verify_volumetric_visualisation, VISUALISATION_SUCCESS
   use directoryUtils_m, only: add_extension, create_file_with_path, create_folder, file_exists, &
                               get_last_component, join_path
   implicit none
   private

   !===========================
   ! Public interface
   !===========================
   public :: init_movie_probe_output
   public :: update_movie_probe_output
   public :: flush_movie_probe_output
   public :: close_movie_probe_output
   public :: configure_movie_probe_publication
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

    subroutine init_movie_probe_output(this, lowerBound, upperBound, field, domain, control, problemInfo, outputTypeExtension)
      type(movie_probe_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in)     :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in)        :: field
      type(domain_t), intent(in)              :: domain
      type(sim_control_t), intent(in)         :: control
      type(problem_info_t), intent(in)        :: problemInfo
      character(len=BUFSIZE), intent(in)      :: outputTypeExtension

      integer :: error
      character(len=BUFSIZE) :: filename
      real(RKIND), pointer :: xsteps(:), ysteps(:), zsteps(:)

      this%mainCoords = lowerBound
      this%auxCoords = upperBound
      this%component = field
      this%domain = domain

      xsteps => problemInfo%xSteps(lowerBound%x:upperBound%x)
      ysteps => problemInfo%ySteps(lowerBound%y:upperBound%y)
      zsteps => problemInfo%zSteps(lowerBound%z:upperBound%z)

      call find_and_store_important_coords(this%mainCoords, this%auxCoords, this%component, problemInfo, this%nPoints, this%coords)
      call alloc_and_init(this%timeStep, BUFSIZE, 0.0_RKIND_tiempo)

      ! Allocate value arrays based on component type
      call alloc_and_init(this%xValueForTime, BUFSIZE, this%nPoints, 0.0_RKIND)
      call alloc_and_init(this%yValueForTime, BUFSIZE, this%nPoints, 0.0_RKIND)
      call alloc_and_init(this%zValueForTime, BUFSIZE, this%nPoints, 0.0_RKIND)

      this%path = get_output_path(this, outputTypeExtension, field, control%mpidir)
      filename = get_last_component(this%path)
      this%filesPath = join_path(this%path, filename)

       call create_folder(this%path, error)
       call create_bin_file(this%filesPath, error)
       call create_movie_files(this, error, xsteps, ysteps, zsteps)
       call initialise_movie_metadata(this, error, control%mpidir)
       if (error /= 0) print *, 'error en creacion'
   end subroutine init_movie_probe_output

   subroutine configure_movie_probe_publication(this, publication_mode, local_participates)
      type(movie_probe_output_t), intent(inout) :: this
      integer, intent(in) :: publication_mode
      logical, intent(in) :: local_participates

      this%publication_mode = publication_mode
      this%local_participates = local_participates
   end subroutine configure_movie_probe_publication

     subroutine create_bin_file(filePath, error)
      character(len=*), intent(in) :: filePath
      integer, intent(out) :: error
      call create_file_with_path(add_extension(filePath, binaryExtension), error)
    end subroutine

   subroutine initialise_movie_metadata(this, error, mpidir)
       type(movie_probe_output_t), intent(inout) :: this
       integer, intent(out) :: error
       integer(kind=SINGLE), intent(in) :: mpidir
       character(len=BUFSIZE) :: base_name

       base_name = get_last_component(this%filesPath)
       this%metadata%probe_id = trim(base_name)
        this%metadata%quantity = get_prefix_extension(this%component, mpidir)
       this%metadata%lower_bound = this%mainCoords
       this%metadata%upper_bound = this%auxCoords
       this%metadata%domain_type = TIME_DOMAIN
       if (allocated(this%metadata%artifacts)) deallocate(this%metadata%artifacts)
       allocate(this%metadata%artifacts(3))
       this%metadata%artifacts(1)%kind = OUTPUT_ARTIFACT_BINARY
       this%metadata%artifacts(1)%relative_path = trim(base_name)//binaryExtension
       this%metadata%artifacts(1)%byte_order = BINARY_ENDIAN_LITTLE
       this%metadata%artifacts(1)%numeric_representation = BINARY_NUMERIC_REAL32
       this%metadata%artifacts(1)%record_bytes = 32
       this%metadata%artifacts(2)%kind = OUTPUT_ARTIFACT_VISUALISATION_METADATA
       this%metadata%artifacts(2)%relative_path = trim(base_name)//'.xdmf'
       this%metadata%artifacts(3)%kind = OUTPUT_ARTIFACT_VISUALISATION_DATA
       this%metadata%artifacts(3)%relative_path = trim(base_name)//'.h5'

       call validate_binary_layout(this%metadata%artifacts(1), error)
       if (error /= BINARY_WRITER_SUCCESS) return
       call publish_initial_probe_metadata(add_extension(this%filesPath, '.json'), this%metadata, error)
       if (error /= OUTPUT_METADATA_SUCCESS) return
       this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_ACTIVE
       error = 0
   end subroutine initialise_movie_metadata

   subroutine create_movie_files(this, error, xsteps, ysteps, zsteps)
      type(movie_probe_output_t), intent(inout) :: this
      real(RKIND), pointer, intent(in) :: xsteps(:), ysteps(:), zsteps(:)
      integer, intent(out) :: error

      type(xdmf_options_t) :: options
      type(xdmf_status_t) :: status
      character(len=BUFSIZE) :: attributeBaseName

      error = 0
      allocate(this%writer)
      options%overwrite = .true.
      options%series_kind = XDMF_SERIES_TIME
      call this%writer%create(trim(this%filesPath), options, status)
      if (status%is_error()) then
         error = 1
         print *, trim(status%message())
         return
      end if

      call this%writer%define_rectilinear_grid('movieProbe', &
         real(xsteps, real64), real(ysteps, real64), real(zsteps, real64), &
         this%grid, status)
      if (status%is_error()) then
         error = 1
         print *, trim(status%message())
         return
      end if

      select case(this%component)
      case(iCur, iMEC, iMHC)
         attributeBaseName = 'CurrenDensity'
         if (this%component == iMEC) attributeBaseName = 'ElectricField'
         if (this%component == iMHC) attributeBaseName = 'MagneticField'
         call define_movie_attribute(this, trim(attributeBaseName)//'X', &
            this%xAttribute, status)
         if (.not. status%is_error()) call define_movie_attribute(this, &
            trim(attributeBaseName)//'Y', this%yAttribute, status)
         if (.not. status%is_error()) call define_movie_attribute(this, &
            trim(attributeBaseName)//'Z', this%zAttribute, status)
      case(iCurX, iEXC, iHXC)
         attributeBaseName = 'CurrenDensity'
         if (this%component == iEXC) attributeBaseName = 'ElectricField'
         if (this%component == iHXC) attributeBaseName = 'MagneticField'
         call define_movie_attribute(this, trim(attributeBaseName)//'X', &
            this%xAttribute, status)
      case(iCurY, iEyC, iHyC)
         attributeBaseName = 'CurrenDensity'
         if (this%component == iEyC) attributeBaseName = 'ElectricField'
         if (this%component == iHyC) attributeBaseName = 'MagneticField'
         call define_movie_attribute(this, trim(attributeBaseName)//'Y', &
            this%yAttribute, status)
      case(iCurZ, iEZC, iHzC)
         attributeBaseName = 'CurrenDensity'
         if (this%component == iEZC) attributeBaseName = 'ElectricField'
         if (this%component == iHzC) attributeBaseName = 'MagneticField'
         call define_movie_attribute(this, trim(attributeBaseName)//'Z', &
            this%zAttribute, status)
      end select
      if (status%is_error()) then
         error = 1
         print *, trim(status%message())
      end if
   end subroutine create_movie_files

   subroutine define_movie_attribute(this, name, attribute, status)
      type(movie_probe_output_t), intent(inout) :: this
      character(len=*), intent(in) :: name
      type(xdmf_attribute_id_t), intent(out) :: attribute
      type(xdmf_status_t), intent(out) :: status

      call this%writer%define_attribute(this%grid, name, XDMF_CENTER_NODE, &
         XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., attribute, status)
   end subroutine define_movie_attribute

   subroutine update_movie_probe_output(this, step, fieldsReference, control, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in)       :: step
      type(fields_reference_t), intent(in)      :: fieldsReference
      type(sim_control_t), intent(in)           :: control
      type(problem_info_t), intent(in)          :: problemInfo

      integer(kind=4) :: request
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
      if (this%nTime /= 0) then
          call write_bin_file(this)
          call write_to_external_xdmf(this)
      end if
       call clear_memory_data(this)
    end subroutine flush_movie_probe_output

   subroutine close_movie_probe_output(this)
      type(movie_probe_output_t), intent(inout) :: this
      type(xdmf_status_t) :: writer_status
      integer :: error
      logical :: writer_ok

      if (this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_COMPLETE .or. &
          this%metadata%lifecycle%state == OUTPUT_LIFECYCLE_FAILED) return
      writer_ok = .true.
      if (associated(this%writer)) then
          call this%writer%close(writer_status)
          if (writer_status%is_error()) then
             writer_ok = .false.
             this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_FAILED
             this%metadata%lifecycle%diagnostic = 'Unable to close movie visualisation'
          end if
          deallocate(this%writer)
       end if
       call verify_volumetric_visualisation(this%filesPath, error)
       if (file_exists(add_extension(this%filesPath, binaryExtension)) .and. &
           error == VISUALISATION_SUCCESS .and. writer_ok) then
          this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_COMPLETE
          this%metadata%lifecycle%diagnostic = ''
       else
          this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_FAILED
          this%metadata%lifecycle%diagnostic = 'Required movie artifacts are incomplete'
       end if
       call publish_final_probe_metadata(add_extension(this%filesPath, '.json'), this%metadata, error)
       if (error /= OUTPUT_METADATA_SUCCESS) then
          this%metadata%lifecycle%state = OUTPUT_LIFECYCLE_FAILED
          this%metadata%lifecycle%diagnostic = 'Unable to publish movie metadata'
       end if
   end subroutine close_movie_probe_output

   !===========================
   ! Private routines
   !===========================

    subroutine write_bin_file(this)
       ! Check type definition for binary format
       type(movie_probe_output_t), intent(inout) :: this
       integer :: i, status, t, unit

       call open_binary_append(add_extension(this%filesPath, binaryExtension), this%metadata%artifacts(1), unit, status)
       if (status /= BINARY_WRITER_SUCCESS) return
       do t = 1, this%nTime
      do i = 1, this%nPoints
         write(unit) this%timeStep(t), this%coords(1,i), this%coords(2,i), this%coords(3,i), this%xValueForTime(t,i), this%yValueForTime(t,i), this%zValueForTime(t,i)
      end do
      end do
      flush (unit)
      close (unit)
   end subroutine

   subroutine write_to_external_xdmf(this)
      type(movie_probe_output_t), intent(inout) :: this

      type(xdmf_status_t) :: status
      integer :: time_index

      do time_index = 1, this%nTime
         call this%writer%begin_step(real(this%timeStep(time_index), real64), status)
         if (status%is_error()) then
            print *, trim(status%message())
            return
         end if
         if (any([iCur, iMEC, iMHC, iCurX, iExC, iHxC] == this%component)) then
            call write_external_attribute(this, this%xAttribute, &
               this%xValueForTime, time_index, status)
         end if
         if (.not. status%is_error() .and. &
             any([iCur, iMEC, iMHC, iCurY, iEyC, iHyC] == this%component)) then
            call write_external_attribute(this, this%yAttribute, &
               this%yValueForTime, time_index, status)
         end if
         if (.not. status%is_error() .and. &
             any([iCur, iMEC, iMHC, iCurZ, iEzC, iHzC] == this%component)) then
            call write_external_attribute(this, this%zAttribute, &
               this%zValueForTime, time_index, status)
         end if
         if (.not. status%is_error()) call this%writer%end_step(status)
         if (status%is_error()) then
            print *, trim(status%message())
            return
         end if
      end do
   end subroutine write_to_external_xdmf

   subroutine write_external_attribute(this, attribute, values, time_index, status)
      type(movie_probe_output_t), intent(inout) :: this
      type(xdmf_attribute_id_t), intent(in) :: attribute
      real(RKIND), intent(in) :: values(:, :)
      integer, intent(in) :: time_index
      type(xdmf_status_t), intent(out) :: status

      integer :: i, nx, ny, nz
      real(real64), allocatable :: field(:, :, :)

      nx = this%auxCoords%x - this%mainCoords%x + 1
      ny = this%auxCoords%y - this%mainCoords%y + 1
      nz = this%auxCoords%z - this%mainCoords%z + 1
      allocate(field(nx, ny, nz))
      field = 0.0_real64
      do i = 1, this%nPoints
         field(this%coords(1, i) - this%mainCoords%x + 1, &
               this%coords(2, i) - this%mainCoords%y + 1, &
               this%coords(3, i) - this%mainCoords%z + 1) = &
            real(values(time_index, i), real64)
      end do
      call this%writer%write_attribute(attribute, reshape(field, [size(field)]), status)
      deallocate(field)
   end subroutine write_external_attribute

   function get_output_path(this, outputTypeExtension, field, mpidir) result(path)
      type(movie_probe_output_t), intent(in) :: this
      character(len=*), intent(in)           :: outputTypeExtension
      integer(kind=SINGLE), intent(in)       :: field, mpidir
      character(len=BUFSIZE)                 :: path, probeBoundsExtension, prefixFieldExtension

      probeBoundsExtension = get_coordinates_extension(this%mainCoords, this%auxCoords, mpidir)
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
