module movieProbeOutput_m
   use FDETYPES_m
   use utils_m
   use report_m
   use outputTypes_m
   use outputUtils_m
   use volumicProbeUtils_m
   use HDF5
   use xdmfAPI_m
   implicit none
   private

   !===========================
   ! Public interface
   !===========================
   public :: init_movie_probe_output
   public :: update_movie_probe_output
   public :: flush_movie_probe_output
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
      if (error/=0) print *, 'error en creacion'
   end subroutine init_movie_probe_output

   subroutine create_bin_file(filePath, error)
      character(len=*), intent(in) :: filePath
      integer, intent(out) :: error
      call create_file_with_path(add_extension(filePath, binaryExtension), error)
   end subroutine 

   subroutine create_movie_files(this, error, xsteps, ysteps, zsteps)
      type(movie_probe_output_t), intent(in) :: this
      real(RKIND), pointer, intent(in) :: xsteps(:), ysteps(:), zsteps(:)
      integer, intent(out) :: error

      real(dp), allocatable, dimension(:, :) :: coordsReal
      integer(HID_T) :: file_id
      character(len=BUFSIZE) :: h5_filename
      character(len=BUFSIZE) :: attributeBaseName
      integer(SINGLE), dimension(3) :: topology_size

      h5_filename = add_extension(this%filesPath, ".h5")
      topology_size(1) = this%auxCoords%x - this%mainCoords%x + 1
      topology_size(2) = this%auxCoords%y - this%mainCoords%y + 1
      topology_size(3) = this%auxCoords%z - this%mainCoords%z + 1

      call H5open_f(error)
      call create_h5_file(trim(h5_filename), file_id)

      call h5_create_rectilinear_coords_dataset(file_id, real(xsteps,dp), real(ysteps,dp), real(zsteps,dp))
      call h5_create_times_dataset(file_id, BUFSIZE)
      call create_h5_data_dataset(file_id, this%component, topology_size)

      call H5Fclose_f(file_id, error)
      call H5close_f(error)
   end subroutine 

   subroutine create_h5_data_dataset(file_id, requestedComponent, topology_size)
      integer(HID_T), intent(in) :: file_id
      integer(SINGLE), intent(in) :: requestedComponent
      integer(SINGLE), dimension(3), intent(in) :: topology_size

      character(len=BUFSIZE) :: attributeBaseName

      select case(requestedComponent)
      case(iCur, iCurX, iCurY, iCurZ); attributeBaseName = 'CurrenDensity'
      case(iMEC, iExC, iEyC, iEzC); attributeBaseName = 'ElectricField'
      case(iMHC, iHxC, iHyC, iHzC); attributeBaseName = 'MagneticField'
      end select

      select case(requestedComponent)
      case(iCur, iMEC, iMHC)
         call h5_init_extendable_dataset(file_id, trim(attributeBaseName)//'X', topology_size,  BUFSIZE)  
         call h5_init_extendable_dataset(file_id, trim(attributeBaseName)//'Y', topology_size,  BUFSIZE)  
         call h5_init_extendable_dataset(file_id, trim(attributeBaseName)//'Z', topology_size,  BUFSIZE)  
      case(iCurX, iEXC, iHXC)
         call h5_init_extendable_dataset(file_id, trim(attributeBaseName)//'X', topology_size, BUFSIZE)  
      case(iCurY, iEyC, iHyC) 
         call h5_init_extendable_dataset(file_id, trim(attributeBaseName)//'Y', topology_size, BUFSIZE)  
      case(iCurZ, iEZC, iHzC) 
         call h5_init_extendable_dataset(file_id, trim(attributeBaseName)//'Z', topology_size, BUFSIZE)  
      end select
   end subroutine

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
      integer :: i
      if (this%nTime /= 0) then
         call write_bin_file(this)
         call write_to_h5_file(this)
         call write_to_xdmf_file(this)
      end if
      call clear_memory_data(this)
   end subroutine flush_movie_probe_output

   !===========================
   ! Private routines
   !===========================

   subroutine write_bin_file(this)
      ! Check type definition for binary format
      type(movie_probe_output_t), intent(inout) :: this
      integer :: i, t, unit

      open (unit=unit, file=add_extension(this%filesPath, binaryExtension), &
            status='old', form='unformatted', position='append', access='stream')
      do t = 1, this%nTime
      do i = 1, this%nPoints
         write(unit) this%timeStep(t), this%coords(1,i), this%coords(2,i), this%coords(3,i), this%xValueForTime(t,i), this%yValueForTime(t,i), this%zValueForTime(t,i)
      end do
      end do
      flush (unit)
      close (unit)
   end subroutine

   subroutine write_to_xdmf_file(this)
      type(movie_probe_output_t), intent(inout) :: this

      character(len=256) :: xdmf_filename
      character(len=256) :: h5_filename
      character(len=256) :: attributeBaseName
      integer :: xdmfunit, error
      integer, dimension(3) :: topologyDimensions
      integer, dimension(4) :: h5_dimensions
      xdmf_filename = add_extension(this%filesPath, ".xdmf")
      h5_filename = add_extension(get_last_component(this%filesPath), ".h5")

      topologyDimensions(1) = this%auxCoords%x - this%mainCoords%x + 1
      topologyDimensions(2) = this%auxCoords%y - this%mainCoords%y + 1
      topologyDimensions(3) = this%auxCoords%z - this%mainCoords%z + 1

      h5_dimensions(1) = this%nTime + this%nTimesFlushed
      h5_dimensions(2) = topologyDimensions(3)
      h5_dimensions(3) = topologyDimensions(2)
      h5_dimensions(4) = topologyDimensions(1)

      select case(this%component)
      case(iCur, iCurX, iCurY, iCurZ); attributeBaseName = 'CurrenDensity'
      case(iMEC, iExC, iEyC, iEzC); attributeBaseName = 'ElectricField'
      case(iMHC, iHxC, iHyC, iHzC); attributeBaseName = 'MagneticField'
      end select
      
      open(newunit=xdmfunit, file=trim(xdmf_filename), status='replace', position='append', iostat=error)
      call xdmf_write_header_file(xdmfunit, 'movieProbe')

      call xdmf_write_topology(xdmfunit, topologyDimensions)
      call xdmf_write_geometry(xdmfunit, topologyDimensions, h5_filename)
      call xdmf_write_time_array(xdmfunit, this%nTime + this%nTimesFlushed, h5_filename)
      
      select case(this%component)
      case(iCur, iMEC, iMHC)
         call xdmf_write_scalar_attribute(xdmfunit, h5_dimensions, h5_filename, trim(attributeBaseName)//'X')  
         call xdmf_write_scalar_attribute(xdmfunit, h5_dimensions, h5_filename, trim(attributeBaseName)//'Y')  
         call xdmf_write_scalar_attribute(xdmfunit, h5_dimensions, h5_filename, trim(attributeBaseName)//'Z')  
      case(iCurX, iEXC, iHXC)
         call xdmf_write_scalar_attribute(xdmfunit, h5_dimensions, h5_filename, trim(attributeBaseName)//'X')  
      case(iCurY, iEyC, iHyC) 
         call xdmf_write_scalar_attribute(xdmfunit, h5_dimensions, h5_filename, trim(attributeBaseName)//'Y')  
      case(iCurZ, iEZC, iHzC) 
         call xdmf_write_scalar_attribute(xdmfunit, h5_dimensions, h5_filename, trim(attributeBaseName)//'Z')  
      end select

      call xdmf_write_footer_file(xdmfunit)

      close(xdmfunit)
   end subroutine

   subroutine write_to_h5_file(this)
      type(movie_probe_output_t), intent(inout) :: this

      integer(HID_T) :: file_id
      integer :: i, error, probeDimensions(3)
      real(dp), allocatable, dimension(:,:,:,:) :: h5Table
      character(len=256) :: h5_filename, h5_filepath
      character(len=256) :: attributeBaseName
      h5_filepath = add_extension(this%filesPath, ".h5")
      h5_filename = get_last_component(h5_filepath)

      !Only stores the volume associated to that probe

      probeDimensions(1) = this%auxCoords%x - this%mainCoords%x + 1
      probeDimensions(2) = this%auxCoords%y - this%mainCoords%y + 1
      probeDimensions(3) = this%auxCoords%z - this%mainCoords%z + 1

      select case(this%component)
      case(iCur, iCurX, iCurY, iCurZ); attributeBaseName = 'CurrenDensity'
      case(iMEC, iExC, iEyC, iEzC); attributeBaseName = 'ElectricField'
      case(iMHC, iHxC, iHyC, iHzC); attributeBaseName = 'MagneticField'
      end select

      call H5open_f(error)
      call H5Fopen_f(trim(h5_filepath), H5F_ACC_RDWR_F, file_id, error)

      call h5_append_rows_to_dataset(file_id, 'times', this%timeStep(:this%nTime))

      allocate(h5Table(probeDimensions(1), probeDimensions(2), probeDimensions(3), this%nTime)) !(x,y,z,t)
      if (any([iCur, iMEC, iMHC, iCurX, iExC, iHxC]==this%component)) then
         h5Table = 0_dp
         do i=1, this%nPoints !Readjust idx into hyperslab (x,y,z,t)
            h5Table(this%coords(1,i) - this%mainCoords%x + 1, &
                    this%coords(2,i) - this%mainCoords%y + 1, &
                    this%coords(3,i) - this%mainCoords%z + 1, &
                    : ) = this%xValueForTime(:this%nTime, i) 
         end do
         call h5_append_rows_to_dataset(file_id, trim(attributeBaseName)//'X', h5Table)
      end if
      if (any([iCur, iMEC, iMHC, iCurY, iEyC, iHyC]==this%component)) then
         h5Table = 0_dp
         do i=1, this%nPoints !Readjust idx into hyperslab (x,y,z,t)
            h5Table(this%coords(1,i) - this%mainCoords%x + 1, &
                    this%coords(2,i) - this%mainCoords%y + 1, &
                    this%coords(3,i) - this%mainCoords%z + 1, &
                    : ) = this%yValueForTime(:this%nTime, i) 
         end do
         call h5_append_rows_to_dataset(file_id, trim(attributeBaseName)//'Y', h5Table)
      end if
      if (any([iCur, iMEC, iMHC, iCurZ, iEzC, iHzC]==this%component)) then
         h5Table = 0_dp
         do i=1, this%nPoints !Readjust idx into hyperslab (x,y,z,t)
            h5Table(this%coords(1,i) - this%mainCoords%x + 1, &
                    this%coords(2,i) - this%mainCoords%y + 1, &
                    this%coords(3,i) - this%mainCoords%z + 1, &
                    : ) = this%zValueForTime(:this%nTime, i) 
         end do
         call h5_append_rows_to_dataset(file_id, trim(attributeBaseName)//'Z', h5Table)
      end if
      deallocate(h5Table)
   
      call H5Fclose_f(file_id, error)
      call H5close_f(error)
      if (error/=0) stop
   end subroutine write_to_h5_file

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
