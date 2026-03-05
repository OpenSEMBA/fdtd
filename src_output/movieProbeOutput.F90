module mod_movieProbeOutput
   use FDETYPES
   USE mod_UTILS
   use Report
   use outputTypes
   use mod_outputUtils
   use mod_volumicProbeUtils
   use HDF5
   use mod_xdmfAPI
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

      this%mainCoords = lowerBound
      this%auxCoords = upperBound
      this%component = field
      this%domain = domain

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
      call create_movie_files(this%filesPath, this%coords ,this%nPoints, error)
   end subroutine init_movie_probe_output

   subroutine create_bin_file(filePath, error)
      character(len=*), intent(in) :: filePath
      integer, intent(out) :: error
      call create_file_with_path(add_extension(filePath, binaryExtension), error)
   end subroutine 

   subroutine create_movie_files(filePath, coords, nPoints, error)
      character(len=*), intent(in) :: filePath
      integer(kind=SINGLE), dimension(:, :), intent(in) :: coords
      integer, intent(in) :: nPoints
      integer, intent(out) :: error
      real(dp), allocatable, dimension(:, :) :: coordsReal

      integer(HID_T) :: file_id
      character(len=BUFSIZE) :: h5_filename

      h5_filename = add_extension(filePath, ".h5")

      call H5open_f(error)
      call create_h5_file(trim(h5_filename), file_id)

      allocate (coordsReal(3, nPoints))
      coordsReal = real(coords, dp)
      call write_dataset(file_id, "coords", coordsReal)
      deallocate(coordsReal)

      call init_extendable_2d_dataset(file_id, 'xVal', nPoints, BUFSIZE)
      call init_extendable_2d_dataset(file_id, 'yVal', nPoints, BUFSIZE)
      call init_extendable_2d_dataset(file_id, 'zVal', nPoints, BUFSIZE)
      call H5Fclose_f(file_id, error)
      call H5close_f(error)
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
            call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
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
            call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
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
            call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
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
            call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
         end select
      end if
   end subroutine update_movie_probe_output

   subroutine flush_movie_probe_output(this)
      type(movie_probe_output_t), intent(inout) :: this
      integer :: i

      call write_bin_file(this)
      call write_to_xdmf_h5(this)

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

   subroutine write_to_xdmf_h5(this)

      type(movie_probe_output_t), intent(inout) :: this

      integer(HID_T) :: file_id
      integer :: t, error, xdmfunit
      real(dp), allocatable, dimension(:, :) :: coordsReal
      character(len=256) :: h5_filename, h5_filepath
      character(len=256) :: xdmf_filename
      character(len=10) :: dimension_string
      character(len=10) :: nCoordsString
      character(len=14) :: stepName
      h5_filepath = add_extension(this%filesPath, ".h5")
      h5_filename = get_last_component(h5_filepath)

      call H5open_f(error)
      call H5Fopen_f(trim(h5_filepath), H5F_ACC_RDWR_F, file_id, error)

      
      write(dimension_string, '(I0,1X,I0)') this%nPoints, this%nTimesFlushed + this%nTime
      write(nCoordsString, '(I0, I0)') this%nPoints, 1
      do t = 1, this%nTime
         write(stepName, '((I5.5))') this%nTimesFlushed + t
         call append_rows_dataset(file_id, "xVal", reshape(this%xValueForTime(t, :), [1, this%nPoints]))
         call append_rows_dataset(file_id, "yVal", reshape(this%yValueForTime(t, :), [1, this%nPoints]))
         call append_rows_dataset(file_id, "zVal", reshape(this%zValueForTime(t, :), [1, this%nPoints]))
         if (mod(this%nTimesFlushed + t, 10) == 0) then
            
            xdmf_filename = add_extension(add_extension(this%filesPath, ".ts_"//stepName), ".xdmf")
            open(newunit=xdmfunit, file=trim(xdmf_filename), status='replace', position='append', iostat=error)

            call xdmf_write_header_file(xdmfunit)

            call xdmf_create_grid_step_info(xdmfunit, stepName, real(this%timeStep(t)), trim(h5_filename), this%nPoints)
            
            call xdmf_write_attribute(xdmfunit, 'xVal')
            call xdmf_write_hyperslab_data_item(xdmfunit, nCoordsString)
            call xdmf_write_h5_acces_data_item(xdmfunit, 0, this%nTimesFlushed + t - 1,  this%nPoints, this%nTimesFlushed + this%nTime)
            call xdmf_close_data_item(xdmfunit)
            call xdmf_write_h5_data_item(xdmfunit, trim(h5_filename), '/xVal', dimension_string)
            call xdmf_close_data_item(xdmfunit)
            call xdmf_close_data_item(xdmfunit)
            call xdmf_close_attribute(xdmfunit)

            call xdmf_write_attribute(xdmfunit, 'yVal')
            call xdmf_write_hyperslab_data_item(xdmfunit, nCoordsString)
            call xdmf_write_h5_acces_data_item(xdmfunit, 0, this%nTimesFlushed + t - 1,  this%nPoints, this%nTimesFlushed + this%nTime)
            call xdmf_close_data_item(xdmfunit)
            call xdmf_write_h5_data_item(xdmfunit, trim(h5_filename), '/yVal', dimension_string)
            call xdmf_close_data_item(xdmfunit)
            call xdmf_close_data_item(xdmfunit)
            call xdmf_close_attribute(xdmfunit)

            call xdmf_write_attribute(xdmfunit, 'zVal')
            call xdmf_write_hyperslab_data_item(xdmfunit, nCoordsString)
            call xdmf_write_h5_acces_data_item(xdmfunit, 0, this%nTimesFlushed + t - 1,  this%nPoints, this%nTimesFlushed + this%nTime)
            call xdmf_close_data_item(xdmfunit)
            call xdmf_write_h5_data_item(xdmfunit, trim(h5_filename), '/zVal', dimension_string)
            call xdmf_close_data_item(xdmfunit)
            call xdmf_close_data_item(xdmfunit)
            call xdmf_close_attribute(xdmfunit)

            call xdmf_close_grid(xdmfunit)
         
            call xdmf_write_footer_file(xdmfunit)
            close(xdmfunit)
         endif
      end do

      call H5Fclose_f(file_id, error)
      call H5close_f(error)
   end subroutine write_to_xdmf_h5

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

end module mod_movieProbeOutput
