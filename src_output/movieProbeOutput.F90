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
   private :: write_vtu_timestep
   private :: update_pvd
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
      character(len=BUFSIZE) :: pdvFileName

      this%mainCoords = lowerBound
      this%auxCoords  = upperBound
      this%component  = field
      this%domain     = domain
      this%path       = get_output_path(this, outputTypeExtension, field, control%mpidir)
      
      pdvFileName = add_extension(get_last_component(this%path), pvdExtension)
      this%pvdPath = join_path(this%path, pdvFileName)

      call create_folder(this%path, error)
      call create_file_with_path(add_extension(this%path, binaryExtension), error)

      call find_and_store_important_coords(this%mainCoords, this%auxCoords, this%component, problemInfo, this%nPoints, this%coords)
      call alloc_and_init(this%timeStep, BuffObse, 0.0_RKIND_tiempo)

      ! Allocate value arrays based on component type
         call alloc_and_init(this%xValueForTime, BuffObse, this%nPoints, 0.0_RKIND)
         call alloc_and_init(this%yValueForTime, BuffObse, this%nPoints, 0.0_RKIND)
         call alloc_and_init(this%zValueForTime, BuffObse, this%nPoints, 0.0_RKIND)
   end subroutine init_movie_probe_output

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
      call write_to_xdmf(this)

      call clear_memory_data(this)
   end subroutine flush_movie_probe_output

   !===========================
   ! Private routines
   !===========================

   subroutine write_bin_file(this)
      ! Check type definition for binary format
      type(movie_probe_output_t), intent(inout) :: this
      integer :: i, t, unit

      open (unit=unit, file=add_extension(this%path, binaryExtension), &
            status='old', form='unformatted', position='append', access='stream')
      do t = 1, this%nTime
      do i = 1, this%nPoints
         write(unit) this%timeStep(t), this%coords(1,i), this%coords(2,i), this%coords(3,i), this%xValueForTime(t,i), this%yValueForTime(t,i), this%zValueForTime(t,i)
      end do
      end do
      close (unit)
   end subroutine


   subroutine read_bin_file(this)
      ! Check type definition for binary format
      type(movie_probe_output_t), intent(inout) :: this
      integer :: unit
      integer :: iostat
      real(kind=RKIND_tiempo) timeStamp
      integer(kind=SINGLE) :: x, y, z
      real(kind=RKIND) :: xVal, yVal, zVal
      integer(kind=4) :: dataSize

      open (unit=unit, file=add_extension(this%path, binaryExtension), &
            status='old', form='unformatted', access='stream', iostat=iostat)
      if (iostat /= 0) then
         print *, 'Error opening file!'
         return
      end if

      ! Read until end-of-file
      do
         read (unit, iostat=iostat) timeStamp, x, y, z, xVal, yVal, zVal
         if (iostat /= 0) exit  ! EOF or error
         print *, timeStamp, x, y, z, xVal, yVal, zVal
      end do

      close (unit)
   end subroutine

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

   subroutine write_vtu_timestep(this, stepIndex, filename)
      type(movie_probe_output_t), intent(in) :: this
      integer, intent(in)                    :: stepIndex
      character(len=*), intent(in)           :: filename

      integer :: npts, i, ierr
      real(kind=RKIND), allocatable :: x(:), y(:), z(:)
      real(kind=RKIND), allocatable :: Componentx(:), Componenty(:), Componentz(:)
      logical :: writeX, writeY, writeZ
      character(len=BUFSIZE) :: requestName

   end subroutine write_vtu_timestep

   subroutine update_pvd(this, stepIndex, PVDfilePath)
      implicit none
      type(movie_probe_output_t), intent(in) :: this
      integer, intent(in) :: stepIndex
      character(len=*), intent(in) :: PVDfilePath
      character(len=64) :: ts
      character(len=BUFSIZE) :: newVTUfilename
      integer :: unit

   end subroutine update_pvd

   subroutine clear_memory_data(this)
      type(movie_probe_output_t), intent(inout) :: this

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
