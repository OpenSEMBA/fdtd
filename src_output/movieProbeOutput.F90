module mod_movieProbeOutput
   use FDETYPES
   use Report
   use outputTypes
   use mod_outputUtils
   implicit none
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: init_movie_probe_output
   public :: update_movie_probe_output
   public :: flush_movie_probe_output
   !===========================

   !===========================
   !  Private interface summary
   !===========================
   ! Data Extraction & Processing
   private :: count_required_coords
   private :: save_current_module
   private :: save_current_component
   private :: save_current
   private :: save_field_module
   private :: save_field_component
   private :: save_field

   ! Output & File Management
   private :: write_vtu_timestep
   private :: update_pvd

   ! Validation Logic (Functions)
   private :: isValidPointForCurrent
   private :: isValidPointForField
   private :: volumicCurrentRequest
   private :: volumicElectricRequest
   private :: volumicMagneticRequest
   private :: componentCurrentRequest
   private :: componentFieldRequest
   !===========================

   abstract interface
      logical function logical_func(component, i, j, k, problemInfo)
         import :: problem_info_t
         type(problem_info_t), intent(in) :: problemInfo
         integer, intent(in) :: component, i, j, k
      end function logical_func
   end interface

contains

   subroutine init_movie_probe_output(this, lowerBound, upperBound, field, domain, control, problemInfo, outputTypeExtension)
      type(movie_probe_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in) :: field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      type(sim_control_t), intent(in) :: control
      type(problem_info_t), intent(in) :: problemInfo

      type(domain_t), intent(in) :: domain

      this%mainCoords = lowerBound
      this%auxCoords = upperBound
      this%component = field !This can refer to electric, magnetic or currentDensity
      this%domain = domain
      this%path = get_output_path()

      call find_and_store_important_coords(this, problemInfo)

      call alloc_and_init(this%timeStep, BuffObse, 0.0_RKIND_tiempo)

      if (any(VOLUMIC_M_MEASURE == this%component)) then
         call alloc_and_init(this%xValueForTime, BuffObse, this%nPoints, 0.0_RKIND)
         call alloc_and_init(this%yValueForTime, BuffObse, this%nPoints, 0.0_RKIND)
         call alloc_and_init(this%zValueForTime, BuffObse, this%nPoints, 0.0_RKIND)
      else
         if (any(VOLUMIC_X_MEASURE == this%component)) then
            call alloc_and_init(this%xValueForTime, BuffObse, this%nPoints, 0.0_RKIND)
         elseif (any(VOLUMIC_Y_MEASURE == this%component)) then
            call alloc_and_init(this%yValueForTime, BuffObse, this%nPoints, 0.0_RKIND)
         elseif (any(VOLUMIC_Z_MEASURE == this%component)) then
            call alloc_and_init(this%zValueForTime, BuffObse, this%nPoints, 0.0_RKIND)
         else
            call StopOnError(control%layoutnumber, control%size, "Unexpected output type for movie probe")
         end if
      end if

   contains
      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_coordinates_extension(this%mainCoords, this%auxCoords, control%mpidir)
         prefixFieldExtension = get_prefix_extension(field, control%mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

   end subroutine init_movie_probe_output

   subroutine update_movie_probe_output(this, step, fieldsReference, control, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(sim_control_t), intent(in) :: control
      type(problem_info_t), intent(in) :: problemInfo
      type(fields_reference_t), intent(in) :: fieldsReference

      integer(kind=4) :: request
      request = this%component

      this%nTime = this%nTime + 1

      if (any(VOLUMIC_M_MEASURE == request)) then
         select case (request)
         case (iCur); call save_current_module(this, fieldsReference, step, problemInfo)
         case (iMEC); call save_field_module(this, fieldsReference%E, request, step, problemInfo)
         case (iMHC); call save_field_module(this, fieldsReference%H, request, step, problemInfo)
         case default; call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_X_MEASURE == request)) then
         select case (request)
         case (iCurX); call save_current_component(this, this%xValueForTime, fieldsReference, step, problemInfo, iEx)
         case (iExC); call save_field_component(this, this%xValueForTime, fieldsReference%E%x, step, problemInfo, iEx)
         case (iHxC); call save_field_component(this, this%xValueForTime, fieldsReference%H%x, step, problemInfo, iHx)
         case default; call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_Y_MEASURE == request)) then
         select case (request)
         case (iCurY); call save_current_component(this, this%yValueForTime, fieldsReference, step, problemInfo, iEy)
         case (iEyC); call save_field_component(this, this%yValueForTime, fieldsReference%E%y, step, problemInfo, iEy)
         case (iHyC); call save_field_component(this, this%yValueForTime, fieldsReference%H%y, step, problemInfo, iHy)
         case default; call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_Z_MEASURE == request)) then
         select case (request)
         case (iCurZ); call save_current_component(this, this%zValueForTime, fieldsReference, step, problemInfo, iEz)
         case (iEzC); call save_field_component(this, this%zValueForTime, fieldsReference%E%z, step, problemInfo, iEz)
         case (iHzC); call save_field_component(this, this%zValueForTime, fieldsReference%H%z, step, problemInfo, iHz)
         case default; call StopOnError(control%layoutnumber, control%size, "Volumic measure not supported")
         end select
      end if
   end subroutine update_movie_probe_output

   subroutine save_current_module(this, fieldsReference, simTime, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      type(fields_reference_t), intent(in) :: fieldsReference
      real(kind=RKIND_tiempo), intent(in) :: simTime
      type(problem_info_t), intent(in) :: problemInfo

      integer :: i, j, k, coordIdx

      this%timeStep(this%nTime) = simTime

      coordIdx = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (isValidPointForCurrent(iCur, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_current(this%xValueForTime, this%nTime, coordIdx, iEx, i, j, k, fieldsReference)
            call save_current(this%yValueForTime, this%nTime, coordIdx, iEy, i, j, k, fieldsReference)
            call save_current(this%zValueForTime, this%nTime, coordIdx, iEz, i, j, k, fieldsReference)
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_current_component(this, currentData, fieldsReference, simTime, problemInfo, fieldDir)
      type(movie_probe_output_t), intent(inout) :: this
      real(kind=RKIND), intent(inout) :: currentData(:, :)
      type(fields_reference_t), intent(in) :: fieldsReference
      real(kind=RKIND_tiempo), intent(in) :: simTime
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: fieldDir

      integer :: i, j, k, coordIdx

      this%timeStep(this%nTime) = simTime

      coordIdx = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (isValidPointForCurrent(fieldDir, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_current(currentData, this%nTime, coordIdx, fieldDir, i, j, k, fieldsReference)
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_current(currentData, timeIdx, coordIdx, field, i, j, k, fieldsReference)
      real(kind=RKIND), intent(inout) :: currentData(:, :)
      integer(kind=SINGLE), intent(in) :: timeIdx, coordIdx, field, i, j, k
      type(fields_reference_t), intent(in) :: fieldsReference

      real(kind=RKIND) :: jdir
      jdir = computeJ(field, i, j, k, fieldsReference)
      currentData(timeIdx, coordIdx) = jdir
   end subroutine

   subroutine save_field_module(this, field, request, simTime, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      type(field_data_t), intent(in) :: field
      real(kind=RKIND_tiempo), intent(in) :: simTime
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: request

      integer :: i, j, k, coordIdx

      this%timeStep(this%nTime) = simTime

      coordIdx = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (isValidPointForField(request, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            call save_field(this%xValueForTime, this%nTime, coordIdx, field%x(i, j, k))
            call save_field(this%yValueForTime, this%nTime, coordIdx, field%y(i, j, k))
            call save_field(this%zValueForTime, this%nTime, coordIdx, field%z(i, j, k))
         end if
      end do
      end do
      end do

   end subroutine

   subroutine save_field_component(this, fieldData, fieldComponent, simTime, problemInfo, fieldDir)
      type(movie_probe_output_t), intent(inout) :: this
      real(kind=RKIND), intent(inout) :: fieldData(:, :)
      real(kind=RKIND), intent(in) :: fieldComponent(:, :, :)
      real(kind=RKIND_tiempo), intent(in) :: simTime
      type(problem_info_t), intent(in) :: problemInfo
      integer, intent(in) :: fieldDir

      integer :: i, j, k, coordIdx

      this%timeStep(this%nTime) = simTime

      coordIdx = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (isValidPointForField(fieldDir, i, j, k, problemInfo)) then
            coordIdx = coordIdx + 1
            coordIdx = coordIdx + 1
            call save_field(fieldData, this%nTime, coordIdx, fieldComponent(i, j, k))
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_field(fieldData, timeIdx, coordIdx, fieldValue)
      real(kind=RKIND), intent(inout) :: fieldData(:, :)
      integer(kind=SINGLE), intent(in) :: timeIdx, coordIdx
      real(kind=RKIND), intent(in) :: fieldValue
      fieldData(timeIdx, coordIdx) = fieldValue
   end subroutine

   subroutine flush_movie_probe_output(this)
      type(movie_probe_output_t), intent(inout) :: this
      integer :: status, i

      do i = 1, this%nTime
         call update_pvd(this, i, this%fileUnitTime)
      end do
      call clear_memory_data()

   contains
      subroutine clear_memory_data()
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

   end subroutine flush_movie_probe_output

   subroutine write_vtu_timestep(this, stepIndex, filename)
      use vtk_fortran
      implicit none

      type(movie_probe_output_t), intent(in) :: this
      integer, intent(in) :: stepIndex
      character(len=*), intent(in) :: filename

      character(len=BUFSIZE) :: requestName
      type(vtk_file) :: vtkOutput
      integer :: ierr, npts, i
      real(kind=RKIND), allocatable :: x(:), y(:), z(:)
      real(kind=RKIND), allocatable :: Componentx(:), Componenty(:), Componentz(:)
      logical :: writeX, writeY, writeZ

      !================= Determine the measure type =================

      if (any(CURRENT_MEASURE == this%component)) then
         requestName = 'Current'
      else if (any(ELECTRIC_FIELD_MEASURE == this%component)) then
         requestName = 'Electric'
      else if (any(MAGNETIC_FIELD_MEASURE == this%component)) then
         requestName = 'Magnetic'
      else
         requestName = 'Unknown'
      end if

      !================= Determine which components to write =================
      writeX = any(VOLUMIC_M_MEASURE == this%component) .or. any(VOLUMIC_X_MEASURE == this%component)
      writeY = any(VOLUMIC_M_MEASURE == this%component) .or. any(VOLUMIC_Y_MEASURE == this%component)
      writeZ = any(VOLUMIC_M_MEASURE == this%component) .or. any(VOLUMIC_Z_MEASURE == this%component)

      !================= Allocate and fill coordinates =================
      npts = this%nPoints
      allocate (x(npts), y(npts), z(npts))
      do i = 1, npts
         x(i) = this%coords(1, i)
         y(i) = this%coords(2, i)
         z(i) = this%coords(3, i)
      end do

      ierr = vtkOutput%initialize(format='ASCII', filename=trim(filename), mesh_topology='UnstructuredGrid')
      ierr = vtkOutput%xml_writer%write_geo(n=npts, x=x, y=y, z=z)

      !================= Allocate and fill component arrays =================
      if (writeX) then
         allocate (Componentx(npts))
         do i = 1, npts
            Componentx(i) = this%xValueForTime(stepIndex, i)
         end do
      end if

      if (writeY) then
         allocate (Componenty(npts))
         do i = 1, npts
            Componenty(i) = this%yValueForTime(stepIndex, i)
         end do
      end if

      if (writeZ) then
         allocate (Componentz(npts))
         do i = 1, npts
            Componentz(i) = this%zValueForTime(stepIndex, i)
         end do
      end if

      !================= Write arrays to VTK =================
      if (writeX) then
         ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='open')
         ierr = vtkOutput%xml_writer%write_dataarray(data_name=trim(adjustl(requestName))//'X', x=Componentx)
         ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='close')
         deallocate (Componentx)
      end if

      if (writeY) then
         ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='open')
         ierr = vtkOutput%xml_writer%write_dataarray(data_name=trim(adjustl(requestName))//'Y', x=Componenty)
         ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='close')
         deallocate (Componenty)
      end if

      if (writeZ) then
         ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='open')
         ierr = vtkOutput%xml_writer%write_dataarray(data_name=trim(adjustl(requestName))//'Z', x=Componentz)
         ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='close')
         deallocate (Componentz)
      end if

      ierr = vtkOutput%xml_writer%finalize()
      deallocate (x, y, z)

   end subroutine write_vtu_timestep

   subroutine update_pvd(this, stepIndex, unitPVD)
      implicit none
      type(movie_probe_output_t), intent(in) :: this
      integer, intent(in) :: stepIndex
      integer, intent(in) :: unitPVD
      character(len=64) :: ts
      character(len=256) :: filename

      ! Generamos nombre del archivo VTU para este timestep
      write (filename, '(A,A,I4.4,A)') trim(this%path), '_ts', stepIndex, '.vtu'

      ! Escribimos el VTU correspondiente
      call write_vtu_timestep(this, stepIndex, filename)

      ! Añadimos entrada en el PVD
      write (ts, '(ES16.8)') this%timeStep(stepIndex)
      write (unitPVD, '(A)') '    <DataSet timestep="'//trim(ts)// &
         '" group="" part="0" file="'//trim(filename)//'"/>'
   end subroutine update_pvd

   subroutine find_and_store_important_coords(this, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      type(problem_info_t), intent(in) :: problemInfo

      call count_required_coords(this, problemInfo)
      call alloc_and_init(this%coords, 3, this%nPoints, 0_SINGLE)
      call store_required_coords(this, problemInfo)
   end subroutine

   subroutine count_required_coords(this, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      type(problem_info_t), intent(in) :: problemInfo

      integer :: i, j, k

      procedure(logical_func), pointer :: checker => null()  ! Pointer to logical function
      integer :: component, count
      call get_checker_and_component(this, checker, component)

      count = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (checker(component, i, j, k, problemInfo)) count = count + 1
      end do
      end do
      end do

      this%nPoints = count

   end subroutine

   subroutine store_required_coords(this, problemInfo)
      type(movie_probe_output_t), intent(inout) :: this
      type(problem_info_t), intent(in) :: problemInfo

      integer :: i, j, k

      procedure(logical_func), pointer :: checker => null()  ! Pointer to logical function
      integer :: component, count
      call get_checker_and_component(this, checker, component)

      count = 0
      do i = this%mainCoords%x, this%auxCoords%x
      do j = this%mainCoords%y, this%auxCoords%y
      do k = this%mainCoords%z, this%auxCoords%z
         if (checker(component, i, j, k, problemInfo)) then
            count = count + 1
            this%coords(1, count) = i
            this%coords(2, count) = j
            this%coords(3, count) = k
         end if
      end do
      end do
      end do
   end subroutine

   subroutine get_checker_and_component(this, checker, component)
      type(movie_probe_output_t), intent(in) :: this
      procedure(logical_func), pointer, intent(out) :: checker
      integer, intent(out) :: component

      select case (this%component)
      case (iCur)
         checker => volumicCurrentRequest
         component = iCur
      case (iMEC)
         checker => volumicElectricRequest
         component = iMEC
      case (iMHC)
         checker => volumicMagneticRequest
         component = iMHC
      case (iCurx)
         checker => componentCurrentRequest
         component = iEx
      case (iExC)
         checker => componentFieldRequest
         component = iEx
      case (iHxC)
         checker => componentFieldRequest
         component = iHx
      case (iCurY)
         checker => componentCurrentRequest
         component = iEy
      case (iEyC)
         checker => componentFieldRequest
         component = iEy
      case (iHyC)
         checker => componentFieldRequest
         component = iHy
      case (iCurZ)
         checker => componentCurrentRequest
         component = iEz
      case (iEzC)
         checker => componentFieldRequest
         component = iEz
      case (iHzC)
         checker => componentFieldRequest
         component = iHz
      end select
   end subroutine

   logical function isValidPointForCurrent(request, i, j, k, problemInfo)
      integer, intent(in) :: i, j, k, request
      type(problem_info_t) :: problemInfo
      select case (request)
      case (iCur)
         isValidPointForCurrent = volumicCurrentRequest(request, i, j, k, problemInfo)
      case (iEx, iEy, iEz)
         isValidPointForCurrent = componentCurrentRequest(request, i, j, k, problemInfo)
      case default
         isValidPointForCurrent = .false.
      end select
   end function

   logical function isValidPointForField(request, i, j, k, problemInfo)
      integer, intent(in) :: i, j, k, request
      type(problem_info_t) :: problemInfo
      select case (request)
      case (iMEC)
         isValidPointForField = volumicElectricRequest(request, i, j, k, problemInfo)
      case (iMHC)
         isValidPointForField = volumicMagneticRequest(request, i, j, k, problemInfo)
      case (iEx, iEy, iEz, iHx, iHy, iHz)
         isValidPointForField = componentFieldRequest(request, i, j, k, problemInfo)
      case default
         isValidPointForField = .false.
      end select
   end function

   logical function volumicCurrentRequest(request, i, j, k, problemInfo)
      integer, intent(in) :: i, j, k, request
      type(problem_info_t), intent(in) :: problemInfo
      volumicCurrentRequest = componentCurrentRequest(iEx, i, j, k, problemInfo) &
                              .or. componentCurrentRequest(iEy, i, j, k, problemInfo) &
                              .or. componentCurrentRequest(iEz, i, j, k, problemInfo)
   end function
   logical function volumicElectricRequest(request, i, j, k, problemInfo)
      integer, intent(in) :: i, j, k, request
      type(problem_info_t), intent(in) :: problemInfo
      volumicElectricRequest = componentFieldRequest(iEx, i, j, k, problemInfo) &
                               .or. componentFieldRequest(iEy, i, j, k, problemInfo) &
                               .or. componentFieldRequest(iEz, i, j, k, problemInfo)
   end function
   logical function volumicMagneticRequest(request, i, j, k, problemInfo)
      integer, intent(in) :: i, j, k, request
      type(problem_info_t), intent(in) :: problemInfo
      volumicMagneticRequest = componentFieldRequest(iHx, i, j, k, problemInfo) &
                               .or. componentFieldRequest(iHy, i, j, k, problemInfo) &
                               .or. componentFieldRequest(iHz, i, j, k, problemInfo)
   end function
   logical function componentCurrentRequest(fieldDir, i, j, k, problemInfo)
      integer, intent(in) :: i, j, k, fieldDir
      type(problem_info_t), intent(in) :: problemInfo
      componentCurrentRequest = isWithinBounds(fieldDir, i, j, k, problemInfo)
      if (componentCurrentRequest) then
         componentCurrentRequest = isPEC(fieldDir, i, j, k, problemInfo) &
                                   .or. isThinWire(fieldDir, i, j, k, problemInfo)
      end if
   end function
   logical function componentFieldRequest(fieldDir, i, j, k, problemInfo)
      integer, intent(in) :: i, j, k, fieldDir
      type(problem_info_t), intent(in) :: problemInfo
      componentFieldRequest = isWithinBounds(fieldDir, i, j, k, problemInfo)
   end function

end module mod_movieProbeOutput
