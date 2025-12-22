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
   private :: get_measurements_coords
   private :: save_current_data
   private :: write_vtu_timestep
   private :: update_pvd
   !===========================

contains

   subroutine init_movie_probe_output(this, lowerBound, upperBound, field, domain, geometryMedia, registeredMedia, sinpml_fullsize, outputTypeExtension, mpidir)
      type(movie_probe_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      type(MediaData_t), dimension(:), intent(in) :: registeredMedia
      type(media_matrices_t), intent(in) :: geometryMedia
      type(limit_t), dimension(:), intent(in)  :: sinpml_fullsize

      type(domain_t), intent(in) :: domain

      this%lowerBound = lowerBound
      this%upperBound = upperBound
      this%fieldComponent = field !This can refer to field or currentDensity
      this%domain = domain
      this%path = get_output_path()

      call get_measurements_coords(this, geometryMedia, registeredMedia, sinpml_fullsize)

      call alloc_and_init(this%timeStep, BuffObse, 0.0_RKIND_tiempo)

      if (any(VOLUMIC_M_MEASURE == this%fieldComponent)) then
         call alloc_and_init(this%xValueForTime, BuffObse, this%nMeasuredElements, 0.0_RKIND)
         call alloc_and_init(this%yValueForTime, BuffObse, this%nMeasuredElements, 0.0_RKIND)
         call alloc_and_init(this%zValueForTime, BuffObse, this%nMeasuredElements, 0.0_RKIND)
      else
         if (any(VOLUMIC_X_MEASURE == this%fieldComponent)) then
            call alloc_and_init(this%xValueForTime, BuffObse, this%nMeasuredElements, 0.0_RKIND)
         elseif (any(VOLUMIC_Y_MEASURE == this%fieldComponent)) then
            call alloc_and_init(this%yValueForTime, BuffObse, this%nMeasuredElements, 0.0_RKIND)
         elseif (any(VOLUMIC_Z_MEASURE == this%fieldComponent)) then
            call alloc_and_init(this%zValueForTime, BuffObse, this%nMeasuredElements, 0.0_RKIND)
         else
            call StopOnError(0, 0, "Unexpected output type for movie probe")
         end if
      end if

   contains
      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_coordinates_extension(this%lowerBound, this%upperBound, mpidir)
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

   end subroutine init_movie_probe_output

   subroutine update_movie_probe_output(this, step, geometryMedia, registeredMedia, sinpml_fullsize, fieldsReference)
      type(movie_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step

      type(media_matrices_t), intent(in) :: geometryMedia
      type(MediaData_t), dimension(:), intent(in) :: registeredMedia
      type(limit_t), dimension(:), intent(in)  :: sinpml_fullsize
      type(fields_reference_t), intent(in) :: fieldsReference

      integer(kind=4) :: request
      request = this%fieldComponent

      this%serializedTimeSize = this%serializedTimeSize + 1

      if (any(VOLUMIC_M_MEASURE == request)) then
         select case (request)
         case (iCur); call save_current_module(this, fieldsReference, step)
         case (iMEC); call save_field_module(this, fieldsReference, request, step)
         case (iMHC); call save_field_module(this, fieldsReference, request, step)
         case default; StopOnError(0, 0, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_X_MEASURE == request)) then
         select case (request)
         case (iCurX); call save_current_component(this%xValueForTime, fieldsReference, step, iEx)
         case (iExC); call save_field_component(this, fieldsReference, request, step)
         case (iHxC); call save_field_component(this, fieldsReference, request, step)
         case default; StopOnError(0, 0, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_Y_MEASURE == request)) then
         select case (request)
         case (iCurY); call save_current_component(this%yValueForTime, fieldsReference, step, iEy)
         case (iEyC); call save_field_component(this, fieldsReference, request, step)
         case (iHyC); call save_field_component(this, fieldsReference, request, step)
         case default; StopOnError(0, 0, "Volumic measure not supported")
         end select

      else if (any(VOLUMIC_Z_MEASURE == request)) then
         select case (request)
         case (iCurZ); call save_current_component(this%zValueForTime, fieldsReference, step, iEz)
         case (iEzC); call save_field_component(this, fieldsReference, request, step)
         case (iHzC); call save_field_component(this, fieldsReference, request, step)
         case default; StopOnError(0, 0, "Volumic measure not supported")
         end select
      end if
   end subroutine update_movie_probe_output

   subroutine save_current_module(this, fieldsReference, simTime)
      type(movie_probe_output_t), intent(inout) :: this
      type(fields_reference_t), intent(in) :: fieldsReference
      real(kind=RKIND_tiempo), intent(in) :: simTime

      integer :: i, j, k, coordIdx

      this%timeStep(this%serializedTimeSize) = simTime

      coordIdx = 0
      do i = this%lowerBound%x, this%upperBound%x
      do j = this%lowerBound%y, this%upperBound%y
      do k = this%lowerBound%z, this%upperBound%z
         if (saveCurrentFrom(i, j, k)) then
            coordIdx = coordIdx + 1
            call save_current(this%xValueForTime, timeIdx, coordIdx, iEx, i, j, k, fieldsReference)
            call save_current(this%yValueForTime, timeIdx, coordIdx, iEy, i, j, k, fieldsReference)
            call save_current(this%zValueForTime, timeIdx, coordIdx, iEz, i, j, k, fieldsReference)
         end if
      end do
      end do
      end do
   end subroutine

   subroutine save_current_component(currentData, fieldsReference, simTime, fieldDir)
      real(kind=RKIND), intent(inout) :: currentData(:, :)
      type(fields_reference_t), intent(in) :: fieldsReference
      real(kind=RKIND_tiempo), intent(in) :: simTime
      integer, intent(in) :: fieldDir

      integer :: i, j, k, coordIdx

      this%timeStep(this%serializedTimeSize) = simTime

      coordIdx = 0
      do i = this%lowerBound%x, this%upperBound%x
      do j = this%lowerBound%y, this%upperBound%y
      do k = this%lowerBound%z, this%upperBound%z
         if (saveCurrentFrom(i, j, k)) then
            coordIdx = coordIdx + 1
            call save_current(currentData, timeIdx, coordIdx, fieldDir, i, j, k, fieldsReference)
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

   subroutine save_field_module(this, fieldsReference, request, simTime)
      type(movie_probe_output_t), intent(inout) :: this
      type(fields_reference_t), intent(in) :: fieldsReference
      integer, intent(in) :: request
      real(kind=RKIND_tiempo), intent(in) :: simTime

      type(field_data_t), pointer :: field
      integer :: i, j, k, coordIdx

      if (request == iMEC) then
         field => fieldsReference%E
      else if (request == iMHC) then
         field => fieldsReference%H
      end if

      this%timeStep(this%serializedTimeSize) = simTime

      coordIdx = 0
      do i = this%lowerBound%x, this%upperBound%x
      do j = this%lowerBound%y, this%upperBound%y
      do k = this%lowerBound%z, this%upperBound%z
         if (saveFieldFrom(i, j, k)) then
            coordIdx = coordIdx + 1
            this%xValueForTime(timeIdx, coordIdx) = field%x(i, j, k)
            this%yValueForTime(timeIdx, coordIdx) = field%y(i, j, k)
            this%zValueForTime(timeIdx, coordIdx) = field%z(i, j, k)
         end if
      end do
      end do
      end do

   end subroutine

   subroutine save_field_component(fieldData, fieldsReference, request, simTime)
      real(kind=RKIND), intent(in) :: fieldData
      type(fields_reference_t), intent(in) :: fieldsReference
      integer, intent(in) :: request
      real(kind=RKIND_tiempo), intent(in) :: simTime

      real(kind=RKIND), pointer :: fieldComponent(:,:,:)
      integer :: i, j, k, coordIdx

      fieldComponent = get_field_component()


      this%timeStep(this%serializedTimeSize) = simTime

      coordIdx = 0
      do i = this%lowerBound%x, this%upperBound%x
      do j = this%lowerBound%y, this%upperBound%y
      do k = this%lowerBound%z, this%upperBound%z
         if (saveFieldFrom(i, j, k)) then
            coordIdx = coordIdx + 1
            this%xValueForTime(timeIdx, coordIdx) = field%x(i, j, k)
            this%yValueForTime(timeIdx, coordIdx) = field%y(i, j, k)
            this%zValueForTime(timeIdx, coordIdx) = field%z(i, j, k)
         end if
      end do
      end do
      end do

   end subroutine


   subroutine flush_movie_probe_output(this)
      type(movie_probe_output_t), intent(inout) :: this
      integer :: status, i

      do i = 1, this%serializedTimeSize
         call update_pvd(this, i, this%PDVUnit)
      end do
      call clear_memory_data()

   contains
      subroutine clear_memory_data()
         this%serializedTimeSize = 0
         this%timeStep = 0.0_RKIND
         this%xValueForTime = 0.0_RKIND
         this%yValueForTime = 0.0_RKIND
         this%zValueForTime = 0.0_RKIND
      end subroutine clear_memory_data

   end subroutine flush_movie_probe_output

   subroutine get_measurements_coords(this, geometryMedia, registeredMedia, sinpml_fullsize)
      type(movie_probe_output_t), intent(inout) :: this
      type(media_matrices_t), intent(in) :: geometryMedia
      type(MediaData_t), dimension(:), intent(in) :: registeredMedia
      type(limit_t), dimension(:), intent(in) :: sinpml_fullsize

      integer(kind=SINGLE) :: i, j, k, field
      integer(kind=SINGLE) :: istart, jstart, kstart, iend, jend, kend
      integer(kind=SINGLE) :: count
      ! Limites de la región de interés
      istart = this%lowerBound%x
      jstart = this%lowerBound%y
      kstart = this%lowerBound%z

      iend = this%upperBound%x
      jend = this%upperBound%y
      kend = this%upperBound%z

      ! Primer barrido para contar cuántos puntos válidos
      count = 0
      select case (this%fieldComponent)
      case (iCur)
         do i = istart, iend
         do j = jstart, jend
         do k = kstart, kend
            do field = iEx, iEz
               if (isWithinBounds(field, i, j, k, sinpml_fullsize)) then
                  if (isPEC(field, i, j, k, geometryMedia, registeredMedia)) then
                     count = count + 1
                  end if
               end if
            end do
         end do
         end do
         end do
      end select

      this%nMeasuredElements = count

      allocate (this%coords(3, this%nMeasuredElements))

      count = 0
      select case (this%fieldComponent)
      case (iCur)
         do i = istart, iend
         do j = jstart, jend
         do k = kstart, kend
            do field = iEx, iEz
               if (isWithinBounds(field, i, j, k, sinpml_fullsize)) then
                  if (isPEC(field, i, j, k, geometryMedia, registeredMedia)) then
                     count = count + 1
                     this%coords(:, count) = [i, j, k]
                  end if
               end if
            end do
         end do
         end do
         end do
      end select

   end subroutine get_measurements_coords

   subroutine save_current_data(this, step, fieldsReference, geometryMedia, registeredMedia, sinpml_fullsize)
      type(movie_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(fields_reference_t), intent(in) :: fieldsReference

      type(media_matrices_t), intent(in) :: geometryMedia
      type(MediaData_t), dimension(:) :: registeredMedia
      type(limit_t), dimension(:), intent(in)  :: sinpml_fullsize

      integer(kind=SINGLE) :: i, j, k, field
      integer(kind=SINGLE) :: istart, jstart, kstart, iend, jend, kend
      integer(kind=SINGLE) :: n

      istart = this%lowerBound%x
      jstart = this%lowerBound%y
      kstart = this%lowerBound%z

      iend = this%upperBound%x
      jend = this%upperBound%y
      kend = this%upperBound%z

      n = 0
      do i = istart, iend
      do j = jstart, jend
      do k = kstart, kend
         do field = iEx, iEz
            if (isWithinBounds(field, i, j, k, SINPML_fullsize)) then
            if (isPEC(field, i, j, k, geometryMedia, registeredMedia)) then
               n = n + 1
               call save_current_component()
            end if
            end if
         end do
      end do
      end do
      end do

      if (n < this%nMeasuredElements) call StopOnError(0, 0, "Missing measurment to update at movie probe")
   end subroutine save_current_data

   subroutine write_vtu_timestep(this, stepIndex, filename)
      use vtk_fortran
      implicit none

      type(movie_probe_output_t), intent(in) :: this
      integer, intent(in) :: stepIndex
      character(len=*), intent(in) :: filename

      type(vtk_file) :: vtkOutput
      integer :: ierr, npts, i
      real(kind=RKIND), allocatable :: x(:), y(:), z(:)
      real(kind=RKIND), allocatable :: Jx(:), Jy(:), Jz(:)

      npts = this%nMeasuredElements

      allocate (x(npts), y(npts), z(npts))
      do i = 1, npts
         x(i) = this%coords(1, i)
         y(i) = this%coords(2, i)
         z(i) = this%coords(3, i)
      end do

      allocate (Jx(npts), Jy(npts), Jz(npts))
      do i = 1, npts
         Jx(i) = this%xValueForTime(stepIndex, i)
         Jy(i) = this%yValueForTime(stepIndex, i)
         Jz(i) = this%zValueForTime(stepIndex, i)
      end do
      ierr = vtkOutput%initialize(format='ASCII', filename=trim(filename), mesh_topology='UnstructuredGrid')
      ierr = vtkOutput%xml_writer%write_geo(n=npts, x=x, y=y, z=z)
      ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='open')
      ierr = vtkOutput%xml_writer%write_dataarray(data_name='CurrentX', x=Jx)
      ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='close')
      ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='open')
      ierr = vtkOutput%xml_writer%write_dataarray(data_name='CurrentY', x=Jy)
      ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='close')
      ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='open')
      ierr = vtkOutput%xml_writer%write_dataarray(data_name='CurrentZ', x=Jz)
      ierr = vtkOutput%xml_writer%write_dataarray(location='node', action='close')
      ierr = vtkOutput%xml_writer%finalize()

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

end module mod_movieProbeOutput
