module mod_frequencySliceProbeOutput
   use FDETYPES
   use Report
   use outputTypes
   use mod_outputUtils
   implicit none
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: init_frequency_slice_probe_output
   public :: update_frequency_slice_probe_output
   public :: flush_frequency_slice_probe_output
   !===========================

   !===========================
   !  Private interface summary
   !===========================
   private :: get_measurements_coords
   private :: save_current_data
   private :: write_vtu_frequency_slice
   private :: update_pvd
   !===========================

contains

   subroutine init_frequency_slice_probe_output(this, lowerBound, upperBound, field, domain, geometryMedia, registeredMedia, sinpml_fullsize, outputTypeExtension, mpidir, timeInterval)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      type(MediaData_t), pointer, dimension(:), intent(in) :: registeredMedia
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(limit_t), pointer, dimension(:), intent(in)  :: sinpml_fullsize

      type(domain_t), intent(in) :: domain
      real(kind=RKIND_tiempo), intent(in) :: timeInterval
      integer :: i

      if (domain%domainType /= FREQUENCY_DOMAIN) call StopOnError(0, 0, "Unexpected domain type for frequency_slice probe")

      this%lowerBound = lowerBound
      this%upperBound = upperBound
      this%fieldComponent = field !This can refer to field or currentDensity
      this%domain = domain
      this%path = get_output_path()
      call get_measurements_coords(this, geometryMedia, registeredMedia, sinpml_fullsize)

      this%nFreq = domain%fnum
      allocate (this%frequencySlice(this%nFreq))
      allocate (this%xValueForFreq(this%nFreq, this%nMeasuredElements))
      allocate (this%yValueForFreq(this%nFreq, this%nMeasuredElements))
      allocate (this%zValueForFreq(this%nFreq, this%nMeasuredElements))
      do i = 1, this%nFreq
         call init_frequency_slice(this%frequencySlice, this%domain)
      end do
      this%xValueForFreq = (0.0_RKIND, 0.0_RKIND)
      this%yValueForFreq = (0.0_RKIND, 0.0_RKIND)
      this%zValueForFreq = (0.0_RKIND, 0.0_RKIND)

      allocate (this%auxExp_E(this%nFreq))
      allocate (this%auxExp_H(this%nFreq))
      do i = 1, this%nFreq
         this%auxExp_E(i) = timeInterval*(1.0E0_RKIND, 0.0E0_RKIND)*Exp(mcpi2*this%frequencySlice(i))   !el dt deberia ser algun tipo de promedio
         this%auxExp_H(i) = this%auxExp_E(i)*Exp(mcpi2*this%frequencySlice(i)*timeInterval*0.5_RKIND)
      end do

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

   end subroutine init_frequency_slice_probe_output

   subroutine update_frequency_slice_probe_output(this, step, geometryMedia, registeredMedia, sinpml_fullsize, fieldsReference)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step

      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t), pointer, dimension(:), intent(in) :: registeredMedia
      type(limit_t), pointer, dimension(:), intent(in)  :: sinpml_fullsize
      type(fields_reference_t), pointer, intent(in) :: fieldsReference

      select case (this%fieldComponent)
      case (iCur)
         call save_current_data(this, step, fieldsReference, geometryMedia, registeredMedia, sinpml_fullsize)
      end select
   end subroutine update_frequency_slice_probe_output

   subroutine flush_frequency_slice_probe_output(this)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      integer :: status, i

      do i = 1, this%nFreq
         call update_pvd(this, i, this%PDVUnit)
      end do

   end subroutine flush_frequency_slice_probe_output

   subroutine get_measurements_coords(this, geometryMedia, registeredMedia, sinpml_fullsize)
      type(frequency_slice_probe_output_t), intent(inout) :: this
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t), pointer, dimension(:) :: registeredMedia
      type(limit_t), pointer, dimension(:), intent(in) :: sinpml_fullsize

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
      type(frequency_slice_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(fields_reference_t), pointer, intent(in) :: fieldsReference

      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t), pointer, dimension(:) :: registeredMedia
      type(limit_t), pointer, dimension(:), intent(in)  :: sinpml_fullsize

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

      if (n < this%nMeasuredElements) call StopOnError(0, 0, "Missing measurment to update at frequency_slice probe")
   contains

      subroutine save_current_component()
         real(kind=RKIND) :: jdir
         integer :: freqIdx
         jdir = computeJ(field, i, j, k, fieldsReference)

         do freqIdx = 1, this%nFreq
            call updateComplexComponent(iEx, field, this%xValueForFreq(freqIdx, n), jdir, this%auxExp_E(freqIdx)**step)
            call updateComplexComponent(iEy, field, this%yValueForFreq(freqIdx, n), jdir, this%auxExp_E(freqIdx)**step)
            call updateComplexComponent(iEz, field, this%zValueForFreq(freqIdx, n), jdir, this%auxExp_E(freqIdx)**step)
         end do
      end subroutine save_current_component

      subroutine updateComplexComponent(direction, fieldIndex, valorComplex, jdir, auxExp)
         integer, intent(in) :: direction, fieldIndex
         complex(kind=CKIND), intent(inout) :: valorComplex
         complex(kind=CKIND), intent(in) :: auxExp
         real(kind=RKIND), intent(in) :: jdir

         complex(kind=CKIND) :: z_cplx = (0.0_RKIND, 0.0_RKIND)

         valorComplex = merge(valorComplex + auxExp*jdir, z_cplx, fieldIndex == direction)
      end subroutine updateComplexComponent
   end subroutine save_current_data

   subroutine write_vtu_frequency_slice(this, freq, filename)
      use vtk_fortran
      implicit none

      type(frequency_slice_probe_output_t), intent(in) :: this
      integer, intent(in) :: freq
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
         Jx(i) = this%xValueForFreq(freq, i)
         Jy(i) = this%yValueForFreq(freq, i)
         Jz(i) = this%zValueForFreq(freq, i)
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

   end subroutine write_vtu_frequency_slice

   subroutine update_pvd(this, freq, unitPVD)
      implicit none
      type(frequency_slice_probe_output_t), intent(in) :: this
      integer, intent(in) :: freq
      integer, intent(in) :: unitPVD
      character(len=64) :: ts
      character(len=256) :: filename

      ! Generamos nombre del archivo VTU para este timestep
      write (filename, '(A,A,I4.4,A)') trim(this%path), '_fq', freq, '.vtu'

      ! Escribimos el VTU correspondiente
      call write_vtu_frequency_slice(this, freq, filename)

      ! Añadimos entrada en el PVD
      write (ts, '(ES16.8)') this%frequencySlice(freq)
      write (unitPVD, '(A)') '    <DataSet timestep="'//trim(ts)// &
         '" group="" part="0" file="'//trim(filename)//'"/>'
   end subroutine update_pvd

end module mod_frequencySliceProbeOutput
