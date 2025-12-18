module mod_volumicProbeOutput
   use FDETYPES
   use mod_domain
   use mod_outputUtils
   implicit none
   private

   !===========================
   !  Public interface summary
   !===========================
   public :: init_volumic_probe_output
   public :: update_volumic_probe_output
   public :: flush_volumic_probe_output
   !===========================

   !===========================
   !  Private interface summary
   !===========================
   private :: isRelevantCell
   private :: isRelevantSurfaceCell
   private :: updateComplexComponent
   private :: count_relevant_geometries
   !===========================

contains

  subroutine init_volumic_probe_output(this, lowerBound, upperBound, field, domain, geometryMedia, registeredMedia, sinpml_fullsize, outputTypeExtension, mpidir, timeInterval)
      type(volumic_current_probe_t), intent(inout) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      type(MediaData_t), pointer, dimension(:) :: registeredMedia
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(limit_t), pointer, dimension(:), intent(in)  :: sinpml_fullsize

      real(kind=RKIND_tiempo), intent(in) :: timeInterval

      type(domain_t), intent(in) :: domain

      integer(kind=SINGLE) :: i, relevantGeometriesCount

      this%lowerBound = lowerBound
      this%upperBound = upperBound
      this%fieldComponent = field
      this%domain = domain
      this%path = get_output_path()

      relevantGeometriesCount = count_relevant_geometries(this, geometryMedia, registeredMedia, sinpml_fullsize)

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         allocate (this%timeStep(BuffObse))
         allocate (this%xValueForTime(BuffObse, relevantGeometriesCount))
         allocate (this%yValueForTime(BuffObse, relevantGeometriesCount))
         allocate (this%zValueForTime(BuffObse, relevantGeometriesCount))
         this%xValueForTime = 0.0_RKIND
         this%yValueForTime = 0.0_RKIND
         this%zValueForTime = 0.0_RKIND
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         this%nFreq = this%nFreq
         allocate (this%frequencySlice(this%nFreq))
         allocate (this%xValueForFreq(this%nFreq, relevantGeometriesCount))
         allocate (this%yValueForFreq(this%nFreq, relevantGeometriesCount))
         allocate (this%zValueForFreq(this%nFreq, relevantGeometriesCount))
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

   end subroutine init_volumic_probe_output

   function count_relevant_geometries(this, geometryMedia, registeredMedia, sinpml_fullsize) result(n)
      type(volumic_current_probe_t), intent(in) :: this
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t), pointer, dimension(:), intent(in) :: registeredMedia
      type(limit_t), pointer, dimension(:), intent(in) :: sinpml_fullsize
      integer(kind=SINGLE) :: i, j, k, field
      integer(kind=SINGLE) :: n

      n = 0_SINGLE
      do i = this%lowerBound%x, this%upperBound%x
      do j = this%lowerBound%y, this%upperBound%y
      do k = this%lowerBound%z, this%upperBound%z
         do field = iEx, iEz
            if (isRelevantCell(field, i, j, k, geometryMedia, registeredMedia, sinpml_fullsize)) then
               n = n + 1
            end if
         end do
         do field = iHx, iHz
            if (isRelevantSurfaceCell(field, i, j, k, this%fieldComponent, geometryMedia, registeredMedia, sinpml_fullsize)) then
               n = n + 1
            end if
         end do
      end do
      end do
      end do
   end function

   subroutine update_volumic_probe_output(this, step, geometryMedia, registeredMedia, sinpml_fullsize, fieldsReference)
      type(volumic_current_probe_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step

      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t), pointer, dimension(:) :: registeredMedia
      type(limit_t), pointer, dimension(:), intent(in)  :: sinpml_fullsize
      type(fields_reference_t), intent(in) :: fieldsReference

      integer(kind=SINGLE) :: Efield, Hfield, i, j, k, conta
      integer(kind=SINGLE) :: i1, i2, j1, j2, k1, k2

      i1 = this%lowerBound%x
      j1 = this%lowerBound%y
      k1 = this%lowerBound%z

      i2 = this%upperBound%x
      j2 = this%upperBound%y
      k2 = this%upperBound%z

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         conta = 0
         this%serializedTimeSize = this%serializedTimeSize + 1
         do i = i1, i2
         do j = j1, j2
         do k = k1, k2
         do Efield = iEx, iEz
            if (isRelevantCell(Efield, i, j, k, geometryMedia, registeredMedia, sinpml_fullsize)) then
               conta = conta + 1
               call save_current(this, Efield, i, j, k, conta, fieldsReference)
            end if
         end do
         do Hfield = iHx, iHz
            if (isRelevantSurfaceCell(Hfield, i, j, k, this%fieldComponent, geometryMedia, registeredMedia, sinpml_fullsize)) then
               conta = conta + 1
               call save_current_surfaces(this, Hfield, i, j, k, conta, fieldsReference)
            end if
         end do
         end do
         end do
         end do
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         conta = 0
         do i = i1, i2
         do j = j1, j2
         do k = k1, k2
         do Efield = iEx, iEz
            if (isRelevantCell(Efield, i, j, k, geometryMedia, registeredMedia, sinpml_fullsize)) then
               conta = conta + 1
               call update_current(this, Efield, i, j, k, conta, fieldsReference, step)
            end if
         end do
         do Hfield = iHx, iHz
            if (isRelevantSurfaceCell(Hfield, i, j, k, this%fieldComponent, geometryMedia, registeredMedia, sinpml_fullsize)) then
               conta = conta + 1
               call update_current_surfaces(this, Hfield, i, j, k, conta, fieldsReference, step)
            end if
         end do
         end do
         end do
         end do
      end if
   contains
      subroutine save_current(this, Efield, i, j, k, conta, field_reference)
         type(fields_reference_t), intent(in) :: field_reference
         type(volumic_current_probe_t), intent(inout) :: this
         integer(kind=SINGLE), intent(in) :: Efield, i, j, k, conta

         real(kind=RKIND) :: jdir

         jdir = computeJ(EField, i, j, k, field_reference)
         this%xValueForTime(this%serializedTimeSize, conta) = merge(jdir, 0.0_RKIND, Efield == iEx)
         this%yValueForTime(this%serializedTimeSize, conta) = merge(jdir, 0.0_RKIND, Efield == iEy)
         this%zValueForTime(this%serializedTimeSize, conta) = merge(jdir, 0.0_RKIND, Efield == iEz)
      end subroutine save_current

      subroutine save_current_surfaces(this, Hfield, i, j, k, conta, field_reference)
         implicit none
         type(fields_reference_t), intent(in) :: field_reference
         type(volumic_current_probe_t), intent(inout) :: this
         integer(kind=SINGLE), intent(in) :: Hfield, i, j, k, conta

         real(kind=RKIND) :: jdir1, jdir2
         jdir1 = computeJ1(HField, i, j, k, field_reference)
         jdir2 = computeJ2(HField, i, j, k, field_reference)

         this%xValueForTime(this%serializedTimeSize, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHz), Hfield == iHx)
         this%yValueForTime(this%serializedTimeSize, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHx), Hfield == iHy)
         this%zValueForTime(this%serializedTimeSize, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHy), Hfield == iHz)
      end subroutine save_current_surfaces

      subroutine update_current(this, Efield, i, j, k, conta, field_reference, step)
         integer(kind=SINGLE), intent(in) :: Efield, i, j, k, conta
         type(volumic_current_probe_t), intent(inout) :: this
         type(fields_reference_t), intent(in) :: field_reference
         real(kind=RKIND_tiempo), intent(in) :: step

         integer(kind=SINGLE) :: freqIdx
         real(kind=RKIND) :: jdir

         jdir = computeJ(Efield, i, j, k, field_reference)
         do freqIdx = 1, this%nFreq
            call updateComplexComponent(iEx, EField, this%xValueForFreq(freqIdx, conta), jdir, this%auxExp_E(freqIdx)**step)
            call updateComplexComponent(iEy, EField, this%yValueForFreq(freqIdx, conta), jdir, this%auxExp_E(freqIdx)**step)
            call updateComplexComponent(iEz, EField, this%zValueForFreq(freqIdx, conta), jdir, this%auxExp_E(freqIdx)**step)
         end do
      end subroutine update_current

      subroutine update_current_surfaces(this, Hfield, i, j, k, conta, field_reference, step)
         integer(kind=SINGLE), intent(in) :: Hfield, i, j, k, conta
         type(volumic_current_probe_t), intent(inout) :: this
         type(fields_reference_t), intent(in) :: field_reference
         real(kind=RKIND_tiempo), intent(in) :: step

         integer(kind=SINGLE) :: freqIdx
         real(kind=RKIND) :: jdir, jdir1, jdir2

         jdir1 = computeJ1(HField, i, j, k, field_reference)
         jdir2 = computeJ2(HField, i, j, k, field_reference)
         do freqIdx = 1, this%nFreq
            jdir = merge(jdir1, jdir2, HField == iHz)
            call updateComplexComponent(iHx, Hfield, this%xValueForFreq(freqIdx, conta), jdir, this%auxExp_H(freqIdx)**step)

            jdir = merge(jdir1, jdir2, HField == iHx)
            call updateComplexComponent(iHy, Hfield, this%yValueForFreq(freqIdx, conta), jdir, this%auxExp_H(freqIdx)**step)

            jdir = merge(jdir1, jdir2, HField == iHy)
            call updateComplexComponent(iHz, Hfield, this%zValueForFreq(freqIdx, conta), jdir, this%auxExp_H(freqIdx)**step)
         end do
      end subroutine update_current_surfaces

   end subroutine update_volumic_probe_output

   subroutine flush_volumic_probe_output
      !!TODO
   end subroutine flush_volumic_probe_output

   logical function isRelevantCell(Efield, I, J, K, geometryMedia, registeredMedia, sinpml_fullsize)
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t), pointer, dimension(:), intent(in) :: registeredMedia
      type(limit_t), pointer, dimension(:), intent(in) :: sinpml_fullsize
      integer(kind=SINGLE), intent(in) :: Efield, I, J, K
      isRelevantCell = .false.

      if (isWithinBounds(Efield, I, J, K, sinpml_fullsize)) then
         if (isThinWire(Efield, I, J, K, geometryMedia, registeredMedia)) then
            isRelevantCell = .true.
         end if
         if (.NOT. isMediaVacuum(Efield, I, J, K, geometryMedia)) then
            if (.NOT. isSplitOrAdvanced(Efield, I, J, K, geometryMedia, registeredMedia)) then
               isRelevantCell = .true.
            end if
         end if
      end if

   end function

   logical function isRelevantSurfaceCell(field, i, j, k, outputType, geometryMedia, registeredMedia, sinpml_fullsize)
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t), pointer, dimension(:), intent(in) :: registeredMedia
      type(limit_t), pointer, dimension(:), intent(in) :: sinpml_fullsize
      integer(kind=SINGLE), intent(in) :: field, i, j, k, outputType

      isRelevantSurfaceCell = .false.
      if (isWithinBounds(field, i, j, k, sinpml_fullsize)) then
         isRelevantSurfaceCell = isPEC(field, i, j, k, geometryMedia, registeredMedia)
      end if

   end function

   subroutine updateComplexComponent(direction, fieldIndex, valorComplex, jdir, auxExp)
      integer, intent(in) :: direction, fieldIndex
      complex(kind=CKIND), intent(inout) :: valorComplex
      complex(kind=CKIND), intent(in) :: auxExp
      real(kind=RKIND), intent(in) :: jdir

      complex(kind=CKIND) :: z_cplx = (0.0_RKIND, 0.0_RKIND)

      valorComplex = merge(valorComplex + auxExp*jdir, z_cplx, fieldIndex == direction)
   end subroutine updateComplexComponent

end module mod_volumicProbeOutput
