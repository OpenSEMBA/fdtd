module mod_volumicProbe
   use FDETYPES
   use mod_domain
   use mod_outputUtils

   implicit none
   private :: isRelevantCell, isRelevantSurfaceCell

   type volumic_current_probe_t
      integer(kind=SINGLE) :: columnas = 4_SINGLE !reference and current components
      type(domain_t) :: domain
      integer(kind=SINGLE) :: xCoord, yCoord, zCoord
      integer(kind=SINGLE) :: x2Coord, y2Coord, z2Coord
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent

      !Intent storage order:
      !(:) == (timeinstance) => timeValue
      !(:,:) == (timeInstance, componentId) => escalar

      !Time Domain (requires first allocation)
      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(:), allocatable :: timeStep
      real(kind=RKIND), dimension(:, :), allocatable :: xValueForTime
      real(kind=RKIND), dimension(:, :), allocatable :: yValueForTime
      real(kind=RKIND), dimension(:, :), allocatable :: zValueForTime

      !Intent storage order:
      !(:) == (frquencyinstance) => timeValue
      !(:,:) == (frquencyinstance, componentId) => escalar

      !Frequency Domain (requires first allocation)
      integer(kind=SINGLE) :: nFreq = 0_SINGLE
      real(kind=RKIND), dimension(:), allocatable :: frequencySlice
      real(kind=CKIND), dimension(:, :), allocatable :: xValueForFreq
      real(kind=CKIND), dimension(:, :), allocatable :: yValueForFreq
      real(kind=CKIND), dimension(:, :), allocatable :: zValueForFreq
   end type volumic_current_probe_t

contains

  subroutine init_volumic_probe_output(this, iCoord, jCoord, kCoord, i2Coord, j2Coord, k2Coord, field, domain, geometryMedia, registeredMedia, sinpml_fullsize, outputTypeExtension, mpidir)
      type(volumic_current_probe_t), intent(inout) :: this
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord
      integer(kind=SINGLE), intent(in) :: i2Coord, j2Coord, k2Coord
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      type(MediaData_t), pointer, dimension(:) :: registeredMedia
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(limit_t), pointer, dimension(:), intent(in)  :: sinpml_fullsize

      type(domain_t), intent(in) :: domain

      integer(kind=SINGLE) :: i, totalPecSurfaces

      this%xCoord = iCoord
      this%yCoord = jCoord
      this%zCoord = kCoord

      this%x2Coord = i2Coord
      this%y2Coord = j2Coord
      this%z2Coord = k2Coord

      this%fieldComponent = field

      this%domain = domain
      this%path = get_output_path()

      totalPecSurfaces = count_pec_surfaces(this, geometryMedia, registeredMedia, sinpml_fullsize)

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         allocate (this%timeStep(BuffObse))
         allocate (this%xValueForTime(BuffObse, totalPecSurfaces))
         allocate (this%yValueForTime(BuffObse, totalPecSurfaces))
         allocate (this%zValueForTime(BuffObse, totalPecSurfaces))
         this%xValueForTime = 0.0_RKIND
         this%yValueForTime = 0.0_RKIND
         this%zValueForTime = 0.0_RKIND
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         this%nFreq = this%domain%fnum
         allocate (this%frequencySlice(this%domain%fnum))
         allocate (this%xValueForFreq(this%domain%fnum, totalPecSurfaces))
         allocate (this%yValueForFreq(this%domain%fnum, totalPecSurfaces))
         allocate (this%zValueForFreq(this%domain%fnum, totalPecSurfaces))
         do i = 1, this%nFreq
            call init_frequency_slice(this%frequencySlice, this%domain)
         end do
         this%xValueForFreq = (0.0_RKIND, 0.0_RKIND)
         this%yValueForFreq = (0.0_RKIND, 0.0_RKIND)
         this%zValueForFreq = (0.0_RKIND, 0.0_RKIND)
      end if

   contains
      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_probe_bounds_coords_extension(iCoord, jCoord, kCoord, i2Coord, j2Coord, k2Coord, mpidir)
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

   end subroutine init_volumic_probe_output

   function count_pec_surfaces(this, geometryMedia, registeredMedia, sinpml_fullsize) result(n)
      type(volumic_current_probe_t), intent(in) :: this
      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t), pointer, dimension(:), intent(in) :: registeredMedia
      type(limit_t), pointer, dimension(:), intent(in) :: sinpml_fullsize
      integer(kind=SINGLE) :: i, j, k, field
      integer(kind=SINGLE) :: n
      n = 0_SINGLE
      do i = this%xCoord, this%x2Coord
      do j = this%yCoord, this%y2Coord
      do k = this%zCoord, this%z2Coord
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
   end function count_pec_surfaces

   subroutine update_volumic_probe_output(this, step, geometryMedia, registeredMedia, sinpml_fullsize, fieldsReference)
      type(volumic_current_probe_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step

      type(media_matrices_t), pointer, intent(in) :: geometryMedia
      type(MediaData_t), pointer, dimension(:) :: registeredMedia
      type(limit_t), pointer, dimension(:), intent(in)  :: sinpml_fullsize
      type(fields_reference_t), pointer, intent(in) :: fieldsReference

      integer(kind=SINGLE) :: Efield, Hfield, i, j, k, conta
      integer(kind=SINGLE) :: i1, i2, j1, j2, k1, k2

      conta = 0

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         this%serializedTimeSize = this%serializedTimeSize + 1
         do k = k1, k2
         do j = j1, j2
         do i = i1, i2
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
      end if
   contains
      subroutine save_current(this, Efield, i, j, k, conta, field_reference)
         type(fields_reference_t), pointer, intent(in) :: field_reference
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
         type(fields_reference_t), pointer, intent(in) :: field_reference
         type(volumic_current_probe_t), intent(inout) :: this
         integer(kind=SINGLE), intent(in) :: Hfield, i, j, k, conta

         real(kind=RKIND) :: jdir1, jdir2
         jdir1 = computeJ1(HField, i, j, k, field_reference)
         jdir2 = computeJ2(HField, i, j, k, field_reference)

         this%xValueForTime(this%serializedTimeSize, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHz), Hfield == iHx)
         this%yValueForTime(this%serializedTimeSize, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHx), Hfield == iHy)
         this%zValueForTime(this%serializedTimeSize, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHy), Hfield == iHz)
      end subroutine save_current_surfaces
   end subroutine update_volumic_probe_output

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

end module mod_volumicProbe
