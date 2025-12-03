module mod_volumicProbe
   use FDETYPES
   use mod_domain
   use mod_outputUtils

   implicit none
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

  subroutine init_volumic_probe_output(this, iCoord, jCoord, kCoord, i2Coord, j2Coord, k2Coord, field, domain, media, simulationMedia, sinpml_fullsize, outputTypeExtension, mpidir)
      type(volumic_current_probe_t), intent(out) :: this
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord
      integer(kind=SINGLE), intent(in) :: i2Coord, j2Coord, k2Coord
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      type(MediaData_t), pointer, dimension(:) :: simulationMedia
      type(media_matrices_t), pointer, intent(in) :: media
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

      totalPecSurfaces = count_pec_surfaces()

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
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

      function count_pec_surfaces() result(n)
         integer(kind=SINGLE) :: i, j, k, field
         integer(kind=SINGLE) :: n, iii, jjj, kkk
         n = 0_SINGLE
         do i = icoord,i2coord
         do j = jcoord,j2coord
         do k = kcoord,k2coord
         do field = iEx,iEz
            if (isWithinBounds(field, iii, jjj, kkk, sinpml_fullsize)) then
               if (isThinWire(field, iii, jjj, kkk, simulationMedia, media)) then
                  n = n + 1
               end if
               if (.not. isMediaVacuum(field, iii, jjj, kkk, media) .and. .not. isSplitOrAdvanced(field, iii, jjj, kkk, media, simulationMedia)) then
                  n = n + 1
               end if
               if (isPECorSurface(field, iii, jjj, kkk, media, simulationMedia) .or. field == blockCurrent(field)) then
                  n = n + 1
               end if
            end if
         end do
         end do
         end do
         end do

      end function count_pec_surfaces
   end subroutine init_volumic_probe_output

   subroutine update_volumic_probe_output(this, step, media, simulationMedia, sinpml_fullsize, fieldsReference)
      type(volumic_current_probe_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step

      type(media_matrices_t), pointer, intent(in) :: media
      type(MediaData_t), pointer, dimension(:) :: simulationMedia
      type(limit_t), pointer, dimension(:), intent(in)  :: sinpml_fullsize
      type(fields_reference_t), pointer, intent(in) :: fieldsReference

      integer(kind=SINGLE) :: Efield, Hfield, iii, jjj, kkk
      integer(kind=SINGLE) :: i1, i2, j1, j2, k1, k2, conta

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         this%serializedTimeSize = this%serializedTimeSize + 1
         conta = 0
         do KKK = k1, k2
         do JJJ = j1, j2
         do III = i1, i2
         do Efield = iEx, iEz
            if (isRelevantCell(Efield, iii, jjj, kkk)) then
               conta = conta + 1
               call save_current(this, Efield, iii, jjj, kkk, conta, fieldsReference)

            end if
         end do
         do Hfield = iHx, iHz
            if (isRelevantSurfaceCell(Hfield, iii, jjj, kkk, this%fieldComponent)) then
               conta = conta + 1
               call save_current_surfaces(this, Hfield, iii, jjj, kkk, conta, fieldsReference)
            end if
         end do
         end do
         end do
         end do
      end if
      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
      end if
   contains
      logical function isRelevantCell(Efield, I, J, K)
         integer(kind=SINGLE), intent(in) :: Efield, I, J, K

         if (isWithinBounds(Efield, I, J, K, sinpml_fullsize)) then
            isRelevantCell = isThinWire(Efield, I, J, K, simulationMedia, media) .OR. &
                (.NOT. isMediaVacuum(Efield, I, J, K, media) .AND. .NOT. isSplitOrAdvanced(Efield, I, J, K, media, simulationMedia))
         else
            isRelevantCell = .false.
         end if

      END FUNCTION isRelevantCell

      logical function isRelevantSurfaceCell(Hfield, I, J, K, outputType)
         integer(kind=SINGLE), intent(in) :: Hfield, I, J, K, outputType

         if (isWithinBounds(Hfield, I, J, K, sinpml_fullsize)) then
       isRelevantSurfaceCell = isPECorSurface(Hfield, iii, jjj, kkk, media, simulationMedia) .or. outputType == blockCurrent(Hfield)
         else
            isRelevantSurfaceCell = .false.
         end if
      end function

      subroutine save_current(this, Efield, iii, jjj, kkk, conta, field_reference)
         type(fields_reference_t), pointer, intent(in) :: field_reference
         type(volumic_current_probe_t), intent(inout) :: this
         integer(kind=SINGLE), intent(in) :: Efield, iii, jjj, kkk, conta

         real(kind=RKIND) :: jdir

         jdir = computeJ(EField, iii, jjj, kkk, field_reference)
         this%xValueForTime(this%serializedTimeSize, conta) = merge(jdir, 0.0_RKIND, Efield == iEx)
         this%yValueForTime(this%serializedTimeSize, conta) = merge(jdir, 0.0_RKIND, Efield == iEy)
         this%zValueForTime(this%serializedTimeSize, conta) = merge(jdir, 0.0_RKIND, Efield == iEz)
      end subroutine save_current

      subroutine save_current_surfaces(this, Hfield, iii, jjj, kkk, conta, field_reference)
         implicit none
         type(fields_reference_t), pointer, intent(in) :: field_reference
         type(volumic_current_probe_t), intent(inout) :: this
         integer(kind=SINGLE), intent(in) :: Hfield, iii, jjj, kkk, conta

         real(kind=RKIND) :: jdir1, jdir2
         jdir1 = computeJ1(HField, iii, jjj, kkk, field_reference)
         jdir2 = computeJ2(HField, iii, jjj, kkk, field_reference)

         this%xValueForTime(this%serializedTimeSize, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHz), Hfield == iHx)
         this%yValueForTime(this%serializedTimeSize, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHx), Hfield == iHy)
         this%zValueForTime(this%serializedTimeSize, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHy), Hfield == iHz)
      end subroutine save_current_surfaces
   end subroutine update_volumic_probe_output

end module mod_volumicProbe
