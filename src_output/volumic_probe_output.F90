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
      real(kind=RKIND_tiempo), dimension(:) :: timeStep
      real(kind=RKIND), dimension(:, :) :: xValueForTime
      real(kind=RKIND), dimension(:, :) :: yValueForTime
      real(kind=RKIND), dimension(:, :) :: zValueForTime

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
      type(media_matrices_t), intent(in) :: media
      type(limit_t), dimension(1:6), intent(in)  :: sinpml_fullsize

      type(domain_t), intent(in) :: domain

      integer(kind=SINGLE) :: i

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
         allocate (timeStep(BuffObse, totalPecSurfaces) &
                   xValueForTime(BuffObse, totalPecSurfaces) &
                   yValueForTime(BuffObse, totalPecSurfaces) &
                   zValueForTime(BuffObse, totalPecSurfaces))
         xValueForTime = 0.0_RKIND
         yValueForTime = 0.0_RKIND
         zValueForTime = 0.0_RKIND
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         this%nFreq = this%domain%fnum
         allocate (this%frequencySlice(this%domain%fnum))
         allocate (this%xValueForFreq(this%domain%fnum, totalPecSurfaces) &
                   this%yValueForFreq(this%domain%fnum, totalPecSurfaces) &
                   this%zValueForFreq(this%domain%fnum, totalPecSurfaces))
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
         integer(kind=SINGLE) :: n = 0_SINGLE
         do concurrent(i=icoord:i2coord, j=jcoord:j2coord, k=kcoord:k2coord, field= iEx:iEz) !Ejecuta todas las combinaciones de (i, j, k, field)
            if (isWithinBounds(field, iii, jjj, kkk)) then
               if (isThinWire(field, iii, jjj, kkk)) then
                  n = n + 1
               end if
               if (.not. isMediaVacuum(field, iii, jjj, kkk, media) .and. .not. isSplitOrAdvanced(field, iii, jjj, kkk)) then
                  n = n + 1
               end if
               if (isPECorSurface(field, iii, jjj, kkk, media, simulationMedia) .or. field == blockCurrent(field)) then
                  n = n + 1
               end if
            end if
         end do

      end function count_pec_surface
   end subroutine init_volumic_probe_output

   subroutine update_volumic_probe_output(this, step)
      type(volumic_current_probe_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step

      integer(kind=SINGLE) :: Efield, iii, jjj, kkk
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
               call save_current(this, Efield, iii, jjj, kkk, conta)
            end if
         end do
         end do
         end do
         end do
      end if
      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
      end if
   contains
      LOGICAL FUNCTION isRelevantCell(Efield, I, J, K)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: Efield, I, J, K

         isRelevantCell = isWithinBounds(Efield, I, J, K) .AND. &
                          (isThinWire(Efield, I, J, K) .OR. &
                           (.NOT. isMediaVacuum(Efield, I, J, K) .AND. &
                            .NOT. isSplitOrAdvanced(Efield, I, J, K)))

      END FUNCTION isRelevantCell

      subroutine save_current(this, Efield, iii, jjj, kkk, conta)
         type(volumic_current_probe_t), intent(inout) :: this
         integer(kind=SINGLE), intent(in) :: Efield, iii, jjj, kkk, conta
         jdir = computeJ(EField, iii, jjj, kkk)
         this%xValueForTime(this%serializedTimeSize, conta) = merge(jdir, 0.0_RKIND, Efield == iEx)
         this%yValueForTime(this%serializedTimeSize, conta) = merge(jdir, 0.0_RKIND, Efield == iEy)
         this%zValueForTime(this%serializedTimeSize, conta) = merge(jdir, 0.0_RKIND, Efield == iEz)
      end subroutine save_current
   end subroutine update_volumic_probe_output

end module mod_volumicProbe
