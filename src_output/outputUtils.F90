module mod_outputUtils
   use FDETYPES
   use mod_domain
   use report
   implicit none
   character(len=4), parameter :: datFileExtension = '.dat', timeExtension = 'tm', frequencyExtension = 'fq'
   integer(kind=SINGLE), parameter :: FILE_UNIT = 400

   type field_data_t
      real(kind=RKIND), pointer, dimension(:, :, :) :: x, y, z
      real(kind=RKIND), pointer, dimension(:) :: deltaX, deltaY, deltaZ
   end type field_data_t

contains

   function get_probe_coords_extension(iCoord, jCoord, kCoord, mpidir) result(ext)
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord, mpidir
      character(len=BUFSIZE) :: ext
      character(len=BUFSIZE)  ::  chari, charj, chark

      write (chari, '(i7)') iCoord
      write (charj, '(i7)') jCoord
      write (chark, '(i7)') kCoord

#if CompileWithMPI
      if (mpidir == 3) then
         ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
      elseif (mpidir == 2) then
         ext = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))
      elseif (mpidir == 1) then
         ext = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))
      else
         call stoponerror('Buggy error in mpidir. ')
      end if
#else
      ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
#endif

      return
   end function get_probe_coords_extension

   function get_probe_bounds_coords_extension(iCoord, jCoord, kCoord, i2Coord, j2Coord, k2Coord, mpidir) result(ext)
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord, i2Coord, j2Coord, k2Coord, mpidir
      character(len=BUFSIZE) :: ext
      character(len=BUFSIZE)  ::  chari, charj, chark, chari2, charj2, chark2

      write (chari, '(i7)') iCoord
      write (charj, '(i7)') jCoord
      write (chark, '(i7)') kCoord

      write (chari2, '(i7)') i2Coord
      write (charj2, '(i7)') j2Coord
      write (chark2, '(i7)') k2Coord

#if CompileWithMPI
      if (mpidir == 3) then
         ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'// &
               trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
      elseif (mpidir == 2) then
         ext = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))//'__'// &
               trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
      elseif (mpidir == 1) then
         ext = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))//'__'// &
               trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
      else
         call stoponerror('Buggy error in mpidir. ')
      end if
#else
      ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'// &
            trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
#endif

      return
   end function get_probe_bounds_coords_extension

   function get_prefix_extension(field, mpidir) result(prefixExtension)
      integer(kind=SINGLE), intent(in)  ::  field, mpidir
      character(len=BUFSIZE)  ::  prefixExtension

#if CompileWithMPI
      prefixExtension = get_rotated_prefix(field, mpidir)
#else
      prefixExtension = prefix(field)
#endif
   end function get_prefix_extension

   function get_rotated_prefix(field, mpidir) result(prefixExtension)
      integer(kind=SINGLE), intent(in)  ::  field, mpidir
      character(len=BUFSIZE)  ::  prefixExtension
      if (mpidir == 3) then
         select case (field)
         case (iEx); prefixExtension = prefix(iEx)
         case (iEy); prefixExtension = prefix(iEy)
         case (iEz); prefixExtension = prefix(iEz)
         case (iJx); prefixExtension = prefix(iJx)
         case (iJy); prefixExtension = prefix(iJy)
         case (iJz); prefixExtension = prefix(iJz)
         case (iQx); prefixExtension = prefix(iQx)
         case (iQy); prefixExtension = prefix(iQy)
         case (iQz); prefixExtension = prefix(iQz)
         case (iVx); prefixExtension = prefix(iVx)
         case (iVy); prefixExtension = prefix(iVy)
         case (iVz); prefixExtension = prefix(iVz)
         case (iHx); prefixExtension = prefix(iHx)
         case (iHy); prefixExtension = prefix(iHy)
         case (iHz); prefixExtension = prefix(iHz)
         case (iBloqueJx); prefix_field = prefix(iBloqueJx)
         case (iBloqueJy); prefix_field = prefix(iBloqueJy)
         case (iBloqueJz); prefix_field = prefix(iBloqueJz)
         case (iBloqueMx); prefix_field = prefix(iBloqueMx)
         case (iBloqueMy); prefix_field = prefix(iBloqueMy)
         case (iBloqueMz); prefix_field = prefix(iBloqueMz)
         case default; prefixExtension = prefix(field)
         end select
      elseif (mpidir == 2) then
         select case (field)
         case (iEx); prefixExtension = prefix(iEz)
         case (iEy); prefixExtension = prefix(iEx)
         case (iEz); prefixExtension = prefix(iEy)
         case (iJx); prefixExtension = prefix(iJz)
         case (iJy); prefixExtension = prefix(iJx)
         case (iJz); prefixExtension = prefix(iJy)
         case (iQx); prefixExtension = prefix(iQz)
         case (iQy); prefixExtension = prefix(iQx)
         case (iQz); prefixExtension = prefix(iQy)
         case (iVx); prefixExtension = prefix(iVz)
         case (iVy); prefixExtension = prefix(iVx)
         case (iVz); prefixExtension = prefix(iVy)
         case (iHx); prefixExtension = prefix(iHz)
         case (iHy); prefixExtension = prefix(iHx)
         case (iHz); prefixExtension = prefix(iHy)
         case (iBloqueJx); prefix_field = prefix(iBloqueJz)
         case (iBloqueJy); prefix_field = prefix(iBloqueJx)
         case (iBloqueJz); prefix_field = prefix(iBloqueJy)
         case (iBloqueMx); prefix_field = prefix(iBloqueMz)
         case (iBloqueMy); prefix_field = prefix(iBloqueMx)
         case (iBloqueMz); prefix_field = prefix(iBloqueMy)
         case default; prefixExtension = prefix(field)
         end select
      elseif (mpidir == 1) then
         select case (field)
         case (iEx); prefixExtension = prefix(iEy)
         case (iEy); prefixExtension = prefix(iEz)
         case (iEz); prefixExtension = prefix(iEx)
         case (iJx); prefixExtension = prefix(iJy)
         case (iJy); prefixExtension = prefix(iJz)
         case (iJz); prefixExtension = prefix(iJx)
         case (iQx); prefixExtension = prefix(iQy)
         case (iQy); prefixExtension = prefix(iQz)
         case (iQz); prefixExtension = prefix(iQx)
         case (iVx); prefixExtension = prefix(iVy)
         case (iVy); prefixExtension = prefix(iVz)
         case (iVz); prefixExtension = prefix(iVx)
         case (iHx); prefixExtension = prefix(iHy)
         case (iHy); prefixExtension = prefix(iHz)
         case (iHz); prefixExtension = prefix(iHx)
         case (iBloqueJx); prefix_field = prefix(iBloqueJy)
         case (iBloqueJy); prefix_field = prefix(iBloqueJz)
         case (iBloqueJz); prefix_field = prefix(iBloqueJx)
         case (iBloqueMx); prefix_field = prefix(iBloqueMy)
         case (iBloqueMy); prefix_field = prefix(iBloqueMz)
         case (iBloqueMz); prefix_field = prefix(iBloqueMx)
         case default; prefixExtension = prefix(field)
         end select
      else
         call stoponerror(0, 0, "Buggy error in mpidir.")
      end if
      return
   end function get_rotated_prefix

   function prefix(campo) result(ext)
      integer(kind=SINGLE), intent(in)  ::  campo
      character(len=BUFSIZE)  ::  ext

      select case (campo)
      case (iEx); ext = 'Ex'
      case (iEy); ext = 'Ey'
      case (iEz); ext = 'Ez'
      case (iVx); ext = 'Vx'
      case (iVy); ext = 'Vy'
      case (iVz); ext = 'Vz'
      case (iHx); ext = 'Hx'
      case (iHy); ext = 'Hy'
      case (iHz); ext = 'Hz'
      case (iBloqueJx); ext = 'Jx'
      case (iBloqueJy); ext = 'Jy'
      case (iBloqueJz); ext = 'Jz'
      case (iBloqueMx); ext = 'Mx'
      case (iBloqueMy); ext = 'My'
      case (iBloqueMz); ext = 'Mz'
      case (iJx); ext = 'Wx'
      case (iJy); ext = 'Wy'
      case (iJz); ext = 'Wz'
      case (iQx); ext = 'Qx'
      case (iQy); ext = 'Qy'
      case (iQz); ext = 'Qz'
      case (iExC); ext = 'ExC'
      case (iEyC); ext = 'EyC'
      case (iEzC); ext = 'EzC'
      case (iHxC); ext = 'HxC'
      case (iHyC); ext = 'HyC'
      case (iHzC); ext = 'HzC'
      case (iMEC); ext = 'ME'
      case (iMHC); ext = 'MH'
      case (iCur); ext = 'BC'
      case (mapvtk); ext = 'MAP'
      case (iCurX); ext = 'BCX'
      case (iCurY); ext = 'BCY'
      case (iCurZ); ext = 'BCZ'
      case (farfield); ext = 'FF'
      case (lineIntegral); ext = 'LI'
      end select
      return
   end function prefix

   function open_file(fileUnit, fileName) result(iostat)
      character(len=*), intent(in) :: fileName
      integer(kind=SINGLE), intent(in) :: fileUnit
      integer(kind=SINGLE) :: iostat

      open (unit=fileUnit, file=fileName, status='OLD', action='WRITE', position='APPEND', iostat=iostat)
      if (iostat /= 0) then
         open (unit=fileUnit, file=fileName, status='NEW', action='WRITE', iostat=iostat)
      end if
      return
   end function open_file

   function close_file(fileUnit) result(iostat)
      integer(kind=SINGLE), intent(in) :: fileUnit
      integer(kind=SINGLE) :: iostat

      close (fileUnit, iostat=iostat)
   end function close_file

   subroutine init_frequency_slice(frequencySlice, domain)
      real(kind=RKIND), dimension(:), intent(out) :: frequencySlice
      type(domain_t), intent(in) :: domain

      integer(kind=SINGLE) :: i

      if (domain%logarithmicSpacing) then
         do i = 1, domain%fnum
            frequencySlice(i) = 10.0_RKIND**(domain%fstart + (i - 1)*domain%fstep)
         end do
      else
         do i = 1, domain%fnum
            frequencySlice(i) = domain%fstart + (i - 1)*domain%fstep
         end do
      end if
   end subroutine init_frequency_slice

   integer function blockCurrent(field)
      integer(kind=4) :: field
      select case (field)
      case (iHx); blockCurrent = iCurX
      case (iHy); blockCurrent = iCurY
      case (iHz); blockCurrent = iCurZ
      case default; call StopOnError(layoutnumber, size, 'field is not H field')
      end select
   end function

   logical function isPECorSurface(field, i, j, k, media, simulationMedia)
      type(MediaData_t), pointer, dimension(:), intent(in) :: simulationMedia
      type(media_matrices_t), intent(in) :: media
      integer(kind=4), intent(in) :: field, i, j, k
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: mediaIndex
      mediaIndex = getMedia(field, i, j, k, media)
      isPECorSurface = simulationMedia(mediaIndex)%is%PEC .or. simulationMedia(mediaIndex)%is%Surface
   end function

   function getMedia(field, i, j, k, media) result(res)
      TYPE(media_matrices_t), INTENT(IN) :: media
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: res
      integer(kind=4) :: field, i, j, k
      select case (field)
      case (iEx); res = media%sggMiEx(i, j, k)
      case (iEy); res = media%sggMiEy(i, j, k)
      case (iEz); res = media%sggMiEz(i, j, k)
      case (iHx); res = media%sggMiHx(i, j, k)
      case (iHy); res = media%sggMiHy(i, j, k)
      case (iHz); res = media%sggMiHz(i, j, k)
      case default; call StopOnError(layoutnumber, size, 'Unrecognized field')
      end select
   end function

   logical function isWithinBounds(field, i, j, k, SINPML_fullsize)
      TYPE(limit_t), DIMENSION(:), INTENT(IN) :: SINPML_fullsize
      integer(kind=4) :: field, i, j, k
      isWithinBounds = (i <= SINPML_fullsize(field)%XE) .and. &
                       (j <= SINPML_fullsize(field)%YE) .and. &
                       (k <= SINPML_fullsize(field)%ZE)
   end function

   logical function isMediaVacuum(field, i, j, k, media)
      TYPE(media_matrices_t), INTENT(IN) :: media
      integer(kind=4) :: field, i, j, k
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: mediaIndex, vacuum = 1
      mediaIndex = getMedia(field, i, j, k, media)
      isMediaVacuum = (mediaIndex == vacuum)
   end function

   logical function isSplitOrAdvanced(field, i, j, k, media, simulationMedia)
      type(MediaData_t), pointer, dimension(:), intent(in) :: simulationMedia
      type(media_matrices_t), intent(in) :: media
      integer(kind=4) :: field, i, j, k
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: mediaIndex
      mediaIndex = getMedia(field, i, j, k, media)
      isSplitOrAdvanced = sgg%med(mediaIndex)%is%split_and_useless .or. &
                          sgg%med(mediaIndex)%is%already_YEEadvanced_byconformal

   end function
end module mod_outputUtils
