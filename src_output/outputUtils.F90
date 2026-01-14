module mod_outputUtils
   use FDETYPES
   use outputTypes
   use mod_domain
   use report
   implicit none
   integer(kind=SINGLE), parameter :: FILE_UNIT = 400

   private

   !===========================
   !  Public interface summary
   !===========================
   public :: new_cell_coordinate
   public :: get_coordinates_extension
   public :: get_prefix_extension
   public :: get_field_component
   public :: get_field_reference
   public :: open_file
   public :: close_file
   public :: create_or_clear_file
   public :: init_frequency_slice
   public :: getBlockCurrentDirection
   public :: isPEC
   public :: isSplitOrAdvanced
   public :: isThinWire
   public :: isMediaVacuum
   public :: isWithinBounds
   public :: isSurface
   public :: isFlush
   public :: computej
   public :: computeJ1
   public :: computeJ2
   public :: fieldo
   !===========================

   !===========================
   !  Private interface summary
   !===========================
   private :: get_rotated_prefix
   private :: prefix
   private :: get_probe_coords_extension
   private :: get_probe_bounds_coords_extension
   private :: get_delta
   !===========================

   interface get_coordinates_extension
      module procedure get_probe_coords_extension, get_probe_bounds_coords_extension
   end interface get_coordinates_extension

contains
   function new_cell_coordinate(x, y, z) result(cell)
      integer(kind=SINGLE), intent(in) :: x, y, z
      type(cell_coordinate_t) :: cell

      cell%x = x
      cell%y = y
      cell%z = z
   end function new_cell_coordinate

   function getMediaIndex(field, i, j, k, CoordToMaterial) result(res)
      integer, intent(in) :: field, i, j, k
      type(media_matrices_t), pointer, intent(in) :: CoordToMaterial

      integer :: res

      select case (field)
      case (iEx); res = CoordToMaterial%sggMiEx(i, j, k)
      case (iEy); res = CoordToMaterial%sggMiEy(i, j, k)
      case (iEz); res = CoordToMaterial%sggMiEz(i, j, k)
      case (iHx); res = CoordToMaterial%sggMiHx(i, j, k)
      case (iHy); res = CoordToMaterial%sggMiHy(i, j, k)
      case (iHz); res = CoordToMaterial%sggMiHz(i, j, k)
      end select

   end function

   function get_probe_coords_extension(coordinates, mpidir) result(ext)
      type(cell_coordinate_t) :: coordinates
      integer(kind=SINGLE), intent(in) ::  mpidir
      character(len=BUFSIZE) :: ext
      character(len=BUFSIZE)  ::  chari, charj, chark

      write (chari, '(i7)') coordinates%x
      write (charj, '(i7)') coordinates%y
      write (chark, '(i7)') coordinates%z

#if CompileWithMPI
      if (mpidir == 3) then
         ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
      elseif (mpidir == 2) then
         ext = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))
      elseif (mpidir == 1) then
         ext = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))
      else
         call stoponerror(0, 0, 'Buggy error in mpidir. ')
      end if
#else
      ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
#endif

      return
   end function get_probe_coords_extension

   function get_probe_bounds_coords_extension(lowerCoordinates, upperCoordinates, mpidir) result(ext)
      type(cell_coordinate_t) :: lowerCoordinates, upperCoordinates
      integer(kind=SINGLE), intent(in) :: mpidir
      character(len=BUFSIZE) :: ext
      character(len=BUFSIZE)  ::  chari, charj, chark, chari2, charj2, chark2

      write (chari, '(i7)') lowerCoordinates%x
      write (charj, '(i7)') lowerCoordinates%y
      write (chark, '(i7)') lowerCoordinates%z

      write (chari2, '(i7)') upperCoordinates%x
      write (charj2, '(i7)') upperCoordinates%y
      write (chark2, '(i7)') upperCoordinates%z

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
         call stoponerror(0, 0, 'Buggy error in mpidir. ')
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
         case (iBloqueJx); prefixExtension = prefix(iBloqueJx)
         case (iBloqueJy); prefixExtension = prefix(iBloqueJy)
         case (iBloqueJz); prefixExtension = prefix(iBloqueJz)
         case (iBloqueMx); prefixExtension = prefix(iBloqueMx)
         case (iBloqueMy); prefixExtension = prefix(iBloqueMy)
         case (iBloqueMz); prefixExtension = prefix(iBloqueMz)
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
         case (iBloqueJx); prefixExtension = prefix(iBloqueJz)
         case (iBloqueJy); prefixExtension = prefix(iBloqueJx)
         case (iBloqueJz); prefixExtension = prefix(iBloqueJy)
         case (iBloqueMx); prefixExtension = prefix(iBloqueMz)
         case (iBloqueMy); prefixExtension = prefix(iBloqueMx)
         case (iBloqueMz); prefixExtension = prefix(iBloqueMy)
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
         case (iBloqueJx); prefixExtension = prefix(iBloqueJy)
         case (iBloqueJy); prefixExtension = prefix(iBloqueJz)
         case (iBloqueJz); prefixExtension = prefix(iBloqueJx)
         case (iBloqueMx); prefixExtension = prefix(iBloqueMy)
         case (iBloqueMy); prefixExtension = prefix(iBloqueMz)
         case (iBloqueMz); prefixExtension = prefix(iBloqueMx)
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

   function fieldo(field, dir) result(fieldo2)
      integer  ::  fieldo2, field
      character(len=1) :: dir
      fieldo2 = -1
      select case (field)
      case (iEx, iEy, iEz, iHx, iHy, iHz); fieldo2 = field
      case (iJx, iVx, iBloqueJx, iExC, iQx); fieldo2 = iEx
      case (iJy, iVy, iBloqueJy, iEyC, iQy); fieldo2 = iEy
      case (iJz, iVz, iBloqueJz, iEzC, iQz); fieldo2 = iEz
      case (iBloqueMx, iHxC); fieldo2 = iHx
      case (iBloqueMy, iHyC); fieldo2 = iHy
      case (iBloqueMz, iHzC); fieldo2 = iHz
      case (iMEC)
         select case (dir)
         CASE ('X', 'x'); fieldo2 = iEx
         CASE ('Y', 'y'); fieldo2 = iEY
         CASE ('Z', 'z'); fieldo2 = iEz
         END SELECT
      case (iMHC)
         select case (dir)
         CASE ('X', 'x'); fieldo2 = ihx
         CASE ('Y', 'y'); fieldo2 = iHY
         CASE ('Z', 'z'); fieldo2 = iHz
         END SELECT
      case (iCur, iCurX, icurY, icurZ, mapvtk)  !los pongo en efield para evitar problemas con el MPI
         select case (dir)
         CASE ('X', 'x'); fieldo2 = iEx
         CASE ('Y', 'y'); fieldo2 = iEY
         CASE ('Z', 'z'); fieldo2 = iEz
         END SELECT
      end select
   end function

   function get_field_component(fieldId, fieldReference) result(component)
      type(fields_reference_t), intent(in) :: fieldReference
      integer(kind=SINGLE), intent(in) :: fieldId
      real(kind=RKIND), pointer, dimension(:, :, :) :: component
      select case (fieldId)
      case (iEx); component => fieldReference%E%x
      case (iEy); component => fieldReference%E%y
      case (iEz); component => fieldReference%E%z
      case (iHx); component => fieldReference%H%x
      case (iHy); component => fieldReference%H%y
      case (iHz); component => fieldReference%H%z
      end select
   end function

   function get_field_reference(fieldId, fieldReference) result(field)
      type(fields_reference_t), intent(in) :: fieldReference
      integer(kind=SINGLE), intent(in) :: fieldId
      type(field_data_t) :: field
      select case (fieldId)
      case (iBloqueJx, iBloqueJy, iBloqueJz)
         field%x => fieldReference%E%x
         field%y => fieldReference%E%y
         field%z => fieldReference%E%z

         field%deltaX => fieldReference%E%deltax
         field%deltaY => fieldReference%E%deltay
         field%deltaZ => fieldReference%E%deltaz
      case (iBloqueMx, iBloqueMy, iBloqueMz)
         field%x => fieldReference%H%x
         field%y => fieldReference%H%y
         field%z => fieldReference%H%z

         field%deltaX => fieldReference%H%deltax
         field%deltaY => fieldReference%H%deltay
         field%deltaZ => fieldReference%H%deltaz
      end select
   end function get_field_reference

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

   integer function getBlockCurrentDirection(field)
      integer(kind=4) :: field
      select case (field)
      case (iHx); getBlockCurrentDirection = iCurX
      case (iHy); getBlockCurrentDirection = iCurY
      case (iHz); getBlockCurrentDirection = iCurZ
      case default; call StopOnError(0, 0, 'field is not H field')
      end select
   end function

   logical function isThinWire(field, i, j, k, problem)
      integer(kind=4), intent(in) :: field, i, j, k
      type(problem_info_t), intent(in) :: problem

      integer(kind=SINGLE) :: mediaIndex

      mediaIndex = getMediaIndex(field, i, j, k, problem%geometryToMaterialData)
      isThinWire = problem%materialList(mediaIndex)%is%ThinWire
   end function

   logical function isPEC(field, i, j, k, problem)
      integer(kind=4), intent(in) :: field, i, j, k
      type(problem_info_t), intent(in) :: problem

      integer(kind=SINGLE) :: mediaIndex

      mediaIndex = getMediaIndex(field, i, j, k, problem%geometryToMaterialData)
      isPEC = problem%materialList(mediaIndex)%is%PEC
   end function

   logical function isSurface(field, i, j, k, problem)
      integer(kind=4), intent(in) :: field, i, j, k
      type(problem_info_t), intent(in) :: problem

      integer(kind=SINGLE) :: mediaIndex

      mediaIndex = getMediaIndex(field, i, j, k, problem%geometryToMaterialData)
      isSurface = problem%materialList(mediaIndex)%is%Surface
   end function

   logical function isWithinBounds(field, i, j, k, problem)
      integer(kind=4), intent(in) :: field, i, j, k
      type(problem_info_t), intent(in) :: problem

      isWithinBounds = (i <= problem%problemDimension(field)%XE) .and. &
                       (j <= problem%problemDimension(field)%YE) .and. &
                       (k <= problem%problemDimension(field)%ZE) .and. &
                       (i >= problem%problemDimension(field)%XI) .and. &
                       (j >= problem%problemDimension(field)%YI) .and. &
                       (k >= problem%problemDimension(field)%ZI)
   end function

   logical function isMediaVacuum(field, i, j, k, problem)
      integer(kind=4), intent(in) :: field, i, j, k
      type(problem_info_t), intent(in) :: problem

      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: mediaIndex
      integer(kind=INTEGERSIZEOFMEDIAMATRICES), parameter :: vacuum = 1

      mediaIndex = getMediaIndex(field, i, j, k, problem%geometryToMaterialData)
      isMediaVacuum = (mediaIndex == vacuum)
   end function

   logical function isSplitOrAdvanced(field, i, j, k, problem)
      integer(kind=4), intent(in) :: field, i, j, k
      type(problem_info_t), intent(in) :: problem

      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: mediaIndex
      mediaIndex = getMediaIndex(field, i, j, k, problem%geometryToMaterialData)

      isSplitOrAdvanced = problem%materialList(mediaIndex)%is%split_and_useless .or. &
                          problem%materialList(mediaIndex)%is%already_YEEadvanced_byconformal
   end function

   function computej(field, i, j, k, fields_reference) result(res)
      implicit none

      ! Input Arguments
      integer(kind=single), intent(in) :: field, i, j, k
      type(fields_reference_t), intent(in) :: fields_reference

      ! Local Variables
      integer(kind=single) :: i_shift_a, j_shift_a, k_shift_a  ! Shift for Term A (Offset for H/M field)
      integer(kind=single) :: i_shift_b, j_shift_b, k_shift_b  ! Shift for Term B (Offset for H/M field)

      integer(kind=single) :: curl_component_a                 ! H/M field component for Term A
      integer(kind=single) :: curl_component_b                 ! H/M field component for Term B

      real(kind=rkind) :: res

      ! -----------------------------------------------------------
      ! 1. Determine Curl Components
      !    The MOD 3 operation cyclically maps the E-field to the two required H-field components.
      ! -----------------------------------------------------------

      ! Component A (The 'next' component in the sequence)
      curl_component_a = 1 + mod(field + 1, 3)

      ! Component B (The 'current' component in the sequence)
      curl_component_b = 1 + mod(field, 3)

      ! -----------------------------------------------------------
      ! 2. Calculate Spatial Shifts (Yee Cell Staggering)
      !    We use MERGE to apply the (i-1) shift only in the relevant direction.
      ! -----------------------------------------------------------

      ! Shift for Term A
      i_shift_a = i - merge(1, 0, curl_component_a == iex)
      j_shift_a = j - merge(1, 0, curl_component_a == iey)
      k_shift_a = k - merge(1, 0, curl_component_a == iez)

      ! Shift for Term B
      i_shift_b = i - merge(1, 0, curl_component_b == iex)
      j_shift_b = j - merge(1, 0, curl_component_b == iey)
      k_shift_b = k - merge(1, 0, curl_component_b == iez)

      ! -----------------------------------------------------------
      ! 3. Calculate J (Curl Difference)
      !    The H/M fields are accessed using an offset (+3) from the E-field index.
      ! -----------------------------------------------------------

      res = &
         ! TERM B: (Positive term in the difference)
         (get_delta(curl_component_b, i, j, k, fields_reference)* &
      ( get_field(curl_component_b + 3, i, j, k, fields_reference) - get_field(curl_component_b + 3, i_shift_b, j_shift_b, k_shift_b, fields_reference) ) &
          ) - &
         ! TERM A: (Negative term in the difference)
         (get_delta(curl_component_a, i, j, k, fields_reference)* &
      ( get_field(curl_component_a + 3, i, j, k, fields_reference) - get_field(curl_component_a + 3, i_shift_a, j_shift_a, k_shift_a, fields_reference) ) &
          )

   end function computej

   function computeJ1(f, i, j, k, fields_reference) result(res)
      implicit none
      integer(kind=4), intent(in) :: f, i, j, k
      type(fields_reference_t), intent(in) :: fields_reference
      integer(kind=4) :: c       ! Complementary H-field index (Hy/Hz)
      real(kind=rkind) :: res
      real(kind=rkind) :: curl_h_term_a, curl_h_term_b, field_diff_term

      ! Calculate complementary H-field index (e.g., if f=1 (Ex), c=5 (Hy) and c+1=6 (Hz) or vice versa depending on definitions)
      ! For f=1 (Ex), c = mod(1-2, 3)+4 = mod(-1, 3)+4 = 2+4 = 6 (Hz).

      c = mod(f - 2, 3) + 4 ! This typically corresponds to H_z for J_x, or H_x for J_y, etc.

      ! First set of H-field terms
      curl_h_term_a = get_delta(c, i, j, k, fields_reference)*get_field(c, i, j, k, fields_reference) + &
                    get_delta(c, i+u(f,iHy), j+u(f,iHz), k+u(f,iHx), fields_reference) * get_field(c, i+u(f,iHy), j+u(f,iHz), k+u(f,iHx), fields_reference)

      ! Second set of H-field terms
    curl_h_term_b = get_delta(c, i, j, k, fields_reference) * get_field(c, i-u(f,iHx), j-u(f,iHy), k-u(f,iHz), fields_reference) + &
                    get_delta(c, i+u(f,iHy), j+u(f,iHz), k+u(f,iHx), fields_reference) * get_field(c, i-u(f,iHx)+u(f,iHy), j-u(f,iHy)+u(f,iHz), k-u(f,iHz)+u(f,iHx), fields_reference)

      ! E-field term (approximates the change in E-field at the J-node)
      field_diff_term = get_delta(f, i, j, k, fields_reference)*( &
                        get_field(f, i - u(f, iHy), j - u(f, iHz), k - u(f, iHx), fields_reference) - &
                        get_field(f, i + u(f, iHy), j + u(f, iHz), k + u(f, iHx), fields_reference))

      ! Final computation: J1 = - ((Curl_H_A) - (Curl_H_B) + (E_diff))
      res = -((curl_h_term_a - curl_h_term_b) + field_diff_term)

   end function computeJ1

   function computeJ2(f, i, j, k, fields_reference) result(res)
      implicit none
      integer(kind=4), intent(in) :: f, i, j, k
      type(fields_reference_t), intent(in) :: fields_reference
      integer(kind=4) :: c       ! Complementary H-field index (Hx/Hy/Hz)
      real(kind=rkind) :: res
      real(kind=rkind) :: curl_h_term_a, curl_h_term_b, field_diff_term

      ! Calculate complementary H-field index (e.g., if f=1 (Ex), c=4 (Hx) or c=5 (Hy))
      ! For f=1 (Ex), c = mod(1-3, 3)+4 = mod(-2, 3)+4 = 1+4 = 5 (Hy). This is the second H-field curl component.
      c = mod(f - 3, 3) + 4

      ! First set of H-field terms
      curl_h_term_a = get_delta(c, i, j, k, fields_reference)*get_field(c, i, j, k, fields_reference) + &
                    get_delta(c, i+u(f,iHz), j+u(f,iHx), k+u(f,iHy), fields_reference) * get_field(c, i+u(f,iHz), j+u(f,iHx), k+u(f,iHy), fields_reference)

      ! Second set of H-field terms
    curl_h_term_b = get_delta(c, i, j, k, fields_reference) * get_field(c, i-u(f,iHx), j-u(f,iHy), k-u(f,iHz), fields_reference) + &
                    get_delta(c, i+u(f,iHz), j+u(f,iHx), k+u(f,iHy), fields_reference) * get_field(c, i-u(f,iHx)+u(f,iHz), j-u(f,iHy)+u(f,iHx), k-u(f,iHz)+u(f,iHy), fields_reference)

      ! E-field term (approximates the change in E-field at the J-node)
      field_diff_term = get_delta(f, i, j, k, fields_reference)*( &
                        get_field(f, i - u(f, iHz), j - u(f, iHx), k - u(f, iHy), fields_reference) - &
                        get_field(f, i + u(f, iHz), j + u(f, iHx), k + u(f, iHy), fields_reference))

      ! Final computation: J2 = (Curl_H_A) - (Curl_H_B) + (E_diff)
      res = (curl_h_term_a - curl_h_term_b) + field_diff_term

   end function computeJ2

   integer function u(field1, field2)
      integer(kind=4) :: field1, field2
      if (field1 == field2) then
         u = 1
      else
         u = 0
      end if
   end function

   function get_field(field, i, j, k, fields_reference) result(res)
      implicit none
      real(kind=rkind) :: res
      integer(kind=4), intent(in) :: field, i, j, k
      type(fields_reference_t), intent(in) :: fields_reference

      ! Retrieves the field value based on the field index (1-3 for E, 4-6 for H)
      select case (field)
      case (iex); res = fields_reference%e%x(i, j, k)
      case (iey); res = fields_reference%e%y(i, j, k)
      case (iez); res = fields_reference%e%z(i, j, k)
      case (ihx); res = fields_reference%h%x(i, j, k)
      case (ihy); res = fields_reference%h%y(i, j, k)
      case (ihz); res = fields_reference%h%z(i, j, k)
      end select
   end function get_field

   function get_delta(field, i, j, k, fields_reference) result(res)
      implicit none
      real(kind=rkind) :: res
      integer(kind=4), intent(in) :: field, i, j, k
      type(fields_reference_t), intent(in) :: fields_reference

      ! Retrieves the spatial step size (delta) corresponding to the field direction
      ! Note: i, j, k are used to select the correct array index if the grid is non-uniform.
      select case (field)
      case (iex); res = fields_reference%e%deltax(i)
      case (iey); res = fields_reference%e%deltay(j)
      case (iez); res = fields_reference%e%deltaz(k)
      case (ihx); res = fields_reference%h%deltax(i)
      case (ihy); res = fields_reference%h%deltay(j)
      case (ihz); res = fields_reference%h%deltaz(k)
      end select
   end function get_delta

   subroutine create_or_clear_file(path, unit_out, err)
      implicit none
      character(len=*), intent(in)  :: path
      integer, intent(out) :: unit_out
      integer, intent(out) :: err
      integer :: unit, ios
      logical :: opened
      character(len=BUFSIZE) :: fname
      integer, parameter :: unit_min = 10, unit_max = 99

      err = 0
      unit_out = -1

      ! --- Find a free unit ---
      do unit = unit_min, unit_max
         inquire (unit=unit, opened=opened, name=fname)
         if (.not. opened) exit       ! Found free unit
         if (trim(fname) == trim(path)) then
            ! Unit is already associated with the same file -> safe to clear
            close (unit)
            exit
         end if
      end do

      ! Check if no free unit was found
      inquire (unit=unit, opened=opened)
      if (opened) then
         err = 1
         return
      end if

      ! --- Open the file, replacing it if it exists ---
      open (unit=unit, file=path, status="replace", action="write", iostat=ios)
      if (ios /= 0) then
         err = 2
         return
      end if

      close (unit)

      ! --- Success ---
      unit_out = unit
   end subroutine create_or_clear_file

end module mod_outputUtils
