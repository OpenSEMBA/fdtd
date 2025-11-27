module mod_outputUtils
   use FDETYPES
   implicit none

contains

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
         case default; prefixExtension = prefix(field)
         end select
      else
         call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
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

      open (unit=fileUnit, file=fileName status='OLD', action='WRITE', possition='APPEND', iostat=iostat)
      if (iostat /= 0) then
         open (unit=fileUnit, file=fileName status='NEW', action='WRITE', iostat=iostat)
      end if
      return
   end function open_file

   function close_file(fileUnit) result(iostat)
      integer(kind=SINGLE), intent(in) :: fileUnit
      integer(kind=SINGLE) :: iostat

      close (fileUnit, iostat=iostat)
   end function close_file
end module mod_outputUtils
