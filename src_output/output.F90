module output
  use FDETYPES
  implicit none
  character(len=4) :: datFileExtension = '.dat'

  type solver_output_t
    type(point_probe_output_t        ), allocatable :: pointProbe
    type(wire_probe_output_t         ), allocatable :: wireProbe
    type(bulk_current_probe_output_t ), allocatable :: bulkCurrentProbe
    type(far_field_t                 ), allocatable :: farField
    type(time_movie_output_t         ), allocatable :: timeMovie
    type(frequency_slice_output_t    ), allocatable :: frequencySlice
  end type solver_output_t

  type point_probe_output_t
    integer(kind=SINGLE) :: columnas = 2_SINGLE
    character(len=BUFSIZE) :: path

  end type point_probe_output_t


  interface init_solver_output
    module procedure &
      init_point_probe_output , &
      init_wire_probe_output , &
      init_bulk_current_probe_output , &
      init_far_field , &
      initime_movie_output , &
      init_frequency_slice_output
  end interface

  interface update_solver_output
    module procedure &
      update_point_probe_output , &
      update_wire_probe_output , &
      update_bulk_current_probe_output , &
      update_far_field , &
      updateime_movie_output , &
      update_frequency_slice_output
  end interface

  interface flush_solver_output
    module procedure &
      flush_point_probe_output , &
      flush_wire_probe_output , &
      flush_bulk_current_probe_output , &
      flush_far_field , &
      flushime_movie_output , &
      flush_frequency_slice_output
  end interface

  interface delete_solver_output
    module procedure &
      delete_point_probe_output , &
      delete_wire_probe_output , &
      delete_bulk_current_probe_output , &
      delete_far_field , &
      deleteime_movie_output , &
      delete_frequency_slice_output
  end interface
contains

    subroutine init_point_probe_output(probeOutput, iCoord, jCoord, kCoord, field, outputTypeExtension, mpidir)
      type(point_probe_output_t), intent(out) :: probeOutput
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension

      character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension

      probeBoundsExtension = get_probe_bounds_extension(iCoord, jCoord, kCoord, mpidir)
      prefixFieldExtension = get_prefix_extension(field, mpidir)
      probeOutput%path = &
        trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//trim(adjustl(probeBoundsExtension))//trim(adjustl(datFileExtension))

    end subroutine init_point_probe_output

    function get_probe_bounds_extension(iCoord, jCoord, kCoord) result(probeBoundsExtension)
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord
      character(len=BUFSIZE) :: probeBoundsExtension
      
      character(len=BUFSIZE)  ::  chari, charj, chark

      write (chari, '(i7)') iCoord
      write (charj, '(i7)') jCoord
      write (chark, '(i7)') kCoord

    #if CompileWithMPI
      if (mpidir == 3) then 
        probeBoundsExtension = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
      elseif (mpidir == 2) then 
        probeBoundsExtension = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))
      elseif (mpidir == 1) then 
        probeBoundsExtension = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))
      else 
        call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
      end if
    #else
      probeBoundsExtension = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
    #endif

      return
    end function get_probe_bounds_extension

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


end module output