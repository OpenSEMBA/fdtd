module xdmf_h5_m
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use FDETYPES_m
   use xdmf_hdf5_m, only: xdmf_writer_t, xdmf_options_t, xdmf_status_t, &
      xdmf_grid_id_t, xdmf_attribute_id_t, XDMF_SERIES_TIME, &
      XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL32, &
      XDMF_NUMERIC_REAL64
   implicit none
   private
   public openh5file, writeh5file, closeh5file, createh5filefromsinglebin

#ifdef CompileWithReal8
   integer, parameter :: FIELD_NUMERIC_TYPE = XDMF_NUMERIC_REAL64
#else
   integer, parameter :: FIELD_NUMERIC_TYPE = XDMF_NUMERIC_REAL32
#endif

   type(xdmf_writer_t) :: writer
   type(xdmf_grid_id_t) :: grid
   type(xdmf_attribute_id_t) :: attribute
   integer(int64) :: field_dimensions(3)
   logical :: writer_defined = .false.

contains

   subroutine openh5file(filename, finalstep, minXabs, maxXabs, minYabs, &
      maxYabs, minZabs, maxZabs)
      character(len=*), intent(in) :: filename
      integer(kind=4), intent(in) :: finalstep
      integer(kind=4), intent(in) :: minXabs, maxXabs, minYabs, maxYabs
      integer(kind=4), intent(in) :: minZabs, maxZabs

      type(xdmf_options_t) :: options
      type(xdmf_status_t) :: status

      options%overwrite = .true.
      options%series_kind = XDMF_SERIES_TIME
      call writer%create(trim(filename), options, status)
      call check_status(status)
      field_dimensions = [int(maxXabs - minXabs + 1, int64), &
         int(maxYabs - minYabs + 1, int64), &
         int(maxZabs - minZabs + 1, int64)]
      writer_defined = .false.
   end subroutine openh5file

   subroutine writeh5file(filename, valor3d, indi, attindi, minXabs, maxXabs, &
      minYabs, maxYabs, minZabs, maxZabs, linez_minZabs_primero, &
      liney_minYabs_primero, linex_minXabs_primero, dz_minZabs, dy_minYabs, &
      dx_minXabs, minZabs_primero, minYabs_primero, minXabs_primero, &
      finalstep, vtkindex)
      character(len=*), intent(in) :: filename
      real(kind=RKIND), intent(in) :: valor3d(:, :, :, :)
      integer(kind=4), intent(in) :: indi
      real(kind=RKIND_tiempo), intent(in) :: attindi
      integer(kind=4), intent(in) :: minXabs, maxXabs, minYabs, maxYabs
      integer(kind=4), intent(in) :: minZabs, maxZabs
      real(kind=RKIND), intent(in) :: linez_minZabs_primero, &
         liney_minYabs_primero, linex_minXabs_primero
      real(kind=RKIND), intent(in) :: dz_minZabs, dy_minYabs, dx_minXabs
      integer(kind=4), intent(in) :: minZabs_primero, minYabs_primero, &
         minXabs_primero, finalstep
      logical, intent(in) :: vtkindex

      type(xdmf_status_t) :: status
      real(real64) :: origin(3), spacing(3)

      if (.not. writer_defined) then
         origin = [real(linex_minXabs_primero, real64), &
            real(liney_minYabs_primero, real64), &
            real(linez_minZabs_primero, real64)]
         spacing = [real(dx_minXabs, real64), real(dy_minYabs, real64), &
            real(dz_minZabs, real64)]
         call writer%define_uniform_grid('grid', field_dimensions, origin, &
            spacing, grid, status)
         call check_status(status)
         call writer%define_attribute(grid, 'values', XDMF_CENTER_NODE, &
            XDMF_ATTRIBUTE_SCALAR, FIELD_NUMERIC_TYPE, .true., attribute, &
            status)
         call check_status(status)
         writer_defined = .true.
      end if

      call writer%begin_step(real(attindi, real64), status)
      call check_status(status)
      if (FIELD_NUMERIC_TYPE == XDMF_NUMERIC_REAL64) then
         call writer%write_attribute(attribute, &
            reshape(real(valor3d(:, :, :, 1), real64), &
            [size(valor3d(:, :, :, 1))]), status)
      else
         call writer%write_attribute(attribute, &
            reshape(real(valor3d(:, :, :, 1), kind=kind(1.0)), &
            [size(valor3d(:, :, :, 1))]), status)
      end if
      call check_status(status)
      call writer%end_step(status)
      call check_status(status)
   end subroutine writeh5file

   subroutine closeh5file(finalstep, att)
      integer(kind=4), intent(in) :: finalstep
      real(kind=RKIND_tiempo), intent(in) :: att(:)

      type(xdmf_status_t) :: status

      call writer%close(status)
      call check_status(status)
      writer_defined = .false.
   end subroutine closeh5file

   subroutine createh5filefromsinglebin(filename, vtkindex)
      character(len=*), intent(in) :: filename
      logical, intent(in) :: vtkindex

      integer(kind=4) :: myunit, fieldob, pass, total_passes
      real(kind=RKIND_tiempo), allocatable :: att(:)
      real(kind=RKIND), allocatable :: valor3d(:, :, :, :)
      logical :: time_domain
      integer(kind=4) :: minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs
      integer(kind=4) :: minZabs_primero, minYabs_primero, minXabs_primero
      integer(kind=4) :: finalstep, step_index, i1, j1, k1
      real(kind=RKIND) :: linez_minZabs_primero, liney_minYabs_primero
      real(kind=RKIND) :: linex_minXabs_primero, dz_minZabs, dy_minYabs
      real(kind=RKIND) :: dx_minXabs
      character(len=BUFSIZE) :: input_name, message

      input_name = filename(1:index(filename, '.h5bin') - 1)
      input_name = trim(adjustl(input_name))
      open(newunit=myunit, file=trim(input_name)//'.h5bin', form='unformatted')
      read(myunit) finalstep, minXabs, maxXabs, minYabs, maxYabs, minZabs, &
         maxZabs, fieldob, time_domain, total_passes
      allocate(valor3d(minXabs:maxXabs, minYabs:maxYabs, minZabs:maxZabs, 1))
      allocate(att(1:finalstep))

      do pass = 1, total_passes
         if (time_domain) then
            if (pass == 1) then
               input_name = trim(adjustl(filename(1:index(filename, '.h5bin') - 1)))//'_time'
            else
               print *, 'Buggy error in valor3d.'
               stop 1
            end if
         else
            if (pass == 1) then
               input_name = trim(adjustl(filename(1:index(filename, '.h5bin') - 1)))//'_mod'
            else if (pass == 2) then
               input_name = trim(adjustl(filename(1:index(filename, '.h5bin') - 1)))//'_phase'
            else
               print *, 'Buggy error in valor3d.'
               stop 1
            end if
         end if

         if (.not. (((fieldob == iMEC) .or. (fieldob == iMHC)) .and. &
            (pass == 2))) then
            call openh5file(input_name, finalstep, minXabs, maxXabs, minYabs, &
               maxYabs, minZabs, maxZabs)
         end if

         valor3d = 0.0_RKIND
         att = 0.0_RKIND_tiempo
         do step_index = 1, finalstep
            read(myunit) minZabs_primero, minYabs_primero, minXabs_primero
            read(myunit) linez_minZabs_primero, liney_minYabs_primero, &
               linex_minXabs_primero
            read(myunit) dz_minZabs, dy_minYabs, dx_minXabs
            read(myunit) att(step_index)
            write(message, *) ' ----> .xdmf file ', att(step_index), '(', &
               step_index, '/', &
               finalstep, ')'
            print *, trim(adjustl(message))
            do k1 = minZabs, maxZabs
               do j1 = minYabs, maxYabs
                  read(myunit) (valor3d(i1, j1, k1, 1), i1 = minXabs, maxXabs)
               end do
            end do
            if (.not. (((fieldob == iMEC) .or. (fieldob == iMHC)) .and. &
               (pass == 2))) then
               call writeh5file(input_name, valor3d, step_index, &
                  att(step_index), &
                  minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs, &
                  linez_minZabs_primero, liney_minYabs_primero, &
                  linex_minXabs_primero, dz_minZabs, dy_minYabs, dx_minXabs, &
                  minZabs_primero, minYabs_primero, minXabs_primero, &
                  finalstep, vtkindex)
            end if
         end do
         if (.not. (((fieldob == iMEC) .or. (fieldob == iMHC)) .and. &
            (pass == 2))) call closeh5file(finalstep, att)
      end do
      close(myunit)
      deallocate(valor3d, att)
   end subroutine createh5filefromsinglebin

   subroutine check_status(result)
      type(xdmf_status_t), intent(in) :: result

      if (result%is_error()) then
         print *, trim(result%message())
         stop 1
      end if
   end subroutine check_status
end module xdmf_h5_m
