module snapxdmf_m
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use FDETYPES_m
   use xdmf_hdf5_m, only: xdmf_writer_t, xdmf_options_t, xdmf_status_t, &
      xdmf_grid_id_t, xdmf_attribute_id_t, XDMF_SERIES_TIME, &
      XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL32
   implicit none
   private
   public WRITE_XDMFSNAP

contains

   subroutine write_xdmfsnap(ninstant, filename, minXabs, maxXabs, minYabs, &
      maxYabs, minZabs, maxZabs, valor3D)
      integer(kind=4), intent(in) :: ninstant
      character(len=BUFSIZE), intent(in) :: filename
      integer(kind=4), intent(in) :: minXabs, maxXabs, minYabs, maxYabs
      integer(kind=4), intent(in) :: minZabs, maxZabs
      real(kind=4), intent(in) :: valor3D(minXabs:maxXabs, minYabs:maxYabs, &
         minZabs:maxZabs, 1:1)

      type(xdmf_writer_t) :: writer
      type(xdmf_options_t) :: options
      type(xdmf_status_t) :: status
      type(xdmf_grid_id_t) :: grid
      type(xdmf_attribute_id_t) :: attribute
      integer(int64) :: dimensions(3)
      real(real64) :: origin(3), spacing(3)

      options%overwrite = .true.
      options%series_kind = XDMF_SERIES_TIME
      call writer%create(trim(filename), options, status)
      call check_status(status)

      dimensions = [int(maxXabs - minXabs + 1, int64), &
         int(maxYabs - minYabs + 1, int64), &
         int(maxZabs - minZabs + 1, int64)]
      origin = [real(minXabs, real64), real(minYabs, real64), &
         real(minZabs, real64)]
      spacing = 1.0_real64
      call writer%define_uniform_grid('snapshot', dimensions, origin, spacing, &
         grid, status)
      call check_status(status)
      call writer%define_attribute(grid, 'values', XDMF_CENTER_NODE, &
         XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL32, .true., attribute, status)
      call check_status(status)
      call writer%begin_step(real(ninstant, real64), status)
      call check_status(status)
      call writer%write_attribute(attribute, &
         reshape(valor3D(:, :, :, 1), [size(valor3D(:, :, :, 1))]), status)
      call check_status(status)
      call writer%end_step(status)
      call check_status(status)
      call writer%close(status)
      call check_status(status)

   contains

      subroutine check_status(result)
         type(xdmf_status_t), intent(in) :: result

         if (result%is_error()) error stop result%message()
      end subroutine check_status
   end subroutine write_xdmfsnap
end module snapxdmf_m
