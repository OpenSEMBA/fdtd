module outputVisualisation_m
   use, intrinsic :: iso_fortran_env, only: int64, real64
   use FDETYPES_m, only: BUFSIZE
   use directoryUtils_m, only: file_exists
   use xdmf_hdf5_m, only: xdmf_writer_t, xdmf_options_t, xdmf_status_t, &
                          xdmf_grid_id_t, xdmf_attribute_id_t, XDMF_SERIES_FREQUENCY, XDMF_SERIES_TIME, &
                          XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, &
                          XDMF_TOPOLOGY_POLYVERTEX
   implicit none(type, external)
   private

   integer, parameter, public :: VISUALISATION_SUCCESS = 0
   integer, parameter, public :: VISUALISATION_IO_ERROR = 1
   integer, parameter, public :: VISUALISATION_WRITER_ERROR = 2

   integer, parameter, public :: VISUALISATION_ATTRIBUTE_X = 1
   integer, parameter, public :: VISUALISATION_ATTRIBUTE_Y = 2
   integer, parameter, public :: VISUALISATION_ATTRIBUTE_Z = 3
   integer, parameter, public :: VISUALISATION_ATTRIBUTE_TAG = 4
   integer, parameter, public :: VISUALISATION_ATTRIBUTE_X_PHASE = 4
   integer, parameter, public :: VISUALISATION_ATTRIBUTE_Y_PHASE = 5
   integer, parameter, public :: VISUALISATION_ATTRIBUTE_Z_PHASE = 6

   type, public :: visualisation_writer_t
      private
      type(xdmf_writer_t), pointer :: writer => null()
      type(xdmf_grid_id_t) :: grid
      type(xdmf_attribute_id_t) :: attributes(6)
   end type visualisation_writer_t

   public :: initialise_movie_visualisation
   public :: initialise_frequency_slice_visualisation
   public :: begin_visualisation_step
   public :: write_visualisation_attribute
   public :: write_visualisation_attribute_hyperslab
   public :: end_visualisation_step
   public :: flush_visualisation
   public :: close_visualisation
   public :: visualisation_is_open
   public :: verify_volumetric_visualisation

contains

   subroutine initialise_movie_visualisation(visualisation, path, x_coordinates, y_coordinates, z_coordinates, &
                                             attribute_names, attribute_enabled, collective_io, communicator, &
                                             status, diagnostic)
      type(visualisation_writer_t), intent(inout) :: visualisation
      character(len=*), intent(in) :: path
      real(real64), intent(in) :: x_coordinates(:), y_coordinates(:), z_coordinates(:)
      character(len=*), intent(in) :: attribute_names(4)
      logical, intent(in) :: attribute_enabled(4), collective_io
      integer, intent(in) :: communicator
      integer, intent(out) :: status
      character(len=BUFSIZE), intent(out) :: diagnostic

      type(xdmf_options_t) :: options
      type(xdmf_status_t) :: writer_status
      integer :: attribute_index

      allocate (visualisation%writer)
      options%overwrite = .true.
      options%series_kind = XDMF_SERIES_TIME
      options%collective_io = collective_io
      if (collective_io) then
         options%communicator = communicator
         options%root_rank = 0
      end if
      call visualisation%writer%create(trim(path), options, writer_status)
      if (.not. writer_status%is_error()) then
         call visualisation%writer%define_rectilinear_grid('movieProbe', x_coordinates, y_coordinates, &
                                                           z_coordinates, visualisation%grid, writer_status)
      end if
      do attribute_index = 1, size(attribute_names)
         if (writer_status%is_error()) exit
         if (attribute_enabled(attribute_index)) then
            call define_attribute(visualisation, attribute_index, trim(attribute_names(attribute_index)), &
                                  writer_status)
         end if
      end do
      call convert_status(writer_status, status, diagnostic)
   end subroutine initialise_movie_visualisation

   subroutine initialise_frequency_slice_visualisation(visualisation, path, points, collective_io, communicator, &
                                                       status, diagnostic)
      type(visualisation_writer_t), intent(inout) :: visualisation
      character(len=*), intent(in) :: path
      real(real64), intent(in) :: points(:, :)
      logical, intent(in) :: collective_io
      integer, intent(in) :: communicator
      integer, intent(out) :: status
      character(len=BUFSIZE), intent(out) :: diagnostic

      type(xdmf_options_t) :: options
      type(xdmf_status_t) :: writer_status
      integer(int64), allocatable :: connectivity(:, :)
      character(len=10), parameter :: attribute_names(6) = [character(len=10) :: &
                                                            'xMagnitude', 'yMagnitude', 'zMagnitude', 'xPhase', 'yPhase', 'zPhase']
      integer :: attribute_index, point_index

      allocate (connectivity(1, size(points, 2)))
      do point_index = 1, size(points, 2)
         connectivity(1, point_index) = int(point_index, int64)
      end do

      allocate (visualisation%writer)
      options%overwrite = .true.
      options%series_kind = XDMF_SERIES_FREQUENCY
      options%collective_io = collective_io
      if (collective_io) then
         options%communicator = communicator
         options%root_rank = 0
      end if
      call visualisation%writer%create(trim(path), options, writer_status)
      if (.not. writer_status%is_error()) then
         call visualisation%writer%define_unstructured_grid('frequencySlice', XDMF_TOPOLOGY_POLYVERTEX, &
                                                            points, connectivity, visualisation%grid, writer_status)
      end if
      do attribute_index = 1, size(attribute_names)
         if (writer_status%is_error()) exit
         call define_attribute(visualisation, attribute_index, trim(attribute_names(attribute_index)), writer_status)
      end do
      call convert_status(writer_status, status, diagnostic)
   end subroutine initialise_frequency_slice_visualisation

   subroutine define_attribute(visualisation, attribute_index, name, status)
      type(visualisation_writer_t), intent(inout) :: visualisation
      integer, intent(in) :: attribute_index
      character(len=*), intent(in) :: name
      type(xdmf_status_t), intent(out) :: status

      call visualisation%writer%define_attribute(visualisation%grid, name, XDMF_CENTER_NODE, &
                              XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., visualisation%attributes(attribute_index), status)
   end subroutine define_attribute

   subroutine begin_visualisation_step(visualisation, coordinate, status, diagnostic)
      type(visualisation_writer_t), intent(inout) :: visualisation
      real(real64), intent(in) :: coordinate
      integer, intent(out) :: status
      character(len=BUFSIZE), intent(out) :: diagnostic
      type(xdmf_status_t) :: writer_status

      call visualisation%writer%begin_step(coordinate, writer_status)
      call convert_status(writer_status, status, diagnostic)
   end subroutine begin_visualisation_step

   subroutine write_visualisation_attribute(visualisation, attribute_index, values, status, diagnostic)
      type(visualisation_writer_t), intent(inout) :: visualisation
      integer, intent(in) :: attribute_index
      real(real64), intent(in) :: values(:)
      integer, intent(out) :: status
      character(len=BUFSIZE), intent(out) :: diagnostic
      type(xdmf_status_t) :: writer_status

      call visualisation%writer%write_attribute(visualisation%attributes(attribute_index), values, writer_status)
      call convert_status(writer_status, status, diagnostic)
   end subroutine write_visualisation_attribute

   subroutine write_visualisation_attribute_hyperslab(visualisation, attribute_index, values, offset, shape, &
                                                      status, diagnostic)
      type(visualisation_writer_t), intent(inout) :: visualisation
      integer, intent(in) :: attribute_index
      real(real64), intent(in) :: values(:)
      integer(int64), intent(in) :: offset(:), shape(:)
      integer, intent(out) :: status
      character(len=BUFSIZE), intent(out) :: diagnostic
      type(xdmf_status_t) :: writer_status

      call visualisation%writer%write_attribute_hyperslab(visualisation%attributes(attribute_index), values, &
                                                          offset, shape, writer_status)
      call convert_status(writer_status, status, diagnostic)
   end subroutine write_visualisation_attribute_hyperslab

   subroutine end_visualisation_step(visualisation, status, diagnostic)
      type(visualisation_writer_t), intent(inout) :: visualisation
      integer, intent(out) :: status
      character(len=BUFSIZE), intent(out) :: diagnostic
      type(xdmf_status_t) :: writer_status

      call visualisation%writer%end_step(writer_status)
      call convert_status(writer_status, status, diagnostic)
   end subroutine end_visualisation_step

   subroutine flush_visualisation(visualisation, status, diagnostic)
      type(visualisation_writer_t), intent(inout) :: visualisation
      integer, intent(out) :: status
      character(len=BUFSIZE), intent(out) :: diagnostic
      type(xdmf_status_t) :: writer_status

      if (.not. associated(visualisation%writer)) then
         status = VISUALISATION_SUCCESS
         diagnostic = ''
         return
      end if
      call visualisation%writer%flush(writer_status)
      call convert_status(writer_status, status, diagnostic)
   end subroutine flush_visualisation

   subroutine close_visualisation(visualisation, status, diagnostic)
      type(visualisation_writer_t), intent(inout) :: visualisation
      integer, intent(out) :: status
      character(len=BUFSIZE), intent(out) :: diagnostic
      type(xdmf_status_t) :: writer_status

      if (.not. associated(visualisation%writer)) then
         status = VISUALISATION_SUCCESS
         diagnostic = ''
         return
      end if
      call visualisation%writer%close(writer_status)
      deallocate (visualisation%writer)
      call convert_status(writer_status, status, diagnostic)
   end subroutine close_visualisation

   logical function visualisation_is_open(visualisation)
      type(visualisation_writer_t), intent(in) :: visualisation

      visualisation_is_open = associated(visualisation%writer)
   end function visualisation_is_open

   subroutine convert_status(writer_status, status, diagnostic)
      type(xdmf_status_t), intent(in) :: writer_status
      integer, intent(out) :: status
      character(len=BUFSIZE), intent(out) :: diagnostic

      if (writer_status%is_error()) then
         status = VISUALISATION_WRITER_ERROR
         diagnostic = writer_status%message()
      else
         status = VISUALISATION_SUCCESS
         diagnostic = ''
      end if
   end subroutine convert_status

   subroutine verify_volumetric_visualisation(path, status)
      character(len=*), intent(in) :: path
      integer, intent(out) :: status

      if (file_exists(trim(path)//'.xdmf') .and. file_exists(trim(path)//'.h5')) then
         status = VISUALISATION_SUCCESS
      else
         status = VISUALISATION_IO_ERROR
      end if
   end subroutine verify_volumetric_visualisation

end module outputVisualisation_m
