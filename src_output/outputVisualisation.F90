module outputVisualisation_m
   use, intrinsic :: iso_fortran_env, only: int64, real64
    use directoryUtils_m, only: create_folder, file_exists
   use xdmf_hdf5_m, only: xdmf_writer_t, xdmf_options_t, xdmf_status_t, &
                           xdmf_grid_id_t, xdmf_attribute_id_t, XDMF_SERIES_FREQUENCY, XDMF_SERIES_TIME, &
                           XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64
   implicit none
   private

   integer, parameter, public :: VISUALISATION_SUCCESS = 0
   integer, parameter, public :: VISUALISATION_IO_ERROR = 1
   integer, parameter, public :: VISUALISATION_WRITER_ERROR = 2

     public :: publish_volumetric_visualisation, publish_frequency_slice_visualisation, verify_volumetric_visualisation

contains

    subroutine publish_volumetric_visualisation(path, grid_name, attribute_name, dimensions, origin, spacing, &
                                               time, values, status)
      character(len=*), intent(in) :: path, grid_name, attribute_name
      integer(int64), intent(in) :: dimensions(3)
      real(real64), intent(in) :: origin(3), spacing(3), time, values(:)
      integer, intent(out) :: status

      type(xdmf_writer_t) :: writer
      type(xdmf_options_t) :: options
      type(xdmf_status_t) :: writer_status
      type(xdmf_grid_id_t) :: grid
      type(xdmf_attribute_id_t) :: attribute
      integer :: directory_end, ios

      status = VISUALISATION_IO_ERROR
      directory_end = scan(trim(path), '/\\', back=.true.)
      if (directory_end > 0) then
         call create_folder(path(:directory_end - 1), ios)
         if (ios /= 0) return
      end if

      options%overwrite = .true.
      options%series_kind = XDMF_SERIES_TIME
      call writer%create(trim(path), options, writer_status)
      if (writer_status%is_error()) then
         status = VISUALISATION_WRITER_ERROR
         return
      end if
      call writer%define_uniform_grid(trim(grid_name), dimensions, origin, spacing, grid, writer_status)
      if (.not. writer_status%is_error()) then
         call writer%define_attribute(grid, trim(attribute_name), XDMF_CENTER_NODE, XDMF_ATTRIBUTE_SCALAR, &
                                      XDMF_NUMERIC_REAL64, .true., attribute, writer_status)
      end if
      if (.not. writer_status%is_error()) call writer%begin_step(time, writer_status)
      if (.not. writer_status%is_error()) call writer%write_attribute(attribute, values, writer_status)
      if (.not. writer_status%is_error()) call writer%end_step(writer_status)
      if (writer_status%is_error()) then
         call writer%close(writer_status)
         status = VISUALISATION_WRITER_ERROR
         return
      end if
      call writer%close(writer_status)
      if (writer_status%is_error()) then
         status = VISUALISATION_WRITER_ERROR
         return
      end if
       status = VISUALISATION_SUCCESS
    end subroutine publish_volumetric_visualisation

    subroutine publish_frequency_slice_visualisation(path, grid_name, dimensions, x_coordinates, y_coordinates, &
                                                     z_coordinates, frequencies, x_values, y_values, z_values, status)
       character(len=*), intent(in) :: path, grid_name
       integer(int64), intent(in) :: dimensions(3)
       real(real64), intent(in) :: x_coordinates(:), y_coordinates(:), z_coordinates(:)
       real(real64), intent(in) :: frequencies(:)
       complex(real64), intent(in) :: x_values(:, :), y_values(:, :), z_values(:, :)
       integer, intent(out) :: status

       type(xdmf_writer_t) :: writer
       type(xdmf_options_t) :: options
       type(xdmf_status_t) :: writer_status
       type(xdmf_grid_id_t) :: grid
       type(xdmf_attribute_id_t) :: x_real, x_imag, y_real, y_imag, z_real, z_imag
       integer :: directory_end, frequency_index, ios

        status = VISUALISATION_IO_ERROR
        if (size(frequencies) /= size(x_values, 1) .or. size(x_values, 2) /= product(dimensions)) return
        if (any(shape(y_values) /= shape(x_values)) .or. any(shape(z_values) /= shape(x_values))) return
        if (any([size(x_coordinates), size(y_coordinates), size(z_coordinates)] /= dimensions)) return
       directory_end = scan(trim(path), '/\\', back=.true.)
       if (directory_end > 0) then
          call create_folder(path(:directory_end - 1), ios)
          if (ios /= 0) return
       end if

       options%overwrite = .true.
        options%series_kind = XDMF_SERIES_FREQUENCY
       call writer%create(trim(path), options, writer_status)
       if (writer_status%is_error()) then
          status = VISUALISATION_WRITER_ERROR
          return
       end if
        call writer%define_rectilinear_grid(trim(grid_name), x_coordinates, y_coordinates, z_coordinates, &
                                            grid, writer_status)
       if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'X_real', XDMF_CENTER_NODE, &
          XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., x_real, writer_status)
       if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'X_imag', XDMF_CENTER_NODE, &
          XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., x_imag, writer_status)
       if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'Y_real', XDMF_CENTER_NODE, &
          XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., y_real, writer_status)
       if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'Y_imag', XDMF_CENTER_NODE, &
          XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., y_imag, writer_status)
       if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'Z_real', XDMF_CENTER_NODE, &
          XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., z_real, writer_status)
       if (.not. writer_status%is_error()) call writer%define_attribute(grid, 'Z_imag', XDMF_CENTER_NODE, &
          XDMF_ATTRIBUTE_SCALAR, XDMF_NUMERIC_REAL64, .true., z_imag, writer_status)
       do frequency_index = 1, size(frequencies)
          if (.not. writer_status%is_error()) call writer%begin_step(frequencies(frequency_index), writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(x_real, real(x_values(frequency_index, :), real64), writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(x_imag, aimag(x_values(frequency_index, :)), writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(y_real, real(y_values(frequency_index, :), real64), writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(y_imag, aimag(y_values(frequency_index, :)), writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(z_real, real(z_values(frequency_index, :), real64), writer_status)
          if (.not. writer_status%is_error()) call writer%write_attribute(z_imag, aimag(z_values(frequency_index, :)), writer_status)
          if (.not. writer_status%is_error()) call writer%end_step(writer_status)
       end do
       call writer%close(writer_status)
       if (writer_status%is_error()) then
          status = VISUALISATION_WRITER_ERROR
       else
          status = VISUALISATION_SUCCESS
       end if
    end subroutine publish_frequency_slice_visualisation

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
