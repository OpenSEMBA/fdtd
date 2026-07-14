module xdmf_hdf5_backend_m
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  use hdf5
  use xdmf_model_m, only: xdmf_status_t, xdmf_options_t, &
    XDMF_ERROR_HDF5, XDMF_ERROR_CONSISTENCY, &
    XDMF_NUMERIC_REAL32, XDMF_NUMERIC_REAL64, &
    XDMF_NUMERIC_INT32, XDMF_NUMERIC_INT64, &
    numeric_type_size, product_int64, set_status_success, set_status_error

  implicit none

  private

  type, public :: hdf5_file_t
    private
    integer(HID_T) :: id = -1_HID_T
  end type hdf5_file_t

  public :: hdf_create_file
  public :: hdf_file_is_open
  public :: hdf_close_file
  public :: hdf_flush_file
  public :: hdf_create_group
  public :: hdf_write_dataset
  public :: hdf_create_series_dataset
  public :: hdf_append_series
  public :: hdf_truncate_series

  interface hdf_write_dataset
    module procedure hdf_write_dataset_r4
    module procedure hdf_write_dataset_r8
    module procedure hdf_write_dataset_i4
    module procedure hdf_write_dataset_i8
  end interface hdf_write_dataset

  interface hdf_append_series
    module procedure hdf_append_series_r4
    module procedure hdf_append_series_r8
    module procedure hdf_append_series_i4
    module procedure hdf_append_series_i8
  end interface hdf_append_series

contains

  subroutine hdf_create_file(file, path, options, status)
    type(hdf5_file_t), intent(out) :: file
    character(len=*), intent(in) :: path
    type(xdmf_options_t), intent(in) :: options
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: access_property
    integer :: access_flag, close_error, error

    call set_status_success(status)
    call h5open_f(error)
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not initialize the HDF5 Fortran interface')
      return
    end if
    access_property = -1_HID_T
    call h5pcreate_f(H5P_FILE_ACCESS_F, access_property, error)
    if (error == 0) then
      call h5pset_fclose_degree_f(access_property, H5F_CLOSE_STRONG_F, error)
    end if
    if (error /= 0) then
      if (access_property >= 0_HID_T) then
        call h5pclose_f(access_property, close_error)
      end if
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not configure HDF5 file ownership')
      return
    end if
    if (options%overwrite) then
      access_flag = H5F_ACC_TRUNC_F
    else
      access_flag = H5F_ACC_EXCL_F
    end if
    call h5fcreate_f(trim(path), access_flag, file%id, error, &
      H5P_DEFAULT_F, access_property)
    call h5pclose_f(access_property, close_error)
    if (error /= 0 .or. file%id < 0_HID_T) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not create HDF5 file: '//trim(path))
      return
    end if
    if (close_error /= 0) then
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'Could not close the HDF5 file access property list')
      call close_after_create_error(file, status)
      return
    end if

    call write_string_attribute(file%id, 'schema_name', &
      'OpenSEMBA XDMF-HDF5', status)
    if (status%is_error()) then
      call close_after_create_error(file, status)
      return
    end if
    call write_string_attribute(file%id, 'schema_version', '1.0', status)
    if (status%is_error()) then
      call close_after_create_error(file, status)
      return
    end if
  end subroutine hdf_create_file

  logical function hdf_file_is_open(file)
    type(hdf5_file_t), intent(in) :: file

    hdf_file_is_open = file%id >= 0_HID_T
  end function hdf_file_is_open

  subroutine hdf_close_file(file, status)
    type(hdf5_file_t), intent(inout) :: file
    type(xdmf_status_t), intent(out) :: status

    integer :: error

    call set_status_success(status)
    if (file%id >= 0_HID_T) then
      call h5fclose_f(file%id, error)
      if (error /= 0) then
        call set_status_error(status, XDMF_ERROR_HDF5, &
          'Could not close the HDF5 file')
      else
        file%id = -1_HID_T
      end if
    end if
  end subroutine hdf_close_file

  subroutine hdf_flush_file(file, status)
    type(hdf5_file_t), intent(in) :: file
    type(xdmf_status_t), intent(out) :: status

    integer :: error

    call set_status_success(status)
    if (file%id < 0_HID_T) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Cannot flush a closed HDF5 file')
      return
    end if

    call h5fflush_f(file%id, H5F_SCOPE_GLOBAL_F, error)
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not flush the HDF5 file')
    end if
  end subroutine hdf_flush_file

  subroutine hdf_create_group(file, path, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: group_id
    integer :: error

    call set_status_success(status)
    call h5gcreate_f(file%id, trim(path), group_id, error)
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not create HDF5 group: '//trim(path))
      return
    end if

    call h5gclose_f(group_id, error)
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not close HDF5 group: '//trim(path))
    end if
  end subroutine hdf_create_group

  subroutine hdf_write_dataset_r4(file, path, data, shape, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(real32), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataset_id
    integer(HSIZE_T) :: buffer_dims(1)
    integer :: error

    call create_fixed_dataset(file, path, shape, &
      hdf_datatype(XDMF_NUMERIC_REAL32), &
      dataset_id, status)
    if (status%is_error()) return

    buffer_dims(1) = int(size(data), HSIZE_T)
    call h5dwrite_f(dataset_id, hdf_datatype(XDMF_NUMERIC_REAL32), &
      data, buffer_dims, error)
    call finish_fixed_write(file, dataset_id, path, error, status)
  end subroutine hdf_write_dataset_r4

  subroutine hdf_write_dataset_r8(file, path, data, shape, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(real64), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataset_id
    integer(HSIZE_T) :: buffer_dims(1)
    integer :: error

    call create_fixed_dataset(file, path, shape, &
      hdf_datatype(XDMF_NUMERIC_REAL64), &
      dataset_id, status)
    if (status%is_error()) return

    buffer_dims(1) = int(size(data), HSIZE_T)
    call h5dwrite_f(dataset_id, hdf_datatype(XDMF_NUMERIC_REAL64), &
      data, buffer_dims, error)
    call finish_fixed_write(file, dataset_id, path, error, status)
  end subroutine hdf_write_dataset_r8

  subroutine hdf_write_dataset_i4(file, path, data, shape, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(int32), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataset_id
    integer(HSIZE_T) :: buffer_dims(1)
    integer :: error

    call create_fixed_dataset(file, path, shape, &
      hdf_datatype(XDMF_NUMERIC_INT32), &
      dataset_id, status)
    if (status%is_error()) return

    buffer_dims(1) = int(size(data), HSIZE_T)
    call h5dwrite_f(dataset_id, hdf_datatype(XDMF_NUMERIC_INT32), &
      data, buffer_dims, error)
    call finish_fixed_write(file, dataset_id, path, error, status)
  end subroutine hdf_write_dataset_i4

  subroutine hdf_write_dataset_i8(file, path, data, shape, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(int64), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataset_id
    integer(HSIZE_T) :: buffer_dims(1)
    integer :: error

    call create_fixed_dataset(file, path, shape, &
      hdf_datatype(XDMF_NUMERIC_INT64), &
      dataset_id, status)
    if (status%is_error()) return

    buffer_dims(1) = int(size(data), HSIZE_T)
    call h5dwrite_f(dataset_id, hdf_datatype(XDMF_NUMERIC_INT64), &
      data, buffer_dims, error)
    call finish_fixed_write(file, dataset_id, path, error, status)
  end subroutine hdf_write_dataset_i8

  subroutine hdf_create_series_dataset(file, path, numeric_type, shape, &
      options, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer, intent(in) :: numeric_type
    integer(int64), intent(in) :: shape(:)
    type(xdmf_options_t), intent(in) :: options
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataspace_id, dataset_id, property_id, datatype
    integer(HSIZE_T), allocatable :: dims(:), maxdims(:), chunks(:)
    integer :: close_error, error, rank

    call set_status_success(status)
    dataspace_id = -1_HID_T
    dataset_id = -1_HID_T
    property_id = -1_HID_T
    rank = size(shape) + 1
    allocate(dims(rank), maxdims(rank), chunks(rank))
    if (size(shape) > 0) then
      dims(:rank - 1) = int(shape, HSIZE_T)
      maxdims(:rank - 1) = int(shape, HSIZE_T)
      call choose_chunk_shape(shape, numeric_type, &
        options%chunk_target_bytes, chunks(:rank - 1))
      chunks(rank) = 1_HSIZE_T
    else
      chunks(rank) = 256_HSIZE_T
    end if
    dims(rank) = 0_HSIZE_T
    maxdims(rank) = H5S_UNLIMITED_F

    call h5screate_simple_f(rank, dims, dataspace_id, error, maxdims)
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not create series dataspace: '//trim(path))
      return
    end if

    call h5pcreate_f(H5P_DATASET_CREATE_F, property_id, error)
    if (error == 0) call h5pset_chunk_f(property_id, rank, chunks, error)
    if (error == 0 .and. options%compression_level > 0) then
      call h5pset_deflate_f(property_id, options%compression_level, error)
    end if
    if (error /= 0) then
      if (property_id >= 0_HID_T) then
        call h5pclose_f(property_id, close_error)
      end if
      call h5sclose_f(dataspace_id, close_error)
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not configure series dataset: '//trim(path))
      return
    end if

    datatype = hdf_datatype(numeric_type)
    call h5dcreate_f(file%id, trim(path), datatype, dataspace_id, &
      dataset_id, error, property_id)
    if (error /= 0) then
      call h5pclose_f(property_id, close_error)
      call h5sclose_f(dataspace_id, close_error)
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not create series dataset: '//trim(path))
      return
    end if

    call h5dclose_f(dataset_id, error)
    call h5pclose_f(property_id, close_error)
    if (error == 0) error = close_error
    call h5sclose_f(dataspace_id, close_error)
    if (error == 0) error = close_error
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'Could not close series dataset resources: '//trim(path))
    end if
  end subroutine hdf_create_series_dataset

  subroutine hdf_append_series_r4(file, path, data, shape, committed_steps, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(real32), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    type(xdmf_status_t), intent(out) :: status

    call append_series_r4_impl(file, path, data, shape, committed_steps, &
      hdf_datatype(XDMF_NUMERIC_REAL32), status)
  end subroutine hdf_append_series_r4

  subroutine hdf_append_series_r8(file, path, data, shape, committed_steps, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(real64), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    type(xdmf_status_t), intent(out) :: status

    call append_series_r8_impl(file, path, data, shape, committed_steps, &
      hdf_datatype(XDMF_NUMERIC_REAL64), status)
  end subroutine hdf_append_series_r8

  subroutine hdf_append_series_i4(file, path, data, shape, committed_steps, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(int32), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    type(xdmf_status_t), intent(out) :: status

    call append_series_i4_impl(file, path, data, shape, committed_steps, &
      hdf_datatype(XDMF_NUMERIC_INT32), status)
  end subroutine hdf_append_series_i4

  subroutine hdf_append_series_i8(file, path, data, shape, committed_steps, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(int64), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    type(xdmf_status_t), intent(out) :: status

    call append_series_i8_impl(file, path, data, shape, committed_steps, &
      hdf_datatype(XDMF_NUMERIC_INT64), status)
  end subroutine hdf_append_series_i8

  subroutine append_series_r4_impl(file, path, data, shape, committed_steps, &
      datatype, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(real32), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    integer(HID_T), intent(in) :: datatype
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataset_id, filespace, memspace
    integer(HSIZE_T), allocatable :: new_dims(:), offset(:), count(:)
    integer(HSIZE_T) :: mem_dims(1)
    integer :: error

    call prepare_append(file, path, shape, committed_steps, dataset_id, &
      filespace, memspace, new_dims, offset, count, mem_dims, status)
    if (status%is_error()) return
    call h5dwrite_f(dataset_id, datatype, data, mem_dims, error, &
      memspace, filespace)
    call finish_append(dataset_id, filespace, memspace, path, shape, &
      committed_steps, error, status)
  end subroutine append_series_r4_impl

  subroutine append_series_r8_impl(file, path, data, shape, committed_steps, &
      datatype, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    real(real64), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    integer(HID_T), intent(in) :: datatype
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataset_id, filespace, memspace
    integer(HSIZE_T), allocatable :: new_dims(:), offset(:), count(:)
    integer(HSIZE_T) :: mem_dims(1)
    integer :: error

    call prepare_append(file, path, shape, committed_steps, dataset_id, &
      filespace, memspace, new_dims, offset, count, mem_dims, status)
    if (status%is_error()) return
    call h5dwrite_f(dataset_id, datatype, data, mem_dims, error, &
      memspace, filespace)
    call finish_append(dataset_id, filespace, memspace, path, shape, &
      committed_steps, error, status)
  end subroutine append_series_r8_impl

  subroutine append_series_i4_impl(file, path, data, shape, committed_steps, &
      datatype, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(int32), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    integer(HID_T), intent(in) :: datatype
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataset_id, filespace, memspace
    integer(HSIZE_T), allocatable :: new_dims(:), offset(:), count(:)
    integer(HSIZE_T) :: mem_dims(1)
    integer :: error

    call prepare_append(file, path, shape, committed_steps, dataset_id, &
      filespace, memspace, new_dims, offset, count, mem_dims, status)
    if (status%is_error()) return
    call h5dwrite_f(dataset_id, datatype, data, mem_dims, error, &
      memspace, filespace)
    call finish_append(dataset_id, filespace, memspace, path, shape, &
      committed_steps, error, status)
  end subroutine append_series_i4_impl

  subroutine append_series_i8_impl(file, path, data, shape, committed_steps, &
      datatype, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(int64), intent(in) :: data(:)
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    integer(HID_T), intent(in) :: datatype
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataset_id, filespace, memspace
    integer(HSIZE_T), allocatable :: new_dims(:), offset(:), count(:)
    integer(HSIZE_T) :: mem_dims(1)
    integer :: error

    call prepare_append(file, path, shape, committed_steps, dataset_id, &
      filespace, memspace, new_dims, offset, count, mem_dims, status)
    if (status%is_error()) return
    call h5dwrite_f(dataset_id, datatype, data, mem_dims, error, &
      memspace, filespace)
    call finish_append(dataset_id, filespace, memspace, path, shape, &
      committed_steps, error, status)
  end subroutine append_series_i8_impl

  subroutine hdf_truncate_series(file, path, shape, committed_steps, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataset_id
    integer(HSIZE_T), allocatable :: dims(:)
    integer :: error, close_error, rank

    call set_status_success(status)
    dataset_id = -1_HID_T
    close_error = 0
    rank = size(shape) + 1
    allocate(dims(rank))
    if (size(shape) > 0) dims(:rank - 1) = int(shape, HSIZE_T)
    dims(rank) = int(committed_steps, HSIZE_T)

    call h5dopen_f(file%id, trim(path), dataset_id, error)
    if (error == 0) call h5dset_extent_f(dataset_id, dims, error)
    if (dataset_id >= 0_HID_T) call h5dclose_f(dataset_id, close_error)
    if (close_error /= 0) then
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'Rollback succeeded but dataset close failed: '//trim(path))
    else if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not roll back series dataset: '//trim(path))
    end if
  end subroutine hdf_truncate_series

  subroutine create_fixed_dataset(file, path, shape, datatype, dataset_id, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(int64), intent(in) :: shape(:)
    integer(HID_T), intent(in) :: datatype
    integer(HID_T), intent(out) :: dataset_id
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataspace_id
    integer(HSIZE_T), allocatable :: dims(:)
    integer :: dataset_close_error, delete_error, error, close_error

    call set_status_success(status)
    dataset_id = -1_HID_T
    dataspace_id = -1_HID_T
    close_error = 0
    if (size(shape) == 0 .or. any(shape <= 0_int64)) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'A fixed HDF5 dataset must have positive dimensions: '//trim(path))
      return
    end if

    allocate(dims(size(shape)))
    dims = int(shape, HSIZE_T)
    call h5screate_simple_f(size(shape), dims, dataspace_id, error)
    if (error == 0) then
      call h5dcreate_f(file%id, trim(path), datatype, dataspace_id, &
        dataset_id, error)
    end if
    if (dataspace_id >= 0_HID_T) call h5sclose_f(dataspace_id, close_error)
    if (close_error /= 0) then
      dataset_close_error = 0
      delete_error = 0
      if (dataset_id >= 0_HID_T) then
        call h5dclose_f(dataset_id, dataset_close_error)
        if (dataset_close_error == 0) then
          dataset_id = -1_HID_T
          call h5ldelete_f(file%id, trim(path), delete_error)
        end if
      end if
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'Could not close fixed dataset resources: '//trim(path))
      return
    end if
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not create HDF5 dataset: '//trim(path))
    end if
  end subroutine create_fixed_dataset

  subroutine finish_fixed_write(file, dataset_id, path, write_error, status)
    type(hdf5_file_t), intent(in) :: file
    integer(HID_T), intent(inout) :: dataset_id
    character(len=*), intent(in) :: path
    integer, intent(in) :: write_error
    type(xdmf_status_t), intent(inout) :: status

    integer :: close_error, delete_error

    delete_error = 0
    call h5dclose_f(dataset_id, close_error)
    dataset_id = -1_HID_T
    if (write_error /= 0 .or. close_error /= 0) then
      call h5ldelete_f(file%id, trim(path), delete_error)
    end if
    if (write_error /= 0) then
      if (delete_error /= 0) then
        call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
          'Dataset write and cleanup failed: '//trim(path))
      else
        call set_status_error(status, XDMF_ERROR_HDF5, &
          'Could not write HDF5 dataset: '//trim(path))
      end if
    else if (close_error /= 0) then
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'Could not close HDF5 dataset: '//trim(path))
    end if
  end subroutine finish_fixed_write

  subroutine prepare_append(file, path, shape, committed_steps, dataset_id, &
      filespace, memspace, new_dims, offset, count, mem_dims, status)
    type(hdf5_file_t), intent(in) :: file
    character(len=*), intent(in) :: path
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps
    integer(HID_T), intent(out) :: dataset_id, filespace, memspace
    integer(HSIZE_T), allocatable, intent(out) :: new_dims(:), offset(:), count(:)
    integer(HSIZE_T), intent(out) :: mem_dims(1)
    type(xdmf_status_t), intent(out) :: status

    integer(HSIZE_T), allocatable :: dims(:), maxdims(:)
    integer :: error, close_error, rank
    logical :: close_failed

    call set_status_success(status)
    dataset_id = -1_HID_T
    filespace = -1_HID_T
    memspace = -1_HID_T
    rank = size(shape) + 1
    allocate(dims(rank), maxdims(rank), new_dims(rank), offset(rank), count(rank))

    call h5dopen_f(file%id, trim(path), dataset_id, error)
    if (error == 0) call h5dget_space_f(dataset_id, filespace, error)
    if (error == 0) then
      call h5sget_simple_extent_dims_f(filespace, dims, maxdims, error)
      if (error >= 0) error = 0
    end if
    if (error /= 0) then
      call close_append_handles(dataset_id, filespace, memspace, close_failed)
      if (close_failed) then
        call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
          'Series open and cleanup failed: '//trim(path))
      else
        call set_status_error(status, XDMF_ERROR_HDF5, &
          'Could not open series dataset: '//trim(path))
      end if
      return
    end if

    if (size(shape) > 0) then
      if (any(dims(:rank - 1) /= int(shape, HSIZE_T))) then
        call close_append_handles(dataset_id, filespace, memspace, close_failed)
        call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
          'Series dataset shape changed: '//trim(path))
        return
      end if
    end if
    if (dims(rank) /= int(committed_steps, HSIZE_T)) then
      call close_append_handles(dataset_id, filespace, memspace, close_failed)
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'Series dataset extent is not synchronized: '//trim(path))
      return
    end if

    new_dims = dims
    new_dims(rank) = dims(rank) + 1_HSIZE_T
    call h5dset_extent_f(dataset_id, new_dims, error)
    call h5sclose_f(filespace, close_error)
    if (close_error /= 0) then
      call close_append_handles(dataset_id, filespace, memspace, close_failed)
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'Could not close the original series dataspace: '//trim(path))
      return
    end if
    filespace = -1_HID_T
    if (error == 0) call h5dget_space_f(dataset_id, filespace, error)

    offset = 0_HSIZE_T
    offset(rank) = dims(rank)
    count = new_dims
    count(rank) = 1_HSIZE_T
    if (error == 0) then
      call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, &
        offset, count, error)
    end if

    mem_dims(1) = int(max(1_int64, product_int64(shape)), HSIZE_T)
    if (error == 0) call h5screate_simple_f(1, mem_dims, memspace, error)
    if (error /= 0) then
      call h5dset_extent_f(dataset_id, dims, close_error)
      call close_append_handles(dataset_id, filespace, memspace, close_failed)
      if (close_error /= 0 .or. close_failed) then
        call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
          'Series setup and rollback failed: '//trim(path))
      else
        call set_status_error(status, XDMF_ERROR_HDF5, &
          'Could not extend series dataset: '//trim(path))
      end if
    end if
  end subroutine prepare_append

  subroutine finish_append(dataset_id, filespace, memspace, path, shape, &
      committed_steps, write_error, status)
    integer(HID_T), intent(inout) :: dataset_id, filespace, memspace
    character(len=*), intent(in) :: path
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: committed_steps, write_error
    type(xdmf_status_t), intent(inout) :: status

    integer(HSIZE_T), allocatable :: rollback_dims(:)
    integer :: error, rank
    logical :: close_failed

    rank = size(shape) + 1
    if (write_error /= 0) then
      allocate(rollback_dims(rank))
      if (size(shape) > 0) rollback_dims(:rank - 1) = int(shape, HSIZE_T)
      rollback_dims(rank) = int(committed_steps, HSIZE_T)
      call h5dset_extent_f(dataset_id, rollback_dims, error)
      if (error /= 0) then
        call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
          'Series append and rollback failed: '//trim(path))
      else
        call set_status_error(status, XDMF_ERROR_HDF5, &
          'Could not append series dataset: '//trim(path))
      end if
    end if
    call close_append_handles(dataset_id, filespace, memspace, close_failed)
    if (close_failed) then
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'Could not close series dataset resources: '//trim(path))
    end if
  end subroutine finish_append

  subroutine close_append_handles(dataset_id, filespace, memspace, close_failed)
    integer(HID_T), intent(inout) :: dataset_id, filespace, memspace
    logical, intent(out) :: close_failed
    integer :: error

    close_failed = .false.
    if (memspace >= 0_HID_T) then
      call h5sclose_f(memspace, error)
      if (error /= 0) close_failed = .true.
    end if
    if (filespace >= 0_HID_T) then
      call h5sclose_f(filespace, error)
      if (error /= 0) close_failed = .true.
    end if
    if (dataset_id >= 0_HID_T) then
      call h5dclose_f(dataset_id, error)
      if (error /= 0) close_failed = .true.
    end if
    memspace = -1_HID_T
    filespace = -1_HID_T
    dataset_id = -1_HID_T
  end subroutine close_append_handles

  subroutine choose_chunk_shape(shape, numeric_type, target_bytes, chunks)
    integer(int64), intent(in) :: shape(:)
    integer, intent(in) :: numeric_type
    integer(int64), intent(in) :: target_bytes
    integer(HSIZE_T), intent(out) :: chunks(:)

    integer :: largest(1)
    integer(int64) :: chunk_elements, target_elements

    chunks = int(max(shape, 1_int64), HSIZE_T)
    target_elements = max(1_int64, target_bytes / &
      max(1_int64, numeric_type_size(numeric_type)))
    do
      chunk_elements = product_int64(int(chunks, int64))
      if (chunk_elements >= 0_int64 .and. &
          chunk_elements <= target_elements) exit
      largest = maxloc(chunks)
      chunks(largest(1)) = max(1_HSIZE_T, (chunks(largest(1)) + 1_HSIZE_T) / 2_HSIZE_T)
    end do
  end subroutine choose_chunk_shape

  integer(HID_T) function hdf_datatype(numeric_type)
    integer, intent(in) :: numeric_type

    select case (numeric_type)
    case (XDMF_NUMERIC_REAL32)
      hdf_datatype = h5kind_to_type(real32, H5_REAL_KIND)
    case (XDMF_NUMERIC_REAL64)
      hdf_datatype = h5kind_to_type(real64, H5_REAL_KIND)
    case (XDMF_NUMERIC_INT32)
      hdf_datatype = h5kind_to_type(int32, H5_INTEGER_KIND)
    case (XDMF_NUMERIC_INT64)
      hdf_datatype = h5kind_to_type(int64, H5_INTEGER_KIND)
    case default; hdf_datatype = -1_HID_T
    end select
  end function hdf_datatype

  subroutine write_string_attribute(location_id, name, value, status)
    integer(HID_T), intent(in) :: location_id
    character(len=*), intent(in) :: name, value
    type(xdmf_status_t), intent(out) :: status

    integer(HID_T) :: dataspace_id, datatype_id, attribute_id
    integer(HSIZE_T) :: dims(1)
    integer(SIZE_T) :: string_size
    integer :: error, close_error

    call set_status_success(status)
    dataspace_id = -1_HID_T
    datatype_id = -1_HID_T
    attribute_id = -1_HID_T
    dims(1) = 1_HSIZE_T
    string_size = int(max(1, len_trim(value)), SIZE_T)
    call h5screate_f(H5S_SCALAR_F, dataspace_id, error)
    if (error == 0) call h5tcopy_f(H5T_FORTRAN_S1, datatype_id, error)
    if (error == 0) call h5tset_size_f(datatype_id, string_size, error)
    if (error == 0) then
      call h5acreate_f(location_id, trim(name), datatype_id, dataspace_id, &
        attribute_id, error)
    end if
    if (error == 0) then
      call h5awrite_f(attribute_id, datatype_id, trim(value), dims, error)
    end if
    if (attribute_id >= 0_HID_T) call h5aclose_f(attribute_id, close_error)
    if (datatype_id >= 0_HID_T) call h5tclose_f(datatype_id, close_error)
    if (dataspace_id >= 0_HID_T) call h5sclose_f(dataspace_id, close_error)
    if (error /= 0) then
      call set_status_error(status, XDMF_ERROR_HDF5, &
        'Could not write HDF5 schema attribute: '//trim(name))
    end if
  end subroutine write_string_attribute

  subroutine close_after_create_error(file, status)
    type(hdf5_file_t), intent(inout) :: file
    type(xdmf_status_t), intent(inout) :: status
    integer :: error

    error = 0
    if (file%id >= 0_HID_T) call h5fclose_f(file%id, error)
    if (error == 0) then
      file%id = -1_HID_T
    else
      call set_status_error(status, XDMF_ERROR_CONSISTENCY, &
        'HDF5 file creation failed and its handle could not be closed')
    end if
  end subroutine close_after_create_error

end module xdmf_hdf5_backend_m
