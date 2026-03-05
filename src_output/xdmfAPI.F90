module mod_xdmfAPI
   use HDF5
   implicit none

   ! HDF5 constants

   integer, parameter :: dp = kind(1.0d0)


   interface write_dataset
      module procedure write_1d_dataset
      module procedure write_2d_dataset
      module procedure write_3d_dataset
   end interface

contains
   subroutine xdmf_write_header_file(unit)
      integer :: unit
      write (unit, '(A)') '<?xml version="1.0" ?>'
      write (unit, '(A)') '<Xdmf Version="3.0">'
      write (unit, '(A)') '  <Domain>'
      write (unit, '(A)') '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
   end subroutine 

   subroutine xdmf_write_footer_file(unit)
      integer :: unit
      write (unit, '(A)') '</Grid>'
      write (unit, '(A)') '</Domain>'
      write (unit, '(A)') '</Xdmf>'
   end subroutine 

   subroutine xdmf_create_grid_step_info(unit, stepName, stepValue, h5_filename, ncoords)
      !Requires file already open
      !h5_file must contain a coords field
      integer, intent(in) :: unit
      character(len=*), intent(in) :: stepName
      real, intent(in) :: stepValue
      character(len=*), intent(in) :: h5_filename
      integer, intent(in) :: ncoords

      write (unit, '(A,A,A)'       ) '<Grid Name="', trim(stepName), '" GridType="Uniform">'
      write (unit, '(A,G0,A)'      ) '<Time Value="',stepValue,'"/>'

      write (unit, '(A,I0,A)'      ) '<Topology TopologyType="Polyvertex" NumberOfElements="',ncoords,'"/>'
      write (unit, '(A)'           ) '<Geometry GeometryType="XYZ">'
      write (unit, '(A,I0,1X,I0,A)') '<DataItem Format="HDF" Dimensions="',ncoords, 3,'" NumberType="Float">'
      write (unit, '(A,A)'         ) trim(h5_filename),':/coords'
      write (unit, '(A)'           ) '</DataItem>'
      write (unit, '(A)'           ) '</Geometry>'
   end subroutine xdmf_create_grid_step_info

   subroutine xdmf_close_grid(unit)
      !Requires file already open
      integer, intent(in) :: unit
      write (unit, '(A)'     ) '</Grid>'
   end subroutine xdmf_close_grid

   subroutine xdmf_write_attribute(unit, attributeName)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: attributeName
      write (unit, '(A,A,A)'  ) '<Attribute Name="',trim(attributeName),'" AttributeType="Scalar" Center="Node">'
   end subroutine xdmf_write_attribute

   subroutine xdmf_close_attribute(unit)
      integer, intent(in) :: unit
      write (unit, '(A)'      ) '</Attribute>'
   end subroutine xdmf_close_attribute

   subroutine xdmf_write_hyperslab_data_item(unit, dimension_string)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: dimension_string !int writen in space separated format Ex: '2 20 3'
      write (unit, '(A,A,A)'  ) '<DataItem ItemType="HyperSlab" Dimensions="',trim(dimension_string),'" Type="HyperSlab">'
   end subroutine xdmf_write_hyperslab_data_item

   subroutine xdmf_write_h5_acces_data_item(unit, row_offset, column_offset, row_count, column_count)
      !Used on cases where we acces parts of an h5 array
      integer, intent(in) :: unit, row_offset, column_offset, row_count, column_count
      write (unit, '(A,A,A)'  ) '<DataItem Dimensions="3 2" Format="XML" NumberType="Int">'
      write (unit, '(I0,1X,I0)') row_offset, column_offset
      write (unit, '(I0,1X,I0)') 1, 1
      write (unit, '(I0,1X,I0)') row_count, column_count
   end subroutine xdmf_write_h5_acces_data_item

   subroutine xdmf_write_h5_data_item(unit, h5_filename, h5_data_path, dimension_string)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: h5_filename
      character(len=*), intent(in) :: h5_data_path
      character(len=*), intent(in) :: dimension_string !int writen in space separated format Ex: '2 20 3'
      write (unit, '(A,A,A)') '<DataItem Format="HDF" Dimensions="',trim(dimension_string),'" NumberType="Float">'
      write (unit, '(A,A,A)') h5_filename,':', h5_data_path 
   end subroutine xdmf_write_h5_data_item

   subroutine xdmf_close_data_item(unit)
      integer, intent(in) :: unit
      write (unit, '(A)') '</DataItem>'
   end subroutine xdmf_close_data_item

   subroutine create_h5_file(filename, file_id)
      character(len=*), intent(in) :: filename
      integer(HID_T), intent(out) :: file_id
      integer :: error

      call H5Fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)
      if (file_id < 0) then
         print *, "Error creating HDF5 file: ", trim(filename)
      end if
      if (error /= 0) then
         print *, "Error raised creating HDF5"
      end if
   end subroutine create_h5_file

   subroutine write_1d_dataset(file_id, dataset_name, data)
      integer(HID_T), intent(in) :: file_id
      character(len=*), intent(in) :: dataset_name
      real(dp), dimension(:), intent(in) :: data

      integer(HID_T) :: dataspace_id, dataset_id
      integer(HSIZE_T), dimension(1) :: dims
      integer :: error

      dims(1) = size(data)

      call H5Screate_simple_f(1, dims, dataspace_id, error)
      if (error /= 0) then
         print *, "Error creating dataspace"
      end if

      call H5Dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, &
                       dataspace_id, dataset_id, error)
      if (error /= 0) then
         print *, "Error creating dataset: ", trim(dataset_name)
      end if

      call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dims, error)

      call H5Dclose_f(dataset_id, error)
      call H5Sclose_f(dataspace_id, error)
   end subroutine write_1d_dataset

   subroutine write_2d_dataset(file_id, dataset_name, data)
      integer(HID_T), intent(in) :: file_id
      character(len=*), intent(in) :: dataset_name
      real(dp), dimension(:, :), intent(in) :: data

      integer(HID_T) :: dataspace_id = 0, dataset_id = 0
      integer(HSIZE_T), dimension(2) :: dims
      integer :: error = 0

      dims(1) = size(data, 1)
      dims(2) = size(data, 2)

      call H5Screate_simple_f(2, dims, dataspace_id, error)
      if (error /= 0) then
         print *, "Error creating dataspace"
      end if

      call H5Dcreate_f(file_id, trim(adjustl(dataset_name)), H5T_NATIVE_DOUBLE, &
                       dataspace_id, dataset_id, error)
      if (error /= 0) then
         print *, "Error creating dataset: ", trim(dataset_name)
      end if

      call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dims, error)

      call H5Dclose_f(dataset_id, error)
      call H5Sclose_f(dataspace_id, error)
   end subroutine write_2d_dataset

   subroutine write_3d_dataset(file_id, dataset_name, data)
      integer(HID_T), intent(in) :: file_id
      character(len=*), intent(in) :: dataset_name
      real(dp), dimension(:, :, :), intent(in) :: data

      integer(HID_T) :: dataspace_id, dataset_id
      integer(HSIZE_T), dimension(3) :: dims
      integer :: error

      dims(1) = size(data, 1)
      dims(2) = size(data, 2)
      dims(3) = size(data, 3)

      call H5Screate_simple_f(3, dims, dataspace_id, error)
      if (error /= 0) then
         print *, "Error creating dataspace"
      end if

      call H5Dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, &
                       dataspace_id, dataset_id, error)
      if (error /= 0) then
         print *, "Error creating dataset: ", trim(dataset_name)
      end if

      call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dims, error)

      call H5Dclose_f(dataset_id, error)
      call H5Sclose_f(dataspace_id, error)
   end subroutine write_3d_dataset

   subroutine init_extendable_2d_dataset(file_id, dataset_name, fixed_dim, chunk_rows)

      integer(HID_T), intent(in) :: file_id
      character(len=*), intent(in) :: dataset_name
      integer, intent(in) :: fixed_dim
      integer, intent(in) :: chunk_rows

      integer(HID_T) :: dataspace_id, dataset_id, plist_id
      integer(HSIZE_T), dimension(2) :: dims
      integer(HSIZE_T), dimension(2) :: maxdims
      integer(HSIZE_T), dimension(2) :: chunk_dims
      integer :: error

      ! Initial size
      dims(1) = 0
      dims(2) = fixed_dim

      maxdims(1) = H5S_UNLIMITED_F
      maxdims(2) = fixed_dim

      call H5Screate_simple_f(2, dims, dataspace_id, error, maxdims)
      call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

      chunk_dims(1) = chunk_rows
      chunk_dims(2) = fixed_dim

      call H5Pset_chunk_f(plist_id, 2, chunk_dims, error)

      call H5Dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, &
                       dataspace_id, dataset_id, error, plist_id)

      if (error /= 0) then
         print *, "Error creating dataset: ", trim(dataset_name)
      end if

      call H5Dclose_f(dataset_id, error)
      call H5Sclose_f(dataspace_id, error)
      call H5Pclose_f(plist_id, error)

   end subroutine init_extendable_2d_dataset

   subroutine append_rows_dataset(file_id, dataset_name, data)
      integer(HID_T), intent(in) :: file_id
      character(len=*), intent(in) :: dataset_name
      real(dp), dimension(:, :), intent(in) :: data  ! shape: (nrows, fixed_size)

      integer(HID_T) :: dataset_id, filespace, memspace
      integer(HSIZE_T), dimension(2) :: dims, new_dims
      integer(HSIZE_T), dimension(2) :: start, dataCount
      integer :: error
      integer :: nrows, fixed_size

      nrows = size(data, 1)
      fixed_size = size(data, 2)

      call H5Dopen_f(file_id, trim(dataset_name), dataset_id, error)
      call H5Dget_space_f(dataset_id, filespace, error)

      ! Get current dimensions
      call H5Sget_simple_extent_dims_f(filespace, dims, new_dims, error)

      ! Extend dataset by nrows
      new_dims(1) = dims(1) + nrows
      new_dims(2) = dims(2)
      call H5Dset_extent_f(dataset_id, new_dims, error)

      ! Get updated dataspace
      call H5Sclose_f(filespace, error)
      call H5Dget_space_f(dataset_id, filespace, error)

      ! Select hyperslab for the new rows
      start(1) = dims(1)        ! start at first new row
      start(2) = 0              ! start at first column
      dataCount(1) = nrows
      dataCount(2) = fixed_size

      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, start, dataCount, error)

      ! Memory dataspace
      call H5Screate_simple_f(2, dataCount, memspace, error)

      ! Write data
      call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dataCount, error, &
                      memspace, filespace)

      ! Close resources
      call H5Sclose_f(memspace, error)
      call H5Sclose_f(filespace, error)
      call H5Dclose_f(dataset_id, error)

   end subroutine append_rows_dataset

end module mod_xdmfAPI
