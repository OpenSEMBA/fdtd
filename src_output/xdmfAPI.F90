module xdmfAPI_m
   use HDF5
   implicit none

   ! HDF5 constants

   integer, parameter :: dp = kind(1.0d0)


   interface h5_write_dataset
      module procedure write_1d_dataset
      module procedure write_2d_dataset
      module procedure write_3d_dataset
   end interface

   interface h5_init_extendable_dataset
      module procedure init_extendable_1D_dataset
      module procedure init_extendable_2D_dataset
      module procedure init_extendable_4D_dataset
   end interface

   interface h5_append_rows_to_dataset
      module procedure append_rows_to_1d_dataset
      module procedure append_rows_to_2d_dataset
      module procedure append_rows_to_4d_dataset
   end interface

contains
   subroutine xdmf_write_header_file(unit, gridName)
      integer, intent(in) :: unit
      character (len=*), intent(in) :: gridName
      write (unit, '(A)'    ) '<?xml version="1.0" ?>'
      write (unit, '(A)'    ) '<Xdmf Version="3.0">'
      write (unit, '(A)'    ) '<Domain>'
      write (unit, '(A,A,A)') '<Grid Name="',gridName,'" GridType="Uniform">'
   end subroutine

   subroutine xdmf_write_footer_file(unit)
      integer :: unit
      write (unit, '(A)') '</Grid>'
      write (unit, '(A)') '</Domain>'
      write (unit, '(A)') '</Xdmf>'
   end subroutine

   subroutine xdmf_write_topology(unit, dimensions)
      integer, intent(in) :: unit
      integer, intent(in), dimension(3) :: dimensions
      write (unit, '(A,I0,1X,I0,1X,I0,A)') '<Topology TopologyType="3DRectMesh" Dimensions="',dimensions(3), dimensions(2), dimensions(1),'"/>' !(z,y,x)
   end subroutine

   subroutine xdmf_write_geometry(unit, dimensions, h5_filename)
      integer, intent(in) :: unit
      integer, intent(in), dimension(3) :: dimensions
      character (len=*), intent(in) :: h5_filename

      write(unit, '(A)') '<Geometry GeometryType="VXVYVZ">'
      ! X coordinates
      write(unit, '(A,I0,A)') '<DataItem Dimensions="', dimensions(1), '" Format="HDF" NumberType="Float" Precision="8">'
      write(unit, '(A,A)'   ) trim(h5_filename),':/coordsX'
      write(unit, '(A)'     ) '</DataItem>'
      ! Y coordinates
      write(unit, '(A,I0,A)') '<DataItem Dimensions="', dimensions(2), '" Format="HDF" NumberType="Float" Precision="8">'
      write(unit, '(A,A)'   ) trim(h5_filename),':/coordsY'
      write(unit, '(A)'     ) '</DataItem>'
      ! Z coordinates
      write(unit, '(A,I0,A)') '<DataItem Dimensions="', dimensions(3), '" Format="HDF" NumberType="Float" Precision="8">'
      write(unit, '(A,A)'   )  trim(h5_filename),':/coordsZ'
      write(unit, '(A)'     ) '</DataItem>'
      ! Close JOIN and Geometry
      write(unit, '(A)'     ) '</Geometry>'
   end subroutine

   subroutine xdmf_write_time_array(unit, nTime, h5_filename)
      integer, intent(in) :: unit
      integer, intent(in) :: nTime
      character (len=*), intent(in) :: h5_filename
      write(unit, '(A)'     ) '<Time TimeType="List">'
      write(unit, '(A,I0,A)') '<DataItem Dimensions="',nTime,'" NumberType="Float" Format="HDF">'
      write(unit, '(A,A)'   ) trim(h5_filename),':/times'
      write(unit, '(A)'     ) '</DataItem>'
      write(unit, '(A)'     ) '</Time>'
   end subroutine

   subroutine xdmf_write_scalar_attribute(unit, dimensions, h5_filename, attributeName)
      integer, intent(in) :: unit
      integer, intent(in), dimension(4) :: dimensions
      character (len=*), intent(in) :: h5_filename
      character (len=*), intent(in) :: attributeName

      write(unit, '(A,A,A)' ) '<Attribute Name="',trim(attributeName),'" Center="Node" AttributeType="Scalar">'
      write(unit, '(A,I0,1X,I0,1X,I0,1X,I0,A)') '<DataItem Dimensions="',dimensions(1), dimensions(2), dimensions(3), dimensions(4),'" NumberType="Float" Format="HDF">'
      write(unit, '(A,A,A)' )  trim(h5_filename), ':/' ,trim(attributeName)
      write(unit, '(A)'     ) '</DataItem>'
      write(unit, '(A)'     ) '</Attribute>'
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

   subroutine xdmf_write_hyperslab_data_item(unit, offset, stride, selection, h5_dimension, h5_file_path, h5_data_path)
      integer, intent(in) :: unit, offset(4), stride(4), selection(4), h5_dimension(4)
      character(len=*), intent(in) :: h5_file_path, h5_data_path
      write (unit, '(A,I0,1X,I0,1X,I0,1X,I0,A)') '<DataItem ItemType="HyperSlab" Dimensions="',selection(1), selection(2), selection(3), selection(4),'Format="XML">'
      call xdmf_write_h5_access_data_item(unit, offset, stride, selection)
      call xdmf_write_h5_data_path(unit, h5_dimension, h5_file_path, h5_data_path)
      call xdmf_close_data_item(unit)
   end subroutine xdmf_write_hyperslab_data_item

   subroutine xdmf_write_h5_access_data_item(unit, offset, stride, selection)
      !Used on cases where we acces parts of an h5 array
      integer, intent(in) :: unit, offset(4), stride(4), selection(4)
      write (unit, '(A)'                   ) '<DataItem Dimensions="3 4" Format="XML">'
      write (unit, '(I0,1X,I0,1X,I0,1X,I0)') offset(1), offset(2), offset(3), offset(4) 
      write (unit, '(I0,1X,I0,1X,I0,1X,I0)') stride(1), stride(2), stride(3), stride(4) 
      write (unit, '(I0,1X,I0,1X,I0,1X,I0)') selection(1), selection(2), selection(3), selection(4) 
      call xdmf_close_data_item(unit)
   end subroutine xdmf_write_h5_access_data_item

   subroutine xdmf_write_h5_data_path(unit, h5_dimension, h5_file_path, h5_data_path)
      integer, intent(in) :: unit, h5_dimension(4)
      character(len=*), intent(in) :: h5_file_path, h5_data_path
      write (unit, '(A,A,A)') '<DataItem Format="HDF" " NumberType="Float" Dimensions="',h5_dimension(1), h5_dimension(2), h5_dimension(3), h5_dimension(4),'>'
      write (unit, '(A,A,A)') h5_file_path,':', h5_data_path
      call xdmf_close_data_item(unit)
   end subroutine

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
   end subroutine

   subroutine h5_create_rectilinear_coords_dataset(file_id, xsteps, ysteps, zsteps)
      integer(HID_T), intent(in) :: file_id
      real(dp), dimension(:), intent(in) :: xsteps, ysteps, zsteps
      call h5_write_dataset(file_id, 'coordsX', xsteps)
      call h5_write_dataset(file_id, 'coordsY', ysteps)
      call h5_write_dataset(file_id, 'coordsZ', zsteps)
   end subroutine

   subroutine h5_create_times_dataset(file_id, chunk_size)
      integer(HID_T), intent(in) :: file_id
      integer, intent(in) :: chunk_size
      call h5_init_extendable_dataset(file_id, 'times', chunk_size)
   end subroutine

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
   end subroutine

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
   end subroutine

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
   end subroutine


   subroutine init_extendable_1d_dataset(file_id, dataset_name, chunk_rows)

      integer(HID_T), intent(in) :: file_id
      character(len=*), intent(in) :: dataset_name
      integer, intent(in) :: chunk_rows

      integer(HID_T) :: dataspace_id, dataset_id, plist_id
      integer(HSIZE_T), dimension(2) :: dims
      integer(HSIZE_T), dimension(2) :: maxdims
      integer(HSIZE_T), dimension(2) :: chunk_dims
      integer :: error

      ! Initial size
      dims(1) = 0

      maxdims(1) = H5S_UNLIMITED_F

      call H5Screate_simple_f(1, dims, dataspace_id, error, maxdims)
      call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

      chunk_dims(1) = chunk_rows

      call H5Pset_chunk_f(plist_id, 1, chunk_dims, error)

      call H5Dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, &
         dataspace_id, dataset_id, error, plist_id)

      if (error /= 0) then
         print *, "Error creating dataset: ", trim(dataset_name)
      end if

      call H5Dclose_f(dataset_id, error)
      call H5Sclose_f(dataspace_id, error)
      call H5Pclose_f(plist_id, error)

   end subroutine

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
      dims(1) = fixed_dim
      dims(2) = 0

      maxdims(1) = fixed_dim
      maxdims(2) = H5S_UNLIMITED_F

      call H5Screate_simple_f(2, dims, dataspace_id, error, maxdims)
      call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

      chunk_dims(1) = fixed_dim
      chunk_dims(2) = chunk_rows

      call H5Pset_chunk_f(plist_id, 2, chunk_dims, error)

      call H5Dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, &
         dataspace_id, dataset_id, error, plist_id)

      if (error /= 0) then
         print *, "Error creating dataset: ", trim(dataset_name)
      end if

      call H5Dclose_f(dataset_id, error)
      call H5Sclose_f(dataspace_id, error)
      call H5Pclose_f(plist_id, error)

   end subroutine

   subroutine init_extendable_4d_dataset(file_id, dataset_name, fixed_dim, chunk_rows)

      integer(HID_T), intent(in) :: file_id
      character(len=*), intent(in) :: dataset_name
      integer, intent(in), dimension(3) :: fixed_dim
      integer, intent(in) :: chunk_rows

      integer(HID_T) :: dataspace_id, dataset_id, plist_id
      integer(HSIZE_T), dimension(4) :: dims
      integer(HSIZE_T), dimension(4) :: maxdims
      integer(HSIZE_T), dimension(4) :: chunk_dims
      integer :: error

      ! Initial size
      dims(1) = fixed_dim(1)
      dims(2) = fixed_dim(2)
      dims(3) = fixed_dim(3)
      dims(4) = 0

      maxdims(1) = fixed_dim(1)
      maxdims(2) = fixed_dim(2)
      maxdims(3) = fixed_dim(3)
      maxdims(4) = H5S_UNLIMITED_F

      call H5Screate_simple_f(4, dims, dataspace_id, error, maxdims)
      call H5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)

      chunk_dims(1) = fixed_dim(1)
      chunk_dims(2) = fixed_dim(2)
      chunk_dims(3) = fixed_dim(3)
      chunk_dims(4) = chunk_rows

      call H5Pset_chunk_f(plist_id, 4, chunk_dims, error)

      call H5Dcreate_f(file_id, trim(dataset_name), H5T_NATIVE_DOUBLE, &
         dataspace_id, dataset_id, error, plist_id)

      if (error /= 0) then
         print *, "Error creating dataset: ", trim(dataset_name)
      end if

      call H5Dclose_f(dataset_id, error)
      call H5Sclose_f(dataspace_id, error)
      call H5Pclose_f(plist_id, error)

   end subroutine

   subroutine append_rows_to_1d_dataset(file_id, dataset_name, data)

      integer(HID_T), intent(in) :: file_id
      character(len=*), intent(in) :: dataset_name
      real(dp), dimension(:), intent(in) :: data

      integer(HID_T) :: dataset_id, filespace, memspace
      integer(HSIZE_T) :: dims(1), maxdims(1)
      integer(HSIZE_T) :: new_dims(1)
      integer(HSIZE_T) :: offsetHyperSlab(1), dimsHyperSlab(1)
      integer :: error
      integer :: nrows
      integer(HSIZE_T) :: nt

      nt = size(data)

      call H5Dopen_f(file_id, trim(dataset_name), dataset_id, error)

      call H5Dget_space_f(dataset_id, filespace, error)

      call H5Sget_simple_extent_dims_f(filespace, dims, maxdims, error)

      new_dims(1) = dims(1) + nt

      call H5Dset_extent_f(dataset_id, new_dims, error)

      call H5Sclose_f(filespace, error)

      call H5Dget_space_f(dataset_id, filespace, error)

      offsetHyperSlab(1) = dims(1)

      dimsHyperSlab(1) = nt

      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offsetHyperSlab, dimsHyperSlab, error)

      call H5Screate_simple_f(1, dimsHyperSlab, memspace, error)

      call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dimsHyperSlab, error, memspace, filespace)

      call H5Sclose_f(memspace, error)
      call H5Sclose_f(filespace, error)
      call H5Dclose_f(dataset_id, error)

   end subroutine

   subroutine append_rows_to_2d_dataset(file_id, dataset_name, data)
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

   end subroutine

   subroutine append_rows_to_4d_dataset(file_id, dataset_name, data)

      integer(HID_T), intent(in) :: file_id
      character(len=*), intent(in) :: dataset_name
      real(dp), dimension(:, :, :, :), intent(in) :: data

      integer(HID_T) :: dataset_id, filespace, memspace
      integer(HSIZE_T), dimension(4) :: dims, maxdims
      integer(HSIZE_T), dimension(4) :: new_dims
      integer(HSIZE_T), dimension(4) :: offsetHyperSlab, dimsHyperSlab
      integer :: error
      integer :: nrows
      integer(HSIZE_T) :: nx, ny, nz, nt
      integer, dimension(3) :: fixed_size

      nx = size(data, 1)
      ny = size(data, 2)
      nz = size(data, 3)
      nt = size(data, 4)

      call H5Dopen_f(file_id, trim(dataset_name), dataset_id, error)

      call H5Dget_space_f(dataset_id, filespace, error)

      call H5Sget_simple_extent_dims_f(filespace, dims, maxdims, error)

      new_dims = dims
      new_dims(4) = dims(4) + nt

      if (dims(1) /= nx .or. dims(2) /= ny .or. dims(3) /= nz) then
        print *, "Error: dataset spatial dimensions do not match input data."
        stop
      end if

      call H5Dset_extent_f(dataset_id, new_dims, error)

      call H5Sclose_f(filespace, error)

      call H5Dget_space_f(dataset_id, filespace, error)

      offsetHyperSlab = 0_HSIZE_T
      offsetHyperSlab(4) = dims(4)

      dimsHyperSlab = (/ nx, ny, nz, nt /)

      call H5Sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, offsetHyperSlab, dimsHyperSlab, error)

      call H5Screate_simple_f(4, dimsHyperSlab, memspace, error)

      call H5Dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, data, dimsHyperSlab, error, memspace, filespace)

      call H5Sclose_f(memspace, error)
      call H5Sclose_f(filespace, error)
      call H5Dclose_f(dataset_id, error)

   end subroutine

end module xdmfAPI_m
