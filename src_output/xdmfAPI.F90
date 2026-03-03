module mod_xdmfAPI
   use HDF5
   implicit none

   ! HDF5 constants
   private
   integer, parameter :: dp = kind(1.0d0)
   public :: dp
   public :: xdmf_create_file
   public :: xdmf_write_timestep
   public :: xdmf_write_step
   public :: xdmf_close_file
   public :: create_h5_file
   public :: h5_close_file
   public :: write_dataset
   public :: xdmf_write_step_header
   public :: xdmf_write_attribute
   public :: xdmf_write_step_footer
   public :: init_extendable_2d_dataset
   public :: append_rows_dataset

   interface write_dataset
      module procedure write_1d_dataset
      module procedure write_2d_dataset
      module procedure write_3d_dataset
   end interface

contains
   subroutine xdmf_create_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit, ierr

      open (newunit=unit, file=filename, status='replace', action='write', &
            form='formatted', iostat=ierr)
      if (ierr /= 0) stop "Cannot create XDMF file"

      write (unit, '(A)') '<?xml version="1.0" ?>'
      write (unit, '(A)') '<Xdmf Version="3.0">'
      write (unit, '(A)') '  <Domain>'
      write (unit, '(A)') '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
      close (unit)
   end subroutine xdmf_create_file

   subroutine xdmf_write_timestep(filename, time, dataset_path, grid_name, dims)
      implicit none
      character(len=*), intent(in) :: filename
      real(dp), intent(in) :: time
      character(len=*), intent(in) :: dataset_path
      character(len=*), intent(in) :: grid_name
      integer, dimension(:), intent(in) :: dims
      integer :: unit
      character(len=50) :: dim_str1, dim_str2

      ! Convert dimensions to string
      write (dim_str1, '(I0," ",I0," ",I0)') dims(3), dims(2), dims(1)
      write (dim_str2, '(I0," ",I0," ",I0)') dims(1), dims(2), dims(3)

      open (newunit=unit, file=filename, status='old', action='write', position='append')

      write (unit, '(A)') '      <Grid Name="'//trim(grid_name)//'" GridType="Uniform">'
      write (unit, '(A,F8.4)') '        <Time Value="', time, '"/>'
      write (unit, '(A)') '        <Topology TopologyType="3DRectMesh" Dimensions="'//trim(dim_str1)//'"/>'
      write (unit, '(A)') '        <Geometry GeometryType="Origin_DxDyDz">'
      write (unit, '(A)') '          <DataItem Dimensions="3" NumberType="Float" Precision="8" Format="XML">0 0 0</DataItem>'
      write (unit, '(A)') '          <DataItem Dimensions="3" NumberType="Float" Precision="8" Format="XML">1 1 1</DataItem>'
      write (unit, '(A)') '        </Geometry>'
      write (unit, '(A)') '        <Attribute Name="Efield" AttributeType="Scalar" Center="Node">'
      write (unit, '(A)'     ) '          <DataItem Dimensions="'//trim(dim_str2)//'" NumberType="Float" Precision="8" Format="HDF">'//trim(dataset_path)//'</DataItem>'
      write (unit, '(A)') '        </Attribute>'
      write (unit, '(A)') '      </Grid>'

      close (unit)
   end subroutine xdmf_write_timestep

   subroutine xdmf_write_step_header(filename, t, time, npoints)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: t
      real(8), intent(in) :: time
      integer, intent(in) :: npoints
      integer :: unit
      character(len=20) :: ts

      write (ts, '(I0)') t

      open (newunit=unit, file=filename, position="append")

      write (unit, '(A)') '  <Grid Name="Step'//trim(ts)//'" GridType="Uniform">'
      write (unit, '(A,F12.6,A)') '    <Time Value="', time, '"/>'
      write (unit, '(A,I0,A)') '    <Topology TopologyType="Polyvertex" NumberOfElements="', npoints, '"/>'

      write (unit, '(A)') '    <Geometry GeometryType="XYZ">'
      write (unit, '(A)') '      <DataItem Dimensions="'//trim(ts)//'" Format="HDF">data.h5:/coords</DataItem>'
      write (unit, '(A)') '    </Geometry>'

      close (unit)
   end subroutine xdmf_write_step_header

   subroutine xdmf_write_attribute(filename, t, attr_name, h5_dataset)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: t
      character(len=*), intent(in) :: attr_name
      character(len=*), intent(in) :: h5_dataset
      integer :: unit
      character(len=20) :: ts

      write (ts, '(I0)') t
      open (newunit=unit, file=filename, position="append")

      write (unit, '(A)') '    <Attribute Name="'//trim(attr_name)//'" Center="Node">'
   write(unit, '(A,I0,A)') '      <DataItem Dimensions="'//trim(ts)//'" Format="HDF">data.h5:/'//trim(h5_dataset)//'['//trim(ts)//',:]</DataItem>'
      write (unit, '(A)') '    </Attribute>'

      close (unit)
   end subroutine xdmf_write_attribute

   subroutine xdmf_write_step_footer(filename)
      character(len=*), intent(in) :: filename
      integer :: unit

      open (newunit=unit, file=filename, position="append")
      write (unit, '(A)') '  </Grid>'
      close (unit)
   end subroutine xdmf_write_step_footer

   subroutine xdmf_write_step(filename, t, time, npoints)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: t
      real(8), intent(in) :: time
      integer, intent(in) :: npoints
      integer :: unit
      character(len=20) :: ts

      write (ts, '(I0)') t

      open (newunit=unit, file=filename, position="append")

      write (unit, '(A)') '  <Grid Name="Step'//trim(ts)//'" GridType="Uniform">'
      write (unit, '(A,F12.6,A)') '    <Time Value="', time, '"/>'

      write (unit, '(A,I0,A)') '    <Topology TopologyType="Polyvertex" NumberOfElements="', npoints, '"/>'

      write (unit, '(A)') '    <Geometry GeometryType="XYZ">'
      write (unit, '(A)') '      <DataItem Dimensions="'//trim(ts)//'" Format="HDF">data.h5:/coords</DataItem>'
      write (unit, '(A)') '    </Geometry>'

      write (unit, '(A)') '    <Attribute Name="Ex" Center="Node">'
      write (unit, '(A,I0,A)') '      <DataItem Dimensions="'//trim(ts)//'" Format="HDF">data.h5:/Ex['//trim(ts)//',:]</DataItem>'
      write (unit, '(A)') '    </Attribute>'

      write (unit, '(A)') '    <Attribute Name="Ey" Center="Node">'
      write (unit, '(A,I0,A)') '      <DataItem Dimensions="'//trim(ts)//'" Format="HDF">data.h5:/Ey['//trim(ts)//',:]</DataItem>'
      write (unit, '(A)') '    </Attribute>'

      write (unit, '(A)') '    <Attribute Name="Ez" Center="Node">'
      write (unit, '(A,I0,A)') '      <DataItem Dimensions="'//trim(ts)//'" Format="HDF">data.h5:/Ez['//trim(ts)//',:]</DataItem>'
      write (unit, '(A)') '    </Attribute>'

      write (unit, '(A)') '  </Grid>'

      close (unit)

   end subroutine

   subroutine xdmf_close_file(filename)
      character(len=*), intent(in) :: filename
      integer :: unit

      open (newunit=unit, file=filename, status='old', action='write', position='append')
      write (unit, '(A)') '    </Grid>'
      write (unit, '(A)') '  </Domain>'
      write (unit, '(A)') '</Xdmf>'
      close (unit)
   end subroutine xdmf_close_file

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
