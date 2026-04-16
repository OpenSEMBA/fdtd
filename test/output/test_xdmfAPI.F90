integer function test_create_h5_file() bind(c) result(err)

   use HDF5
   use mod_xdmfAPI
   use mod_assertionTools
   use directoryUtils_m
   implicit none

   character(len=14), parameter :: folder='testing_folder'
   character(len=1024) :: file
   integer(HID_T) :: file_id
   integer :: error
   logical :: exists
   integer :: test_err = 0

   call create_folder(folder, error)
   file = join_path(folder, "test_api.h5")

   call H5open_f(error)

   call create_h5_file(trim(adjustl(file)), file_id)

   test_err = test_err + assert_true(file_id > 0, "create_h5_file returned invalid id")

   call H5Fclose_f(file_id, error)

   inquire(file=trim(adjustl(file)), exist=exists)
   test_err = test_err + assert_true(exists, "HDF5 file was not created")

   call H5close_f(error)

   err = test_err
   call remove_folder(folder, error)
end function 

integer function test_write_1d_dataset() bind(c) result(err)

   use HDF5
   use mod_xdmfAPI
   use mod_assertionTools
   use directoryUtils_m
   implicit none

   character(len=14), parameter :: folder='testing_folder'
   character(len=1024) :: file
   integer(HID_T) :: file_id
   real(8), allocatable :: data(:)
   integer :: error
   integer :: test_err = 0
   integer :: n, i

   call create_folder(folder, error)
   file = join_path(folder, "test_api.h5")

   call H5open_f(error)

   n = 10
   allocate(data(n))
   data = [(real(i,8), i=1,n)]

   call create_h5_file(trim(adjustl(file)), file_id)

   call write_dataset(file_id, "/data1d", data)

   test_err = test_err + assert_integer_equal(size(data), 10, "1D dataset size mismatch")

   call H5Fclose_f(file_id, error)

   deallocate(data)

   call H5close_f(error)

   err = test_err
   call remove_folder(folder, error)
end function 

integer function test_write_2d_dataset() bind(c) result(err)

   use HDF5
   use mod_xdmfAPI
   use mod_assertionTools
   use directoryUtils_m
   implicit none

   character(len=14), parameter :: folder='testing_folder'
   character(len=1024) :: file
   integer(HID_T) :: file_id
   real(8), allocatable :: data(:,:)
   integer :: error
   integer :: test_err = 0

   call create_folder(folder, error)
   file = join_path(folder, "test_api.h5")

   call H5open_f(error)

   allocate(data(4,5))
   data = 1.0d0

   call create_h5_file(trim(adjustl(file)), file_id)

   call write_dataset(file_id, "/data2d", data)

   test_err = test_err + assert_integer_equal(size(data,1),4,"2D dim1 mismatch")
   test_err = test_err + assert_integer_equal(size(data,2),5,"2D dim2 mismatch")

   call H5Fclose_f(file_id, error)

   deallocate(data)

   call H5close_f(error)

   err = test_err
   call remove_folder(folder, error)
end function 

integer function test_write_3d_dataset() bind(c) result(err)

   use HDF5
   use mod_xdmfAPI
   use mod_assertionTools
   use directoryUtils_m
   implicit none

   character(len=14), parameter :: folder='testing_folder'
   character(len=1024) :: file
   integer(HID_T) :: file_id
   real(8), allocatable :: data(:,:,:)
   integer :: error
   integer :: test_err = 0

   call create_folder(folder, error)
   file = join_path(folder, "test_api.h5")

   call H5open_f(error)

   allocate(data(3,3,3))
   data = 2.0d0

   call create_h5_file(trim(adjustl(file)), file_id)

   call write_dataset(file_id, "/data3d", data)

   test_err = test_err + assert_integer_equal(size(data,1),3,"3D dim1 mismatch")
   test_err = test_err + assert_integer_equal(size(data,2),3,"3D dim2 mismatch")
   test_err = test_err + assert_integer_equal(size(data,3),3,"3D dim3 mismatch")

   call H5Fclose_f(file_id, error)

   deallocate(data)

   call H5close_f(error)

   err = test_err

   call remove_folder(folder, error)
end function 

integer function test_xdmf_file_creation() bind(c) result(err)

   use mod_xdmfAPI
   use mod_assertionTools
   use directoryUtils_m
   implicit none

   character(len=14), parameter :: folder='testing_folder'
   character(len=1024) :: file
   integer :: test_err = 0
   integer :: error, unit
   logical :: exists
   integer :: dims(3)

   call create_folder(folder, error)
   file = join_path(folder, "test_api.xdmf")

   dims = [4,4,4]

   
   open(newunit=unit, file=trim(file), position='append')
   call xdmf_write_header_file(unit, 'movieProbe')

   call xdmf_create_grid_step_info(unit,"step0",0.0,"data.h5",dims(1)*dims(2)*dims(3))

   call xdmf_write_attribute(unit,"Efield")

   call xdmf_write_h5_data_item(unit,"data.h5","/Efield","4 4 4")

   call xdmf_close_data_item(unit)
   call xdmf_close_attribute(unit)
   call xdmf_close_grid(unit)

   call xdmf_write_footer_file(unit)
   close(unit)


   inquire(file=trim(file), exist=exists)

   test_err = test_err + assert_true(exists, "XDMF file not created")

   err = test_err
   call remove_folder(folder, error)

end function

integer function test_xdmf_file_with_h5data() bind(c) result(err)
    use HDF5
    use mod_xdmfAPI
    use mod_assertionTools
    use directoryUtils_m
    implicit none

    character(len=20), parameter :: folder = "testing_folder"
    character(len=1024) :: xdmf_file, h5_file
    integer :: test_err = 0
    integer :: error, t, unit
    logical :: exists
    integer :: dims(3)
    integer(HID_T) :: file_id
    real(dp), allocatable :: Efield(:,:), coords(:,:)
    real(dp) :: time
    character(len=20) :: ts
    integer :: i,j, npoints

    call create_folder(folder, error)

    xdmf_file = join_path(folder, "test_api.xdmf")
    h5_file   = join_path(folder, "data.h5")

    call H5open_f(error)

    ! Create HDF5 file
    call create_h5_file(trim(h5_file), file_id)

    dims = [3,3,1]      ! 2D grid stored as 3D with depth 1
    npoints = dims(1)*dims(2)

    ! Allocate and write coords data: shape (npoints,3)
    allocate(coords(npoints,3))
    do j = 1, dims(2)
       do i = 1, dims(1)
          coords((j-1)*dims(1)+i,1) = real(i-1, dp)  ! X
          coords((j-1)*dims(1)+i,2) = real(j-1, dp)  ! Y
          coords((j-1)*dims(1)+i,3) = 0.0_dp         ! Z
       end do
    end do
    call write_dataset(file_id,"coords",coords)
    deallocate(coords)

    ! Create XDMF file
    open(newunit=unit,file=trim(xdmf_file),position='append')
    call xdmf_write_header_file(unit)

    do t = 1, 5
        time = real(t-1,dp)*0.1_dp

        allocate(Efield(dims(1),dims(2)))
        do j=1,dims(2)
           do i=1,dims(1)
              Efield(i,j) = i + j + t - 1
           end do
        end do

        write(ts,'("Efield_",I0)') t

        ! Write timestep data
        call write_dataset(file_id,trim(ts),Efield)

        ! XDMF grid
        call xdmf_create_grid_step_info(unit,trim(ts),real(time),trim(h5_file),npoints)
        call xdmf_write_attribute(unit,"Efield")
        call xdmf_write_h5_data_item(unit,trim(h5_file),"/"//trim(ts),"3 3 1")
        call xdmf_close_data_item(unit)
        call xdmf_close_attribute(unit)
        call xdmf_close_grid(unit)

        deallocate(Efield)
    end do
    call xdmf_write_footer_file(unit)
    close(unit)
    call H5Fclose_f(file_id,error)

    inquire(file=trim(xdmf_file),exist=exists)
    test_err = test_err + assert_true(exists,"XDMF file not created")

    call remove_folder(folder,error)
    call H5close_f(error)

    err = test_err
end function