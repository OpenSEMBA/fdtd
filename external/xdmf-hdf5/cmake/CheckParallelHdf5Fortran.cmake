include(CMakePushCheckState)
include(CheckFortranSourceCompiles)

function(opensemba_check_parallel_hdf5_fortran result)
  set(${result} FALSE PARENT_SCOPE)

  if(NOT HDF5_IS_PARALLEL)
    return()
  endif()

  cmake_push_check_state(RESET)
  set(CMAKE_REQUIRED_INCLUDES
    ${HDF5_Fortran_INCLUDE_DIRS}
    ${MPI_Fortran_INCLUDE_DIRS}
  )
  set(CMAKE_REQUIRED_LIBRARIES
    ${HDF5_Fortran_LIBRARIES}
    ${HDF5_LIBRARIES}
    MPI::MPI_Fortran
  )
  check_fortran_source_compiles([=[
program parallel_hdf5_fortran_probe
  use hdf5
  use mpi
  implicit none
  integer :: error
  integer(HID_T) :: file_access_plist, transfer_plist

  call h5pcreate_f(H5P_FILE_ACCESS_F, file_access_plist, error)
  call h5pset_fapl_mpio_f(file_access_plist, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  call h5pcreate_f(H5P_DATASET_XFER_F, transfer_plist, error)
  call h5pset_dxpl_mpio_f(transfer_plist, H5FD_MPIO_COLLECTIVE_F, error)
end program parallel_hdf5_fortran_probe
]=] OPENSEMBA_PARALLEL_HDF5_FORTRAN_PROBE_WORKS SRC_EXT F90)
  cmake_pop_check_state()

  if(OPENSEMBA_PARALLEL_HDF5_FORTRAN_PROBE_WORKS)
    set(${result} TRUE PARENT_SCOPE)
  endif()
endfunction()
