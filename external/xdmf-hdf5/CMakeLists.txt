cmake_minimum_required(VERSION 3.15)

project(
  opensemba-xdmf-hdf5
  VERSION 0.1.0
  LANGUAGES Fortran
)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

if(DEFINED SEMBA_FDTD_ENABLE_TEST)
  set(_xdmf_hdf5_test_default ${SEMBA_FDTD_ENABLE_TEST})
else()
  if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    set(_xdmf_hdf5_test_default ON)
  else()
    set(_xdmf_hdf5_test_default OFF)
  endif()
endif()

option(
  XDMF_HDF5_BUILD_TESTING
  "Build the XDMF/HDF5 conformance suite"
  ${_xdmf_hdf5_test_default}
)

find_package(HDF5 REQUIRED COMPONENTS Fortran)

if(DEFINED SEMBA_FDTD_ENABLE_MPI)
  set(_xdmf_hdf5_mpi_default ${SEMBA_FDTD_ENABLE_MPI})
else()
  set(_xdmf_hdf5_mpi_default OFF)
endif()
option(
  XDMF_HDF5_ENABLE_MPI
  "Build the external writer MPI API"
  ${_xdmf_hdf5_mpi_default}
)

set(XDMF_HDF5_MPI_ENABLED FALSE)
if(XDMF_HDF5_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  set(XDMF_HDF5_MPI_ENABLED TRUE)
endif()

set(XDMF_HDF5_PARALLEL_AVAILABLE FALSE)
if(XDMF_HDF5_MPI_ENABLED AND SEMBA_FDTD_HAVE_PARALLEL_HDF5_FORTRAN)
  set(XDMF_HDF5_PARALLEL_AVAILABLE TRUE)
elseif(XDMF_HDF5_MPI_ENABLED AND HDF5_IS_PARALLEL)
  include(CMakePushCheckState)
  include(CheckFortranSourceCompiles)
  cmake_push_check_state(RESET)
  set(CMAKE_REQUIRED_INCLUDES
    ${HDF5_Fortran_INCLUDE_DIRS}
    ${MPI_Fortran_INCLUDE_DIRS}
  )
  set(CMAKE_REQUIRED_LIBRARIES
    ${HDF5_Fortran_LIBRARIES}
    MPI::MPI_Fortran
  )
  check_fortran_source_compiles([=[
program xdmf_hdf5_mpi_probe
  use hdf5
  use mpi
  implicit none
  integer :: error
  integer(HID_T) :: file_access_plist, transfer_plist

  call h5pcreate_f(H5P_FILE_ACCESS_F, file_access_plist, error)
  call h5pset_fapl_mpio_f(file_access_plist, MPI_COMM_WORLD, MPI_INFO_NULL, error)
  call h5pcreate_f(H5P_DATASET_XFER_F, transfer_plist, error)
  call h5pset_dxpl_mpio_f(transfer_plist, H5FD_MPIO_COLLECTIVE_F, error)
end program xdmf_hdf5_mpi_probe
]=] XDMF_HDF5_MPI_PROBE_WORKS SRC_EXT F90)
  cmake_pop_check_state()
  if(XDMF_HDF5_MPI_PROBE_WORKS)
    set(XDMF_HDF5_PARALLEL_AVAILABLE TRUE)
  endif()
endif()

add_library(opensemba_xdmf_hdf5
  src/xdmf_model.F90
  src/xdmf_hdf5_backend.F90
  src/xdmf_xml.F90
  src/xdmf_hdf5.F90
)
add_library(OpenSEMBA::xdmf-hdf5 ALIAS opensemba_xdmf_hdf5)

set_target_properties(opensemba_xdmf_hdf5 PROPERTIES
  EXPORT_NAME xdmf-hdf5
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  Fortran_STANDARD 2008
  Fortran_STANDARD_REQUIRED YES
)

target_include_directories(opensemba_xdmf_hdf5
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/opensemba-xdmf-hdf5>
)
if(NOT TARGET HDF5::HDF5)
  add_library(HDF5::HDF5 INTERFACE IMPORTED)
  set_target_properties(HDF5::HDF5 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${HDF5_Fortran_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${HDF5_Fortran_LIBRARIES}"
  )
endif()
target_link_libraries(opensemba_xdmf_hdf5 PUBLIC HDF5::HDF5)
if(XDMF_HDF5_MPI_ENABLED)
  target_link_libraries(opensemba_xdmf_hdf5 PUBLIC MPI::MPI_Fortran)
  target_compile_definitions(opensemba_xdmf_hdf5 PRIVATE XDMF_HDF5_WITH_MPI)
endif()
if(XDMF_HDF5_PARALLEL_AVAILABLE)
  target_compile_definitions(opensemba_xdmf_hdf5 PRIVATE XDMF_HDF5_WITH_PARALLEL_HDF5)
endif()

if(XDMF_HDF5_BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

install(
  TARGETS opensemba_xdmf_hdf5
  EXPORT OpenSEMBAXdmfHdf5Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(
  DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/opensemba-xdmf-hdf5"
  FILES_MATCHING PATTERN "*.mod"
)

configure_package_config_file(
  cmake/OpenSEMBAXdmfHdf5Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/OpenSEMBAXdmfHdf5Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/OpenSEMBAXdmfHdf5"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/OpenSEMBAXdmfHdf5ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  EXPORT OpenSEMBAXdmfHdf5Targets
  FILE OpenSEMBAXdmfHdf5Targets.cmake
  NAMESPACE OpenSEMBA::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/OpenSEMBAXdmfHdf5"
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/OpenSEMBAXdmfHdf5Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenSEMBAXdmfHdf5ConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/OpenSEMBAXdmfHdf5"
)
install(
  FILES README.md LICENSE
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/doc/opensemba-xdmf-hdf5"
)
