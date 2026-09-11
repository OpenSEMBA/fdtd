cmake_minimum_required(VERSION 3.15)

project(
  xdmf-hdf5
  VERSION 0.1.0
  LANGUAGES Fortran
)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

if(CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
  set(_xdmf_hdf5_test_default ON)
else()
  set(_xdmf_hdf5_test_default OFF)
endif()

option(
  XDMF_HDF5_BUILD_TESTING
  "Build the XDMF/HDF5 conformance suite"
  ${_xdmf_hdf5_test_default}
)

find_package(HDF5 REQUIRED COMPONENTS Fortran)

option(
  XDMF_HDF5_ENABLE_MPI
  "Build the external writer MPI API"
  OFF
)

set(XDMF_HDF5_MPI_ENABLED FALSE)
if(XDMF_HDF5_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS Fortran)
  set(XDMF_HDF5_MPI_ENABLED TRUE)
endif()

set(XDMF_HDF5_PARALLEL_AVAILABLE FALSE)
if(XDMF_HDF5_MPI_ENABLED AND HDF5_IS_PARALLEL)
  include("${CMAKE_CURRENT_LIST_DIR}/cmake/CheckParallelHdf5Fortran.cmake")
  xdmf_check_parallel_hdf5_fortran(XDMF_HDF5_PARALLEL_AVAILABLE)
endif()

add_library(xdmf_hdf5
  src/xdmf_model.F90
  src/xdmf_hdf5_backend.F90
  src/xdmf_xml.F90
  src/xdmf_hdf5.F90
)
add_library(XDMF::HDF5 ALIAS xdmf_hdf5)

set_target_properties(xdmf_hdf5 PROPERTIES
  EXPORT_NAME HDF5
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  Fortran_STANDARD 2018
  Fortran_STANDARD_REQUIRED YES
  Fortran_EXTENSIONS OFF
)

target_include_directories(xdmf_hdf5
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/xdmf-hdf5>
)
if(NOT TARGET HDF5::HDF5)
  add_library(HDF5::HDF5 INTERFACE IMPORTED)
  set_target_properties(HDF5::HDF5 PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${HDF5_Fortran_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${HDF5_Fortran_LIBRARIES}"
  )
endif()
target_link_libraries(xdmf_hdf5 PUBLIC HDF5::HDF5)
if(XDMF_HDF5_MPI_ENABLED)
  target_link_libraries(xdmf_hdf5 PUBLIC MPI::MPI_Fortran)
  target_compile_definitions(xdmf_hdf5 PRIVATE XDMF_HDF5_WITH_MPI)
endif()
if(XDMF_HDF5_PARALLEL_AVAILABLE)
  target_compile_definitions(xdmf_hdf5 PRIVATE XDMF_HDF5_WITH_PARALLEL_HDF5)
  target_compile_definitions(xdmf_hdf5 PUBLIC XDMF_HDF5_PARALLEL_AVAILABLE)
endif()

if(XDMF_HDF5_BUILD_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

install(
  TARGETS xdmf_hdf5
  EXPORT XDMFHdf5Targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(
  DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/xdmf-hdf5"
  FILES_MATCHING PATTERN "*.mod"
)

configure_package_config_file(
  cmake/XDMFHdf5Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/XDMFHdf5Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/XDMFHdf5"
)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/XDMFHdf5ConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  EXPORT XDMFHdf5Targets
  FILE XDMFHdf5Targets.cmake
  NAMESPACE XDMF::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/XDMFHdf5"
)
install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/XDMFHdf5Config.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/XDMFHdf5ConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/XDMFHdf5"
)
install(
  FILES README.md LICENSE
  DESTINATION "${CMAKE_INSTALL_DATAROOTDIR}/doc/xdmf-hdf5"
)
