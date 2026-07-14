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
