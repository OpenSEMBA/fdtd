cmake_minimum_required(VERSION 3.15)

project(xdmf-hdf5-example LANGUAGES Fortran)

find_package(OpenSEMBAXdmfHdf5 CONFIG REQUIRED)

add_executable(xdmf-hdf5-example basic_writer.F90)
target_link_libraries(xdmf-hdf5-example PRIVATE OpenSEMBA::xdmf-hdf5)
set_target_properties(xdmf-hdf5-example PROPERTIES
  Fortran_STANDARD 2008
  Fortran_STANDARD_REQUIRED YES
)
