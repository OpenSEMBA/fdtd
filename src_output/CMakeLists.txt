message(STATUS "Creating build system for output")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_library(fdtd-output
  "output.F90"
  "outputTypes.F90"
  "outputBinary.F90"
  "outputVisualisation.F90"
  "outputDecomposition.F90"
  "outputCollective.F90"
  "outputTransport.F90"
  "domain.F90"
  "outputUtils.F90"
  "volumicProbeUtils.F90"
  "lineProbeOutput.F90"
  "pointProbeOutput.F90"
  "wireProbeOutput.F90"
  "bulkProbeOutput.F90"
  "movieProbeOutput.F90"
  "frequencySliceProbeOutput.F90"
  "farFieldProbeOutput.F90"
  "mapVTKOutput.F90"
)

target_link_libraries(fdtd-output PRIVATE
  semba-components
  semba-types
  fdtd-utils
  XDMF::HDF5
)

if(SEMBA_FDTD_ENABLE_MPI)
  target_link_libraries(fdtd-output PRIVATE MPI::MPI_Fortran)
endif()


add_library(vtkAPI
  "vtkAPI.F90"
)

target_link_libraries(vtkAPI
 fdtd-utils
)
