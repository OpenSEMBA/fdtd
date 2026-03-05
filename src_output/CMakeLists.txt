add_library(fdtd-output 
  "output.F90"
  "outputTypes.F90"
  "domain.F90"
  "outputUtils.F90"
  "volumicProbeUtils.F90"
  "pointProbeOutput.F90"
  "wireProbeOutput.F90"
  "bulkProbeOutput.F90"
  "movieProbeOutput.F90"
  "frequencySliceProbeOutput.F90"
  "farFieldProbeOutput.F90"
  "mapVTKOutput.F90"
)
add_library(vtkAPI
  "vtkAPI.F90"
)

add_library(xdmfAPI
  "xdmfAPI.F90"
)

target_include_directories(xdmfAPI PRIVATE ${HDF5_INCLUDE_DIRS})

target_link_libraries(vtkAPI
 fdtd-utils
)

