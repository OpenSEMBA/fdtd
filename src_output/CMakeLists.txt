message(STATUS "Creating build system for output")

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

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

	target_link_libraries(fdtd-output PRIVATE
    semba-components
    semba-types
    fdtd-utils
    ${HDF5_LIBRARIES} 
    ${HDF5_HL_LIBRARIES}
	)


add_library(vtkAPI
  "vtkAPI.F90"
)

add_library(xdmfAPI
  "xdmfAPI.F90"
)

	target_link_libraries(xdmfAPI PRIVATE
    ${HDF5_LIBRARIES} 
    ${HDF5_HL_LIBRARIES}
	)

target_link_libraries(vtkAPI
 fdtd-utils
)

