message(STATUS "Creating build system for mtlnsolver")

project(mtln)
enable_language(C CXX Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_subdirectory(interface)

# Precompiled XSPICE code models, loaded at run time from the executable folder.
# They land in bin/, so they are shipped with the release archives.
if (WIN32)
	set(NGSPICE_CODEMODELS_OS windows)
else()
	set(NGSPICE_CODEMODELS_OS linux)
endif()
file(GLOB NGSPICE_CODEMODELS "${CMAKE_CURRENT_SOURCE_DIR}/../testData/spinit/${NGSPICE_CODEMODELS_OS}/*.cm")
file(COPY ${NGSPICE_CODEMODELS} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})

add_library(mtlnsolver
	"mtl.F90"
	"mtl_bundle.F90"
	"mtln_solver.F90"
	"multipolar_expansion.F90"
	"network.F90"
	"network_manager.F90"
	"utils.F90"
	"probes.F90"
	"generators.F90"
	"dispersive.F90"
	"rational_approximation.F90"
	"circuit.F90"
	"ngspice_interface.F90"
	"preprocess.F90"
)

target_link_libraries(mtlnsolver PRIVATE
	${LAPACK_LIB} 
	${BLAS_LIB} 
	fhash
	ngspice_interface
	semba-types
	semba-reports
	fdtd-utils
	${MPI_Fortran_LIBRARIES}
)

if(TARGET jsonfortran)
	target_link_libraries(mtlnsolver PRIVATE jsonfortran)
endif()
