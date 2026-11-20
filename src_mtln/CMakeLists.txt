message(STATUS "Creating build system for mtlnsolver")

project(mtln)
enable_language(C CXX Fortran)

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

add_subdirectory(interface)

# XSPICE code models are loaded at run time from the executable folder and
# shipped from bin/ in the release archives. When they are built from the
# ngspice submodule (cmake/NgspiceCodemodels.cmake) that build places them in
# bin/. The precompiled ones come from another ngspice build and are not ABI
# compatible with the linked library, so they are only a fallback.
if (NOT SEMBA_FDTD_BUILD_NGSPICE_CODEMODELS)
	if (WIN32)
		set(NGSPICE_CODEMODELS_OS windows)
	else()
		set(NGSPICE_CODEMODELS_OS linux)
	endif()
	file(GLOB NGSPICE_CODEMODELS "${CMAKE_CURRENT_SOURCE_DIR}/../testData/spinit/${NGSPICE_CODEMODELS_OS}/*.cm")
	file(COPY ${NGSPICE_CODEMODELS} DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
endif()

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
